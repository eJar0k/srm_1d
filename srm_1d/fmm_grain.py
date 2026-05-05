"""
fmm_grain.py — Fast Marching Method grain support (openMotor bridge)
=====================================================================

Bridges srm_1d to openMotor's FmmGrain machinery (skfmm-based regression
maps) for arbitrary cross-sections: Finocyl, Star, Moonburner, X, D, C,
RodTube, Custom.

WHAT THIS MODULE DOES
    1. Lazily adds the local openMotor checkout to sys.path so
       `import motorlib` and `import mathlib` work.
    2. Runs openMotor's FmmGrain pipeline (initGeometry → generateCoreMap →
       generateRegressionMap) on a populated grain instance.
    3. Samples the resulting regression map at fine resolution and packs
       it into an `FmmTable` of (regression_depth, perimeter, port_area)
       arrays suitable for srm_1d's @njit hot loop.

WHAT THIS MODULE DOES NOT
    - Reimplement skfmm or openMotor's FmmGrain — we call upstream.
    - Run openMotor's solver. We use the cross-section processor only.
    - Modify the time loop or geometry compilation directly. That's
      Phase 2, in grain_geometry.py.

DESIGN NOTE
    FMM tables are per-grain (one set per FmmTable). Per-cell regression
    in the hot loop is preserved — the Ma erosive model produces axially
    varying burn rates, and per-cell regression captures that variation
    while the heavy data (the regressionMap and sampled tables) stays
    shared at the grain level.
"""

import sys
import types
from dataclasses import dataclass
from pathlib import Path
import numpy as np

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


# Cached so we only do path manipulation once per process
_OPENMOTOR_PATH_CONFIGURED = False


# ================================================================
# Pure-Python (Numba-JIT) replacement for openMotor's Cython
# `mathlib._find_perimeter_cy._get_perimeter`. Avoids requiring the
# user to compile a Cython extension (MSVC build tools on Windows,
# etc.). The algorithm is a verbatim port of openMotor's marching-
# squares perimeter computation; only the language differs.
# ================================================================

@njit(cache=True)
def _marching_squares_perimeter(arr, level):
    """
    Marching-squares perimeter of the iso-contour at `level` in a 2D
    float64 array. Mirrors openMotor's
    `mathlib._find_perimeter_cy._get_perimeter`, including the
    "skip squares outside the motor radius" optimization.

    Returns the perimeter in array-index distance units (one unit =
    one pixel-width); the caller (openMotor's getCorePerimeter)
    multiplies by `mapToLength` to get meters.
    """
    perimeter = 0.0
    rows, cols = arr.shape[0], arr.shape[1]
    half_r = rows / 2.0
    half_c = cols / 2.0
    radius_thresh = half_r - 3.0

    for r0 in range(3, rows - 4):
        for c0 in range(3, cols - 4):
            r1 = r0 + 1
            c1 = c0 + 1

            ul = arr[r0, c0]
            ur = arr[r0, c1]
            ll = arr[r1, c0]
            lr = arr[r1, c1]

            sc = 0
            if ul > level:
                sc += 1
            if ur > level:
                sc += 2
            if ll > level:
                sc += 4
            if lr > level:
                sc += 8

            if sc == 0 or sc == 15:
                continue

            # Skip squares outside the motor radius (cylindrical mask).
            dist_r = r0 + 0.5 - half_r
            dist_c = c0 + 0.5 - half_c
            if radius_thresh < (dist_r * dist_r + dist_c * dist_c) ** 0.5:
                continue

            # _get_fraction (linear interp): zero-safe when both ends equal
            top = 0.0 if ul == ur else (level - ul) / (ur - ul)
            bottom = 0.0 if ll == lr else (level - ll) / (lr - ll)
            left = 0.0 if ll == ul else (level - ll) / (ul - ll)
            right = 0.0 if lr == ur else (level - lr) / (ur - lr)

            add = 0.0
            if sc == 1 or sc == 14:
                add = (top * top + (1.0 - left) ** 2) ** 0.5
            elif sc == 2 or sc == 13:
                add = ((1.0 - top) ** 2 + (1.0 - right) ** 2) ** 0.5
            elif sc == 3 or sc == 12:
                add = ((right - left) ** 2 + 1.0) ** 0.5
            elif sc == 4 or sc == 11:
                add = (left * left + bottom * bottom) ** 0.5
            elif sc == 5 or sc == 10:
                add = ((top - bottom) ** 2 + 1.0) ** 0.5
            elif sc == 6:
                add = (((1.0 - top) ** 2 + (1.0 - right) ** 2) ** 0.5
                       + (left * left + bottom * bottom) ** 0.5)
            elif sc == 7 or sc == 8:
                add = ((1.0 - bottom) ** 2 + right * right) ** 0.5
            elif sc == 9:
                add = ((top * top + (1.0 - left) ** 2) ** 0.5
                       + ((1.0 - bottom) ** 2 + right * right) ** 0.5)
            perimeter += add

    return perimeter


def _get_perimeter_replacement(image, level, vertex_connect_high, returning_contours):
    """
    Python-level wrapper matching the Cython
    `_get_perimeter(image, level, vertex_connect_high, returning_contours)`
    signature. Returns (perimeter, segments) where segments is always
    an empty list — openMotor's `getCorePerimeter` only uses [0]
    (the perimeter), and we don't need contour reassembly for the
    solver path. If `returning_contours` is requested (e.g. for
    plotting via openMotor's `getRegressionData`), this raises.
    """
    if returning_contours:
        raise NotImplementedError(
            "srm_1d's pure-Python find_perimeter shim doesn't support "
            "returning_contours=True. Build openMotor's Cython extension "
            "(`python setup.py build_ext --inplace` in the openMotor dir) "
            "if you need contour visualizations."
        )
    arr = np.ascontiguousarray(np.asarray(image, dtype=np.float64))
    return float(_marching_squares_perimeter(arr, float(level))), []


def _install_fake_find_perimeter_cy():
    """
    Inject a fake `mathlib._find_perimeter_cy` module exposing our
    Numba `_get_perimeter` so openMotor's mathlib package can import
    without requiring the Cython build. Idempotent.
    """
    mod_name = 'mathlib._find_perimeter_cy'
    if mod_name in sys.modules:
        return
    fake = types.ModuleType(mod_name)
    fake._get_perimeter = _get_perimeter_replacement
    sys.modules[mod_name] = fake


def _setup_openmotor_path():
    """
    Add the local openMotor checkout to sys.path so `motorlib` and
    `mathlib` are importable. Inject the find_perimeter shim before
    `mathlib` imports its missing Cython dependency. Idempotent.

    Expected layout (per reference_openmotor_source memory):
        Erosive Burning Solver/
            openMotor/
                openMotor/      ← inner package dir, contains motorlib/
                    motorlib/
                    mathlib/
            Claude Testing/
                v0.4 Adapter Layer/
                    srm_1d/
                        fmm_grain.py    ← __file__
    """
    global _OPENMOTOR_PATH_CONFIGURED
    if _OPENMOTOR_PATH_CONFIGURED:
        return

    here = Path(__file__).resolve().parent             # srm_1d/
    # ../../ → "v0.4 Adapter Layer/"
    # ../../../ → "Claude Testing/"
    # ../../../../ → "Erosive Burning Solver/"
    om_root = here.parent.parent.parent / "openMotor" / "openMotor"

    if not om_root.is_dir():
        raise ImportError(
            f"openMotor checkout not found at expected path:\n  {om_root}\n"
            "Clone it next to this project:\n"
            "    cd 'Erosive Burning Solver'\n"
            "    git clone https://github.com/reilleya/openMotor"
        )

    om_str = str(om_root)
    if om_str not in sys.path:
        sys.path.insert(0, om_str)

    # Inject the Cython-replacement BEFORE mathlib tries to import it.
    _install_fake_find_perimeter_cy()

    # Verify both required packages import.
    try:
        import motorlib  # noqa: F401
        import mathlib   # noqa: F401
    except ImportError as e:
        raise ImportError(
            f"Found openMotor at {om_root}, but `import motorlib`/`import mathlib` "
            f"failed: {e}\n"
            "Make sure scikit-fmm is installed: pip install scikit-fmm"
        ) from e

    _OPENMOTOR_PATH_CONFIGURED = True


# ================================================================
# FmmTable: tabulated grain shape data ready for the hot loop
# ================================================================

@dataclass
class FmmTable:
    """
    Pre-tabulated FMM grain shape data — perimeter and port area as
    1D functions of regression depth. Built once at simulation start
    by sampling an openMotor FmmGrain's regression map.

    Attributes
    ----------
    reg_depth : ndarray (T,)
        Regression depths sampled, [m]. Uniform grid from 0 to wall_web.
    perimeter : ndarray (T,)
        Burning core perimeter at each regression depth, [m].
        Used as `C_burn` in the srm_1d solver.
    port_area : ndarray (T,)
        Open cross-sectional area at each regression depth, [m²].
        Used as `A_port`. Includes the casting-tube fully-open value
        once regression > wall_web (after burnout).
    wall_web : float
        Maximum regression depth for this grain, [m]. Beyond this
        the grain is fully consumed (port_area = casting-tube area).
    grain_outer_diameter : float
        Casting tube ID, [m]. (openMotor: `props['diameter']`.)
    grain_length : float
        Axial length, [m]. (openMotor: `props['length']`.)
    inhibited_fwd : bool
    inhibited_aft : bool
        End-face inhibition flags. (Translated from openMotor's
        `inhibitedEnds` enum: Top→fwd, Bottom→aft.)
    geom_name : str
        openMotor `geomName` of the source grain (e.g. 'Finocyl').
        Informational only.
    """
    reg_depth: np.ndarray
    perimeter: np.ndarray
    port_area: np.ndarray
    wall_web: float
    grain_outer_diameter: float
    grain_length: float
    inhibited_fwd: bool
    inhibited_aft: bool
    geom_name: str

    @property
    def n_samples(self) -> int:
        return len(self.reg_depth)

    @property
    def initial_perimeter(self) -> float:
        """C_burn at regress=0."""
        return float(self.perimeter[0])

    @property
    def initial_port_area(self) -> float:
        """A_port at regress=0."""
        return float(self.port_area[0])


# Maps openMotor's 'inhibitedEnds' enum to (inhibit_fwd, inhibit_aft)
# (Top = forward / head end, Bottom = aft / nozzle end.)
_INHIBIT_MAP = {
    'Neither': (False, False),
    'Top':     (True,  False),
    'Bottom':  (False, True),
    'Both':    (True,  True),
}


def from_openmotor(om_grain, map_dim: int = 1001) -> FmmTable:
    """
    Run openMotor's FmmGrain pipeline on the supplied grain instance and
    extract its perimeter / port-area tables as numpy arrays.

    Parameters
    ----------
    om_grain : openMotor FmmGrain (e.g. Finocyl, Star, Moonburner)
        Properties (diameter, length, inhibitedEnds, and grain-specific
        ones like coreDiameter/numFins/finWidth/...) must already be set.
    map_dim : int
        FMM regression-map resolution. openMotor default is 1001. Higher
        is more accurate but quadratically slower (skfmm.distance is
        O(mapDim²)). Sim-level config knob.

    Returns
    -------
    FmmTable

    Notes
    -----
    The regression map is sampled at one entry per FMM pixel-width in
    regression depth (matches openMotor's internal sampling density),
    plus a final sample slightly beyond burnout to make the table cover
    [0, wall_web] inclusively for safe interpolation at the boundary.
    """
    _setup_openmotor_path()

    # Late import — only available after _setup_openmotor_path()
    from motorlib.grain import FmmGrain
    if not isinstance(om_grain, FmmGrain):
        raise TypeError(
            f"Expected an openMotor FmmGrain subclass, got {type(om_grain).__name__}"
        )

    # The full simulationSetup wants a config object; do its three steps directly.
    om_grain.initGeometry(map_dim)
    om_grain.generateCoreMap()
    om_grain.generateRegressionMap()

    wall_web = float(om_grain.wallWeb)
    if wall_web <= 0.0:
        raise ValueError(
            f"FmmGrain {type(om_grain).__name__} returned wallWeb={wall_web}; "
            "check grain properties (likely zero core or zero web)."
        )

    # Sample at ~one point per FMM pixel-width of regression. openMotor's
    # internal sampling uses int(maxDist*mapDim)+2 entries; we mirror.
    n_samples = max(int(wall_web * map_dim / om_grain.props['diameter'].getValue()) + 2,
                    64)
    reg_depth = np.linspace(0.0, wall_web, n_samples)

    # Sample perimeter and port area directly from openMotor's API.
    perimeter = np.empty(n_samples)
    port_area = np.empty(n_samples)
    casting_area = np.pi / 4.0 * om_grain.props['diameter'].getValue() ** 2
    for i, r in enumerate(reg_depth):
        if r >= wall_web:
            # At/past burnout, perimeter and face_area both go to 0; the
            # cell is fully consumed and the port equals the casting tube.
            perimeter[i] = 0.0
            port_area[i] = casting_area
        else:
            perimeter[i] = float(om_grain.getCorePerimeter(r))
            port_area[i] = float(om_grain.getPortArea(r))

    inh_str = om_grain.props['inhibitedEnds'].getValue()
    if inh_str not in _INHIBIT_MAP:
        raise ValueError(
            f"Unknown inhibitedEnds value {inh_str!r}. "
            f"Expected one of {list(_INHIBIT_MAP)}."
        )
    inh_fwd, inh_aft = _INHIBIT_MAP[inh_str]

    return FmmTable(
        reg_depth=reg_depth,
        perimeter=perimeter,
        port_area=port_area,
        wall_web=wall_web,
        grain_outer_diameter=float(om_grain.props['diameter'].getValue()),
        grain_length=float(om_grain.props['length'].getValue()),
        inhibited_fwd=inh_fwd,
        inhibited_aft=inh_aft,
        geom_name=getattr(type(om_grain), 'geomName', None) or type(om_grain).__name__,
    )


def from_ric_grain(ric_grain: dict, map_dim: int = 1001) -> FmmTable:
    """
    Convenience wrapper: take a .ric grain dict (with 'type' and
    'properties' keys) and return an FmmTable. Dispatches on the grain
    type to instantiate the right openMotor class.
    """
    _setup_openmotor_path()

    gtype = ric_grain['type']
    om_grain = _instantiate_openmotor_grain(gtype)
    om_grain.setProperties(ric_grain['properties'])
    return from_openmotor(om_grain, map_dim=map_dim)


# Lazy registry of openMotor FmmGrain subclasses by geomName.
_OM_GRAIN_CLASSES = None


def _instantiate_openmotor_grain(geom_name: str):
    """Look up an openMotor grain class by its `geomName` and return a
    fresh instance. Raises ValueError for unknown / non-FMM types."""
    global _OM_GRAIN_CLASSES
    if _OM_GRAIN_CLASSES is None:
        _setup_openmotor_path()
        from motorlib.grain import FmmGrain
        from motorlib.grains.finocyl import Finocyl
        from motorlib.grains.star import StarGrain
        from motorlib.grains.moonBurner import MoonBurner
        from motorlib.grains.cGrain import CGrain
        from motorlib.grains.dGrain import DGrain
        from motorlib.grains.xCore import XCore
        from motorlib.grains.custom import CustomGrain

        # NOTE: RodTubeGrain is a PerforatedGrain (analytic), not FmmGrain;
        # it can be added later by extending from_openmotor to accept
        # PerforatedGrain and skipping the FMM setup steps for analytic types.
        registry = {}
        for cls in (Finocyl, StarGrain, MoonBurner, CGrain, DGrain,
                    XCore, CustomGrain):
            if issubclass(cls, FmmGrain):
                registry[cls.geomName] = cls
        _OM_GRAIN_CLASSES = registry

    if geom_name not in _OM_GRAIN_CLASSES:
        raise ValueError(
            f"Grain type {geom_name!r} is not a supported FMM grain. "
            f"Supported FMM types: {sorted(_OM_GRAIN_CLASSES)}"
        )
    return _OM_GRAIN_CLASSES[geom_name]()


# ================================================================
# Numba-compiled lookup helper
# ================================================================

@njit(cache=True)
def fmm_table_lookup(regress, reg_arr, val_arr, n_samples):
    """
    O(1) linear interpolation into an FMM table sampled on a UNIFORM
    regression-depth grid (which is what `from_openmotor` produces via
    np.linspace). Used inside the @njit hot loop for per-cell
    perimeter/port-area lookups.

    Parameters
    ----------
    regress : float
        Regression depth at this cell, [m].
    reg_arr : ndarray (n_samples,)
        Uniform regression-depth grid from 0 to wall_web.
    val_arr : ndarray (n_samples,)
        Tabulated values to interpolate (perimeter or port_area).
    n_samples : int
        Length of reg_arr / val_arr.

    Returns
    -------
    float
        Linearly interpolated value. Clamped at the table endpoints.
    """
    if regress <= reg_arr[0]:
        return val_arr[0]
    if regress >= reg_arr[n_samples - 1]:
        return val_arr[n_samples - 1]

    # Uniform grid → O(1) index calc
    dr = reg_arr[1] - reg_arr[0]
    idx_f = (regress - reg_arr[0]) / dr
    idx_lo = int(idx_f)
    if idx_lo >= n_samples - 1:
        return val_arr[n_samples - 1]
    f = idx_f - idx_lo
    v0 = val_arr[idx_lo]
    v1 = val_arr[idx_lo + 1]
    return v0 + f * (v1 - v0)
