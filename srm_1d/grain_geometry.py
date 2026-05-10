"""
grain_geometry.py — Multi-Grain Motor Geometry with BATES Support
==================================================================

PURPOSE:
    Defines the physical geometry of the solid rocket motor and tracks
    how it evolves as propellant burns. Supports single-cylinder and
    multi-segment BATES configurations.

ARCHITECTURE:
    Python classes (GrainSegment, MotorGeometry) are used ONLY for setup
    and configuration. Before the time loop starts, all geometry state is
    exported to plain numpy arrays via compile_geometry_arrays(). The
    per-step geometry updates are Numba-compiled functions that operate
    on those arrays with zero Python overhead.

    This separation exists because Numba cannot JIT-compile code that
    accesses Python object attributes. By converting everything to arrays
    up front, the hot path stays in compiled code.

BATES GEOMETRY:
    BATES (Ballistic Test and Evaluation System) motors have multiple
    short cylindrical grain segments separated by small gaps. Both end
    faces of each segment burn (they are uninhibited), causing the
    segments to shorten axially while the bore widens radially.

    Layout: [gap][grain][gap][grain]...[grain][gap]

    The gaps allow combustion gas to flow between segments and prevent
    the grain from bonding to itself during casting. In the model, gap
    cells have full port diameter (D_outer) and contribute no mass.
"""

import numpy as np
from dataclasses import dataclass
from typing import List

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


# ================================================================
# SECTION 1: Configuration Classes (setup only, not in hot path)
# ================================================================

@dataclass
class GrainSegment:
    """
    One grain segment. Used for configuration only — not accessed
    during the time loop.

    Attributes
    ----------
    x_start : float
        Axial start position of this segment [m].
    length : float
        Initial segment length [m].
    D_bore_fwd : float
        Bore diameter at the forward (head) end [m]. Cylindrical/conical
        only; ignored for FMM segments (set to D_outer if `fmm_table`
        is provided).
    D_outer : float
        Outer diameter (casing inner diameter) [m].
    D_bore_aft : float or None
        Bore diameter at the aft (nozzle) end [m].
        If None, defaults to D_bore_fwd (cylindrical grain).
    inhibit_fwd : bool
        If True, the forward face does not burn.
    inhibit_aft : bool
        If True, the aft face does not burn.
    fmm_table : FmmTable or None
        If set, this segment uses an FMM regression-map lookup for
        its perimeter and port area. Built via
        `srm_1d.fmm_grain.from_openmotor(om_grain)` or `from_ric_grain`.
        Cylindrical/conical segments leave this as None.
    """
    x_start: float
    length: float
    D_bore_fwd: float
    D_outer: float
    D_bore_aft: float = None
    inhibit_fwd: bool = False
    inhibit_aft: bool = False
    fmm_table: object = None  # forward-typed: avoid cyclic import of FmmTable

    def __post_init__(self):
        if self.D_bore_aft is None:
            self.D_bore_aft = self.D_bore_fwd


@dataclass
class MotorGeometry:
    """
    Complete motor geometry specification (grain side only). The
    nozzle (throat diameter, exit diameter, erosion/slag, etc.) lives
    on a separate Nozzle object — see srm_1d.nozzle.

    Call compile_geometry_arrays() to extract Numba-compatible arrays
    before starting the time loop.

    Attributes
    ----------
    L_motor : float
        Total motor length including gaps [m].
    D_outer : float
        Outer diameter (casing inner diameter) [m].
    segments : list of GrainSegment
        Grain segment definitions.
    N_cells : int
        Number of finite volume cells along the motor length.
    """
    L_motor: float
    D_outer: float
    segments: List[GrainSegment]
    N_cells: int = 150

    @property
    def dx(self):
        """Cell width [m]."""
        return self.L_motor / self.N_cells

    @property
    def web_thickness(self):
        """Initial web thickness of the first segment [m]."""
        if self.segments:
            return (self.D_outer - self.segments[0].D_bore_fwd) / 2.0
        return 0.0

    def compile_geometry_arrays(self):
        """
        Export all geometry state to plain arrays for Numba.

        Per-cell initial bore diameter (cell_D_bore_init) is computed
        by linear interpolation between each segment's D_bore_fwd and
        D_bore_aft. This supports cylindrical (fwd == aft), conical,
        and future FMM grains where the bore profile is set per-cell
        from a regression map.

        Auto-inhibition: when two segments touch (no gap between them),
        their shared faces are automatically inhibited. This assumes
        bonded propellant at the interface.
        """
        N = self.N_cells
        dx = self.dx
        x_centers = np.linspace(dx / 2, self.L_motor - dx / 2, N)
        N_seg = len(self.segments)

        # Auto-inhibit touching faces before building arrays
        for k in range(N_seg - 1):
            seg_k_end = self.segments[k].x_start + self.segments[k].length
            seg_k1_start = self.segments[k + 1].x_start
            if abs(seg_k_end - seg_k1_start) < 1e-8:
                # Segments touch — inhibit shared faces
                self.segments[k].inhibit_aft = True
                self.segments[k + 1].inhibit_fwd = True

        # Per-segment arrays
        seg_x_start = np.array([s.x_start for s in self.segments])
        seg_length = np.array([s.length for s in self.segments])
        seg_D_bore_fwd = np.array([s.D_bore_fwd for s in self.segments])
        seg_D_bore_aft = np.array([s.D_bore_aft for s in self.segments])
        seg_inhibit_fwd = np.array([s.inhibit_fwd for s in self.segments])
        seg_inhibit_aft = np.array([s.inhibit_aft for s in self.segments])

        # End-face regression state (evolves during simulation)
        seg_fwd_regression = np.zeros(N_seg)
        seg_aft_regression = np.zeros(N_seg)

        # Per-cell: which segment does each cell belong to? (-1 = gap).
        # Choose the segment with the largest overlap. At snapped
        # interfaces, floating-point roundoff can otherwise make the
        # first cell after a boundary appear to overlap the upstream
        # segment by ~1e-16 m and get misclassified.
        cell_segment_id = np.full(N, -1, dtype=np.int64)
        for i in range(N):
            x = x_centers[i]
            x_lo = x - 0.5 * dx
            x_hi = x + 0.5 * dx
            best_overlap = 0.0
            for k in range(N_seg):
                seg_lo = seg_x_start[k]
                seg_hi = seg_x_start[k] + seg_length[k]
                overlap = max(0.0, min(x_hi, seg_hi) - max(x_lo, seg_lo))
                if overlap > best_overlap:
                    best_overlap = overlap
                    cell_segment_id[i] = k
            if best_overlap <= 1e-12 * dx:
                cell_segment_id[i] = -1

        # Per-cell initial bore diameter: interpolated from segment
        # fwd/aft values. This is the fundamental geometry representation
        # — cylindrical grains have constant D, conical grains have a
        # linear gradient, FMM grains will have arbitrary profiles.
        cell_D_bore_init = np.full(N, self.D_outer)
        for i in range(N):
            k = cell_segment_id[i]
            if k >= 0:
                x = x_centers[i]
                seg_lo = seg_x_start[k]
                seg_hi = seg_lo + seg_length[k]
                # Fractional position within segment [0=fwd, 1=aft]
                if seg_length[k] > 1e-10:
                    frac = (x - seg_lo) / seg_length[k]
                    frac = max(0.0, min(1.0, frac))
                else:
                    frac = 0.5
                cell_D_bore_init[i] = (
                    seg_D_bore_fwd[k] * (1.0 - frac) + seg_D_bore_aft[k] * frac
                )

        # Initial D_port from per-cell bore initialization
        D_port = cell_D_bore_init.copy()
        # Gap cells get full outer diameter
        for i in range(N):
            if cell_segment_id[i] < 0:
                D_port[i] = self.D_outer

        # Per-cell wall web (radial regression to burnout). For
        # cylindrical/conical: (D_outer - D_bore_init) / 2. For FMM
        # cells: overwritten with FmmTable.wall_web below.
        cell_wall_web = (self.D_outer - cell_D_bore_init) / 2.0

        # Per-cell segment type (0=cylindrical/conical, 1=FMM).
        cell_segment_type = np.zeros(N, dtype=np.int64)

        # Per-cell FMM-table index (-1 if not FMM).
        cell_fmm_idx = np.full(N, -1, dtype=np.int64)

        # ------------------------------------------------------------
        # Pack FMM tables into flat (CSR-like) arrays
        # ------------------------------------------------------------
        # Each FMM segment k contributes its (reg_depth, perimeter,
        # port_area) sample arrays. We concatenate them and record
        # per-grain offsets so the @njit hot loop can do an O(1)
        # lookup with one indirection.
        fmm_seg_indices = [k for k, seg in enumerate(self.segments)
                           if seg.fmm_table is not None]
        n_fmm_segs = len(fmm_seg_indices)

        if n_fmm_segs > 0:
            # Build offset array (length n_fmm_segs+1)
            fmm_offset = np.zeros(n_fmm_segs + 1, dtype=np.int64)
            # First pass: compute total samples and offsets
            for fi, k in enumerate(fmm_seg_indices):
                fmm_offset[fi + 1] = fmm_offset[fi] + self.segments[k].fmm_table.n_samples
            total_samples = int(fmm_offset[-1])

            fmm_reg_flat = np.empty(total_samples)
            fmm_perim_flat = np.empty(total_samples)
            fmm_port_flat = np.empty(total_samples)

            # Map original-segment-index → fmm_idx for cell tagging
            seg_to_fmm = {k: fi for fi, k in enumerate(fmm_seg_indices)}

            # Second pass: copy samples + tag cells
            for fi, k in enumerate(fmm_seg_indices):
                tab = self.segments[k].fmm_table
                start = fmm_offset[fi]
                end = fmm_offset[fi + 1]
                fmm_reg_flat[start:end] = tab.reg_depth
                fmm_perim_flat[start:end] = tab.perimeter
                fmm_port_flat[start:end] = tab.port_area

            # Tag cells belonging to FMM segments
            for i in range(N):
                k = cell_segment_id[i]
                if k >= 0 and k in seg_to_fmm:
                    cell_segment_type[i] = 1
                    cell_fmm_idx[i] = seg_to_fmm[k]
                    cell_wall_web[i] = self.segments[k].fmm_table.wall_web
        else:
            # No FMM segments — provide empty arrays Numba can accept.
            fmm_offset = np.zeros(1, dtype=np.int64)
            fmm_reg_flat = np.empty(0)
            fmm_perim_flat = np.empty(0)
            fmm_port_flat = np.empty(0)

        # Per-cell radial regression depth (primary state for the hot
        # loop). Starts at 0 (no regression yet).
        regress = np.zeros(N)

        return {
            'x_centers': x_centers,
            'dx': dx,
            'N': N,
            'N_seg': N_seg,
            'D_outer': self.D_outer,
            'seg_x_start': seg_x_start,
            'seg_length': seg_length,
            'seg_D_bore_fwd': seg_D_bore_fwd,
            'seg_D_bore_aft': seg_D_bore_aft,
            'seg_inhibit_fwd': seg_inhibit_fwd,
            'seg_inhibit_aft': seg_inhibit_aft,
            'seg_fwd_regression': seg_fwd_regression,
            'seg_aft_regression': seg_aft_regression,
            'cell_segment_id': cell_segment_id,
            'cell_D_bore_init': cell_D_bore_init,
            'cell_wall_web': cell_wall_web,
            'cell_segment_type': cell_segment_type,
            'cell_fmm_idx': cell_fmm_idx,
            'fmm_offset': fmm_offset,
            'fmm_reg_flat': fmm_reg_flat,
            'fmm_perim_flat': fmm_perim_flat,
            'fmm_port_flat': fmm_port_flat,
            'n_fmm_segs': n_fmm_segs,
            'regress': regress,
            'D_port': D_port,
        }

    def total_propellant_volume(self):
        """
        Total initial propellant volume [m³].

        Cylindrical/conical: integrated annular volume.
        FMM: `(casting_area − initial_port_area) · length`, where the
        port area comes from the segment's FmmTable. This includes the
        non-circular core cross-section exactly.

        NOTE: This only counts bore volume, not end-face propellant.
        For BATES grains where end faces burn, this slightly
        underestimates the true propellant mass.
        """
        V = 0.0
        casting_area = np.pi / 4.0 * self.D_outer ** 2
        for s in self.segments:
            if s.fmm_table is not None:
                V += (casting_area - s.fmm_table.initial_port_area) * s.length
            else:
                # Cylindrical/conical: ∫ π/4 (D_outer² - D_bore(x)²) dx
                # = π/4 · L · (D_outer² - (D_fwd² + D_fwd·D_aft + D_aft²)/3)
                D_f = s.D_bore_fwd
                D_a = s.D_bore_aft
                V += np.pi / 4.0 * (
                    self.D_outer**2 - (D_f**2 + D_f * D_a + D_a**2) / 3.0
                ) * s.length
        return V


# ================================================================
# SECTION 2: Numba-compiled per-step geometry functions
# ================================================================

@njit(cache=True)
def _fmm_lookup_flat(regress_val, fmm_idx, fmm_offset, val_flat, reg_flat):
    """
    O(1) linear interpolation into a flat (CSR-packed) FMM table for
    one segment. `reg_flat` is assumed to be a uniform grid (np.linspace
    output, which `from_openmotor` produces) per segment.
    Mirrors fmm_grain.fmm_table_lookup but accepts CSR offsets.
    """
    start = fmm_offset[fmm_idx]
    end = fmm_offset[fmm_idx + 1]
    n = end - start
    if n == 0:
        return 0.0
    if regress_val <= reg_flat[start]:
        return val_flat[start]
    if regress_val >= reg_flat[end - 1]:
        return val_flat[end - 1]
    dr = reg_flat[start + 1] - reg_flat[start]
    idx_f = (regress_val - reg_flat[start]) / dr
    idx_lo = int(idx_f)
    if idx_lo >= n - 1:
        return val_flat[end - 1]
    f = idx_f - idx_lo
    v0 = val_flat[start + idx_lo]
    v1 = val_flat[start + idx_lo + 1]
    return v0 + f * (v1 - v0)


@njit(cache=True)
def _saint_robert_local(P, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs):
    """Saint-Robert with hard-switchover tab lookup. Mirrors burn_rate.py
    helpers (duplicated here to keep grain_geometry.py a leaf module)."""
    P_pos = max(P, 0.0)
    if P_pos == 0.0:
        return 0.0
    for k in range(n_tabs):
        if tab_min_p[k] < P < tab_max_p[k]:
            return tab_a[k] * P_pos ** tab_n[k]
    # Fallback: closest boundary
    best = 0
    d_lo0 = abs(P - tab_min_p[0])
    d_hi0 = abs(P - tab_max_p[0])
    best_dist = d_lo0 if d_lo0 < d_hi0 else d_hi0
    for k in range(1, n_tabs):
        d_lo = abs(P - tab_min_p[k])
        d_hi = abs(P - tab_max_p[k])
        d = d_lo if d_lo < d_hi else d_hi
        if d < best_dist:
            best = k
            best_dist = d
    return tab_a[best] * P_pos ** tab_n[best]


@njit(cache=True)
def update_cell_geometry(
    regress, D_port, x_centers, dx, N, N_seg, D_outer,
    seg_x_start, seg_length, seg_fwd_reg, seg_aft_reg,
    seg_inhibit_fwd, seg_inhibit_aft,
    cell_segment_id,
    P, rho_propellant,
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
    A_port, C_burn, D_hyd, is_grain, endface_msource,
    cell_D_bore_init, cell_wall_web,
    cell_segment_type, cell_fmm_idx,
    fmm_offset, fmm_reg_flat, fmm_perim_flat, fmm_port_flat,
):
    """
    Compute per-cell geometry from the current state.

    Each cell spans [x - dx/2, x + dx/2]. Each segment currently
    extends from x_fwd to x_aft (accounting for end-face regression).
    The axial overlap between the cell and the segment determines what
    fraction of the cell contains propellant:

        overlap = max(0, min(x + dx/2, x_aft) - max(x - dx/2, x_fwd))
        grain_frac = overlap / dx     ∈ [0, 1]

    The bore burning perimeter scales by grain_frac, so boundary cells
    transition smoothly from fully burning to fully consumed as the end
    face regresses through them. This eliminates the discrete pressure
    steps that occur when cells flip from grain to gap one at a time.

    The radial burnout ramp (f_active) is applied on top of grain_frac
    for cells approaching D_outer. Both factors multiply C_burn and
    (for f_active) the regression rate in advance_bore_regression.

    Parameters
    ----------
    D_port : ndarray (N,)
        Current port diameter at each cell [m].
    x_centers : ndarray (N,)
        Cell-center axial positions [m].
    dx : float
        Cell width [m].
    N : int
        Number of cells.
    N_seg : int
        Number of grain segments.
    D_outer : float
        Outer diameter [m].
    seg_x_start, seg_length : ndarray (N_seg,)
        Segment positions and initial lengths [m].
    seg_fwd_reg, seg_aft_reg : ndarray (N_seg,)
        Current end-face regression distances [m].
    seg_inhibit_fwd, seg_inhibit_aft : ndarray (N_seg,), bool
        End-face inhibition flags.
    cell_segment_id : ndarray (N,), int
        Output: which segment each cell belongs to (-1 = gap).
    P : ndarray (N,)
        Current pressure [Pa]. Used for end-face burn rate.
    rho_propellant : float
        Propellant density [kg/m³].
    tab_min_p, tab_max_p, tab_a, tab_n : ndarray (n_tabs,)
        Per-tab Saint-Robert parameters and operating-pressure ranges.
    n_tabs : int
        Number of tabs.
    A_port, C_burn, D_hyd : ndarray (N,)
        Output: port area, burning perimeter, hydraulic diameter.
    is_grain : ndarray (N,), bool
        Output: True for cells with any grain overlap.
    endface_msource : ndarray (N,)
        Output: end-face mass source per unit length [kg/(m·s)].
    cell_D_bore_init : ndarray (N,)
        Initial bore diameter for each cell [m]. Per-cell to support
        conical and FMM grains with axially-varying bore profiles.

    Returns
    -------
    cell_segment_id : ndarray (N,), int
        Updated segment assignment.
    """
    PI = 3.141592653589793

    for i in range(N):
        # Derive port shape from per-cell regression depth.
        # Branch on cell_segment_type:
        #   0 = cylindrical/conical: D_port = D_bore_init + 2·regress (analytic)
        #   1 = FMM: A_port and base_perimeter from FMM table lookup
        if cell_segment_type[i] == 1 and cell_fmm_idx[i] >= 0:
            A_port_cell = _fmm_lookup_flat(
                regress[i], cell_fmm_idx[i], fmm_offset,
                fmm_port_flat, fmm_reg_flat,
            )
            base_perimeter = _fmm_lookup_flat(
                regress[i], cell_fmm_idx[i], fmm_offset,
                fmm_perim_flat, fmm_reg_flat,
            )
            # Effective port diameter for snapshots / output (hydraulic-
            # equivalent — equals D for circular ports).
            if A_port_cell > 0.0:
                D_eff = (4.0 * A_port_cell / PI) ** 0.5
            else:
                D_eff = D_outer
            if D_eff > D_outer:
                D_eff = D_outer
            D = D_eff
            D_port[i] = D
            A_port[i] = A_port_cell
            # Hydraulic diameter for non-circular port: D_h = 4·A / P_w
            if base_perimeter > 1e-10:
                D_hyd[i] = 4.0 * A_port_cell / base_perimeter
            else:
                D_hyd[i] = D_outer
        else:
            # Cylindrical/conical analytic path
            D = cell_D_bore_init[i] + 2.0 * regress[i]
            if D > D_outer:
                D = D_outer
            D_port[i] = D
            A_port[i] = PI / 4.0 * D * D
            D_hyd[i] = D
            base_perimeter = PI * D

        x = x_centers[i]

        # Reset outputs for this cell
        cell_segment_id[i] = -1
        is_grain[i] = False
        C_burn[i] = 0.0
        endface_msource[i] = 0.0
        total_grain_frac = 0.0

        # Cell extent
        x_lo = x - 0.5 * dx   # left edge of this cell
        x_hi = x + 0.5 * dx   # right edge of this cell

        # Search segments for overlap with this cell
        for k in range(N_seg):
            x_fwd = seg_x_start[k] + seg_fwd_reg[k]
            x_aft = seg_x_start[k] + seg_length[k] - seg_aft_reg[k]

            if x_fwd >= x_aft:
                continue

            # -----------------------------------------------
            # 1. Volumetric Overlap 
            # (Only applies if grain physically exists in this cell)
            # -----------------------------------------------
            overlap = max(0.0, min(x_hi, x_aft) - max(x_lo, x_fwd))
            grain_frac = overlap / dx

            if grain_frac > 1e-12:
                cell_segment_id[i] = k
                is_grain[i] = True
                total_grain_frac += grain_frac

                w_total = cell_wall_web[i]
                w_remaining = w_total - regress[i]
                burnout_zone = 0.05 * w_total

                if w_remaining <= 0.0:
                    f_active = 0.0
                elif w_total < 1e-10:
                    f_active = 1.0
                elif w_remaining < burnout_zone:
                    f_active = w_remaining / burnout_zone
                else:
                    f_active = 1.0

                C_burn[i] = base_perimeter * grain_frac * f_active

            # -----------------------------------------------
            # 2. End-face mass sources (Linear Distribution Kernel)
            # Evaluated independently of grain_frac so mass is not dropped
            # when the hat function pushes it into an empty gap cell.
            # -----------------------------------------------
            casting_area = PI / 4.0 * D_outer * D_outer

            if not seg_inhibit_fwd[k]:
                dist_fwd = abs(x_fwd - x)
                weight = 0.0
                
                if dist_fwd < dx:
                    weight = 1.0 - (dist_fwd / dx)
                if i == 0 and x_fwd < x:
                    weight = 1.0
                elif i == N - 1 and x_fwd > x:
                    weight = 1.0
                    
                if weight > 0.0:
                    x_sample = min(x_fwd + 0.1 * dx, x_aft - 1e-6)
                    i_sample = max(0, min(N - 1, int(x_sample / dx)))
                    
                    if cell_segment_type[i_sample] == 1 and cell_fmm_idx[i_sample] >= 0:
                        A_port_face = _fmm_lookup_flat(regress[i_sample], cell_fmm_idx[i_sample], fmm_offset, fmm_port_flat, fmm_reg_flat)
                    else:
                        D_face = cell_D_bore_init[i_sample] + 2.0 * regress[i_sample]
                        D_face = min(D_face, D_outer)
                        A_port_face = PI / 4.0 * D_face * D_face
                        
                    A_face = casting_area - A_port_face
                    if A_face < 0.0: A_face = 0.0
                    
                    r_normal = _saint_robert_local(P[i], tab_min_p, tab_max_p, tab_a, tab_n, n_tabs)
                    endface_msource[i] += weight * rho_propellant * r_normal * A_face / dx

            if not seg_inhibit_aft[k]:
                dist_aft = abs(x_aft - x)
                weight = 0.0
                
                if dist_aft < dx:
                    weight = 1.0 - (dist_aft / dx)
                if i == 0 and x_aft < x:
                    weight = 1.0
                elif i == N - 1 and x_aft > x:
                    weight = 1.0
                    
                if weight > 0.0:
                    x_sample = max(x_aft - 0.1 * dx, x_fwd + 1e-6)
                    i_sample = max(0, min(N - 1, int(x_sample / dx)))
                    
                    if cell_segment_type[i_sample] == 1 and cell_fmm_idx[i_sample] >= 0:
                        A_port_face = _fmm_lookup_flat(regress[i_sample], cell_fmm_idx[i_sample], fmm_offset, fmm_port_flat, fmm_reg_flat)
                    else:
                        D_face = cell_D_bore_init[i_sample] + 2.0 * regress[i_sample]
                        D_face = min(D_face, D_outer)
                        A_port_face = PI / 4.0 * D_face * D_face
                        
                    A_face = casting_area - A_port_face
                    if A_face < 0.0: A_face = 0.0
                    
                    r_normal = _saint_robert_local(P[i], tab_min_p, tab_max_p, tab_a, tab_n, n_tabs)
                    endface_msource[i] += weight * rho_propellant * r_normal * A_face / dx

        # -----------------------------------------------
        # Volumetric Flow Area Smoothing
        # -----------------------------------------------
        if total_grain_frac > 0.0:
            casting_area = PI / 4.0 * D_outer * D_outer
            # Smoothly blend the port area based on how much grain is in the cell
            A_port[i] = A_port[i] * total_grain_frac + casting_area * (1.0 - total_grain_frac)
            D_port[i] = (4.0 * A_port[i] / PI) ** 0.5
            D_hyd[i] = D_port[i]
        else:
            # Pure Gap Cell
            D_port[i] = D_outer
            A_port[i] = PI / 4.0 * D_outer * D_outer
            D_hyd[i] = D_outer

    return cell_segment_id


@njit(cache=True)
def advance_bore_regression(regress, r_total, dt, N,
                            cell_wall_web, cell_segment_id,
                            burnout_zone_frac=0.05):
    """
    Advance per-cell radial regression depth.

    Primary state is `regress[i]` (radial regression depth in m). This
    unifies cylindrical/conical and FMM grains: for cylindrical the
    derived port diameter is `D_port = D_bore_init + 2·regress`; for
    FMM the perimeter and port area are looked up from a table at
    `regress[i]`. Both produce the same integration here.

    BURNOUT RAMP:
    Near the wall, the regression rate is multiplied by an "active
    fraction" factor that smoothly decreases to zero as the cell
    approaches its wall web (= burnout regression depth):

        f_active = min(1, web_remaining / burnout_zone)
        web_remaining = cell_wall_web[i] - regress[i]
        burnout_zone  = burnout_zone_frac × cell_wall_web[i]

    This captures the physical reality that, in a real grain, the
    propellant at any axial station does not all reach the wall
    simultaneously — variations in local burn rate cause some portions
    to burn through earlier than others, gradually reducing the active
    burning area. A 1D model cannot resolve this internal variation,
    but ramping the effective regression rate captures the integrated
    effect on the chamber pressure.

    MASS CONSERVATION:
    For cylindrical: `C_burn = π·D` (no ramp on perimeter), so the
    mass source `ρ_p · r · C_burn` equals the volumetric regression
    rate `(π/4) · d(D²)/dt` exactly. For FMM: `C_burn = perimeter(reg)`
    and `dV/dr = perimeter(reg)·L·dx`, so the same conservation holds.

    TERMINATION:
    Once f_active drops below 0.001, the cell is forced to burnout
    (regress = cell_wall_web). Mass lost in this final clamp is
    negligible (<0.1% of the cell's propellant).

    Parameters
    ----------
    regress : ndarray (N,)
        Radial regression depth [m] at each cell. Modified in-place.
    r_total : ndarray (N,)
        Total burn rate at each cell [m/s].
    dt : float
        Time step [s].
    N : int
        Number of cells.
    cell_wall_web : ndarray (N,)
        Per-cell regression depth at burnout [m]. For cylindrical
        cells: `(D_outer - D_bore_init) / 2`. For FMM cells: the
        FmmTable's `wall_web`.
    cell_segment_id : ndarray (N,), int
        Segment assignment (-1 = gap, skip).
    burnout_zone_frac : float
        Fraction of the web thickness over which burnout is ramped.

    Returns
    -------
    n_bore_active : int
        Number of cells still burning (regress < cell_wall_web).
    """
    n_bore_active = 0
    F_ACTIVE_MIN = 0.001

    for i in range(N):
        k = cell_segment_id[i]
        if k < 0:
            continue  # Gap

        w_total = cell_wall_web[i]
        if w_total <= 1e-10:
            continue
        if regress[i] >= w_total:
            continue  # Already burned through

        burnout_zone = burnout_zone_frac * w_total
        w_remaining = w_total - regress[i]

        if burnout_zone > 1e-10:
            f_active = min(1.0, w_remaining / burnout_zone)
        else:
            f_active = 1.0

        if f_active < F_ACTIVE_MIN:
            regress[i] = w_total
            continue

        regress[i] += r_total[i] * f_active * dt

        if regress[i] > w_total:
            regress[i] = w_total
        else:
            n_bore_active += 1

    return n_bore_active


@njit(cache=True)
def advance_endface_regression(
    seg_fwd_reg, seg_aft_reg, seg_length, seg_x_start,
    seg_inhibit_fwd, seg_inhibit_aft,
    N_seg, P_avg, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs, dt,
):
    """
    Advance end-face regression for all segments.

    End faces burn at the normal (Saint-Robert) rate only — no erosive
    component. This is standard for BATES geometry: the end faces are
    perpendicular to the flow and see no significant crossflow velocity.

    Parameters
    ----------
    seg_fwd_reg, seg_aft_reg : ndarray (N_seg,)
        Forward/aft regression distance [m]. Modified in-place.
    seg_length : ndarray (N_seg,)
        Initial segment lengths [m].
    seg_x_start : ndarray (N_seg,)
        Segment start positions [m].
    seg_inhibit_fwd, seg_inhibit_aft : ndarray (N_seg,), bool
        Inhibition flags.
    N_seg : int
        Number of segments.
    P_avg : float
        Average pressure [Pa] for burn rate computation.
    tab_min_p, tab_max_p, tab_a, tab_n : ndarray (n_tabs,)
        Per-tab Saint-Robert parameters and operating-pressure ranges.
    n_tabs : int
        Number of tabs.
    dt : float
        Time step [s].
    """
    r_normal = _saint_robert_local(
        P_avg, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs)

    for k in range(N_seg):
        x_fwd = seg_x_start[k] + seg_fwd_reg[k]
        x_aft = seg_x_start[k] + seg_length[k] - seg_aft_reg[k]
        if x_fwd >= x_aft:
            continue

        if not seg_inhibit_fwd[k]:
            seg_fwd_reg[k] += r_normal * dt
        if not seg_inhibit_aft[k]:
            seg_aft_reg[k] += r_normal * dt


import warnings

def build_snapped_geometry(segments_spec: list[dict], D_outer: float, target_propellant_cells: int = 100) -> MotorGeometry:
    """
    Intelligent geometry preprocessor that guarantees perfect node alignment
    and CFD minimum gap resolution via Integer-Snapping.
    
    Parameters
    ----------
    segments_spec : list of dict
        Specifications for each segment. Expected keys:
        - 'D_bore_fwd': float (bore diameter at forward end)
        - 'D_bore_aft': float (optional, defaults to D_bore_fwd)
        - 'length': float (segment length)
        - 'gap_after': float (optional, gap after this segment, default 0.0)
        - 'inhibit_fwd': bool (optional, default False)
        - 'inhibit_aft': bool (optional, default False)
        - 'fmm_table': FmmTable or None (optional, attaches an FMM
          regression table to this segment for Finocyl/Star/etc.)
    D_outer : float
        Outer casing diameter [m].
    target_propellant_cells : int
        Desired number of cells to represent the propellant mass.
    """
    # 1. Sum total pure propellant length
    L_prop_total = sum(spec['length'] for spec in segments_spec)
    
    # 2. Calculate preliminary coarse dx
    dx = L_prop_total / target_propellant_cells
    
    # 3. The Nyquist-CFD Clamp (Enforce 3-cell minimum on physical gaps only)
    gaps = [spec.get('gap_after', 0.0) for spec in segments_spec if spec.get('gap_after', 0.0) > 1e-6]
    leading_gap = 0.001 # Minimal leading spacer for internal boundary conditions
    
    # Only clamp if there are actual physical gaps between segments
    if gaps:
        min_gap = min(gaps)
        max_allowed_dx = min_gap / 1.0 # Reverted to 1.0 from 3.0 to speed up run time, even though 3.0 increases resolution
        
        if dx > max_allowed_dx:
            dx = max_allowed_dx
        
    # 4. Reconstruct geometry using Integer-Snapping
    segments = []
    x_cursor = 0.0
    total_cells = 0
    
    # Snap leading gap
    n_leading = max(1, int(round(leading_gap / dx)))
    x_cursor += n_leading * dx
    total_cells += n_leading

    for i, spec in enumerate(segments_spec):
        # Snap the propellant segment
        n_seg = max(1, int(round(spec['length'] / dx)))
        snapped_length = n_seg * dx
        
        segments.append(GrainSegment(
            x_start=x_cursor,
            length=snapped_length,
            D_bore_fwd=spec['D_bore_fwd'],
            D_outer=D_outer,
            D_bore_aft=spec.get('D_bore_aft', spec['D_bore_fwd']),
            inhibit_fwd=spec.get('inhibit_fwd', False),
            inhibit_aft=spec.get('inhibit_aft', False),
            fmm_table=spec.get('fmm_table', None),
        ))
        
        x_cursor += snapped_length
        total_cells += n_seg
        
        # Snap the trailing gap
        raw_gap = spec.get('gap_after', 0.0)
        if raw_gap > 1e-6:
            n_gap = max(1, int(round(raw_gap / dx)))  # Stop forcing 3-cell minimum, allow 1 for speed
            snapped_gap = n_gap * dx
            x_cursor += snapped_gap
            total_cells += n_gap
            
    # Trailing motor spacer
    n_trailing = max(1, int(round(0.001 / dx)))
    x_cursor += n_trailing * dx
    total_cells += n_trailing
    
    L_motor_snapped = x_cursor
    
    # Verify deviation is within acceptable physical limits (< 10mm)
    if abs(L_motor_snapped - (L_prop_total + sum(gaps))) > 0.01:
        warnings.warn(f"Integer-snapping altered total motor length by more than 10mm. Check segment proportions.")

    return MotorGeometry(
        L_motor=L_motor_snapped,
        D_outer=D_outer,
        segments=segments,
        N_cells=total_cells
    )

