"""
station_viz.py — headless backend for the per-station axial visualization
(v0.8.x; see ``docs/v0_8_0/STATION_VIZ_DESIGN.md``).

srm_1d is a 1D axial solver, so per-cell data is its native resolution.
The GUI per-grain channels are an *aggregation* of per-cell axial snapshots;
a **station** selector exposes the resolution that already exists and lets the
user inspect axial variation *within* a grain (e.g. the fore→aft mass-flux
gradient that drives erosive burning). "Per-grain" becomes the special case
"one station per grain."

This module is the **headless, Qt-free** half of the feature (design phases
1–2): it turns a ``run_simulation`` result into

  * a compact **axial payload** (``build_axial_payload``) — a decimated
    ``[n_frames × n_cells]`` field matrix per plottable quantity, plus the
    time base, cell positions, and the cell→grain map; and
  * a default **station model** (``default_stations``) — fore/mid/aft cells
    per grain with boundary reassignment, ready for the GUI selector.

The GUI panel (design phases 3–5) lives on the openMotor-fork side and
consumes these structures; keeping the logic here makes it unit-testable in
the canonical repo without a display.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Cell→grain gap sentinel (matches grain_geometry.compile_geometry_arrays:
# grain cells carry their grain index k >= 0, inter-segment gap cells carry
# -1).
GAP_SENTINEL = -1

# Snapshot per-cell fields exposed to the station plotter by default. Each
# is carried directly in ``result['snapshots'][s][name]`` (length n_cells).
# Mass flux ``G`` is intentionally omitted until a per-cell density (or R)
# snapshot field exists — we don't fabricate it from incomplete state.
DEFAULT_FIELDS = (
    'P', 'u', 'Mach', 'T', 'r_total', 'r_erosive', 'D_port', 'regress',
)


# ================================================================
# Phase 1 — axial payload
# ================================================================

@dataclass
class AxialPayload:
    """Compact, decimated per-cell axial field carried for the station panel.

    Attributes
    ----------
    snap_times : ndarray (n_frames,)
        Decimated snapshot time base [s].
    x_cell : ndarray (n_cells,)
        Cell-center axial positions [m] (constant across frames).
    cell_segment_id : ndarray (n_cells,), int
        Cell→grain index; gap cells carry ``GAP_SENTINEL`` (-1).
    fields : dict[str, ndarray (n_frames, n_cells)]
        One matrix per plottable quantity (see ``DEFAULT_FIELDS``).
    frame_index : ndarray (n_frames,), int
        Indices into the original snapshot list that were kept (for
        traceability / debugging).
    """
    snap_times: np.ndarray
    x_cell: np.ndarray
    cell_segment_id: np.ndarray
    fields: dict
    frame_index: np.ndarray = field(default_factory=lambda: np.array([], int))

    @property
    def n_frames(self) -> int:
        return int(self.snap_times.shape[0])

    @property
    def n_cells(self) -> int:
        return int(self.x_cell.shape[0])

    def series(self, field_name: str, cell_index: int) -> np.ndarray:
        """Time series of ``field_name`` sliced at ``cell_index``."""
        return self.fields[field_name][:, cell_index]


def _decimate_frame_indices(n_snaps: int, max_frames: int) -> np.ndarray:
    """Uniform-stride frame subsample that always keeps the first and last.

    Returns sorted, unique indices into ``range(n_snaps)`` of length
    <= ``max_frames``. ``max_frames <= 0`` keeps every frame.
    """
    if n_snaps <= 0:
        return np.array([], dtype=np.int64)
    if max_frames <= 0 or n_snaps <= max_frames:
        return np.arange(n_snaps, dtype=np.int64)
    idx = np.linspace(0, n_snaps - 1, max_frames)
    return np.unique(np.round(idx).astype(np.int64))


def build_axial_payload(
    result: dict,
    fields=DEFAULT_FIELDS,
    max_frames: int = 240,
    cell_segment_id: Optional[np.ndarray] = None,
) -> Optional[AxialPayload]:
    """Extract a decimated :class:`AxialPayload` from a ``run_simulation`` result.

    Parameters
    ----------
    result : dict
        Output of :func:`srm_1d.run_simulation` (needs ``'snapshots'``; uses
        ``'cell_segment_id'`` / ``'x_cell'`` when present).
    fields : iterable of str
        Snapshot per-cell field names to carry. Names absent from a snapshot
        are skipped (with no error) so the contract degrades gracefully.
    max_frames : int
        Decimation budget on the time axis (first + last always kept).
        ``<= 0`` keeps every snapshot.
    cell_segment_id : ndarray, optional
        Cell→grain map; falls back to ``result['cell_segment_id']``. Required
        either here or in the result.

    Returns
    -------
    AxialPayload or None
        ``None`` if the result carries no snapshots.
    """
    snapshots = result.get('snapshots')
    if not snapshots:
        return None

    if cell_segment_id is None:
        cell_segment_id = result.get('cell_segment_id')
    if cell_segment_id is None:
        raise ValueError(
            "cell_segment_id not in result and not supplied; cannot build "
            "the station payload (re-run with v0.8.x simulation, or pass it)."
        )
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)

    n_snaps = len(snapshots)
    keep = _decimate_frame_indices(n_snaps, max_frames)
    snap_times = np.array([snapshots[i]['t'] for i in keep], dtype=float)

    x_cell = result.get('x_cell')
    if x_cell is None:
        x_cell = snapshots[0]['x']
    x_cell = np.asarray(x_cell, dtype=float)
    n_cells = x_cell.shape[0]

    # Only carry fields that are actually present in the snapshots.
    present = [f for f in fields if f in snapshots[0]]
    field_mats = {
        name: np.empty((keep.shape[0], n_cells), dtype=float)
        for name in present
    }
    for row, src in enumerate(keep):
        snap = snapshots[src]
        for name in present:
            field_mats[name][row, :] = np.asarray(snap[name])

    return AxialPayload(
        snap_times=snap_times,
        x_cell=x_cell,
        cell_segment_id=cell_segment_id,
        fields=field_mats,
        frame_index=keep,
    )


# ================================================================
# Phase 2 — station model + default population
# ================================================================

@dataclass
class Station:
    """A single axial probe: one cell index, classified to its owning grain.

    Attributes
    ----------
    grain : int
        Owning grain index, or ``GAP_SENTINEL`` for an inter-segment gap.
    cell_index : int
        The cell this station samples (the selection primitive).
    position_m : float
        Cell-center distance from the head [m] (a derived label value).
    active : bool
        Whether the station is currently displayed.
    role : str
        ``'fore'`` / ``'mid'`` / ``'aft'`` for auto-placed stations,
        ``'gap'`` for a gap probe, ``'custom'`` for a user-added one.
    label : str
        Human-readable label (grain + role + cell + mm-from-head).
    """
    grain: int
    cell_index: int
    position_m: float
    active: bool
    role: str
    label: str


def _grain_label(grain: int) -> str:
    return "Gap" if grain == GAP_SENTINEL else f"Grain {grain + 1}"


def _make_label(grain: int, role: str, cell_index: int, position_m: float) -> str:
    return (
        f"{_grain_label(grain)} {role} · cell {cell_index} · "
        f"{position_m * 1e3:.0f} mm"
    )


def grain_cell_spans(cell_segment_id: np.ndarray) -> dict:
    """Map each grain index to its sorted array of owned cell indices.

    Gap cells (``GAP_SENTINEL``) are excluded. Grains are keyed by index in
    ascending order.
    """
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)
    grains = sorted(int(g) for g in np.unique(cell_segment_id) if g >= 0)
    return {g: np.flatnonzero(cell_segment_id == g) for g in grains}


def gap_cell_indices(cell_segment_id: np.ndarray) -> np.ndarray:
    """Indices of inter-segment gap cells (``cell_segment_id == GAP_SENTINEL``)."""
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)
    return np.flatnonzero(cell_segment_id == GAP_SENTINEL)


def default_stations(
    cell_segment_id: np.ndarray,
    x_cell: np.ndarray,
) -> list:
    """Build the default fore/mid/aft station set, one group per grain.

    Mirrors openMotor's "first grain shown" default to avoid clutter:

      * three stations per grain — **fore / mid / aft** cells of its span;
      * **fore default-ON, mid + aft default-OFF**;
      * grains with fewer than three cells collapse (deduped by cell index),
        so no two default stations point at the same cell.

    Because fore/mid/aft are chosen from each grain's *actual* owned cells
    (not a position-derived nominal), they can never land in a gap or outside
    the span — the §7 boundary reassignment is therefore implicit here; it
    only matters for user-driven position picks (handled in the GUI).

    Stations are returned grouped by grain in ascending grain order, and
    within a grain in fore→mid→aft order.
    """
    x_cell = np.asarray(x_cell, dtype=float)
    spans = grain_cell_spans(cell_segment_id)
    stations: list = []
    for grain, cells in spans.items():
        n = cells.shape[0]
        if n == 0:
            continue
        # role -> cell index, deduped while preserving fore/mid/aft order.
        candidates = [
            ('fore', int(cells[0])),
            ('mid', int(cells[n // 2])),
            ('aft', int(cells[-1])),
        ]
        seen: set = set()
        for role, ci in candidates:
            if ci in seen:
                continue
            seen.add(ci)
            pos = float(x_cell[ci])
            stations.append(Station(
                grain=grain,
                cell_index=ci,
                position_m=pos,
                active=(role == 'fore'),
                role=role,
                label=_make_label(grain, role, ci, pos),
            ))
    return stations


def make_station(
    cell_index: int,
    cell_segment_id: np.ndarray,
    x_cell: np.ndarray,
    active: bool = True,
    role: str = 'custom',
) -> Station:
    """Build a single station at ``cell_index``, classifying its owning grain.

    Used for user add-station actions: the grain is read from
    ``cell_segment_id`` (a gap cell yields a ``GAP_SENTINEL`` grain with role
    forced to ``'gap'``).
    """
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)
    x_cell = np.asarray(x_cell, dtype=float)
    if not (0 <= cell_index < cell_segment_id.shape[0]):
        raise IndexError(
            f"cell_index {cell_index} out of range [0, {cell_segment_id.shape[0]})"
        )
    grain = int(cell_segment_id[cell_index])
    if grain == GAP_SENTINEL:
        role = 'gap'
    pos = float(x_cell[cell_index])
    return Station(
        grain=grain,
        cell_index=int(cell_index),
        position_m=pos,
        active=active,
        role=role,
        label=_make_label(grain, role, int(cell_index), pos),
    )


# ================================================================
# Cell classification for the GUI station selector (v0.8.x)
# ================================================================
# A station in the rich selector is just a cell index; its category (head /
# grain / gap / aft) and grain-relative role (fore/mid/aft) are DERIVED from
# the index so editing the index reclassifies it automatically. These helpers
# are the headless, testable core the GUI widget formats labels from.

def cell_categories(cell_segment_id: np.ndarray) -> list:
    """Ordered category descriptors covering every cell, head→aft:

      Head (cells before the first grain), each Grain, each inter-grain Gap
      (numbered 1..), and Aft (cells after the last grain). Only categories
      that actually contain cells are emitted.

    Each entry: ``{'kind','label','short','lo','hi','grain','gap'}`` where
    ``kind`` ∈ {head, grain, gap, aft}, ``lo``/``hi`` are inclusive cell
    bounds, ``grain`` is the grain index (grain kind) and ``gap`` the 1-based
    gap number (gap kind), else -1.
    """
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)
    n = int(cell_segment_id.shape[0])
    spans = grain_cell_spans(cell_segment_id)
    if not spans:
        return [{'kind': 'head', 'label': 'Cells', 'short': '',
                 'lo': 0, 'hi': n - 1, 'grain': -1, 'gap': -1}]
    bounds = [(g, int(cells[0]), int(cells[-1])) for g, cells in spans.items()]
    bounds.sort(key=lambda b: b[1])
    cats = []
    if bounds[0][1] > 0:
        cats.append({'kind': 'head', 'label': 'Head', 'short': 'Head',
                     'lo': 0, 'hi': bounds[0][1] - 1, 'grain': -1, 'gap': -1})
    gap_no = 0
    for i, (g, lo, hi) in enumerate(bounds):
        cats.append({'kind': 'grain', 'label': 'Grain {}'.format(g + 1),
                     'short': 'G{}'.format(g + 1), 'lo': lo, 'hi': hi,
                     'grain': g, 'gap': -1})
        if i < len(bounds) - 1:
            nlo = bounds[i + 1][1]
            if nlo > hi + 1:
                gap_no += 1
                cats.append({'kind': 'gap', 'label': 'Gap {}'.format(gap_no),
                             'short': 'Gap{}'.format(gap_no), 'lo': hi + 1,
                             'hi': nlo - 1, 'grain': -1, 'gap': gap_no})
    if bounds[-1][2] < n - 1:
        cats.append({'kind': 'aft', 'label': 'Aft', 'short': 'Aft',
                     'lo': bounds[-1][2] + 1, 'hi': n - 1, 'grain': -1, 'gap': -1})
    return cats


def grain_role(cell_index: int, grain: int, cell_segment_id: np.ndarray) -> str:
    """``'fore'``/``'mid'``/``'aft'`` if ``cell_index`` is that grain's
    fore/mid/aft cell (same anchors as :func:`default_stations`), else ``''``.
    The role is dynamic: it re-applies whenever a station lands on the cell."""
    if grain < 0:
        return ''
    cells = grain_cell_spans(cell_segment_id).get(grain)
    if cells is None or cells.shape[0] == 0:
        return ''
    n = cells.shape[0]
    if cell_index == int(cells[0]):
        return 'fore'
    if cell_index == int(cells[-1]):
        return 'aft'
    if cell_index == int(cells[n // 2]):
        return 'mid'
    return ''


def classify_cell(cell_index: int, cell_segment_id: np.ndarray,
                  x_cell: Optional[np.ndarray] = None) -> dict:
    """Classify a cell for the station selector: which category it belongs to
    and (for grain cells) its fore/mid/aft role. Returns a dict with
    ``cell_index, kind, grain, gap, category_label, category_short, role,
    n_cells`` and ``position_m`` (when ``x_cell`` is given)."""
    cell_segment_id = np.asarray(cell_segment_id, dtype=np.int64)
    cats = cell_categories(cell_segment_id)
    cat = next((c for c in cats if c['lo'] <= cell_index <= c['hi']), None)
    if cat is None:
        cat = {'kind': 'head', 'label': '?', 'short': '?', 'grain': -1, 'gap': -1}
    role = grain_role(cell_index, cat['grain'], cell_segment_id) if cat['kind'] == 'grain' else ''
    out = {
        'cell_index': int(cell_index),
        'kind': cat['kind'],
        'grain': int(cat['grain']),
        'gap': int(cat['gap']),
        'category_label': cat['label'],
        'category_short': cat['short'],
        'role': role,
        'n_cells': int(cell_segment_id.shape[0]),
    }
    if x_cell is not None:
        out['position_m'] = float(np.asarray(x_cell, dtype=float)[cell_index])
    return out


def station_full_label(classification: dict) -> str:
    """Compact label for the plot legend / grain-tab column header, e.g.
    ``'G1 fore (c3)'``, ``'Head (c1)'``, ``'Gap1 (c47)'``, ``'Aft (c103)'``."""
    parts = [classification['category_short']]
    if classification['role']:
        parts.append(classification['role'])
    base = ' '.join(p for p in parts if p)
    return '{} (c{})'.format(base, classification['cell_index']) if base \
        else '(c{})'.format(classification['cell_index'])
