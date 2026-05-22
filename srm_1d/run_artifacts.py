"""
run_artifacts.py — Per-run artifact directory helper for example scripts.

PURPOSE:
    Example scripts (hasegawa_motor_a.py, run_template.py, etc.) save
    plots and CSV exports to ``artifacts/<motor>/<file>.png`` and
    overwrite previous runs in place. That makes it impossible to keep
    a calibration trail or compare two runs without manually copying
    files between simulations.

    This module gives examples a one-line way to land their artifacts
    in a per-run stamped directory:

        artifacts/<motor>/<YYYY-MM-DDTHH-MM-SS>_<short_sha>[-dirty]/
                          └── pressure.png
                          └── flow.png
                          └── summary.png

USAGE:
    from srm_1d.run_artifacts import artifact_dir
    out = artifact_dir('hasegawa_a')
    plot_pressure(..., save_path=str(out / 'pressure.png'))

DESIGN NOTES:
    - The git SHA is captured at run time (not at module import) so
      a long-running sweep that crosses commits doesn't mis-label.
    - "-dirty" suffix is appended when the working tree has any
      uncommitted changes; this surfaces "I edited code mid-sweep"
      cases that would otherwise be invisible after the fact.
    - If git isn't available (e.g. tarball install), the SHA segment
      is silently dropped and the timestamp alone disambiguates runs.
    - Subdirectories are created lazily on first call; caller is
      responsible for passing the returned path to its save sites.
"""

from __future__ import annotations
import datetime as _dt
import subprocess
from pathlib import Path
from typing import Optional


_REPO_ROOT_CACHE: Optional[Path] = None


def _repo_root() -> Path:
    """Locate the project root (the directory containing ``srm_1d/``)."""
    global _REPO_ROOT_CACHE
    if _REPO_ROOT_CACHE is None:
        # srm_1d/run_artifacts.py → srm_1d/ → repo root
        _REPO_ROOT_CACHE = Path(__file__).resolve().parents[1]
    return _REPO_ROOT_CACHE


def _git_short_sha(repo_root: Path) -> Optional[str]:
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return sha or None
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def _git_is_dirty(repo_root: Path) -> bool:
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return bool(out)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return False


def artifact_dir(motor_name: str, root: Optional[Path] = None) -> Path:
    """
    Return a unique per-run directory for a motor's artifacts.

    Layout:
        <repo_root>/artifacts/<motor_name>/<stamp>/

    where ``<stamp>`` is ``YYYY-MM-DDTHH-MM-SS_<short_sha>[-dirty]``.

    Parameters
    ----------
    motor_name : str
        Short motor identifier used as the subdir name (e.g.
        ``'hasegawa_a'``, ``'BALLSstick'``). Match the existing
        ``artifacts/<motor>/`` convention.
    root : Path, optional
        Override for the repo root. Default: auto-detected as the
        parent of the ``srm_1d/`` package.

    Returns
    -------
    Path
        The freshly-created directory. Caller writes files inside.
    """
    if root is None:
        root = _repo_root()
    root = Path(root)

    stamp = _dt.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    sha = _git_short_sha(root)
    suffix = ""
    if sha:
        suffix = f"_{sha}"
        if _git_is_dirty(root):
            suffix += "-dirty"

    out = root / "artifacts" / motor_name / f"{stamp}{suffix}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_figure(fig, path, *, dpi: int = 150, verbose: bool = True) -> Path:
    """
    Save a matplotlib Figure with srm_1d's standard export defaults.

    This is the lightweight mirror of openMotor's `imageExporter`: one
    place that owns the bbox / dpi / format conventions so individual
    callers don't repeat them. The full openMotor pattern (channel-
    object `SimulationResult` + generic `plot_channels` like
    `GraphWidget`) is deferred until the openMotor frontend
    integration work scopes up (see memory
    `[[srm-1d-long-term-openmotor-integration-goal]]`).

    Replaces the inline ``fig.savefig(save_path, dpi=150)`` calls
    scattered through plotting.py. Adds ``bbox_inches='tight'``
    (matches openMotor) so titles and legends don't get cropped.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save. Caller manages the figure lifecycle (creation,
        closing) — this helper only writes to disk.
    path : str or Path
        Destination. Extension picks the format (matplotlib decides
        backend). Typically built from ``artifact_dir(...)`` for
        per-run organization, e.g.
        ``save_figure(fig, artifact_dir('hasegawa_a') / 'pressure.png')``.
    dpi : int, optional
        Resolution. Default 150 dpi (openMotor's image-export default).
    verbose : bool, optional
        If True (default), print ``"Saved <path>"`` so the caller
        sees where the file landed without managing it themselves.

    Returns
    -------
    Path
        The full path written, for the caller to log / verify.
    """
    out_path = Path(path)
    fig.savefig(str(out_path), bbox_inches="tight", dpi=dpi)
    if verbose:
        print(f"Saved {out_path}")
    return out_path
