"""Refresh motor ``.ric`` files through openMotor itself (igniter provider D7).

    python -m srm_1d.tools.refresh_ric_igniters <paths...> [--dry-run]

``.ric`` files are openMotor's format and are never hand-edited. This tool
regenerates them with openMotor's own code path: ``uilib.fileIO.loadFile``
(which applies openMotor's migrations) -> ``motorlib.motor.Motor.applyDict``
(which applies the igniter provider's exact-match upgrade of superseded
library snapshots, DESIGN D5) -> ``Motor.getDict`` -> ``uilib.fileIO.saveFile``.

Safety check (per file, BEFORE saving): the parsed data before and after the
round trip are compared and every changed key path is printed. Permitted:
the ``igniter.pyrogen`` material (only when it was actually upgraded), the
keys openMotor's save normally adds (``ADDED_BY_OPENMOTOR``), and an ADDED
``/config/<key>`` whose value equals the fork's ``MotorConfig`` default for that
key (``MotorConfig`` fills absent keys with its defaults on load, and srm_1d
reads only ``ambPressure`` from a motor's config, so this is pure
normalisation; the default is read from the fork, never hardcoded). ANY other
added / removed / changed key means the round trip is lossy for that file: it
is NOT saved and the run exits non-zero.

Needs an openMotor checkout (``SRM1D_OPENMOTOR_PATH`` or the usual sibling
lookup), PyQt6 and the srm_1d plugin (this package). The user's real openMotor
log directory is not touched (redirected to a temp dir while ``uilib`` loads).
"""
import argparse
import os
import sys
import tempfile
import warnings

# Keys openMotor's save adds to / normalises in a .ric (path prefixes, as
# '/'-joined key paths). Only an ADDED key (absent before) is permitted here.
ADDED_BY_OPENMOTOR = (
    '/config/taperSlices',        # MotorConfig gained taperSlices
    '/solverConfigs',             # per-motor solver run-config overrides
)


def _normalise(obj):
    """Tuples -> lists so YAML tuple tags do not register as differences."""
    if isinstance(obj, dict):
        return {k: _normalise(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_normalise(v) for v in obj]
    return obj


def diff_paths(a, b, path=''):
    """``[(path, kind)]`` for every difference; kind in changed/added/removed."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            p = '{}/{}'.format(path, k)
            if k not in a:
                out.append((p, 'added'))
            elif k not in b:
                out.append((p, 'removed'))
            else:
                out += diff_paths(a[k], b[k], p)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff_paths(x, y, '{}[{}]'.format(path, i))
    elif a != b:
        out.append((path, 'changed'))
    return out


def motor_config_defaults(motormod):
    """``{key: default}`` of the fork's ``MotorConfig`` (read, not hardcoded)."""
    return motormod.MotorConfig().getProperties()


def _is_permitted(path, kind, material_upgraded, after=None, config_defaults=None):
    """``after`` = the parsed post-round-trip data (to read an added value)."""
    parts = path.split('/')
    if (kind == 'added' and len(parts) == 3 and parts[1] == 'config'
            and config_defaults is not None and after is not None
            and parts[2] in config_defaults
            and after.get('config', {}).get(parts[2]) == config_defaults[parts[2]]):
        return True
    if path == '/igniter/pyrogen' or path.startswith('/igniter/pyrogen/'):
        return material_upgraded
    return kind == 'added' and any(
        path == p or path.startswith(p + '/') for p in ADDED_BY_OPENMOTOR)


def _import_openmotor():
    """Lazily import the fork + this package's plugin; clear error if absent."""
    import platformdirs
    try:
        from srm_1d.fmm_grain import _setup_openmotor_path
        _setup_openmotor_path()
        import PyQt6.QtWidgets  # noqa: F401  (before uilib, as in app.py)
        # uilib.logger opens the user's real log file at import: redirect it.
        platformdirs.user_log_dir = lambda *a, **k: os.path.join(
            tempfile.mkdtemp(prefix='srm1d_refresh_log_'), 'log')
        import srm_1d.srm1d_plugin  # noqa: F401  (registers the provider)
        import motorlib.motor
        import motorlib.solvers
        import uilib.fileIO as fileIO
    except ImportError as exc:
        raise SystemExit(
            'refresh_ric_igniters needs the openMotor checkout (set '
            'SRM1D_OPENMOTOR_PATH to the inner openMotor/openMotor dir) and '
            'PyQt6; import failed: {}'.format(exc))
    if motorlib.solvers.get_igniter_provider() is None:
        raise SystemExit('no igniter provider registered (openMotor checkout '
                         'lacks the igniter hooks, or the srm_1d plugin failed '
                         'to register)')
    return motorlib.motor, fileIO


def refresh(path, motormod, fileIO, dry_run=False):
    """Round-trip one file. Returns True if it is (or would be) saved OK."""
    from srm_1d.openmotor_adapter import (
        IgniterMaterialUpgradedWarning, StaleIgniterMaterialWarning,
        load_igniter)
    name = os.path.basename(path)
    loaded = fileIO.loadFile(path, fileIO.fileTypes.MOTOR)
    motor = motormod.Motor()
    motor.applyDict(loaded)
    new = motor.getDict()

    before, after = _normalise(loaded), _normalise(new)
    ig_before = (before.get('igniter') or {}).get('pyrogen')
    ig_after = (after.get('igniter') or {}).get('pyrogen')
    upgraded = ig_before != ig_after
    material = (ig_after or {}).get('name', '(no igniter)')

    diffs = diff_paths(before, after)
    cfg_defaults = motor_config_defaults(motormod)
    bad = [(p, k) for p, k in diffs
           if not _is_permitted(p, k, upgraded, after, cfg_defaults)]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        if after.get('igniter') is not None:
            load_igniter(after)
    stale = [w for w in caught
             if issubclass(w.category, StaleIgniterMaterialWarning)]
    again = [w for w in caught
             if issubclass(w.category, IgniterMaterialUpgradedWarning)]

    print('{}: material={} upgraded={} stale_after={} re-upgrade_after={}'.format(
        name, material, 'yes' if upgraded else 'no', len(stale), len(again)))
    for p, k in diffs:
        tag = 'ok ' if (p, k) not in bad else 'BAD'
        print('    [{}] {:<8} {}'.format(tag, k, p))
    for w in stale:
        print('    stale-warning left: {}'.format(str(w.message)[:240]))
    if bad:
        print('    STOP: {} unpermitted change(s); {} NOT saved'.format(
            len(bad), name))
        return False
    if again:
        print('    STOP: an upgraded block upgrades again (not idempotent)')
        return False
    if dry_run:
        print('    (dry run: not saved)')
        return True
    fileIO.saveFile(path, new, fileIO.fileTypes.MOTOR)
    check = _normalise(fileIO.loadFile(path, fileIO.fileTypes.MOTOR))
    if check != after:
        print('    STOP: re-reading the saved file differs from what was written')
        return False
    print('    saved')
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('paths', nargs='+', help='.ric files to refresh')
    ap.add_argument('--dry-run', action='store_true',
                    help='report what would change; write nothing')
    args = ap.parse_args(argv)
    motormod, fileIO = _import_openmotor()
    print('permitted additions by openMotor save:', ', '.join(ADDED_BY_OPENMOTOR))
    ok = True
    for path in args.paths:
        ok = refresh(path, motormod, fileIO, args.dry_run) and ok
    print('RESULT:', 'OK' if ok else 'FAILED (see STOP lines)')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
