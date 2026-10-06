"""Offscreen GUI smoke for the plugin-provided igniter (I2 step 10).

Run (all three env vars from docs/igniter_provider/TASKS.md exported):
    QT_QPA_PLATFORM=offscreen python docs/igniter_provider/smoke/gui_smoke.py

Starts the real openMotor App (as main.py does, without entering the event
loop) with ``platformdirs`` pointed at a temp dir (the user's real config/log
dirs are never touched), loads a temp copy of the COMMITTED
motors/machbusterNew.ric (legacy BPNV igniter), drives the motor-config editor
widgets, saves to a temp file and reloads it.
"""
import copy
import os
import subprocess
import sys
import tempfile

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
TMP = tempfile.mkdtemp(prefix='om_gui_smoke_')

import platformdirs  # noqa: E402
platformdirs.user_data_dir = lambda *a, **k: os.path.join(TMP, 'data')
platformdirs.user_log_dir = lambda *a, **k: os.path.join(TMP, 'log')
os.makedirs(os.path.join(TMP, 'data'), exist_ok=True)

OM = os.environ['SRM1D_OPENMOTOR_PATH']
sys.path.insert(0, OM)
os.chdir(OM)                       # app.py resolves resources relative to cwd
sys.argv = [os.path.join(OM, 'main.py')]

import srm_1d.fmm_grain  # noqa: E402
srm_1d.fmm_grain._install_fake_find_perimeter_cy()   # Cython ext not built here

import matplotlib.pyplot  # noqa: E402,F401
from PyQt6.QtWidgets import QApplication  # noqa: E402
from app import App  # noqa: E402
import motorlib.motor  # noqa: E402
import motorlib.solvers  # noqa: E402
import uilib.fileIO as fileIO  # noqa: E402

SOLVER = 'srm_1d-transient'


def check(cond, msg):
    print(('PASS ' if cond else 'FAIL ') + msg)
    if not cond:
        raise SystemExit(1)


def diff_paths(a, b, path=''):
    """Every key path where two parsed structures differ."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(path + '/' + str(k) + ' (added/removed)')
            else:
                out += diff_paths(a[k], b[k], path + '/' + str(k))
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff_paths(x, y, path + '[%d]' % i)
    elif a != b:
        out.append(path)
    return out


src = subprocess.run(
    ['git', '-C', ROOT, 'show', 'HEAD:motors/machbusterNew.ric'],
    capture_output=True, check=True).stdout
ric = os.path.join(TMP, 'machbusterNew.ric')
with open(ric, 'wb') as f:
    f.write(src)

app = App(sys.argv)
# Modal dialogs would block an offscreen run forever: never show them.
app.outputMessage = lambda content, title='openMotor': print('[dialog]', content)


def _no_dialog_exception(exc, text, title='openMotor - Error'):
    raise AssertionError('unexpected error dialog: %s %s' % (text, exc))


app.outputException = _no_dialog_exception
app.promptYesNo = lambda *a, **k: True
check(SOLVER in motorlib.solvers.list_solvers(), 'srm_1d solver discovered')
provider = motorlib.solvers.get_igniter_provider()
check(provider is not None, 'igniter provider registered')
lib = {m['name']: m for m in provider.get_igniter_materials()}
check(sorted(app.pyrogenManager.getNames()) == sorted(lib),
      'GUI pyrogen library holds every srm_1d material: '
      + ', '.join(app.pyrogenManager.getNames()))

committed = fileIO.loadFile(ric, fileIO.fileTypes.MOTOR)
check(committed['igniter']['pyrogen']['T_flame'] == 2800.0,
      'committed file embeds the legacy BPNV (T_flame 2800)')

check(app.fileManager.load(ric), 'motor loaded through FileManager')
cm = app.fileManager.getCurrentMotor()
mat = cm.igniterData['pyrogen']
check(mat == lib['BPNV'] and mat['gas_mass_fraction'] == 0.781,
      'loaded motor carries the CURRENT BPNV (gas_mass_fraction 0.781)')
baseline = copy.deepcopy(cm.getDict())          # upgraded, nothing edited yet

win = app.window
editor = win.ui.motorEditor
app.preferencesManager.setActiveSolver(SOLVER)
names = app.pyrogenManager.getNames()


def open_config():
    editor.loadMotorConfig(app.fileManager.getCurrentMotor(), SOLVER, names)
    check('basket_fill_fraction' in editor.propertyEditors,
          'config page shows the provider chamber fields')


# --- 1. a no-edit Apply must change nothing (D4) -----------------------------
open_config()
editor.apply()
after = app.fileManager.getCurrentMotor().getDict()
d = diff_paths(baseline, after)
print('no-edit Apply diff paths:', d)
check(not [p for p in d if p.startswith('/igniter')],
      'no-edit Apply leaves the igniter block bit-identical (D4)')
# openMotor's own Apply stores the per-solver run-config (solverConfigs); that
# is pre-existing, non-igniter behaviour -> the baseline for later checks.
check(d == ['/solverConfigs/srm_1d-transient (added/removed)'],
      'the only other change a no-edit Apply makes is solverConfigs[srm_1d]')
baseline = copy.deepcopy(after)

# --- 2. edit two chamber fields + re-pick the same material ------------------
open_config()
editor.propertyEditors['basket_fill_fraction'].setValue(0.4)
editor.propertyEditors['mass'].setValue(0.002)
check(editor.pyrogenSelector.currentText() == 'BPNV', 'picker shows BPNV')
editor.apply()
edited = app.fileManager.getCurrentMotor()
out = os.path.join(TMP, 'edited.ric')
app.fileManager.fileName = out
app.fileManager.save()

reloaded = fileIO.loadFile(out, fileIO.fileTypes.MOTOR)
ig = reloaded['igniter']
check(ig['pyrogen']['gas_mass_fraction'] == 0.781 and ig['pyrogen'] == lib['BPNV'],
      'reloaded material is the current BPNV (gas_mass_fraction 0.781)')
check(abs(ig['basket_fill_fraction'] - 0.4) < 1e-9, 'edited basket_fill_fraction persisted: %r'
      % ig['basket_fill_fraction'])
check(abs(ig['mass'] - 0.002) < 1e-9, 'edited mass persisted: %r' % ig['mass'])
expect = copy.deepcopy(baseline)
expect['igniter']['basket_fill_fraction'] = ig['basket_fill_fraction']
expect['igniter']['mass'] = ig['mass']
d = diff_paths(expect, reloaded)
check(d == [], 'no other key changed vs the upgraded baseline: %s' % d)

# --- 3. a real library pick (different material) ------------------------------
app.fileManager.load(out)
editor.loadMotorConfig(app.fileManager.getCurrentMotor(), SOLVER, names)
editor.pyrogenSelector.setCurrentText('MTV')
editor.apply()
out2 = os.path.join(TMP, 'picked.ric')
app.fileManager.fileName = out2
app.fileManager.save()
r2 = fileIO.loadFile(out2, fileIO.fileTypes.MOTOR)['igniter']
check(r2['pyrogen'] == lib['MTV'], 'library pick deep-copied MTV into the block')
check(all(r2[k] == ig[k] for k in ig if k != 'pyrogen'),
      'chamber keys unchanged by the pick')

# --- 4. the pyrogen LIBRARY editor (D4 write-back) ----------------------------
menu = app.pyrogenManager._ensureMenu()
mgr = app.pyrogenManager
mgr.pyrogens[0]['Cp_gas'] = None            # schema-unknown / None key
mgr.pyrogens[0]['extra'] = 'keep'
before = copy.deepcopy(mgr.pyrogens)
menu.setupPyroList()
menu.listWidget.setCurrentRow(0)
menu.editPyro()
menu.pyroEditor.apply()                      # no edit
check(mgr.pyrogens == before, 'library editor: no-edit Apply leaves every entry bit-identical')
menu.listWidget.setCurrentRow(0)
menu.editPyro()
menu.pyroEditor.propertyEditors['gas_mass_fraction'].setValue(0.7)
menu.pyroEditor.apply()
check(diff_paths(before, mgr.pyrogens) == ['[0]/gas_mass_fraction'],
      'library editor: only the edited key changed')
check(abs(mgr.pyrogens[0]['gas_mass_fraction'] - 0.7) < 1e-12 and
      mgr.pyrogens[0]['extra'] == 'keep' and mgr.pyrogens[0]['Cp_gas'] is None,
      'library editor: edited value stored, unknown/None keys kept')

# --- 5. committed motors/ untouched ------------------------------------------
st = subprocess.run(['git', '-C', ROOT, 'status', '--short', 'motors'],
                    capture_output=True, text=True).stdout
check(st.strip() == '', 'srm_1d worktree motors/ untouched by the smoke')
print('GUI SMOKE: PASS')
