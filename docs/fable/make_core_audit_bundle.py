"""Build the blind core-audit bundle for Fable task 2.

The bundle is a folder OUTSIDE the repo holding only:
  srm_1d/      the package at a git ref, minus every .md file, the GUI plugin
               and the dev tools (the memory/time guard stays)
  motors/      three analytic-grain motors (BATES / Conical), committed versions
  references/  the primary sources for the implemented closures
  SPEC.md      the intended physics (docs/fable/core_audit/SPEC.md at the ref)
  TASK.md      the auditor's instructions (docs/fable/core_audit/TASK.md)
  probes/      empty package for the auditor's probes
It is a fresh git repo (branch `bundle`, one commit), so no project history,
docs, DEVNOTES, CLAUDE.md or memory goes with it.

Usage (repo root, pyenv 3.10.5 python):
    python docs/fable/make_core_audit_bundle.py [--ref HEAD] [--out C:/tmp/fable_core_audit]
"""
import argparse
import io
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

MOTORS = ['motors/hasegawa_a.ric', 'motors/machbusterNew.ric',
          'motors/example_bates.ric']

# (source under docs/references, name in the bundle). Copied from the working
# tree: most reference PDFs are not tracked in git.
REFERENCES = [
    ('International Journal of Aerospace Engineering - 2020 - Ma - A New '
     'Erosive Burning Model of Solid Propellant Based on Heat.pdf',
     'ma2020_erosive.pdf'),
    ('gnielinski2013.pdf', 'gnielinski2013.pdf'),
    ('peretz1974.pdf', 'peretz1974.pdf'),   # Princeton AMS Rep. 1100 (1973)
    ('peretz1972.pdf', 'peretz1972.pdf'),   # AIAA 72-1119
    ('MIL-P-46994B.pdf', 'MIL-P-46994B.pdf'),
]

# Removed from the package copy: the openMotor GUI plugin (out of scope) and
# the developers' own diagnostic / verification tools (the auditor writes its
# own). The guard stays: the auditor must run every job through it.
DROP = [
    'srm_1d/srm1d_plugin.py',
    'srm_1d/tools/ignition_diagnostics.py',
    'srm_1d/tools/refresh_ric_igniters.py',
    'srm_1d/tools/sensitivity.py',
    'srm_1d/tools/solid_replay.py',
    'srm_1d/tools/solid_verify.py',
]

GITIGNORE = "__pycache__/\n*.pyc\n*.nbi\n*.nbc\n*.npz\n*.png\n"


def git(*args, **kw):
    return subprocess.run(['git', '-C', str(REPO), *args], check=True,
                          capture_output=True, **kw).stdout


def extract(ref, paths, out):
    data = git('archive', '--format=tar', ref, *paths)
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        tar.extractall(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--ref', default='HEAD')
    ap.add_argument('--out', default='C:/tmp/fable_core_audit')
    a = ap.parse_args()
    out = Path(a.out).resolve()
    if REPO in out.parents or out == REPO:
        sys.exit(f'--out must be outside the repo ({REPO})')
    if out.exists() and any(out.iterdir()):
        sys.exit(f'{out} exists and is not empty; pick another --out or '
                 f'remove it yourself')
    out.mkdir(parents=True, exist_ok=True)
    sha = git('rev-parse', a.ref, text=True).strip()

    extract(sha, ['srm_1d', *MOTORS], out)
    spec_dir = 'docs/fable/core_audit'
    extract(sha, [f'{spec_dir}/SPEC.md', f'{spec_dir}/TASK.md'], out)
    for name in ('SPEC.md', 'TASK.md'):
        shutil.move(str(out / spec_dir / name), str(out / name))
    shutil.rmtree(out / 'docs')

    for md in (out / 'srm_1d').rglob('*.md'):
        md.unlink()
    for rel in DROP:
        p = out / rel
        if p.exists():
            p.unlink()
    for d in (out / 'srm_1d').rglob('__pycache__'):
        shutil.rmtree(d)

    refs = out / 'references'
    refs.mkdir()
    for src, dst in REFERENCES:
        shutil.copy2(REPO / 'docs' / 'references' / src, refs / dst)

    (out / 'probes').mkdir()
    (out / 'probes' / '__init__.py').write_text('')
    (out / '.gitignore').write_text(GITIGNORE)

    env = dict(os.environ, GIT_AUTHOR_NAME='bundle', GIT_COMMITTER_NAME='bundle',
               GIT_AUTHOR_EMAIL='bundle@localhost',
               GIT_COMMITTER_EMAIL='bundle@localhost')
    for cmd in (['init', '-q', '-b', 'bundle'], ['add', '-A'],
                ['commit', '-q', '-m', 'audit bundle']):
        subprocess.run(['git', '-C', str(out), *cmd], check=True, env=env)

    print(f'bundle: {out}  (srm_1d at {sha[:10]})')
    print('next: open a NEW Claude Code session in that folder and say:')
    print('  "Read TASK.md and SPEC.md, then do the task."')


if __name__ == '__main__':
    main()
