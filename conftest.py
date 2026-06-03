"""Pytest bootstrap: ensure the repo root is importable.

Flat layout — the ``srm_1d`` package and the dev-only ``tests`` / ``examples``
packages all live at the repo root. Putting the repo root on ``sys.path``
lets the suite import ``srm_1d`` (also available via ``pip install -e .``)
as well as the dev packages ``tests.*`` and ``examples.*`` regardless of the
editable-install mode.
"""
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
