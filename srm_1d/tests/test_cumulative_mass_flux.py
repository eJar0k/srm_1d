"""
tests/test_cumulative_mass_flux.py — OBSOLETE shim
====================================================

This file originally held Phase B v1 kernel tests for the cumulative-G
+ Dittus-Boelter Re^0.8 augmentation formulation (commit 065d193 +
6360f62). That formulation was found to double-count with PISO's
local-Re tracking and amplify the ignition spike rather than smooth it.

Phase B was reformulated as flame-front-marker augmentation (gate on
recent upstream ignition rather than cumulative-G magnitude). The
replacement kernel `_compute_flame_front_augment` is exercised by
tests/test_flame_front_augment.py.

This file is kept as an empty marker so the obsolete-test commit
history is searchable. Future cleanups may remove it entirely.
"""
# No tests in this file — see test_flame_front_augment.py.
