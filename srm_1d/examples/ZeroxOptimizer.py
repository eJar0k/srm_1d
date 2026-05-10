"""Legacy entry point for Zerox LHS optimization.

v0.7.0 removed the old exponential igniter knobs. Use the pyrogen-based
Zerox LHS driver instead:

    python -m srm_1d.examples.zerox_lhs
"""

from srm_1d.examples.zerox_lhs import main


if __name__ == "__main__":
    main()
