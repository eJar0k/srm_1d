"""
nozzle.py — Nozzle Performance Model
======================================

Computes thrust, Isp, and thrust coefficient from chamber conditions
and nozzle geometry. Throat evolution (erosion and slag buildup) is
pressure-dependent. The Nozzle data class and method names follow
openMotor's `motorlib.nozzle.Nozzle` so the adapter layer stays thin.

Physics from Sutton's Rocket Propulsion Elements (Ch. 3):
    - Ideal thrust coefficient from isentropic expansion
    - Exit pressure from area ratio via Newton iteration
    - Divergence loss from conical half-angle
    - Throat-aspect-ratio loss (RasAero "Departures from Ideal Performance")
    - Skin friction loss (constant 0.99, openMotor convention)
    - User-set efficiency multiplier

Adjusted thrust coefficient (openMotor Nozzle.getAdjustedThrustCoeff):

    CF_adj = divLoss × throatLoss × efficiency
           × (skinLoss × CF_ideal + (1 − skinLoss))

The skin-loss term only attenuates the part of CF above 1 (the
momentum/pressure component above the static-pressure baseline).

Throat model:
    EROSION (throat grows): proportional to chamber pressure.
        erosion_rate [m/s] = erosion_coeff [μm/(s·MPa)] × P [MPa] × 1e-6
    SLAG BUILDUP (throat shrinks): inversely proportional to pressure.
        slag_rate [m/s]    = slag_coeff [(m·MPa)/s] / P [MPa]
    Net diameter change:
        dD/dt = 2 × (erosion_rate − slag_rate)

    The simulation integrates D_throat in-loop (feeding back into the
    nozzle BC) when erosion or slag is nonzero. Post-processing in
    compute_motor_performance replays the history for thrust.

Units:
    srm_1d keeps engineering-readable units internally
    (μm/(s·MPa) for erosion_coeff, (m·MPa)/s for slag_coeff). The
    openMotor adapter converts at the boundary.
"""

import math
import numpy as np
from dataclasses import dataclass

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


# Skin-loss constant from openMotor (motorlib/nozzle.py:Nozzle.getSkinLosses).
# Hardcoded because users typically can't measure this themselves.
_SKIN_LOSS = 0.99


@dataclass
class Nozzle:
    """
    Nozzle configuration. Field names match openMotor's
    motorlib.nozzle.Nozzle (snake_case'd).

    Parameters
    ----------
    D_throat : float
        Initial throat diameter [m]. (openMotor: 'throat')
    D_exit : float
        Exit diameter [m]. Determines expansion ratio. (openMotor: 'exit')
    efficiency : float
        User-set c* multiplier (openMotor 'efficiency'). NOT a roll-up
        of geometric losses — those are computed separately.
        Typical: 0.95–0.98.
    div_angle : float
        Divergence half-angle [deg]. (openMotor: 'divAngle')
        Typical: 12–15. Set 0 for ideal (no divergence loss).
    conv_angle : float
        Convergence half-angle [deg]. (openMotor: 'convAngle')
        Stored for completeness/openMotor-alignment; not yet consumed
        by the thrust calc.
    throat_length : float
        Throat length [m]. (openMotor: 'throatLength')
        Drives the throat-aspect-ratio loss factor.
    erosion_coeff : float
        Throat erosion coefficient [μm/(s·MPa)].
        (openMotor: 'erosionCoeff', stored in m/(s·Pa); adapter converts.)
        Typical: 0 for non-eroding throats; 0.05–0.5 for graphite with
        aluminized propellant. Positive = throat grows.
    slag_coeff : float
        Slag (alumina deposition) coefficient [(m·MPa)/s].
        (openMotor: 'slagCoeff', same units.)
        Typical: 0 for non-aluminized propellant. Positive = throat shrinks.
    """
    D_throat: float
    D_exit: float
    efficiency: float = 0.95
    div_angle: float = 15.0
    conv_angle: float = 30.0
    throat_length: float = 0.0
    erosion_coeff: float = 0.0
    slag_coeff: float = 0.0

    # --------------------------------------------------------------
    # Geometry / area
    # --------------------------------------------------------------
    def throat_area(self, d_throat=0.0):
        """Throat area [m²]. d_throat adds erosion/slag offset to D_throat."""
        D = self.D_throat + d_throat
        return np.pi / 4.0 * D * D

    def exit_area(self):
        """Exit area [m²]."""
        return np.pi / 4.0 * self.D_exit * self.D_exit

    @property
    def expansion_ratio(self):
        """Initial expansion ratio ε = A_exit / A_throat (uses initial D_throat)."""
        return self.exit_area() / self.throat_area(0.0)

    # --------------------------------------------------------------
    # Loss factors (openMotor naming: divergence/throat/skin)
    # --------------------------------------------------------------
    def divergence_losses(self):
        """λ = (1 + cos(α)) / 2 for a conical exhaust."""
        return (1.0 + math.cos(math.radians(self.div_angle))) / 2.0

    def throat_losses(self, d_throat=0.0):
        """
        Loss from throat-aspect ratio (throat_length / D_throat).
        From RasAero "Departures from Ideal Performance":
            aspect > 0.45      → 0.95
            aspect ≤ 0.45      → 0.99 − 0.0333 × aspect
        """
        D = self.D_throat + d_throat
        if D <= 0.0:
            return 0.95
        aspect = self.throat_length / D
        if aspect > 0.45:
            return 0.95
        return 0.99 - 0.0333 * aspect

    def skin_losses(self):
        """Constant 0.99 (openMotor convention; see _SKIN_LOSS)."""
        return _SKIN_LOSS

    # --------------------------------------------------------------
    # Pressure / thrust coefficient (openMotor signature alignment)
    # --------------------------------------------------------------
    def exit_pressure(self, gamma, P_chamber):
        """Solve for exit pressure given chamber pressure and γ."""
        ratio = exit_pressure_from_expansion_ratio(gamma, self.expansion_ratio)
        return ratio * P_chamber

    def ideal_thrust_coeff(self, P_chamber, P_ambient, gamma, d_throat=0.0,
                           exit_pres=None):
        """
        Ideal thrust coefficient (Sutton eq. 3-30, openMotor formulation).
        d_throat is the change in throat diameter due to erosion or slag.
        """
        if P_chamber == 0.0:
            return 0.0
        if exit_pres is None:
            exit_pres = self.exit_pressure(gamma, P_chamber)
        A_e = self.exit_area()
        A_t = self.throat_area(d_throat)

        gm1 = gamma - 1.0
        gp1 = gamma + 1.0
        term1 = (2.0 * gamma * gamma) / gm1
        term2 = (2.0 / gp1) ** (gp1 / gm1)
        term3 = 1.0 - (exit_pres / P_chamber) ** (gm1 / gamma)
        if term3 < 0.0:
            term3 = 0.0
        momentum = (term1 * term2 * term3) ** 0.5
        pressure = ((exit_pres - P_ambient) * A_e) / (A_t * P_chamber)
        return momentum + pressure

    def adjusted_thrust_coeff(self, P_chamber, P_ambient, gamma, d_throat=0.0,
                              exit_pres=None):
        """
        Adjusted CF including geometric and skin losses + user efficiency.
        Mirrors openMotor's Nozzle.getAdjustedThrustCoeff:

            CF_adj = divLoss × throatLoss × efficiency
                   × (skinLoss × CF_ideal + (1 − skinLoss))
        """
        cf_ideal = self.ideal_thrust_coeff(
            P_chamber, P_ambient, gamma, d_throat, exit_pres
        )
        div_loss = self.divergence_losses()
        throat_loss = self.throat_losses(d_throat)
        skin_loss = self.skin_losses()
        return (div_loss * throat_loss * self.efficiency
                * (skin_loss * cf_ideal + (1.0 - skin_loss)))

    # --------------------------------------------------------------
    # Convenience predicates
    # --------------------------------------------------------------
    @property
    def has_throat_change(self):
        """True if either erosion or slag coefficients are nonzero."""
        return self.erosion_coeff != 0.0 or self.slag_coeff != 0.0


# ================================================================
# Numba-compiled scalar functions (used by the time loop and by
# compute_thrust_history; not on the Nozzle class because @njit
# cannot reach into Python dataclass instances)
# ================================================================

@njit(cache=True)
def exit_pressure_from_expansion_ratio(gamma, expansion_ratio):
    """
    Pressure ratio P_e/P_c from area ratio via isentropic Mach Newton iteration.
    """
    if expansion_ratio <= 1.0:
        return 1.0

    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    M = 1.0 + 0.5 * (expansion_ratio - 1.0)

    for _ in range(50):
        term = 1.0 + 0.5 * gm1 * M * M
        A_ratio = (1.0 / M) * ((2.0 / gp1) * term) ** (gp1 / (2.0 * gm1))

        dM = M * 1e-6
        M2 = M + dM
        term2 = 1.0 + 0.5 * gm1 * M2 * M2
        A_ratio2 = (1.0 / M2) * ((2.0 / gp1) * term2) ** (gp1 / (2.0 * gm1))
        dAdM = (A_ratio2 - A_ratio) / dM

        correction = (A_ratio - expansion_ratio) / dAdM
        M = M - correction
        M = max(M, 1.001)

        if abs(A_ratio - expansion_ratio) < 1e-8:
            break

    P_ratio = (1.0 + 0.5 * gm1 * M * M) ** (-gamma / gm1)
    return P_ratio


@njit(cache=True)
def ideal_thrust_coefficient(gamma, expansion_ratio, P_chamber, P_ambient):
    """
    Ideal CF (Sutton eq. 3-30). Returns (CF, P_exit).
    """
    if P_chamber < 1e3:
        return 0.0, P_ambient

    gm1 = gamma - 1.0
    gp1 = gamma + 1.0

    pe_pc = exit_pressure_from_expansion_ratio(gamma, expansion_ratio)
    P_exit = pe_pc * P_chamber

    momentum_term = (2.0 * gamma * gamma / gm1) * (
        (2.0 / gp1) ** (gp1 / gm1)
    ) * (1.0 - pe_pc ** (gm1 / gamma))

    if momentum_term < 0.0:
        momentum_term = 0.0

    C_F = np.sqrt(momentum_term) + (pe_pc - P_ambient / P_chamber) * expansion_ratio

    return max(C_F, 0.0), P_exit


@njit(cache=True)
def _throat_aspect_loss(throat_length, D_throat):
    """Throat-aspect-ratio loss factor (RasAero), Numba scalar version."""
    if D_throat <= 0.0:
        return 0.95
    aspect = throat_length / D_throat
    if aspect > 0.45:
        return 0.95
    return 0.99 - 0.0333 * aspect


@njit(cache=True)
def compute_thrust_isp(
    P_chamber, gamma, expansion_ratio,
    div_loss, efficiency, throat_loss, skin_loss,
    A_throat, P_ambient, g0,
    c_star,
):
    """
    Single-instant thrust, delivered CF, Isp using openMotor's adjusted-CF
    formula:
        CF_adj = divLoss × throatLoss × efficiency
              × (skinLoss × CF_ideal + (1 − skinLoss))
    """
    C_F_ideal, P_exit = ideal_thrust_coefficient(
        gamma, expansion_ratio, P_chamber, P_ambient
    )

    C_F_delivered = (div_loss * throat_loss * efficiency
                     * (skin_loss * C_F_ideal + (1.0 - skin_loss)))
    thrust = C_F_delivered * P_chamber * A_throat

    if c_star > 0.0:
        Isp = C_F_delivered * c_star / g0
    else:
        Isp = 0.0

    return thrust, C_F_delivered, Isp, P_exit


@njit(cache=True)
def compute_thrust_history(
    time_arr, P_head_arr, N_points,
    gamma, D_throat_initial, D_exit,
    erosion_coeff, slag_coeff,
    div_loss, efficiency, throat_length, skin_loss,
    P_ambient, c_star, g0,
):
    """
    Thrust, CF, Isp, P_exit, D_throat over the pressure history.
    Throat diameter is integrated forward via forward Euler:
        e_rate = erosion_coeff [μm/(s·MPa)] × P [MPa] × 1e-6
        s_rate = slag_coeff   [(m·MPa)/s]   / P [MPa]
        dD/dt = 2 × (e_rate − s_rate)
    """
    thrust = np.zeros(N_points)
    C_F = np.zeros(N_points)
    Isp = np.zeros(N_points)
    P_exit = np.zeros(N_points)
    D_throat_arr = np.zeros(N_points)

    PI = 3.141592653589793
    A_exit = PI / 4.0 * D_exit * D_exit

    D_t = D_throat_initial
    D_throat_arr[0] = D_t

    for i in range(N_points):
        P = P_head_arr[i]
        P_MPa = P / 1e6

        if i > 0:
            dt = time_arr[i] - time_arr[i - 1]
            e_rate = erosion_coeff * 1e-6 * P_MPa
            if P_MPa > 0.01:
                s_rate = slag_coeff / P_MPa
            else:
                s_rate = 0.0
            D_t = D_t + 2.0 * (e_rate - s_rate) * dt
            D_t = max(D_t, 1e-6)
            D_throat_arr[i] = D_t

        A_t = PI / 4.0 * D_t * D_t
        eps = A_exit / A_t
        throat_loss = _throat_aspect_loss(throat_length, D_t)

        thrust[i], C_F[i], Isp[i], P_exit[i] = compute_thrust_isp(
            P, gamma, eps,
            div_loss, efficiency, throat_loss, skin_loss,
            A_t, P_ambient, g0, c_star,
        )

    return thrust, C_F, Isp, P_exit, D_throat_arr


# ================================================================
# Post-processing: combine simulation result + Nozzle into perf dict
# ================================================================

def compute_motor_performance(result, nozzle, propellant, P_ambient=None):
    """
    Compute thrust/Isp from simulation results.

    Parameters
    ----------
    result : dict
        Output from run_simulation (keys: 'time', 'P_head', optionally
        'D_throat', 'P_ambient').
    nozzle : Nozzle
    propellant : Propellant
    P_ambient : float or None
        Ambient pressure [Pa]. If None, read from `result['P_ambient']`
        (which run_simulation populates from its kwarg). Falls back to
        sea level (101325) if neither is set.
    """
    from .propellant import critical_flow_function, R_UNIVERSAL

    if P_ambient is None:
        P_ambient = result.get('P_ambient', 101325.0)

    time_arr = result['time']
    P_arr = result['P_head']
    N = len(time_arr)
    g0 = 9.80665

    rep_tab = propellant.representative_tab()
    gamma = rep_tab.gamma
    Gamma = critical_flow_function(gamma)
    R_spec = R_UNIVERSAL / rep_tab.molecular_weight
    c_star = np.sqrt(R_spec * rep_tab.T_flame) / Gamma

    div_loss = nozzle.divergence_losses()
    skin_loss = nozzle.skin_losses()

    thrust, C_F, Isp, P_exit, D_throat_arr = compute_thrust_history(
        time_arr, P_arr, N,
        gamma, nozzle.D_throat, nozzle.D_exit,
        nozzle.erosion_coeff, nozzle.slag_coeff,
        div_loss, nozzle.efficiency, nozzle.throat_length, skin_loss,
        P_ambient, c_star, g0,
    )

    # If the simulation already computed a throat history (because
    # erosion/slag was active during the sim), use that — it's what
    # actually drove the nozzle BC.
    if 'D_throat' in result and len(result['D_throat']) == N:
        D_throat_arr = result['D_throat']
        PI = 3.141592653589793
        A_exit_val = PI / 4.0 * nozzle.D_exit ** 2
        for i in range(N):
            D_t = D_throat_arr[i]
            A_t = PI / 4.0 * D_t * D_t
            eps = A_exit_val / A_t
            throat_loss = _throat_aspect_loss(nozzle.throat_length, D_t)
            thrust[i], C_F[i], Isp[i], P_exit[i] = compute_thrust_isp(
                P_arr[i], gamma, eps,
                div_loss, nozzle.efficiency, throat_loss, skin_loss,
                A_t, P_ambient, g0, c_star,
            )

    total_impulse = 0.0
    for i in range(N - 1):
        total_impulse += 0.5 * (thrust[i] + thrust[i + 1]) * (time_arr[i + 1] - time_arr[i])

    peak_thrust = np.max(thrust) if N > 0 else 0.0
    threshold = 0.05 * peak_thrust
    burning = thrust > threshold
    if np.any(burning):
        burn_indices = np.where(burning)[0]
        burn_time = time_arr[burn_indices[-1]] - time_arr[burn_indices[0]]
        avg_thrust = total_impulse / burn_time if burn_time > 0 else 0.0
    else:
        burn_time = 0.0
        avg_thrust = 0.0

    avg_Isp = float(np.mean(Isp[burning])) if np.any(burning) else 0.0
    designation = _motor_designation(total_impulse, avg_thrust)

    if np.any(burning):
        avg_CF = float(np.mean(C_F[burning]))
        peak_CF = float(np.max(C_F[burning]))
    else:
        avg_CF = 0.0
        peak_CF = 0.0

    D_throat_final = D_throat_arr[-1] if len(D_throat_arr) > 0 else nozzle.D_throat
    throat_change_mm = (D_throat_final - nozzle.D_throat) * 1000

    return {
        'thrust': thrust,
        'C_F': C_F,
        'Isp': Isp,
        'P_exit': P_exit,
        'D_throat': D_throat_arr,
        'total_impulse': total_impulse,
        'average_thrust': avg_thrust,
        'peak_thrust': peak_thrust,
        'average_Isp': avg_Isp,
        'burn_time': burn_time,
        'motor_designation': designation,
        'c_star': c_star,
        'average_C_F': avg_CF,
        'peak_C_F': peak_CF,
        'D_throat_initial': nozzle.D_throat,
        'D_throat_final': D_throat_final,
        'throat_change_mm': throat_change_mm,
        'P_ambient': P_ambient,
    }


def _motor_designation(total_impulse_Ns, avg_thrust_N):
    """NAR/TRA letter designation. A = 1.26-2.50 N·s, doubles each letter."""
    letters = "ABCDEFGHIJKLMNOP"
    if total_impulse_Ns <= 0:
        return "N/A"
    log_impulse = np.log2(total_impulse_Ns / 1.26)
    idx = int(np.floor(log_impulse))
    idx = max(0, min(idx, len(letters) - 1))
    letter = letters[idx]
    return f"{letter}{avg_thrust_N:.0f}"


def print_performance_summary(perf, nozzle):
    """Formatted performance summary."""
    print(f"\n{'='*65}")
    print(f"MOTOR PERFORMANCE")
    print(f"  Designation:    {perf['motor_designation']}")
    print(f"  Total Impulse:  {perf['total_impulse']:.1f} N·s")
    print(f"  Average Thrust: {perf['average_thrust']:.1f} N")
    print(f"  Peak Thrust:    {perf['peak_thrust']:.1f} N")
    print(f"  Delivered Isp:  {perf['average_Isp']:.1f} s")
    print(f"  c*:             {perf['c_star']:.1f} m/s")
    print(f"  Avg C_F:        {perf['average_C_F']:.3f}")
    print(f"  Peak C_F:       {perf['peak_C_F']:.3f}")
    print(f"  Burn Time:      {perf['burn_time']:.2f} s")
    print(f"  Expansion Ratio:{nozzle.expansion_ratio:.2f} (initial)")
    print(f"  Divergence Loss:{nozzle.divergence_losses():.4f}")
    print(f"  Throat Loss:    {nozzle.throat_losses():.4f}  "
          f"(L_t/D_t = {nozzle.throat_length/max(nozzle.D_throat,1e-9):.3f})")
    print(f"  Skin Loss:      {nozzle.skin_losses():.4f}")
    print(f"  Efficiency:     {nozzle.efficiency:.3f}")
    if abs(perf.get('throat_change_mm', 0)) > 0.001:
        D_i = perf['D_throat_initial'] * 1000
        D_f = perf['D_throat_final'] * 1000
        delta = perf['throat_change_mm']
        direction = "eroded" if delta > 0 else "slagged"
        print(f"  Throat:         {D_i:.2f} -> {D_f:.2f} mm "
              f"({direction} {abs(delta):.3f} mm)")
    print(f"{'='*65}")
