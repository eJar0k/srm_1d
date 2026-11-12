# What the code is meant to compute

This is the intended physics of the simulator, stated without reference to how
the code does it. Where the code disagrees with this document, the
disagreement is a finding, unless physics or a cited source shows that this
document is the one that is wrong. That is a finding too.

Units are SI throughout the numerical core.

## 1. Problem

A solid rocket motor: a cylindrical case holding one or more propellant grain
segments, with a nozzle at the aft end. Burning propellant surfaces line the
gas-filled port. The code simulates the internal gas flow and pressure from
igniter firing to burnout:
- the port starts at ambient pressure and temperature, filled with ambient gas;
- an igniter injects hot gas, which heats the propellant surfaces;
- each surface ignites when its temperature reaches an ignition temperature,
  then burns and adds gas;
- the gas leaves through the nozzle.

The outputs that matter most are the head-end pressure history, the thrust
history, and in particular the startup pressure transient.

## 2. Flow model

Quasi-one-dimensional, unsteady, compressible flow along the motor axis x,
from the head end (x = 0) to the nozzle entrance:
- the port cross-section A(x, t) varies along the axis and grows in time as
  the propellant surface regresses;
- the gas is an ideal-gas mixture of several species: ambient air, the
  combustion products of each propellant, and the igniter products;
- each species k has a constant molar mass, a constant C_p,k and γ_k, and a
  sensible enthalpy h_k = C_p,k T;
- the mixture R, C_p and γ follow from the local mass fractions Y_k (ideal-gas
  mixing);
- equation of state: P = ρ R_mix T.

Conservation, per unit length of the axis:
- **Mass:** the mixture and each species are conserved. Mass enters through
  three routes:
  - propellant burning, at ρ_p · r per unit burning area, where ρ_p is the
    propellant density and r the local burn rate;
  - the burning end faces of grain segments;
  - the igniter.
- **Momentum:** axial momentum, with:
  - the pressure-area force;
  - wall friction (Darcy friction factor, smooth or rough pipe);
  - any axial momentum injected by the igniter jet.
  Gas produced at the burning surface enters normal to the axis and carries no
  axial momentum.
- **Energy:** total energy, internal plus kinetic. Gas leaves a burning surface
  at its species' flame temperature. Wall heat loss and radiation exchange,
  where modelled, are sinks or sources.

The port wall moves outward as the propellant burns. **Derive the correct
conservation form for a control volume whose boundary regresses, and check the
code's source terms against your derivation.** Do not assume the code's form
is right.

Total mass, energy and species inventories must close against the integrated
sources and the nozzle outflow, to within the time-discretization error, and to
round-off for any terms the scheme claims to treat exactly.

## 3. Discretization (intended)

- Uniform cells of width dx along the axis, on a staggered grid: scalars (P,
  T, ρ, Y_k) at cell centres, velocity at faces.
- Grain segment lengths and gaps are snapped to whole cells.
- A pressure-based predictor-corrector (PISO) step:
  1. a momentum predictor;
  2. two pressure corrections from the continuity constraint, each a
     tridiagonal system;
  3. the energy and species equations.
  The density change is linked to the pressure change through the equation of
  state, so the scheme must stay valid for compressible, transonic flow.
- Adaptive time step: an acoustic CFL limit dt = CFL · dx / max(|u| + a), plus
  limits from the source terms.
- The burn rate and the geometry may be updated less often than every flow
  step. Any such sub-cycling must not break conservation.

## 4. Boundaries and end volumes

- **Head end:** a closed wall (zero velocity).
  - There may be a small free volume between the head closure and the first
    propellant face. It is lumped into the first cell's storage volume, not
    resolved with extra cells.
  - The igniter, in its plenum form, injects mass, enthalpy and axial momentum
    here.
- **Aft end:** the nozzle convergent volume is lumped into the last cell's
  storage volume.
- **Nozzle:** the outlet face is the nozzle throat, of area A_t (which can grow
  through throat erosion).
  - Flow through it is quasi-steady and isentropic from the **stagnation
    state** of the gas entering the nozzle.
  - It is choked when the back pressure is below the critical ratio of that
    stagnation pressure, and subsonic isentropic otherwise.
  - Reverse inflow from the ambient is allowed when the chamber pressure falls
    below ambient.
  - The pressure equation needs the derivative of the boundary mass flow with
    respect to cell pressure, and it must be consistent with the flow function
    actually used.
- Thrust uses a nozzle thrust-coefficient model. It is post-processing and has
  low priority in this audit.

## 5. Burn rate

- **Base rate:** Saint-Robert, r_0 = a P^n. The propellant may have several
  (a, n) tabs over pressure ranges; the lookup selects the tab containing the
  local pressure.
- **Erosive increment:** Ma et al. 2020 (`references/ma2020_erosive.pdf`), a
  heat-transfer model:
  - a friction factor (Haaland) feeds a Nusselt correlation (Gnielinski 2013,
    `references/gnielinski2013.pdf`) for the convective coefficient h_0;
  - h_0 carries a transpiration correction for surface blowing;
  - a surface energy balance gives the erosive increment;
  - the total rate r = r_0 + r_e is implicit in r and is solved iteratively.
  Check every equation, variable definition and unit against the paper. Note
  where the paper is ambiguous and which reading the code takes.
- Inputs are the surface roughness (default 32 µm), a constant κ (default 0.44),
  and the gas transport properties (μ, k, C_p) per propellant tab.
- Only ignited cells burn. Burnout ends the burning when the web is consumed,
  with no loss of mass conservation.

## 6. Ignition and solid heating

- Each grain cell has a propellant surface temperature T_s, starting at the
  ambient temperature.
- **Solid:** transient 1-D heat conduction into the propellant, treated as
  semi-infinite on the ignition timescale. Properties: k_s, ρ_p and C_ps.
  Boundary condition: the convective flux h (T_drive − T_s).
- **The cell ignites when T_s reaches T_ignition** (default 756 K), and
  thereafter burns at the rate in section 5.
- **Driver temperature:** the gas recovery temperature,
  T_aw = T + r_f u² / (2 C_p), with recovery factor r_f = Pr^(1/3) for
  turbulent flow and Pr^(1/2) for laminar flow.
- **Bore heat-transfer coefficient h:** the Peretz correlation for igniter-gas
  convective heating of the port. The source is Peretz, Caveny, Kuo &
  Summerfield, Princeton AMS Report 1100 (1973), `references/peretz1974.pdf`,
  Eqs. III-24, III-25 and III-31. The companion AIAA paper 72-1119 is
  `references/peretz1972.pdf`.
  - Nu = 0.0256 (x/d)^(−0.1) Re^0.8 Pr^0.4 (A_p/A_t)^0.4;
  - A_p is the local port area and A_t the throat area.
- **Near the head end, with the plenum igniter:** h is the larger of the bore
  value and the Peretz igniter-jet correlations (Eqs. III-28 to III-30 in the
  same report).
  - The jet Reynolds number is built from the igniter mass flow and the port
  diameter.
  - The jet correlations apply only over the first ten port diameters.
- Check that each correlation is transcribed correctly: definitions of Re, d,
  x, the properties and the temperature at which they are evaluated, and units.
- Optional heat paths (radiation from burning neighbours or igniter pellets,
  igniter surface flux) are off by default. They are low priority.

## 7. Igniter

- The igniter material (the "pyrogen") burns by its own Saint-Robert law and
  produces a gas species with its own flame temperature and properties. Only
  a gas fraction of its mass becomes gas.
- **Topologies:**
  - *forward plenum:* a 0-D chamber at the head end that vents through an
    orifice, choked or subsonic;
  - *head basket* or *aft basket:* uncontained pellets in the port, burning at
    the local port pressure of the cells they occupy.
- **Pellets:** MIL-P-46994B Type I-A cylinders (`references/MIL-P-46994B.pdf`).
  The burning area follows the exact burn-back of a cylinder.
- The automatic charge-sizing and retry logic is a policy wrapper around the
  simulation. It is out of scope except for crashes or conservation effects.

## 8. Grain geometry

- **BATES:** a cylindrical bore with burning end faces unless they are
  inhibited.
- **Conical:** a linearly tapered bore.
- More complex cross-sections (fast-marching regression tables) need an
  external openMotor checkout that this bundle does not include. Out of scope.
- Per cell, from the local web burned: the burning perimeter, the port area,
  and the remaining web. Check the analytic perimeter and area against closed
  forms.
- An uninhibited end face burns as an annulus that recedes axially. Its gas
  enters the cells next to the face, with weights that sum to one.

## 9. Not in scope

The openMotor GUI plugin, plotting, channels and visualization, the sensitivity
and sweep tools, file import and export (except where they change what the
solver receives), and calibration of default parameter values.
