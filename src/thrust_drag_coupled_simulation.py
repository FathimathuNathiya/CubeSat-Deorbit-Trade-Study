"""
Coupled low-thrust propulsion and atmospheric drag simulation for CubeSat
mission.

This module extends drag_deorbit_simulation.py by adding a retrograde
thrust term, active over a fixed burn window, on top of gravity and drag.
Propellant mass is tracked as a seventh state variable (in addition to
position and velocity), decreases according to the mass-flow rate determined
by the thrust-exhaust velocity relationship.

This addresses the second extension suggested by Dr. Dario Modenini
(University of Bologna): investigating coupled scenarios in which low-thrust
propulsion and atmospheric drag interact, rather than treating deorbit as 
a purely propulsion-driven process. The electric propulsion scenario below
directly addresses this (low-thrust, long burn duration); the green chemical
propulsion scenario is included as an additional comparison point, using a
comparatively high thrust, short-duration burn, to contrast against the
low-thrust case.

IMPORTANT FRAMING REMARK: the propellant masses used here (0.6223 kg green,
0.1118 kg electric) are taken directly from CubeSat_Deorbit_Analysis.pdf
(Table 3), which sized them for an idealized, instantaneous Hohmann-style
transfer from 500 km to 120 km. This simulation does not reproduce that
transfer - it spends the same propellant budget as a continuous low-thrust
burn instead, with real gravity and drag acting throughout. The two
models are deliberately different (preliminary analytical sizing vs.
higher-fidelity dynamic simulation), and the point of this extension is 
to compare how the same propellant allowance performs under each - not
to confirm the simulation reaches 120 km, which it will not over the 
short window simulated here.

A reentry event (see drag_deorbit_simulation.reentry_event) stops
integration once altitude drops to REENTRY_ALTITUDE_KM, for the same
reason it's needed there: continuing to integrate past reentry has no
physical meaning and produces runaway / negative altitude over long
simulated durations.

Shared constants, rotation matrices, and the geodetic conversion are
imported directly from drag_deorbit_simulation.py rather than redefined
here, so that both extensions remain consistent with a single source of
truth for the underlying orbital mechanics.
"""


import numpy as np
from datetime import timedelta
from pymsis import calculate
from scipy.integrate import solve_ivp

from mission_analysis import DRY_MASS_KG, ISP_GREEN_S, ISP_ELECTRIC_S
from paths import DATA_DIR
from drag_deorbit_simulation import (
    MU_EARTH,
    EARTH_EQUATORIAL_RADIUS,
    SIDEREAL_DAY,
    CROSS_SECTIONAL_AREA,
    DRAG_COEFFICIENT,
    ORBITAL_RADIUS,
    ECCENTRICITY,
    INCLINATION,
    RAAN,
    ARGUMENT_OF_PERIGEE,
    TRUE_ANOMALY,
    REFERENCE_START_DATE,
    REENTRY_ALTITUDE_KM,
    rotation_matrix_x,
    rotation_matrix_z,
    ecef_to_geodetic,
    reentry_event
)


# Propulsion parameters

STANDARD_GRAVITY = 9.80665                                             # Standard gravitational acceleration, used to convert Isp to exhaust velocity (m/s^2)
SIMULATION_DURATION_DAYS = 365   

# Propellant mass, from the analytical trade study (CubeSat_Deorbit_Analysis.pdf,
# Table 3) - see framing not above. Sized for the 500 km -> 120 km transfer
# delta-V of 109.08 m/s at 12 kg dry mass; applied here as a fixed fuel
# budget for a continuous burn instead.

GREEN_PROPELLANT_MASS = 0.6223                                         # report Table 3 (kg)
ELECTRIC_PROPELLANT_MASS = 0.1118                                      # report Table 3 (kg)

# Electric thrust: report Table 5, 30 W operating point (a representative
# mid-range point from the report's five-point power sweep).

ELECTRIC_THRUST = 0.005099                                             # (N)

# Green (chemical) thrust: representative of a small CubeSat monopropellant
# thruster (CU Aerospace MPUC / CMP-X), which demonstrated thrust > 0.5 N at
# Isp > 180 s in thrust-stand testing - closely matching this study's 
# ISP_GREEN_S assumption of 220 s.
# Source: NASA Small Spacecraft Technology State of the Art (In-Space
# Propulsion).

GREEN_THRUST = 0.5                                                     # (N)                                     


# Propulsion scenarios

# Both use the same continuous-retrograde-burn model (see state_derivative below);
# burn duration for each is derived from its propellant mass and mass flow rate, 
# so each scenario burns for exactly as long as its available propellant allows.

def mass_flow_rate(thrust, isp):
    """
    Propellant mass flow rate for a given thrust and specific impulse,
    from the thrust-mass-flow relationship, thrust = exhaust_velocity * (dm/dt).
    
    Parameters
    ----------
    thrust : float
            Thrust force, Newtons.
    isp : float
            Specific impulse, seconds.
            
    Returns
    -------
    float
            Mass flow rate, kg/s.
    """

    exhaust_velocity = isp * STANDARD_GRAVITY                          # (m/s)

    return thrust / exhaust_velocity


# Propulsion scenarios

SCENARIOS = {"electric" : {"isp" : ISP_ELECTRIC_S, "thrust" : ELECTRIC_THRUST, "propellant_mass" : ELECTRIC_PROPELLANT_MASS},
             "green" : {"isp" : ISP_GREEN_S, "thrust" : GREEN_THRUST, "propellant_mass" : GREEN_PROPELLANT_MASS}}

for _, _params in SCENARIOS.items():
    flow_rate = mass_flow_rate(_params["thrust"], _params["isp"])
    _params["burn_duration_s"] = _params["propellant_mass"] / flow_rate


# Initial state: PQW -> ECI, with propellant mass as the 7th state variable

def compute_initial_state(propellant_mass):
    """
    Compute the CubeSat's initial position, velocity, and propellant mass.
    
    Position and velocity are computed identically to 
    drag_deorbit_simulation.compute_initial_state(); this version adds
    initial propellant mass as a seventh state variable.
    
    Parameters
    ----------
    propellant_mass : float
            Initial propellant mass available for this scenario, kg.

    Returns
    -------
    numpy.ndarray
            Seven-element initial state [x, y, z, vx, vy, vz, m_propellant]
            in ECI (km and km/s) plus propellant mass (kg).
    """

    semi_latus_rectum = ORBITAL_RADIUS * (1 - ECCENTRICITY ** 2)

    r_pqw = np.array([ORBITAL_RADIUS * np.cos(np.radians(TRUE_ANOMALY)),
                      ORBITAL_RADIUS * np.sin(np.radians(TRUE_ANOMALY)),
                      0.0])
    v_pqw = np.array([-np.sqrt(MU_EARTH / semi_latus_rectum) * np.sin(np.radians(TRUE_ANOMALY)),
                      np.sqrt(MU_EARTH / semi_latus_rectum) * (ECCENTRICITY + np.cos(np.radians(TRUE_ANOMALY))),
                              0.0])

    pqw_to_eci_matrix = (rotation_matrix_z(RAAN) @ rotation_matrix_x(INCLINATION) @ rotation_matrix_z(ARGUMENT_OF_PERIGEE))

    r_eci = pqw_to_eci_matrix @ r_pqw
    v_eci = pqw_to_eci_matrix @ v_pqw

    return np.concatenate([r_eci, v_eci, [propellant_mass]])


# State-derivative function (gravity + drag + thrust) for solve_ivp

def make_state_derivative(thrust, burn_duration_s, exhaust_velocity):
    """
    Build a state-derivative function for a specific propulsion scenario.
    
    Returns a closure over thrust, burn_duration_s, and exhaust_velocity,
    so the same underlying physics (gravity, drag, retrograde thrust,
    propellant consumption) can be reused across scenarios without
    duplicating code - only the scenario's parameters differ between calls.
    
    Parameters
    ----------
    thrust : float
            Thrust force, Newtons.
    burn_duration_s : float
            Duration the engine fires, from t = 0 seconds.
    exhaust_velocity : float
            Effective exhaust velocity (Isp * standard gravity), m/s.
            
    Returns
    -------
    callable
            Function (t, state) -> rate of change, suitable for solve_ivp.
    """
    def state_derivative(t, state):
        """
        Compute the rate of change of the CubeSat's state under gravity,
        atmospheric drag, and retrograde thrust.
    
        Adds a thrust term to the gravity + drag model in
        drag_deorbit_simulation.state_derivative(). Thrust fires retrograde
        (opposite the current ECI velocity vector) for the first
        BURN_DURATION_S seconds of the simulation, after which the engine is off
        and the CubeSat coasts under gravity and drag alone. Propellant mass is
        tracked as the seventh state component and is consumed according to
        the thrust mass-flow relationship (thrust = exhaust_velocity *
        mass_flow_rate); current total mass (dry + remaining propellant) is
        used for both the thrust and drag acceleration calculations, since
        both scale as force / mass - this is the coupling between the two
        effects that this extension specifically investigates.
    
        J2 perturbation is not included (see drag_deorbit_simulation.py).
    
        Parameters
        ----------
        t : float
                Elapsed simulation time, seconds.
        state : array_like
                Seven-element state [x, y, z, vx, vy, vz, m_propellant] in ECI
                (km and km/s) plus remaining propellant mass (kg).
            
        Returns
        -------
        list of float
                Rate of change [vx, vy, vz, ax, ay, az, dm/dt]
                (km/s, km/s^2, and kg/s).
        """

        x, y, z, vx, vy, vz, m_propellant = state

        current_mass = DRY_MASS_KG + max(m_propellant, 0.0)            # Current total spacecraft mass (kg)

        # Gravity (spherical two_body model)

        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        g = -MU_EARTH / r ** 2
        g_x = g * x / r
        g_y = g * y / r
        g_z = g * z / r

        # ECI -> ECEF (needed only for the density lookup below)

        earth_rotation_angle = (360.0 / SIDEREAL_DAY) * t 
        r_ecef = rotation_matrix_z(earth_rotation_angle) @ np.array([x, y, z])
        x_ecef, y_ecef, z_ecef = r_ecef                                # ECEF position components (km)

        longitude_deg, latitude_deg, altitude_km = ecef_to_geodetic(x_ecef, y_ecef, z_ecef)

        # Atmospheric density (NRLMSISE-00 via pymsis)

        current_date = REFERENCE_START_DATE + timedelta(seconds=t)
        pymsis_output = calculate(dates = current_date, lons = longitude_deg, lats = latitude_deg, alts = altitude_km)
        density = np.asarray(pymsis_output[..., 0]).squeeze().item()   # index 0 = MASS_DENSITY (kg/m^3)

        # Drag, relative to the co-rotating atmosphere

        earth_rotation_rate = 2 * np.pi / SIDEREAL_DAY                 # (radians/s)
        v_atm_x = -earth_rotation_rate * y                             # Velocity of co-rotating atmosphere (km/s)
        v_atm_y = earth_rotation_rate * x
        v_atm_z = 0.0

        relative_velocity_x = vx - v_atm_x
        relative_velocity_y = vy - v_atm_y
        relative_velocity_z = vz - v_atm_z
        relative_velocity = np.sqrt(relative_velocity_x ** 2 + relative_velocity_y ** 2 + relative_velocity_z ** 2)
        relative_velocity_m_s = relative_velocity * 1000.0             # km/s -> m/s, to match density in kg/m^3

        drag_force = -0.5 * density * DRAG_COEFFICIENT * CROSS_SECTIONAL_AREA * relative_velocity_m_s ** 2
        drag_acceleration = drag_force / current_mass
        drag_acceleration_km_s2 = drag_acceleration / 1000.0

        drag_x = drag_acceleration_km_s2 * (relative_velocity_x / relative_velocity)
        drag_y = drag_acceleration_km_s2 * (relative_velocity_y / relative_velocity)
        drag_z = drag_acceleration_km_s2 * (relative_velocity_z / relative_velocity)

        # Retrograde thrust: active only during the burn window while propellant remains

        eci_velocity = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)            # CubeSat's velocity in the ECI frame (km/s)

        engine_is_firing = (t <= burn_duration_s) and (m_propellant > 0.0)

        if engine_is_firing:
            thrust_acceleration_km_s2 = -(thrust / current_mass) / 1000.0

            thrust_x = thrust_acceleration_km_s2 * (vx / eci_velocity) # Retrograde acceleration: opposite to CubeSat's velocity direction
            thrust_y = thrust_acceleration_km_s2 * (vy / eci_velocity)
            thrust_z = thrust_acceleration_km_s2 * (vz / eci_velocity)

            mass_flow_rate = thrust / exhaust_velocity                 # Rocket thrust relation: thrust = v_e * (dm/dt) (kg/s)
            dm_dt = -mass_flow_rate
        else:
            thrust_x = thrust_y = thrust_z = 0.0
            dm_dt = 0.0                                                # Rate of change of remaining propellant mass (kg)

        # Combined acceleration

        total_acceleration_x = g_x + drag_x + thrust_x
        total_acceleration_y = g_y + drag_y + thrust_y
        total_acceleration_z = g_z + drag_z + thrust_z

        return [vx, vy, vz, total_acceleration_x, total_acceleration_y, total_acceleration_z, dm_dt]

    return state_derivative


# Run a single scenario

def run_scenario(scenario_name):
    """
    Propagate the CubeSat's orbit under gravity, drag, and a fixed duration
    retrograde thrust burn for one named propulsion scenario, until either
    the configured duration elapses or reentry occurs.
    
    Parameters
    ----------
    scenario_name : str
            Key into SCENARIOS ("electric" or "green").

    Returns
    -------
    scipy.integrate.Oderesult
            Solution object with .t (time points) and .y (state history,
            including remaining propellant mass as the seventh row).
    """

    params = SCENARIOS[scenario_name]
    exhaust_velocity = params["isp"] * STANDARD_GRAVITY

    y0 = compute_initial_state(params["propellant_mass"])
    t_span = (0, SIMULATION_DURATION_DAYS * 24 * 3600)
    t_eval = np.linspace(*t_span, 1000)

    derivative_fn = make_state_derivative(thrust = params["thrust"], burn_duration_s = params["burn_duration_s"], exhaust_velocity = exhaust_velocity)

    ode_result = solve_ivp(fun = derivative_fn, t_span = t_span, y0 = y0, t_eval = t_eval, events = reentry_event, rtol = 1e-8, atol = 1e-8)

    return ode_result


def save_scenario_results(scenario_name, solution):
    """
    Extract altitude and propellant history from a solution and save both,
    along with the time array, using scenario-prefixed filenames. Trims to
    only the portion solve_ivp actually integrated, in case a reentry
    event ended the run early.
    
    Parameters
    ----------
    scenario_name : str
            Key into SCENARIOS ("electric" or "green"); used as a filename prefix.
    solution : scipy.integrate.Oderesult
            Result from run_scenario().
            
    Returns
    -------
    numpy.ndarray
            Altitude history, km (returned for immediate use).
    """

    n_valid = solution.y.shape[1]
    x_vals, y_vals, z_vals = solution.y[0][:n_valid], solution.y[1][:n_valid], solution.y[2][:n_valid]
    propellant_mass_vals = solution.y[6][:n_valid]

    r_vals = np.sqrt(x_vals ** 2 + y_vals ** 2 + z_vals ** 2)
    altitude_vals = r_vals - EARTH_EQUATORIAL_RADIUS

    np.save(DATA_DIR / f"thrust_drag_{scenario_name}_time.npy", solution.t)
    np.save(DATA_DIR / f"thrust_drag_{scenario_name}_altitude.npy", altitude_vals)
    np.save(DATA_DIR / f"thrust_drag_{scenario_name}_propellant.npy", propellant_mass_vals)

    return altitude_vals


if __name__ == "__main__":
    for scenario_name in SCENARIOS:
        params = SCENARIOS[scenario_name]
        solution = run_scenario(scenario_name)
        altitude_vals = save_scenario_results(scenario_name, solution)
        propellant_mass_vals = solution.y[6]
        propellant_remaining = max(0.0, propellant_mass_vals[-1])

        print(f"--- {scenario_name} ---")
        print(f"Isp = {params["isp"]} s, thrust = {params["thrust"]} N, "
              f"burn_duration = {params["burn_duration_s"] / 3600:.2f} hours "
              f"(derived from {params["propellant_mass"]} kg propellant budget, per report Table 3)")

        if solution.t_events[0].size > 0:
            reentry_time_days = solution.t_events[0][0] / 86400.0
            print(f"Reentry (altitude = {REENTRY_ALTITUDE_KM} km) reached at t = {reentry_time_days:.2f} days")
        else:
            print(f"No reentry within {SIMULATION_DURATION_DAYS} days; simulation ran to completion")

        print(f"Start altitude: {altitude_vals[0]:.3f} km")
        print(f"End altitude: {altitude_vals[-1]:.3f} km")
        print(f"Propellant remaining: {propellant_remaining:.4f} kg "
              f"(started with {params["propellant_mass"]} kg)")
        print()