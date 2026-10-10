"""
Sensitivity of drag-driven orbital decay to solar activity.

Reuse the gravity + drag model from drag_deorbit_simulation.py, run at
four real historical dates spanning Solar Cycle 25's solar minimum,
rising phase, near-peak, and solar maximum, to examine how atmospheric 
density (and therefore decay rate) varies with solar activity - the 
comparison requested in the third extension suggested by Dr. Dario 
Modenini (University of Bologna): extending the sensitivity analysis to 
orbital altitude decay under varying solar activity conditions, which
affects atmospheric density.

Dates and approximate F10.7 (10.7 cm solar radio flux, SFU):
- Solar minimum: December 2019 (smoothed sunspot number = 1.8; F10.7 ~ 68 SFU)
  Source: NOAA / NASA Solar Cycle 25 Prediction Panel,
  https://www.swpc.noaa.gov/products/solar-cycle-progression
- Rising phase: April 2022 (F10.7 ~ 116 SFU, 13-month smoothed)
  Source: NASA / NOAA Solar Cycle 25 progression estimates
- Near-peak: December 2022 (F10.7 ~ 142 SFU, 13-month smoothed)
  Source: NASA / NOAA Solar Cycle 25 progression estimates
- Solar maximum: August 2024 (smoothed sunspot number = 156.8; F10.7 ~ 245.6 SFU)
  Source: https://www.swpc.noaa.gov/products/solar-cycle-progression

Only REFERENCE_START_DATE differs between the two runs; all other
orbital / drag parameters are identical, so any difference in decay is
attributable to solar activity alone.

Important clarification on what each solar activity scenario refer
to: each simulation's reference date advances continuously throughout
its own 365-day propagation (atmospheric density is queried from pymsis
at the satellite's actual simulated date, not frozen at t = 0), so each
run experiences the real, evolving solar activity for a full year
starting from its reference date - not a constant level held fixed
throughout. The scenario labels describe each run's starting condition,
not a condition sustained unchanged for the full simulated year.

This does not undermine the comparison, since each starting regime
stays broadly within its own activity level for the full propagated 
year: the solar minimum run (December 2019 start) remains in the early,
still-low part of the rising phase a year later (December 2020, smoothed
sunspot number only ~ 8 - 10, versus 1.8 at the start); the solar maximum
run (August 2024 start) remains in a broadly elevated, only slowly-declining
phase a year later (August 2025, smoothed sunspot number still ~ 108 - 115, 
versus 156.8 at the start). Both runs therefore stay solidly within the 
activity regime their label describes throughout the comparison window.
"""


import numpy as np
from datetime import datetime, timedelta
from pymsis import calculate
from scipy.integrate import solve_ivp

from paths import DATA_DIR
from drag_deorbit_simulation import (
    MU_EARTH,
    EARTH_EQUATORIAL_RADIUS,
    SIDEREAL_DAY,
    CROSS_SECTIONAL_AREA,
    DRAG_COEFFICIENT,
    DRY_MASS_KG,
    REENTRY_ALTITUDE_KM,
    rotation_matrix_z,
    ecef_to_geodetic,
    compute_initial_state,
    reentry_event
)

SIMULATION_DURATION_DAYS = 365

SOLAR_SCENARIOS = {"solar_min" : datetime(2019, 12, 15),
                   "rising_phase" : datetime(2022, 4, 15),
                   "near_peak" : datetime(2022, 12, 15), 
                   "solar_max" : datetime(2024, 8, 15)}                    


def make_state_derivative(reference_date):
    """
    Build a state-derivative function (gravity + drag only) for a specific
    reference date, so atmospheric density is looked up relative to that
    date's real historical solar activity.
    
    Identical physics to drag_deorbit_simulation.state_derivative; only
    the reference date used for the pymsis density lookup differs between
    calls, via closure.
    
    Parameters
    ----------
    reference_date : datetime.datetime
            Simulation epoch (t = 0); solar activity is looked up relative to
            this date plus elapsed simulation time.
            
    Returns
    -------
    callable
            Function(t, state) -> rate of change, suitable for solve_ivp.
    """

    def state_derivative(t, state):
        x, y, z, vx, vy, vz = state

        # Gravity (spherical two-body model)

        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        g = - MU_EARTH / r ** 2
        g_x = g * x / r
        g_y = g * y / r
        g_z = g * z / r

        # ECI -> ECEF (needed only for the density lookup below)

        earth_rotation_angle = (360.0 / SIDEREAL_DAY) * t                               # (degrees)
        r_ecef = rotation_matrix_z(earth_rotation_angle) @ np.array([x, y, z])
        x_ecef, y_ecef, z_ecef = r_ecef

        longitude, latitude, altitude = ecef_to_geodetic(x_ecef, y_ecef, z_ecef)        # longitude (degrees), latitude (degrees), altitude (km)

        # Atmospheric density (NRLMSISE-00 via pymsis), at this scenario's date

        current_date = reference_date + timedelta(seconds = t)
        pymsis_output = calculate(dates = current_date, lons = longitude, lats = latitude, alts = altitude)
        density = np.asarray(pymsis_output[...,0]).squeeze().item()                     # index 0 = MASS_DENSITY (kg / m^3)

        # Drag, relative to the co-rotating atmosphere

        earth_rotation_rate = 2 * np.pi / SIDEREAL_DAY
        v_atm_x = -earth_rotation_rate * y
        v_atm_y = earth_rotation_rate * x
        v_atm_z = 0.0

        relative_velocity_x = vx - v_atm_x
        relative_velocity_y = vy - v_atm_y
        relative_velocity_z = vz - v_atm_z
        relative_velocity = np.sqrt(relative_velocity_x ** 2 + relative_velocity_y ** 2 + relative_velocity_z ** 2)
        relative_velocity_m_s = relative_velocity * 1000.0                              # km/s -> m/s, to match density in kg/m^3

        drag_force = -0.5 * density * DRAG_COEFFICIENT * CROSS_SECTIONAL_AREA * relative_velocity_m_s ** 2
        drag_acceleration = drag_force / DRY_MASS_KG
        drag_acceleration_km_s2 = drag_acceleration / 1000.0

        drag_acceleration_x = drag_acceleration_km_s2 * (relative_velocity_x / relative_velocity)
        drag_acceleration_y = drag_acceleration_km_s2 * (relative_velocity_y / relative_velocity)
        drag_acceleration_z = drag_acceleration_km_s2 * (relative_velocity_z / relative_velocity)

        total_acceleration_x = g_x + drag_acceleration_x 
        total_acceleration_y = g_y + drag_acceleration_y 
        total_acceleration_z = g_z + drag_acceleration_z 

        return [vx, vy, vz, total_acceleration_x, total_acceleration_y, total_acceleration_z]

    return state_derivative


def run_scenario(reference_date):
    """
    Propagate the CubeSat's orbit under gravity and drag at a specific
    reference date, until either SIMULATION_DURATION_DAYS elapses or
    reentry occurs.
    
    Parameters
    ----------
    reference_date : datetime.datetime
            Simulation epoch for this scenario's density lookups.
            
    Returns
    -------
    scipy.integrate.OdeResult
    """

    y0 = compute_initial_state()
    t_span = (0, SIMULATION_DURATION_DAYS * 24 * 3600)
    t_eval = np.linspace(*t_span, 1000)

    derivative_fn = make_state_derivative(reference_date)

    return solve_ivp(fun = derivative_fn, t_span = t_span, y0 = y0, t_eval = t_eval, events = reentry_event, rtol = 1e-10, atol = 1e-12)


def altitude_from_solution(solution):
    """
    Extract altitude history (km) from a solve_ivp solution, trimmed to what was integrated.
    """


    n_valid = solution.y.shape[1]
    x_vals, y_vals, z_vals = solution.y[0][:n_valid], solution.y[1][:n_valid], solution.y[2][:n_valid]
    r_vals = np.sqrt(x_vals ** 2 + y_vals ** 2 + z_vals ** 2)

    return r_vals - EARTH_EQUATORIAL_RADIUS


if __name__ == "__main__":
    decay_summary = []

    for scenario_name, reference_date in SOLAR_SCENARIOS.items():
        print(f"Running {scenario_name} ({reference_date.date()})...")
        solution = run_scenario(reference_date)
        altitude_vals = altitude_from_solution(solution)
        n_valid = len(altitude_vals)

        decay_km = altitude_vals[0] - altitude_vals[-1]
        elapsed_days = solution.t[n_valid - 1] / 86400.0
        reentry_detected = solution.t_events[0].size > 0

        status = f"reentry at {elapsed_days:.2f} days" if reentry_detected else f"no reentry within {SIMULATION_DURATION_DAYS} days"
        print(f"{status} | Start: {altitude_vals[0]:.3f} km | End: {altitude_vals[-1]:.3f} km | "
              f"Decay: {decay_km:.3f} km \n")

        decay_summary.append((scenario_name, reference_date, elapsed_days, decay_km, reentry_detected))

        np.save(DATA_DIR / f"{scenario_name}_time.npy", solution.t[:n_valid])
        np.save(DATA_DIR / f"{scenario_name}_altitude.npy", altitude_vals)

    np.save(DATA_DIR / "solar_scenario_name.npy", np.array([s[0] for s in decay_summary]))
    np.save(DATA_DIR / "solar_scenario_decay_km.npy", np.array([s[3] for s in decay_summary]))
    np.save(DATA_DIR / "solar_scenario_elapsed_days.npy", np.array([s[2] for s in decay_summary]))
    np.save(DATA_DIR / "solar_scenario_reentry_detected.npy", np.array([s[4] for s in decay_summary]))

    