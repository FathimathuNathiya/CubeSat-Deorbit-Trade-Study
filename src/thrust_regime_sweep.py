"""
Thrust regime sweep: drag's relative contribution across varying electric
propulsion thrust levels.

Reuses the coupled gravity + drag + thrust model from
thrust_drag_coupled_simulation.py, run at each of the five power / thrust
operating points from CubeSat_Deorbit_Analysis.pdf (Table 5), to examine
how drag's contribution to total altitude change varies with burn
duration and thrust level - the specific comparison requested in the
first bullet of Dr. Dario Modenini's second extension ("examining how
drag's relative contribution varies across different burn durations
and thrust regimes").

As in thrust_drag_coupled_simulation.py, the electric propellant budget
(0.1118 kg, from report's Table 3) is applied at each power level as a
continuous burn - see the framing note in that file for why this is not
a reproduction of the report's idealized 500 km -> 120 km transfer.

Each regime runs until either SIMULATION_DURATION_DAYS elapses or
reentry occurs (altitude reaches REENTRY_ALTITUDE_KM), whichever comes
first (see reentry_event in drag_deorbit_simulation.py). Since regimes
reach reentry at different times, "total decay" and "drag's share" are
each computed over that regime's own elapsed time, not a common window -
this is noted explicitly in the saved regime_elapsed_days array and in
the printed output, since it means regimes are not being compared over
identical durations.
"""


import numpy as np
from scipy.integrate import solve_ivp

from drag_deorbit_simulation import EARTH_EQUATORIAL_RADIUS, REENTRY_ALTITUDE_KM, reentry_event
from paths import DATA_DIR
from thrust_drag_coupled_simulation import (
    STANDARD_GRAVITY,
    SIMULATION_DURATION_DAYS,
    ELECTRIC_PROPELLANT_MASS,
    ISP_ELECTRIC_S,
    compute_initial_state,
    make_state_derivative,
    mass_flow_rate
)


# Thrust regimes: report Table 5 (Available Power (W), Thrust (N))

THRUST_REGIMES = [
    {"power" : 5, "thrust" : 0.000850},
    {"power" : 10, "thrust" : 0.001700},
    {"power" : 20, "thrust" : 0.003399},
    {"power" : 30, "thrust" : 0.005099},
    {"power" : 50, "thrust" : 0.008498}
]

for regime in THRUST_REGIMES:
    regime["isp"] = ISP_ELECTRIC_S
    regime["propellant_mass"] = ELECTRIC_PROPELLANT_MASS
    regime["burn_duration_s"] = ELECTRIC_PROPELLANT_MASS / mass_flow_rate(regime["thrust"], ISP_ELECTRIC_S)
    regime["burn_within_window"] = regime["burn_duration_s"] <= SIMULATION_DURATION_DAYS * 24 * 3600


# Run one thrust regime

def run_regime(regime):
    """
    Propagate the CubeSat's orbit for one thrust / power operating point,
    until either SIMULATION_DURATION_DAYS elapses or reentry occurs.
    
    Parameters
    ----------
    regime : dict
            One entry from THRUST_REGIMES, with thrust, isp, propellant_mass,
            and burn_duration already populated.
            
    Returns
    -------
    scipy.integrate.OdeResult
            Solution object with .t and .y (state history).
    """

    exhaust_velocity = regime["isp"] * STANDARD_GRAVITY
    y0 = compute_initial_state(regime["propellant_mass"])
    t_span = (0, SIMULATION_DURATION_DAYS * 24 * 3600)
    t_eval = np.linspace(*t_span, 1000)

    derivative_fn = make_state_derivative(thrust = regime["thrust"], burn_duration_s = regime["burn_duration_s"], exhaust_velocity = exhaust_velocity)

    return solve_ivp(fun = derivative_fn, t_span = t_span, y0 = y0, t_eval = t_eval, events = reentry_event, rtol = 1e-8, atol = 1e-8)


# Run drag-only baseline once, for comparison at each regime

def run_drag_only_baseline():
    """
    Propagate the same simulation window under gravity and drag alone
    (zero thrust), used as the reference for isolating each regime's 
    thrust-specific contribution.
    
    Returns
    -------
    scipy.integrate.OdeResult
    """

    exhaust_velocity = ISP_ELECTRIC_S * STANDARD_GRAVITY              # unused when thrust = 0, kept for signature
    y0 = compute_initial_state(propellant_mass = 0.0)                 # no propellant carried
    t_span = (0, SIMULATION_DURATION_DAYS * 24 * 3600)
    t_eval = np.linspace(*t_span, 1000)

    derivative_fn = make_state_derivative(thrust = 0.0, burn_duration_s = 0.0, exhaust_velocity = exhaust_velocity)

    return solve_ivp(fun = derivative_fn, t_span = t_span, y0 = y0, t_eval = t_eval, events = reentry_event, rtol = 1e-8, atol = 1e-8)


def altitude_from_solution(solution):
    """
    Extract altitude history (km) from a solve_ivp solution, trimmed to
    only the portion actually integrated (in case reentry ended the run
    early).
    """

    n_valid = solution.y.shape[1]
    x_vals = solution.y[0][:n_valid]
    y_vals = solution.y[1][:n_valid]
    z_vals = solution.y[2][:n_valid]
    r_vals = np.sqrt(x_vals ** 2 + y_vals ** 2 + z_vals ** 2)

    return r_vals - EARTH_EQUATORIAL_RADIUS


if __name__ == "__main__":


    print("Running drag-only baseline...")
    baseline_solution = run_drag_only_baseline()
    baseline_altitude = altitude_from_solution(baseline_solution)
    baseline_time_days = baseline_solution.t[:len(baseline_altitude)] / 86400.0
    baseline_decay = baseline_altitude[0] - baseline_altitude[-1]

    if baseline_solution.t_events[0].size > 0:
        print(f"Baseline reentry at t = {baseline_solution.t_events[0][0] / 86400.0:.2f} days")
        
    print(f"Drag-only decay over {baseline_time_days[-1]:.2f} days: {baseline_decay:.4f} km \n")

    power_levels = []                                                 # (W)
    thrust_levels = []                                                # (N)
    burn_duration_hr = []                                             # (hour)
    burn_within_window_flags = []
    elapsed_days = []
    reentry_reached_flags = []
    total_decay = []                                                  # (km)
    drag_share_percent = []

    for regime in THRUST_REGIMES:
        note = "" if regime["burn_within_window"] else "[burn_exceeds_simulation_window]"
        print(f"Running {regime["power"]} W ({regime["thrust"]} N, "
              f"burn = {regime["burn_duration_s"] / 3600:.2f} hr){note} ... ")
        solution = run_regime(regime)
        altitude_vals = altitude_from_solution(solution)
        regime_time_days = solution.t[:len(altitude_vals)] / 86400.0
        decay = altitude_vals[0] - altitude_vals[-1]
        reentry_reached = solution.t_events[0].size > 0

        # Compare against the drag-only baseline over this regime's own
        # elapsed time (via interpolation), since regime's each reentry
        # at different times and durations are not directly comparable.

        baseline_decay_at_same_time = np.interp(regime_time_days[-1], baseline_time_days, baseline_altitude[0] - baseline_altitude)
        drag_share = (baseline_decay_at_same_time / decay) * 100.0 if decay > 0 else float("nan")

        power_levels.append(regime["power"])
        thrust_levels.append(regime["thrust"])
        burn_duration_hr.append(regime["burn_duration_s"] / 3600.0)
        burn_within_window_flags.append(regime["burn_within_window"])
        elapsed_days.append(regime_time_days[-1])
        reentry_reached_flags.append(reentry_reached)
        total_decay.append(decay)
        drag_share_percent.append(drag_share)

        reentry_note = f"reentry at {regime_time_days[-1]:.2f} days" if reentry_reached else \
                       f"no reentry within {SIMULATION_DURATION_DAYS} days"
        print(f"Elapsed: {regime_time_days[-1]:.2f} days ({reentry_note}) | "
              f"Total decay: {decay:.4f} km | "
              f"Drag's share (over same elapsed time): {drag_share:.2f} % \n")

    np.save(DATA_DIR / "regime_power.npy", np.array(power_levels))
    np.save(DATA_DIR / "regime_thrust.npy", np.array(thrust_levels))
    np.save(DATA_DIR / "regime_burn_duration_hr.npy", np.array(burn_duration_hr))
    np.save(DATA_DIR / "regime_burn_within_window.npy", np.array(burn_within_window_flags))
    np.save(DATA_DIR / "regime_elapsed_days.npy", np.array(elapsed_days))
    np.save(DATA_DIR / "regime_reentry_reached.npy", np.array(reentry_reached_flags))
    np.save(DATA_DIR / "regime_total_decay.npy", np.array(total_decay))
    np.save(DATA_DIR / "regime_drag_share_percent.npy", np.array(drag_share_percent))

    