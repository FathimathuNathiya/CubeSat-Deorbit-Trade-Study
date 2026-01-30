"""
Power Constraint and Deorbit Burn Time Trade Study for CubeSat Deorbit Mission.

This script performs:
Power Constraint Analysis
    Computes thrust and mass flow rate for different available power levels.
Deorbit Time Trade Study
    Computes the required burn time for deorbit using power-limited thrust.
"""

from deorbit_delta_v import compute_deorbit_delta_v
from propellant_mass import propellant_mass, g0
from mission_analysis import (
    INITIAL_ALTITUDE_KM, TARGET_ALTITUDE_KM, DRY_MASS_KG, ISP_ELECTRIC_S)


# Power Constraint Analysis

def power_constraint_analysis(available_power_w):
    """
    Compute thrust and mass flow rate for a given available power level.

    Parameters
    ----------
    available_power_w : float
        Available onboard electrical power in Watts.

    Returns
    -------
    dict
        Dictionary containing thrust and mass flow rate.
    """

    delta_v = compute_deorbit_delta_v(INITIAL_ALTITUDE_KM, TARGET_ALTITUDE_KM)
    propellant_required = propellant_mass(delta_v, DRY_MASS_KG, ISP_ELECTRIC_S)
    exhaust_velocity = ISP_ELECTRIC_S * g0
    thrust = (2 * available_power_w) / exhaust_velocity
    mass_flow_rate = thrust / exhaust_velocity

    return {"power_w": available_power_w, "thrust": thrust, "mass_flow_rate": mass_flow_rate, "propellant_required": propellant_required}


# Deorbit Time Trade Study

def deorbit_burn_time(power_result):
    """
    Compute burn time required for deorbit using power-limited thrust.

    Parameters
    ----------
    power_result : dict
        Dictionary returned from power_constraint_analysis().

    Returns
    -------
    dict
        Dictionary containing burn time in hours.
    """

    propellant_required = power_result["propellant_required"]
    mass_flow_rate = power_result["mass_flow_rate"]
    burn_time_sec = propellant_required / mass_flow_rate
    burn_time_hr = burn_time_sec / 3600

    return {"power_w": power_result["power_w"], "burn_time_hr": burn_time_hr}


# Print Results

def print_results(power_results, time_results):
    """
    Print power constraint analysis and time trade study results in a clean format.
    """
    print("Power Constraint and Deorbit Time Trade Study")
    print("---------------------------------------------")
    for power_results, time_results in zip(power_results, time_results):
        print(f"Available Power: {power_results['power_w']} W")
        print(f"Thrust: {power_results['thrust']:.6f} N")
        print(f"Burn Time: {time_results['burn_time_hr']:.2f} hours \n")


if __name__ == "__main__":
    power_levels = [5, 10, 20, 30, 50]
    power_results = []
    time_results = []

    for power in power_levels:
        p_result = power_constraint_analysis(power)
        t_result = deorbit_burn_time(p_result)
        power_results.append(p_result)
        time_results.append(t_result)

    print_results(power_results, time_results)
