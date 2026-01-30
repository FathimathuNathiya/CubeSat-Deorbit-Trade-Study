"""
Mission-level analysis for CubeSat deorbit feasibility.

This script integrates orbital mechanics, deorbit delta-V,
and propulsion modeling to estimate the propellant mass
required for a controlled end-of-life deorbit maneuver.
"""

from deorbit_delta_v import compute_deorbit_delta_v
from propellant_mass import propellant_mass

# Mission parameters

INITIAL_ALTITUDE_KM = 500.0
TARGET_ALTITUDE_KM = 120.0

DRY_MASS_KG = 12.0

ISP_GREEN_S = 220.0
ISP_ELECTRIC_S = 1200.0


def mission_analysis():
    """
    Perform the mission-level deorbit analysis and compare propulsion options.

    Returns
    -------
    dict
        Dictionary containing delta-V and propellant mass results.
    """

    delta_v = compute_deorbit_delta_v(INITIAL_ALTITUDE_KM, TARGET_ALTITUDE_KM)

    propellant_green = propellant_mass(delta_v, DRY_MASS_KG, ISP_GREEN_S)

    propellant_electric = propellant_mass(delta_v, DRY_MASS_KG, ISP_ELECTRIC_S)

    return {
        "delta_v": delta_v,
        "propellant_green": propellant_green,
        "propellant_electric": propellant_electric
    }


def print_mission_results(results):
    """
    Print the mission analysis results.
    """
    print(f"Required deorbit delta-V: {results['delta_v']:.2f} m/s")
    print(
        f"Green propulsion propellant mass: {results['propellant_green']:.4f} kg")
    print(
        f"Electric propulsion propellant mass: {results['propellant_electric']:.4f} kg")


if __name__ == "__main__":
    results = mission_analysis()
    print_mission_results(results)
