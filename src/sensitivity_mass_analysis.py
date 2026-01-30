"""
Sensitivity analysis of propellant mass with respect to CubeSat dry mass
for Green and Electric propulsion systems.
"""

import numpy as np
from mission_analysis import (
    INITIAL_ALTITUDE_KM, TARGET_ALTITUDE_KM, ISP_GREEN_S, ISP_ELECTRIC_S)
from deorbit_delta_v import compute_deorbit_delta_v
from propellant_mass import propellant_mass


def compute_mass_sensitivity():

    # Range of Cubesat dry masses to study (kg)

    satellite_masses = np.linspace(5, 25, 10)

    green_propellant = []
    electric_propellant = []

    # Compute required delta-V for the mission

    delta_v = compute_deorbit_delta_v(INITIAL_ALTITUDE_KM, TARGET_ALTITUDE_KM)

    # Perform sensitivity analysis over satellite mass

    for mass in satellite_masses:
        green_propellant.append(propellant_mass(delta_v, mass, ISP_GREEN_S))
        electric_propellant.append(
            propellant_mass(delta_v, mass, ISP_ELECTRIC_S))

    print("CubeSat Dry Mass (kg) | Green Propellant Mass (kg) | Electric Propellant Mass (kg)")
    print("----------------------------------------------------------------------------------")
    for m, g, e in zip(satellite_masses, green_propellant, electric_propellant):
        print(f"{m:.1f} \t\t\t {g:.4f} \t\t\t {e:.4f}")

    return satellite_masses, green_propellant, electric_propellant


if __name__ == '__main__':
    compute_mass_sensitivity()
