"""
Propellant mass calculations for CubeSat deorbit maneuvers.

This module computes the propellant mass required to achieve a specified
delta-V using the Tsiolkovsky rocket equation. It is intended for
mission-level propulsion analysis and sustainability studies in
Low Earth Orbit (LEO).
"""

import numpy as np

# Physical constants (SI units)

g0 = 9.80665         # Standard gravity (m/s^2)


def propellant_mass(delta_v, dry_mass, isp):
    """
    Compute the propellant mass required for a given delta-V.

    Parameters
    ----------
    delta_v : float
        Required velocity change in meters per second (m/s).

    dry_mass : float
        Spacecraft dry mass (excluding propellant) in kilograms (kg).

    isp : float
        Specific impulse of the propulsion system in seconds (s).

    Returns
    -------
    float
        Required propellant mass in kilograms (kg).
    """

    mass_ratio = np.exp(delta_v / (isp * g0))
    initial_mass = dry_mass * mass_ratio

    return initial_mass - dry_mass


def print_propellant_mass(delta_v, dry_mass, isp):
    """
    Print the propellant mass required for a deorbit maneuver.
    """
    propellant_mass_kg = propellant_mass(delta_v, dry_mass, isp)
    print(f"Required propellant mass: {propellant_mass_kg:.4f} kg")
