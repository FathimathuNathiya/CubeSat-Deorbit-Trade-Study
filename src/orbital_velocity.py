"""
Orbital mechanics utilities for CubeSat mission analysis.

This module implements basic circular orbit models in Low Earth Orbit (LEO).
It is intended as a building block for deorbit analysis and propulsion trade studies 
in sustainable spacecraft systems research.
"""

import numpy as np

# Physical constants (SI units)

MU_EARTH = 3.986e14        # Earth's gravitational parameter (m^3/s^2)
R_EARTH = 6371e3           # Earth's mean radius (m)


def orbital_velocity(altitude_km):
    """
    Compute the circular orbital velocity (m/s) at a given altitude above Earth's surface.

    Parameters
    ----------
    altitude_km : float
        Altitude above Earth in kilometers

    Returns
    -------
    float
        Orbital velocity in meters per second
    """

    r = R_EARTH + altitude_km * 1e3
    return np.sqrt(MU_EARTH / r)


def print_orbital_velocity(altitude_km):
    """
    Prints the orbital velocity in km/s for a given altitude.
    """
    v = orbital_velocity(altitude_km)
    print(f"Orbital velocity at {altitude_km} km: {v / 1000:.2f} km/s")
