"""
Deorbit maneuver analysis for CubeSat missions.

This module computes the delta-V required to transition a spacecraft
from a circular Low Earth Orbit to an elliptical orbit with a specified
perigee altitude. It builds directly on basic orbital mechanics models 
and is intended for mission-level sustainability studies.
"""

import numpy as np
from orbital_velocity import orbital_velocity, MU_EARTH, R_EARTH


def compute_deorbit_delta_v(initial_altitude_km, perigee_altitude_km):
    """
    Compute the delta_V required to lower the perigee of a circular orbit.

    Parameters
    ----------
    initial_altitude_km : float
        Initial circular orbit altitude above Earth (km).

    perigee_altitude_km : float
        Target perigee altitude after the maneuver (km).

    Returns
    -------
    float
        Required delta-V (m/s).
    """

    r_initial = R_EARTH + initial_altitude_km * 1e3
    r_perigee = R_EARTH + perigee_altitude_km * 1e3

    v_circular = orbital_velocity(initial_altitude_km)

    semi_major_axis = 0.5 * (r_initial + r_perigee)

    v_transfer_apogee = np.sqrt(
        MU_EARTH * (2.0 / r_initial - 1.0 / semi_major_axis))

    return v_circular - v_transfer_apogee


def print_deorbit_delta_v(initial_altitude_km, perigee_altitude_km):
    """
    Print the delta-V required for a deorbit maneuver.
    """
    delta_v = compute_deorbit_delta_v(initial_altitude_km, perigee_altitude_km)
    print(
        f"Delta-V required to lower perigee from {initial_altitude_km:.0f} km to {perigee_altitude_km:.0f} km: {delta_v:.2f} m/s")
