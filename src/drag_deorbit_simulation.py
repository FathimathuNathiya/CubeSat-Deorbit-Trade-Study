"""
Drag-driven orbital decay simulation for Cubesat mission.

This module extends the analytical trade study with a numerical propagator
that models the Cubesat's orbit decaying under atmospheric drag alone (no
propulsion). Gravity and drag accelerations are integrated forward in time
using scipy's solve_ivp, with atmospheric density obtained from the 
NRLMSISE-00 model (via pymsis) at the satellite's true geodetic position
at every timestep.

This addresses the first extension suggested by Dr. Dario Modenini 
(University of Bologna): incorporating aerodynamic drag into the deorbit
model using analytical drag formulas for a circular orbit.

Propagation frame: ECI (Earth-Centered Inertial) - required so that gravity
and drag acceleration can be computed without fictitious rotating-frame
forces. Position is converted to ECEF and then geodetic latitude / longitude /
altitude only where needed (atmospheric density lookup), not for propagation.

A reentry event stops integration once altitude drops to REENTRY_ALTITUDE_KM
(120 km, matching the report's target altitude), since the orbit has no
further physical meaning below that point - this prevents the integrator
from being run past reentry into a regime where the underlying two-body
plus drag model is no longer valid (which otherwise surfaces as runaway or
negative altitude over long simulated durations).
"""


import numpy as np
from datetime import datetime, timedelta
from pymsis import calculate
from scipy.integrate import solve_ivp

from mission_analysis import INITIAL_ALTITUDE_KM, DRY_MASS_KG
from paths import DATA_DIR


# Physical and orbital constants
                                                                
MU_EARTH = 398600.0                                                # Earth's gravitational parameter (km^3/s^2)
EARTH_EQUATORIAL_RADIUS = 6378.137                                 # WGS84 semi-major axis (km)
EARTH_POLAR_RADIUS = 6356.752                                      # WGS84 semi-minor axis (km)
SIDEREAL_DAY = 23 * 3600 + 56 * 60 + 4                             # Earth's rotation period relative to stars (s)

CROSS_SECTIONAL_AREA = 0.02                                        # 6 U Cubesat, smallest face (10 x 20 cm), nose-first assumption (m^2)
DRAG_COEFFICIENT = 2.2                                             # Standard value for Cubesat-class objects in free molecular flow

# Orbital elements (circular, near-polar sun-synchronous orbit)

ORBITAL_RADIUS = INITIAL_ALTITUDE_KM + EARTH_EQUATORIAL_RADIUS     # (km)
ECCENTRICITY = 0.0                                                 # Circular orbit; standard assumption for LEO Cubesats
INCLINATION = 98.0                                                 # Sun-synchronous, consistent with real Cubesat missions (degrees)
RAAN = 0.0                                                         # Arbitrary reference orientation; not tied to a real launch date (degrees)
ARGUMENT_OF_PERIGEE = 0.0                                          # Undefined for circular orbit; set to 0 by convention (degrees)
TRUE_ANOMALY = 0.0                                                 # Arbitrary starting point along the (circular) orbit (degrees)

REFERENCE_START_DATE = datetime(2025, 8, 4)                        # Arbitrary baseline epoch for atmospheric density lookup

SIMULATION_DURATION_DAYS = 365

REENTRY_ALTITUDE_KM = 120.0                                        # Integration stops once altitude reaches this value


# Rotation matrices (PQW -> ECI -> ECEF chain)

def rotation_matrix_x(theta_deg):
    """
    Build a 3 x 3 rotation matrix about the x-axis.

    Parameters
    ----------
    theta_deg : float
            Rotation angle in degrees.

    Returns
    -------
    numpy.ndarray
            3 x 3 rotation matrix.
    """

    theta = np.radians(theta_deg)                                # (radians)
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]
                     ])


def rotation_matrix_z(theta_deg):
    """
    Build a 3 x 3 rotation matrix about the z-axis.

    Parameters
    ----------
    theta_deg : float
            Rotation angle in degrees.

    Returns
    -------
    numpy.ndarray
            3 x 3 rotation matrix.
    """

    theta = np.radians(theta_deg)                               # (radians)
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])


# Initial state: PQW -> ECI

def compute_initial_state():
    """
    Compute the Cubesat's initial position and velocity in the ECI frame.

    Position and velocity are first computed in the perifocal (PQW) frame
    using the classical orvital elements, then rotated into ECI via the 
    standard 3-1-3 Euler sequence (RAAN, inclination, argument of perigee).

    Returns
    -------
    numpy.ndarray
            Six-element initial state [x, y, z, vx, vy, vz] in ECI (km and km/s).
    """

    semi_latus_rectum = ORBITAL_RADIUS * (1 - ECCENTRICITY ** 2) # (km)

    r_p = ORBITAL_RADIUS * np.cos(np.radians(TRUE_ANOMALY))
    r_q = ORBITAL_RADIUS * np.sin(np.radians(TRUE_ANOMALY))
    r_w = 0.0
    r_pqw = np.array([r_p, r_q, r_w])                            # (km)

    v_p = -np.sqrt(MU_EARTH / semi_latus_rectum) * np.sin(np.radians(TRUE_ANOMALY))
    v_q = np.sqrt(MU_EARTH / semi_latus_rectum) * (ECCENTRICITY + np.cos(np.radians(TRUE_ANOMALY)))
    v_w = 0.0
    v_pqw = np.array([v_p, v_q, v_w])                            # (km/s)

    rotation = rotation_matrix_z(RAAN) @ rotation_matrix_x(INCLINATION) @ rotation_matrix_z(ARGUMENT_OF_PERIGEE)

    r_eci = rotation @ r_pqw
    v_eci = rotation @ v_pqw

    return np.concatenate([r_eci, v_eci])


# Geodetic conversion (ECEF -> latitude, longitude, altitude)

ellipsoid_eccentricity = np.sqrt(1 - (EARTH_POLAR_RADIUS ** 2 / EARTH_EQUATORIAL_RADIUS ** 2))

def ecef_to_geodetic(x_ecef, y_ecef, z_ecef):
    """
    Convert ECEF Cartesian coordinates to geodetic latitude, longitude,
    and altitude using the WGS84 ellipsoid (Bowring's iterative method).
    
    Parameters
    ----------
    x_ecef, y_ecef, z_ecef : float
            ECEF position components, km.

    Returns
    -------
    tuple of float
            (longitude_deg, latitude_deg, altitude_km)       
    """

    radial_distance_xy = np.sqrt(x_ecef ** 2 + y_ecef ** 2)         # Distance from rotation (z) axis (km)
    longitude = np.arctan2(y_ecef, x_ecef)                          # (radians) Converted to degrees at return

    # Initial latitude guess (spherical approximation)

    latitude = np.arctan2(z_ecef, radial_distance_xy * 
                          (1 - ellipsoid_eccentricity ** 2))        # (radians) Converted to degrees at return

    # Iterative refinement (Bowrin's method)

    for _ in range(5):
            radius_of_curvature = EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - ellipsoid_eccentricity ** 2 * np.sin(latitude) ** 2)
            latitude = np.arctan2(z_ecef + ellipsoid_eccentricity ** 2 * radius_of_curvature * np.sin(latitude), radial_distance_xy)
    
    altitude = (radial_distance_xy / 
                np.cos(latitude)) - radius_of_curvature            # (km)

    return np.degrees(longitude), np.degrees(latitude), altitude


# Reentry event, for solve_ivp

def reentry_event(t, state):
     """
     Zero-crossing function for solve_ivp: reaches zero when altitude
     drops to REENTRY_ALTITUDE_KM.
     
     Parameters
     ----------
     t : float
            Elapsed simulation time, seconds (unused, required by solve_ivp's
            event function signature).
    state : array_like
            Current state; only the first three components (x, y, z, in ECI, km)
            are used.
    
    Returns
    -------
    float
            Current altitude minus REENTRY_ALTITUDE_KM; solve_ivp stops
            integration when this crosses zero (see .terminal / .direction
            attributes set below).
    """

     x, y, z = state[0], state[1], state[2]
     r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
     altitude = r - EARTH_EQUATORIAL_RADIUS
     return altitude - REENTRY_ALTITUDE_KM


reentry_event.terminal = True                                      # stop integrating when this event fires
reentry_event.direction = -1                                       # only trigger on a decreasing crossing (altitude falling, not rising)


# State-derivative function (gravity + drag), for solve_ivp

def state_derivative(t, state):
    """
    Compute the rate of change of the Cubesat's state under gravity and
    atmospheric drag.
     
    Called repeatedly by solve_ivp. All position / velocity components are in
    the ECI frame; ECEF and geodetic coordinates are computed internally,
    solely to obtain atmospheric density from NRLMSISE-00 at the satellite's
    current true position.

    J2 perturbation is not included; its primary effect is on orbital plane
    orientation (RAAN, argument of perigee drift) rather than altitude
    decay, which is the focus of this study.

    Parameters
    ----------
    t : float
            Elapsed simulation time, seconds.
    state : array_like
            Six-element state [x, y, z, vx, vy, vz] in ECI (km and km/s).

    Returns
    -------
    list of float
            Rate of change [vx, vy, vz, ax, ay, az] (km/s and km/s^2).

    """

    x, y, z, vx, vy, vz = state

    # Gravitational acceleration using the spherical two-body model (J2 perturbation not included)

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)                          # (km)
    g = -MU_EARTH / r ** 2                                         # Gravitational acceleration (km/s^2)
    g_x = g * x / r 
    g_y = g * y / r
    g_z = g * z / r

    # ECI -> ECEF (needed only for the density lookup below)
                                        
    earth_rotation_rate = 2 * np.pi / SIDEREAL_DAY                 # (radians/s)
    earth_rotation_angle = 360.0 / SIDEREAL_DAY * t                # (degrees)
    r_ecef = rotation_matrix_z(earth_rotation_angle) @ np.array([x, y, z])
    x_ecef, y_ecef, z_ecef = r_ecef                                # ECEF position components (km)

    longitude_deg, latitude_deg, altitude_km = ecef_to_geodetic(x_ecef, y_ecef, z_ecef)

    # Atmospheric density (NRLMSISE-00 via pymsis)

    current_date = REFERENCE_START_DATE + timedelta(seconds=t)
    pymsis_output = calculate(dates=current_date, lons=longitude_deg, lats=latitude_deg, alts=altitude_km)

    density = np.asarray(pymsis_output[..., 0]).squeeze().item()          # index 0 = MASS_DENSITY (kg/m^3)

    # Drag, relative to the co-rotating atmosphere

    v_atm_x = -earth_rotation_rate * y
    v_atm_y = earth_rotation_rate * x
    v_atm_z = 0.0

    relative_velocity_x = vx - v_atm_x
    relative_velocity_y = vy - v_atm_y
    relative_velocity_z = vz - v_atm_z
    relative_velocity = np.sqrt(relative_velocity_x ** 2 + relative_velocity_y ** 2 + relative_velocity_z ** 2)
    relative_velocity_m_s = relative_velocity * 1000.0             # km/s -> m/s, to match density in kg/m^3

    drag_force = -0.5 * density * DRAG_COEFFICIENT * CROSS_SECTIONAL_AREA * relative_velocity_m_s ** 2
    drag_acceleration = drag_force / DRY_MASS_KG                   # (m/s^2)
    drag_acceleration_km_s2 = drag_acceleration / 1000.0           # (km/s^2)

    drag_acceleration_x = drag_acceleration_km_s2 * (relative_velocity_x / relative_velocity)
    drag_acceleration_y = drag_acceleration_km_s2 * (relative_velocity_y / relative_velocity)
    drag_acceleration_z = drag_acceleration_km_s2 * (relative_velocity_z / relative_velocity)

    # Combined acceleration

    total_acceleration_x = g_x + drag_acceleration_x
    total_acceleration_y = g_y + drag_acceleration_y
    total_acceleration_z = g_z + drag_acceleration_z

    return [vx, vy, vz, total_acceleration_x, total_acceleration_y, total_acceleration_z]


# Run simulation

def run_simulation():
    """
    Propagate the Cubesat's orbit under gravity and drag over the
    configured simulation duration.
    
    Returns
    -------
    scipy.integrate.OdeResult
            Solution object with .t (time points) and .y (state history), and
            .t_events / .status indicating whether reentry occurred before
            the configured duration elapsed.
    """

    y0 = compute_initial_state()
    t_span = (0, SIMULATION_DURATION_DAYS * 24 * 3600)
    t_eval = np.linspace(*t_span, 1000)

    ode_result = solve_ivp(fun=state_derivative, t_span=t_span, y0=y0, t_eval=t_eval, events = reentry_event, rtol=1e-8, atol=1e-8)

    return ode_result

if __name__ == "__main__":
    solution = run_simulation()

    # t_eval points beyond a terminal event are not populated; trim to what solve_ivp actually integrated

    n_valid = solution.y.shape[1]
    x_vals, y_vals, z_vals = solution.y[0][:n_valid], solution.y[1][:n_valid], solution.y[2][:n_valid]
    r_vals = np.sqrt(x_vals ** 2 + y_vals ** 2 + z_vals ** 2)
    altitude_vals = r_vals - EARTH_EQUATORIAL_RADIUS

    if solution.t_events[0].size > 0:
        reentry_time_days = solution.t_events[0][0] / 86400.0
        print(f"Reentry altitude = {REENTRY_ALTITUDE_KM} km reached at t = {reentry_time_days:.2f} days")
    else:
         print(f"No reentry within {SIMULATION_DURATION_DAYS} days; simulation ran to completion")    

    print(f"Start altitude: {altitude_vals[0]:.3f} km")
    print(f"End altitude: {altitude_vals[-1]:.3f} km")

    np.save(DATA_DIR / "drag_decay_time.npy", solution.t)
    np.save(DATA_DIR / "drag_decay_altitude.npy", altitude_vals)

