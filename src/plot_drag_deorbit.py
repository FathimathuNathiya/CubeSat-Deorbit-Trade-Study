"""
Plot drag_driven altitude decay for the Cubesat mission.

Reads the time and alltitude history saved by drag_deorbit_simulation.py
and produces a plot of altitude versus elapsed time , showing the orbit 
decaying under atmospheric drag alone (no propulsion).

This is the figure corresponding to the first extension suggested by
Dr. Dario Modenini (University of Bologna): aerodynamic drag modeling
for a circular orbit.
"""


import numpy as np
import matplotlib.pyplot as plt

from paths import DATA_DIR, FIGURES_DIR


# Load data

time_s = np.load(DATA_DIR / "drag_decay_time.npy")
altitude_km = np.load(DATA_DIR / "drag_decay_altitude.npy")

time_days = time_s / 86400.0                                          # Seconds to days, for a readable x-axis


# Plot

plt.figure(figsize = (11, 6))
plt.plot(time_days, altitude_km)

plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")
plt.title("Cubesat Altitude Decay Under Atmospheric Drag \n (500 km initial altitude, no propulsion)")
plt.grid(True, alpha = 0.3)

start_altitude = altitude_km[0]
end_altitude = altitude_km[-1]
total_decay = start_altitude - end_altitude

plt.annotate(f"Start: {start_altitude:.3f} km",
             xy = (time_days[0], altitude_km[0]),
             xytext = (2, 3),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {end_altitude:.3f} km \n"
             f"(t: {time_days[-1]} days) \n"
             f"(decay: {total_decay:.3f} km)",
             xy = (time_days[-1], altitude_km[-1]),
             xytext = (-50, 20),
             textcoords = "offset points",
             fontsize = 9)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "drag_deorbit_decay.png", dpi = 300)
plt.show()

print(f"Total altitude decay over {time_days[-1]:.1f} days: {total_decay:.3f} km")
