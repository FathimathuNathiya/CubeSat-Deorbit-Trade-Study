"""
Plot drag-driven altitude decay across four solar cycle phases.

Reads the time / altitude history saved by solar_activity_sensitivity.py
for all four scenarios (solar minimum, rising phase, near-peak, solar
maximum) and plots them together, quantifying how solar activity (via
its effect on atmospheric density) changes deorbit rate - the figure for
Dr. Dario Modenini's third extension.
"""


import numpy as np
import matplotlib.pyplot as plt

from paths import DATA_DIR, FIGURES_DIR


SCENARIO_ORDER = ["solar_min", "rising_phase", "near_peak", "solar_max"]
SCENARIO_LABELS = {"solar_min" : "Solar minimum (December 2019)",
                   "rising_phase" : "Rising phase (April 2022)",
                   "near_peak" : "Near-peak (December 2022)",
                   "solar_max" : "Solar maximum (August 2024)"}
SCENARIO_COLORS = {"solar_min" : "tab:cyan",
                   "rising_phase" : "tab:green",
                   "near_peak" : "tab:orange",
                   "solar_max" : "tab:red"}

def load_scenario(name):
    time = np.load(DATA_DIR / f"{name}_time.npy")
    altitude = np.load(DATA_DIR / f"{name}_altitude.npy")

    return time / 86400.0, altitude


# Plot 1: altitude vs. time, all four phases overlaid

plt.figure(figsize = (10, 7))
decay_values = []

for name in SCENARIO_ORDER:
    time_days, altitude = load_scenario(name)
    plt.plot(time_days, altitude, color = SCENARIO_COLORS[name], linewidth = 1.5, label = SCENARIO_LABELS[name])
    decay_values.append(altitude[0] - altitude[-1])

    start_altitude = altitude[0]
    end_altitude = altitude[-1]
    total_decay = start_altitude - end_altitude

    plt.annotate(f"Start: {start_altitude:.3f} km",
                 xy = (time_days[0], start_altitude),
                 xytext = (5, 5),
                 textcoords = "offset points",
                 fontsize = 8)

    plt.annotate(f"End: {end_altitude:.3f} km \n"
                 f"(Decay: {total_decay:.3f} km)",
                 xy = (time_days[-1], end_altitude),
                 xytext = (20, -20),
                 textcoords = "offset points",
                 ha = "right",
                 fontsize = 8)

plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")
plt.title("Drag-Driven Altitude Decay Across the Solar Cycle \n"
          "(500 km initial altitude, no propulsion)")
plt.legend()
plt.grid(True, alpha = 0.3)
plt.tight_layout()
plt.savefig(FIGURES_DIR / "solar_activity_sensitivity_altitude.png", dpi = 300)


# Plot 2: total decay vs. solar cycle phase, as a bar chart

plt.figure(figsize = (8, 5.5))
bars = plt.bar([SCENARIO_LABELS[name].split(" (")[0] for name in SCENARIO_ORDER],
               decay_values,
               color = [SCENARIO_COLORS[name] for name in SCENARIO_ORDER])
for bar, value in zip(bars, decay_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:.2f}", ha = "center", va = "bottom", fontsize = 10)

plt.ylabel("Total altitude decay (km)")
plt.title("Total Altitude Decay by Solar Cycle Phase")
plt.grid(True, alpha = 0.3, axis = "y")
plt.tight_layout()
plt.savefig(FIGURES_DIR / "solar_activity_sensitivity_bar.png", dpi = 300)

plt.show()

print("Summary:")
for name, decay in zip(SCENARIO_ORDER, decay_values):
    print(f"{SCENARIO_LABELS[name]:<30} : {decay:.3f} km decay")
print(f"\nRatio of altitude decay, solar max / solar min: {decay_values[-1] / decay_values[0]:.2f} times the solar minimum altitude decay")