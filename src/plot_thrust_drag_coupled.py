"""
Plot coupled thrust + drag altitude decay for both propulsion scenarios.

Reads the time, altitude, and propellant history saved by
thrust_drag_coupled_simulation.py for both the electric and green
propulsion scenarios, and plots them together for direct comparison.
Optionally overlays the drag_only baseline from
drag_deorbit_simulation.py, to visualize thrust's contribution 
relative to drag alone.

If either scenario reached reentry before the configured simulation
duration elapsed (see reentry_event in drag_deorbit_simulation.py), the
saved arrays already stop at that point, so no special handling is
needed here - the curves simply end at their respective reentry times.

See the framing note at the top of thrust_drag_coupled_simulation.py:
both scenarios spend the report's analytically sized propellant budget
(Table 3) as a continuous burn, rather than reproducing the report's
idealized instantaneous 500 km -> 120 km tranfer.

This is the figure corresponding to the second extension suggested by
Dr. Dario Modenini (University of Bologna): exploring coupled scenarios
where low-thrust propulsion and atmospheric drag interact. The electric
propulsion scenario below directly addresses this (low-thrust, long 
burn duration); the green chemical propulsion scenario is included as 
an additional comparison point, using a comparatively high thrust, 
short-duration burn, to contrast against the low-thrust case.
"""


import numpy as np
import matplotlib.pyplot as plt

from paths import DATA_DIR, FIGURES_DIR


def load_scenario(name):
    """Load the saved time, altitude, and propellant arrays for one scenario."""

    time = np.load(DATA_DIR / f"thrust_drag_{name}_time.npy")                      # (s)
    altitude = np.load(DATA_DIR / f"thrust_drag_{name}_altitude.npy")              # (km)
    propellant_mass = np.load(DATA_DIR / f"thrust_drag_{name}_propellant.npy")     # (kg)

    return time / 86400.0, altitude, propellant_mass                               # time converted to days


electric_time, electric_altitude, electric_propellant_mass = load_scenario("electric")
green_time, green_altitude, green_propellant_mass = load_scenario("green")


# Optional drag-only baseline (from drag_deorbit_simulation.py)

drag_only_time = np.load(DATA_DIR / "drag_decay_time.npy") / 86400.0               # (days)
drag_only_altitude = np.load(DATA_DIR / "drag_decay_altitude.npy")                 # (km)


# Plot 1: Altitude vs. time, both propulsion scenarios (along with drag-only baseline)

plt.figure(figsize = (9,6))
plt.plot(electric_time, electric_altitude, color = "tab:orange", linewidth = 1.5, label = "Electric propulsion + drag")
plt.plot(green_time, green_altitude, color = "tab:green", linewidth = 1.5, label = "Green propulsion + drag")
plt.plot(drag_only_time, drag_only_altitude, color = "tab:blue", linewidth = 1.2, linestyle = "--", label = "Drag only (no thrust)")

plt.axhline(y = 120.0, color = "gray", linewidth = 1.5, linestyle = ":", label = "Reentry threshold (120 km)")

plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")
plt.title("CubeSat Altitude Decay: Coupled Thrust + Drag \n"
          "(500 km initial altitude; report's propellant budget spent as a continuous burn)")
plt.legend()
plt.grid(True, alpha = 0.3)

electric_decay = electric_altitude[0] - electric_altitude[-1]
green_decay = green_altitude[0] - green_altitude[-1]
drag_only_decay = drag_only_altitude[0] - drag_only_altitude[-1]

plt.annotate(f"Start: {electric_altitude[0]:.3f} km",
             xy = (electric_time[0], electric_altitude[0]),
             xytext = (1, 1),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {electric_altitude[-1]:.3f} km \n"
             f"(t: {electric_time[-1]:.1f} days) \n"
             f"(decay: {electric_decay:.3f} km)",
             xy = (electric_time[-1], electric_altitude[-1]),
             xytext = (1, -2),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"Start: {green_altitude[0]:.3f} km",
             xy = (green_time[0], green_altitude[0]),
             xytext = (1, 1),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {green_altitude[-1]:.3f} km \n"
             f"(t: {green_time[-1]:.1f} days) \n"
             f"(decay: {green_decay:.3f} km)",
             xy = (green_time[-1], green_altitude[-1]),
             xytext = (1, -2),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"Start: {drag_only_altitude[0]:.3f} km",
             xy = (drag_only_time[0], drag_only_altitude[0]),
             xytext = (1, 1),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {drag_only_altitude[-1]:.3f} km \n"
             f"(t: {drag_only_time[-1]:.1f} days) \n"
             f"(decay: {drag_only_decay:.3f} km)",
             xy = (drag_only_time[-1], drag_only_altitude[-1]),
             xytext = (-60, -30),
             textcoords = "offset points",
             fontsize = 9)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "thrust_drag_coupled_altitude.png", dpi = 300)

electric_decay = electric_altitude[0] - electric_altitude[-1]
green_decay = green_altitude[0] - green_altitude[-1]
drag_only_decay = drag_only_altitude[0] - drag_only_altitude[-1]
print(f"Electric propulsion + drag: {electric_decay:.3f} km altitude decay over "
      f"{electric_time[-1]:.1f} days (end altitude: {electric_altitude[-1]:.3f} km)")
print(f"Green propulsion + drag: {green_decay:.3f} km altitude decay over "
      f"{green_time[-1]:.1f} days (end altitude: {green_altitude[-1]:.3f} km)")
print(f"Drag only (baseline): {drag_only_decay:.3f} km altitude decay over "
      f"{drag_only_time[-1]:.1f} days (end altitude: {drag_only_altitude[-1]:.3f} km)")
print(f"\nNote: comparing thrust's contribution against the drag-only baseline is only meaningful "
      f"over the same elapsed time - the baseline and thrust scenarios above may span different "
      f"durations if a scenario reached reentry early.")


# Plot 2: Propellant mass remaining vs. time, both propulsion scenarios

plt.figure(figsize = (9,5))
plt.plot(electric_time, electric_propellant_mass, color = "tab:orange", linewidth = 1.5, label = "Electric propulsion")
plt.plot(green_time, green_propellant_mass, color = "tab:green", linewidth = 1.5, label = "Green propulsion")

plt.xlabel("Time (days)")
plt.ylabel("Propellant remaining (kg)")
plt.title("Propellant Consumption During Burn")
plt.legend()
plt.grid(True, alpha = 0.3)

electric_propellant_remaining = electric_propellant_mass[0] - electric_propellant_mass[-1]
green_propellant_remaining = green_propellant_mass[0] - green_propellant_mass[-1]

plt.annotate(f"Start: {electric_propellant_mass[0]:.4f} kg",
             xy = (electric_time[0], electric_propellant_mass[0]),
             xytext = (0, 1),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {max(electric_propellant_mass[-1], 0.0):.4f} kg \n"
             f"(t: {electric_time[-1]:.1f} days) \n"
             f"(propellant remaining: {electric_propellant_remaining:.4f} kg)",
             xy = (electric_time[-1], electric_propellant_mass[-1]),
             xytext = (-130, 2),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"Start: {green_propellant_mass[0]:.4f} kg",
             xy = (green_time[0], green_propellant_mass[0]),
             xytext = (0, 1),
             textcoords = "offset points",
             fontsize = 9)
plt.annotate(f"End: {max(green_propellant_mass[-1], 0.0):.4f} kg \n"
             f"(t: {green_time[-1]:.1f} days) \n"
             f"(propellant remaining: {green_propellant_remaining:.4f} kg)",
             xy = (green_time[-1], green_propellant_mass[-1]),
             xytext = (-30, 2),
             textcoords = "offset points",
             fontsize = 9)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "thrust_drag_coupled_propellant.png", dpi = 300)

plt.show()
