"""
Plotting Script for Propulsive Burn Time Trade Study.

This script generates a plot of propulsive burn time as a function of available
onboard power for the CubeSat electric propulsion deorbit maneuver.
"""

import matplotlib.pyplot as plt
from power_and_burn_time_analysis import power_constraint_analysis, deorbit_burn_time

power_levels = [5, 10, 20, 30, 50]
burn_time = []

for power in power_levels:
    power_result = power_constraint_analysis(power)
    time_result = deorbit_burn_time(power_result)
    burn_time.append(time_result["burn_time_hr"])

# Plot

plt.figure()
plt.plot(power_levels, burn_time, marker="o")
plt.xlabel("Available Power (W)")
plt.ylabel("Propulsive Burn Time (hours)")
plt.title("Propulsive Burn Time as a Function of Available Onboard Power")
plt.grid(True)

# Save and show plot

plt.savefig("power_vs_burn_time.png")
plt.show()
