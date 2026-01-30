"""
Plotting script for Power Constraint Analysis.

This script generates a plot of thrust as a function of available onboard power
for the electric propulsion system used in the CubeSat deorbit study. 
"""

import matplotlib.pyplot as plt
from power_and_burn_time_analysis import power_constraint_analysis

power_levels = [5, 10, 20, 30, 50]
thrust_values = []

for power in power_levels:
    result = power_constraint_analysis(power)
    thrust_values.append(result["thrust"])

# Plot

plt.figure(figsize=(5, 5))
plt.plot(power_levels, thrust_values, marker="o")
plt.xlabel("Available Power (W)")
plt.ylabel("Thrust (N)")
plt.title("Thrust as a Function of Available Onboard Power")
plt.grid(True)

# Save and show plot

plt.savefig("power_vs_thrust.png", dpi=300, bbox_inches="tight")
plt.show()
