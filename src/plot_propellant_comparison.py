"""
Plot comparison of propellant mass required for Green vs Electric propulsion
using results from the mission analysis model.
"""

import matplotlib.pyplot as plt
from mission_analysis import mission_analysis

results = mission_analysis()

# Extract propellant mass values from results dictionary

propulsion_systems = ["Green Propulsion", "Electric Propulsion"]
propellant_mass = [results["propellant_green"], results["propellant_electric"]]

# Create bar plot

plt.figure(figsize=(6, 6))
plt.bar(propulsion_systems, propellant_mass)
for i, v in enumerate(propellant_mass):
    plt.text(i, v + 0.01, f"{v:.4f}", ha="center")
plt.xlabel("Propulsion System")
plt.ylabel("Required Propellant Mass (kg)")
plt.title("Propellant Mass Comparison for CubeSat Deorbit")
plt.grid(axis="y", linestyle="--", alpha=0.6)

# Save and show plot

plt.savefig("propellant_mass_comparison.png",
            dpi=300, bbox_inches="tight")
plt.show()
