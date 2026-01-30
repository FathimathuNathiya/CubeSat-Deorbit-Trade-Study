"""
Plot of propellant mass sensitivity with respect to CubeSat dry mass
for Green and Electric propulsion systems.
"""

import matplotlib.pyplot as plt
from sensitivity_mass_analysis import compute_mass_sensitivity

satellite_masses, green_propellant, electric_propellant = compute_mass_sensitivity()

# Create plot

plt.figure(figsize=(6, 5))
plt.plot(satellite_masses, green_propellant,
         marker="o", label="Green Propulsion")
plt.plot(satellite_masses, electric_propellant,
         marker="o", label="Electric Propulsion")
plt.xlabel("CubeSat Dry Mass (kg)")
plt.ylabel("Required Propellant Mass (kg)")
plt.title("Sensitivity of Propellant Mass to CubeSat Dry Mass")
plt.grid(True)
plt.legend()

# Save and show plot

plt.savefig("mass_sensitivity_analysis.png", dpi=300, bbox_inches="tight")
plt.show()
