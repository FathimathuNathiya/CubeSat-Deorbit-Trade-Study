# CubeSat Deorbit Trade Study

A systems-level analysis of CubeSat deorbit strategies considering propellant mass, spacecraft mass sensitivity, and power-limited electric propulsion burn time.
This project uses Python-based analytical models to compare green chemical propulsion and electric propulsion for controlled CubeSat deorbit missions.

## Project Objectives

- Compute required deorbit delta-V for a representative LEO mission  
- Compare propellant mass requirements for different propulsion systems  
- Perform sensitivity analysis with respect to CubeSat dry mass  
- Analyze power constraints and resulting propulsive burn time for electric propulsion

## Folder Structure
```
CubeSat_Deorbit_Trade_Study/
│
├── README.md
│
├── src/
|   |── orbital_velocity.py
|   |── deorbit_delta_v.py
|   |── propellant_mass.py
│   ├── mission_analysis.py
│   ├── sensitivity_mass_analysis.py
│   ├── power_and_burn_time_analysis.py
|   |── plot_propellant_comparison.py
│   ├── plot_mass_sensitivity.py
│   ├── plot_power_vs_thrust.py
│   └── plot_power_vs_burn_time.py
│
├── report/
│   └── CubeSat_Deorbit_Analysis.docx
│
└── figures/
    ├── mass_sensitivity.png
    ├── propellant_comparison.png
    ├── power_vs_thrust.png
    └── power_vs_burn_time.png```

## Methodology Summary

- Orbital mechanics used to compute deorbit delta-V between circular orbits  
- Tsiolkovsky rocket equation used to calculate propellant mass  
- CubeSat dry mass varied from 5 kg to 25 kg for sensitivity analysis  
- Electric propulsion thrust limited by available onboard power  
- Burn time computed for fixed delta-V under different power levels

## How to Run the Code

1. Install Python (3.8+ recommended)
2. Install required libraries:
   pip install numpy matplotlib
3. Run analysis scripts from the `src` folder:
   python mission_analysis.py
   python sensitivity_mass_analysis.py
   python power_and_burn_time_analysis.py
4. Generate plots:
   python plot_mass_sensitivity.py
   python plot_propellant_comparison.py
   python plot_power_vs_thrust.py
   python plot_power_vs_burn_time.py

## Key Results

- Electric propulsion significantly reduces propellant mass compared to chemical propulsion  
- Propellant mass increases strongly with CubeSat dry mass  
- Available onboard power is the primary limiting factor for electric propulsion feasibility  
- Low power levels result in long burn durations for deorbit maneuvers

## Tools Used

- Python  
- NumPy  
- Matplotlib  
- Visual Studio Code

## Author

Fathimathu Nathiya
Pre-University Research Project
Year : 2026

## License


This project is intended for academic and educational purposes only.



