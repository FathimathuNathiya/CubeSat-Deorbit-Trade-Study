# CubeSat Deorbit Trade Study

A systems-level analysis of CubeSat deorbit strategies considering propellant mass, propulsion type, spacecraft mass sensitivity, power-limited electric propulsion burn time, atmospheric drag constraints, and coupled propulsion-drag effects.
This project uses Python-based analytical models and simulations to compare green chemical propulsion and electric propulsion for controlled CubeSat deorbit missions.

## Project Objectives

- Compute required deorbit delta-V for a representative LEO mission  
- Compare propellant mass requirements for different propulsion systems  
- Perform sensitivity analysis with respect to CubeSat dry mass  
- Analyze power constraints and resulting propulsive burn time for electric propulsion
- Model atmospheric drag and its contribution to orbital decay
- Compare propulsion-driven deorbit with natural atmospheric drag
- Investigates coupled low-thrust propulsion and atmospheric drag effects
- Track altitude evolution and propellant consumption during coupled simulations

## Methodology

The analysis combines orbital mechanics, propulsion modeling, and numerical simulation.

### Orbital Mechanics

- Circular-orbit velocity calculations
- Deorbit delta-V estimation
- Orbital decay analysis

### Propulsion Analysis

- Tsiolkovsky rocket equation for propellant mass estimation
- Chemical propulsion analysis
- Electric propulsion analysis
- Power-limited electric propulsion burn-time calculations
- Propellant mass evolution during thrusting

### Atmospheric drag

- Atmospheric drag modeling
- Numerical orbital propagation
- Altitude decay analysis
- Comparison between atmospheric drag and active propulsion

### Coupled Propulsion and Drag

The project also investigates the combined effect of continuous thrust and atmospheric drag on spacecraft altitude and propellant consumption.

## Repository Structure
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
|   ├── paths.py
|   ├── drag_deorbit_simulation.py
|   ├── thrust_drag_coupled_simulation.py
|   |── plot_propellant_comparison.py
│   ├── plot_mass_sensitivity.py
│   ├── plot_power_vs_thrust.py
│   ├── plot_power_vs_burn_time.py
|   ├── plot_drag_deorbit.py
|   └── plot_thrust_drag_coupled.py
│
├── figures/
|   |── mass_sensitivity.png
|   |── propellant_comparison.png
|   |── power_vs_thrust.png
│   ├── power_vs_burn_time.png
│   ├── drag_deorbit_decay.png
│   ├── thrust_drag_coupled_altitude.png
|   └── thrust_drag_coupled_propellant.png
|
├── data/
|   |── drag_decay_altitude.npy
|   |── drag_decay_time.npy
|   |── thrust_drag_electric_altitude.npy
│   ├── thrust_drag_electric_propellant.npy
│   ├── thrust_drag_electric_time.npy
│   ├── thrust_drag_green_altitude.npy
|   ├── thrust_drag_green_propellant.npy 
|   └── thrust_drag_green_time.npy
|
├── report/
│   └── CubeSat_Deorbit_Analysis.docx
│

```

## How to Run the Code

1. Install Python (3.8+ recommended)
2. Install required libraries:
   ```
   pip install numpy matplotlib scipy pymsis
3. Run analysis and simulation scripts from the `src` folder:
   ```
   python mission_analysis.py
   python sensitivity_mass_analysis.py
   python power_and_burn_time_analysis.py
   drag_deorbit_simulation.py
   thrust_drag_coupled_simulation.py
4. Generate plots:
   ```
   python plot_mass_sensitivity.py
   python plot_propellant_comparison.py
   python plot_power_vs_thrust.py
   python plot_power_vs_burn_time.py
   plot_drag_deorbit.py
   plot_thrust_drag_coupled.py

## Key Results

- Electric propulsion significantly reduces propellant mass compared to chemical propulsion.  
- Propellant mass increases strongly with CubeSat dry mass.  
- Available onboard power is the primary limiting factor for electric propulsion feasibility.  
- Low power levels result in long burn durations for deorbit maneuvers.
- Atmospheric drag can contribute significantly to long-term orbital decay.
- Active propulsion and atmospheric drag can be modeled together to study coupled deorbit behavior.
- Coupled thrust-drag simulations provide a more representative picture of spacecraft altitude evolution and propellant consumption.

## Tools Used

- Python  
- NumPy  
- Matplotlib
- SciPy
- pymsis
- datetime
- pathlib  
- Visual Studio Code

## Author

**Fathimathu Nathiya**<br>
Pre-University Research Project<br>
Year : 2026

## License

This project is intended for academic and educational purposes only.









