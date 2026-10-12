# CubeSat Deorbit Trade Study

A systems-level trade study of CubeSat deorbit strategies under mass and power constraints, combining an analytical propulsion comparison with a numerical orbit propagation framework modeling atmospheric drag, coupled propulsion-drag dynamics, and solar activity sensitivity.

## Research Objectives

- Compare propellant mass requirements for green chemical and electric propulsion in a representative CubeSat deorbit maneuver (500 km → 120 km), using the 
  Tsiolkovsky rocket equation.

- Quantify how propellant requirements scale with CubeSat dry mass (5–25 kg) and how electric propulsion thrust and burn duration scale with available onboard 
  power (5–50 W).

- Extend this analytical baseline with a numerical propagation framework to examine three physical effects not captured by closed-form delta-V estimates:
  1. **Aerodynamic drag** acting on the orbit independently of propulsion.
  2. **Coupled propulsion-drag interaction**, including how drag's relative contribution to total decay varies across thrust regimes.
  3. **Sensitivity to solar activity**, via its effect on atmospheric density and decay rate.

*This extended scope was developed in response to detailed technical feedback from Dr. Dario Modenini (University of Bologna, u3S Laboratory).*

## Repository Structure

```text
CubeSat_Deorbit_Trade_Study/
│
├── README.md
│
├── src/
|   **Analytical Modeling Scripts**
|   |── orbital_velocity.py
|   |── deorbit_delta_v.py
|   |── propellant_mass.py
│   ├── mission_analysis.py
│   ├── sensitivity_mass_analysis.py
│   ├── power_and_burn_time_analysis.py
│   │
│   │**Shared Configuration** 
│   |── paths.py *(paths.py defines common directory paths used by the numerical simulation scripts to save simulation data and generated figures.)*
│   |
|   │**Numerical Simulation Scripts**
│   |── drag_deorbit_simulation.py
│   |── thrust_drag_coupled_simulation.py
│   |── thrust_regime_sweep.py
│   |── solar_activity_sensitivity.py
│   |
|   │**Plotting Scripts - Analytical Results**
|   |── plot_propellant_comparison.py
│   ├── plot_mass_sensitivity.py
│   ├── plot_power_vs_thrust.py
│   ├── plot_power_vs_burn_time.py
│   |
|   │**Ploting Scripts - numerical Results**
│   ├── plot_drag_deorbit.py
│   ├── plot_thrust_drag_coupled.py
│   ├── plot_thrust_regime_sweep.py
│   └── plot_solar_activity_sensitivity.py
│
├── figures/
│   │**Figures - Analytical Results**
│   ├── propellant_mass_comparison.png
│   ├── mass_sensitivity_analysis.png
│   ├── power_vs_thrust.png
│   ├── power_vs_burn_time.png
│   │
│   │**Figures - Numerical Results**  
│   ├── drag_deorbit_decay.png
│   ├── thrust_drag_coupled_altitude.png
│   ├── thrust_drag_coupled_propellant.png
│   ├── total_decay_by_regime.png
│   ├── drag_relative_contribution.png
│   ├── solar_activity_altitude_decay.png
│   └── solar_phase_decay_comparison.png
│ 
├── data/
│    ├── drag_decay_time.npy
│    ├── drag_decay_altitude.npy
│    ├── thrust_drag_electric_time.npy
│    ├── thrust_drag_electric_altitude.npy
│    ├── thrust_drag_electric_propellant.npy
│    ├── thrust_drag_green_time.npy
│    ├── thrust_drag_green_altitude.npy
│    ├── thrust_drag_green_propellant.npy
│    ├── regime_power.npy
│    ├── regime_thrust.npy
│    ├── regime_burn_duration_hr.npy
│    ├── regime_burn_within_window.npy
│    ├── regime_elapsed_days.npy
│    ├── regime_reentry_reached.npy
│    ├── regime_total_decay.npy
│    ├── regime_drag_share_percent.npy
│    ├── solar_min_time.npy
│    ├── solar_min_altitude.npy
│    ├── rising_phase_time.npy
│    ├── rising_phase_altitude.npy
│    ├── near_peak_time.npy
│    ├── near_peak_altitude.npy
│    ├── solar_max_time.npy
│    ├── solar_max_altitude.npy
│    ├── solar_scenario_name.npy
│    ├── solar_scenario_decay_km.npy
│    ├── solar_scenario_elapsed_days.npy
│    └── solar_scenario_reentry_detected.npy   
│
└── report/
    └── CubeSat_Deorbit_Trade_Study.pdf
```

## Methodology Summary

**Analytical Stage**
- Orbital velocity and deorbit Δv calculated using orbital mechanics and a Hohmann-transfer approximation.
- Propellant mass estimated using the Tsiolkovsky rocket equation.
- Sensitivity analysis performed across spacecraft dry mass and available onboard power to evaluate propellant requirements and electric propulsion performance.
- Electric propulsion burn duration calculated for different power levels under the assumed operating conditions.

**Numerical Stage**
- Equations of motion integrated in the Earth-centred inertial (ECI) frame, accounting for spherical two-body gravity and atmospheric drag.
-	Atmospheric density obtained using the NRLMSISE-00 model through pymsis, with geodetic coordinates calculated using the WGS84 ellipsoid and Bowring’s method.
- Drag-only orbital decay simulated to evaluate natural altitude reduction without propulsion.
- Coupled propulsion–drag simulations performed with continuous retrograde thrust, while tracking propellant consumption as an additional state variable.
- Thrust-regime sensitivity analysed across different electric-propulsion power levels.
- Solar-activity sensitivity investigated using four reference dates spanning Solar Cycle 25.
- Deorbit defined as a terminal event when the spacecraft reaches an altitude of 120 km.

## Key Results

**Analytical Results**
- Electric propulsion requires substantially less propellant than green chemical propulsion for the same deorbit manoeuvre.
- Propellant requirements increase with spacecraft dry mass.
- Available onboard power is a primary constraint on electric propulsion, with lower power levels resulting in longer burn durations.

**Numerical Simulation Results**
- **Drag-only orbital decay:** The simulated spacecraft altitude decreased from 500 km to approximately 497.03 km over one year.
- **Coupled propulsion–drag simulations:** The spacecraft reached the 120 km deorbit threshold after 163.7 days with electric propulsion and 62.1 days with  
  green chemical propulsion.
- **Thrust-regime sensitivity:** Atmospheric drag contributed approximately 0.43–0.48% of the total altitude decay across the tested electric-propulsion power 
  levels of 5 – 50 W.
- **Solar-activity sensitivity:** Over one-year simulations, altitude decay ranged from approximately 0.52 km under solar-minimum conditions to 7.04 km under 
  solar-maximum conditions, a factor of approximately 13.57.

## Tools and Libraries Used

- Python  
- NumPy  
- SciPy
- Matplotlib  
- pymsis
- datetime
- pathlib
- Visual Studio Code

## How to Run the Code

1. Install Python (3.8+ recommended)
2. Install required libraries:
   ```bash
   pip install numpy scipy matplotlib pymsis
   ```
3. Run the following scripts:
   **Analytical Baseline**
   ```bash
   python src/mission_analysis.py
   python src/plot_propellant_comparison.py
   python src/sensitivity_mass_analysis.py
   python src/plot_mass_sensitivity.py
   python src/power_and_burn_time_analysis.py
   python src/plot_power_vs_thrust.py
   python src/plot_power_vs_burn_time.py
   ``` 

   **Extension 1: Drag Only Orbital Decay**
   ```bash
   python src/drag_deorbit_simulation.py
   python src/plot_drag_deorbit.py
   ```

   **Extension 2: Coupled Thrust and Atmospheric Drag**
   ```bash
   python src/thrust_drag_coupled_simulation.py
   python src/plot_thrust_drag_coupled.py
   python src/thrust_regime_sweep.py
   python src/plot_thrust_regime_sweep.py
   ```

   **Extension 3: Solar Activity Sensitivity**
   ```bash
   python src/solar_activity_sensitivity.py
   python src/plot_solar_activity_sensitivity.py
   ```

## Author

Fathimathu Nathiya

Year : 2026
