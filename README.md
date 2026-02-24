# Helmsdale Geothermal Feasibility Study: 2D Thermal Model

## Overview
This repository contains the MATLAB codebase for assessing the geothermal feasibility of the Helmsdale transect. The model implements a **2D variable-coefficient heat transport solver**, accounting for conductive heat transfer, radiogenic heat production, and fault-focused advection. 

This code was developed as part of the EARTH5016 consultancy project to evaluate the accessibility of 120°C geothermal targets within a 2–3 km drilling window.

## Repository Contents
The repository contains two standalone functional user scripts that reproduce all verification benchmarks, validation profiles, and parameter study results discussed in the final consultancy report:

1. **`project1.m`**: 
   - The primary baseline model. 
   - Constructs the geological geometry (He1/He2 granites, sedimentary basin, fault zone, and basement) based on the BGS A-B cross-section data pack.
   - Outputs: 2D thermal structure plots, depth-temperature profile comparisons against He1 borehole validation data, and a CSV table of key isotherm depths.

2. **`run_parameter_study.m`**: 
   - A standalone script that executes a 3x3 parameter sensitivity sweep.
   - Investigates the two primary physical controls: **Fault-zone advective transport strength ($KD$)** and **He2 Granite radiogenic heat production ($Q_r$)**.
   - Outputs: A summary CSV table (`helmsdale_paramstudy_summary.csv`) and three sensitivity plots demonstrating the nonlinear impacts of advection on isotherm uplift.

## Requirements
- **Software**: MATLAB (No special toolboxes are strictly required; base MATLAB is sufficient).
- **Version**: Tested on recent versions of MATLAB (e.g., R2022a or newer).

## Instructions: How to Run the Code

### 1. Running the Baseline Model & Validation
1. Open MATLAB and navigate to the directory containing the scripts.
2. Open `project1_final_report_v2.m`.
3. Click **Run** (or press `F5`).
4. **Expected Output**: The console will display the numerical setup (e.g., grid resolution, CFL time steps, and physical runtime) and track the progress of the explicit solver. Once complete, it will generate two figures (the 2D thermal structure and a geology zoning check) and save them as PNG files in the same directory. A CSV file named `helmsdale_finalreport_v2_isotherm_depths.csv` will also be generated.

### 2. Running the Parameter Sensitivity Study
1. Open `run_parameter_study.m` in MATLAB.
2. Click **Run** (or press `F5`).
3. **Expected Output**: The script will sequentially compute 9 modeling scenarios. The console will output a summary table showing the relationship between $Q_r$, $KD$, and the 100°C/120°C target depths. It will automatically generate and save three sensitivity plots and a summary data file (`helmsdale_paramstudy_summary.csv`).

## Numerical Design & Physical Rationale
To ensure the robustness of the numerical solution, the following designs were implemented:
- **Governing Physics**: Solves the unsteady advection-diffusion equation. The diffusion term is implemented in a conservative variable-coefficient form using harmonic mean averaging at geological interfaces.
- **Stability (CFL)**: Explicit time stepping is dynamically constrained by taking the minimum of the diffusive and advective limits, fortified by a safety factor.
- **Advection**: A first-order upwind scheme is explicitly utilized to suppress non-physical oscillations at the high-velocity fault boundary.
- **Simulation Runtime**: To avoid arbitrary end-times, the physical runtime is determined by the diffusive equilibration timescale across the domain depth ($\tau \approx L^2 / \pi^2 \kappa_{rep}$).