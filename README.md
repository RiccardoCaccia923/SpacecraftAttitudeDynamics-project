# Spacecraft Attitude Determination and Control (ADCS) Simulator
This repository contains a full Guidance, Navigation, and Control (GNC) simulation environment for a Mini-Satellite. 
The project is developed in **MATLAB & Simulink** and models the spacecraft's orbital mechanics, attitude dynamics, sensor measurements, state estimation, and control execution.

## Features
* **Attitude Dynamics & Kinematics**: Rigid body dynamics with DCM-based kinematics.
* **Environmental Models**: Includes orbital propagation and external disturbance torques (perturbations).
* **Sensor Suite Models**: 
  * High-precision **Star Tracker**.
  * Three-axis **Magnetometer**.
* **Estimation (Navigation)**: 
  * A **Multiplicative Extended Kalman Filter (MEKF)** for full state estimation (Attitude and Angular Rates) in a gyro-less configuration.
* **Control Modes**:
  * **Detumbling**: B-Dot control law using Magnetic Torquers to dissipate initial tip-off kinetic energy.
  * **Fine Pointing**: Reaction Wheel control for precise attitude targeting.

## Repository Structure
The project is modularized to separate configurations, models, and post-processing tools:

* `main.m` - The primary entry point script. Initializes the workspace, loads parameters, runs the simulation, and plots the results.
* `mainModel.slx` - The core Simulink model containing the ADCS architecture.
* `/configs/` - Initialization scripts for various subsystems.
* `/plots/` - Scripts and functions for data visualization and post-processing.
* `/miscellaneous/` - Helper functions and utilities.
* `/validations/` - Scripts used for unit testing and subsystem validation.
* `/submodules/` - Additional system dependencies and block configurations.

## Prerequisites
To run this simulation, you will need:
* **MATLAB** (Tested on R2024b)
* **Simulink**
* Aerospace Toolbox (Optional but recommended for certain orbital functions)

## How to Run
Running the simulation is straightforward. 
Everything is handled by the main script.

## Execution Time Optimization
Running the full script can be time-consuming due to the heavy computational load of the plotting phase. 
If you only want to run the core simulation without the visual rendering, 
it is highly recommended to **comment out the plotting section**. 
The relevant visual results and graphs are already documented and available in the provided `SAD_report.pdf` file.