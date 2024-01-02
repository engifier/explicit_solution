# Explicit Saint-Venant Equations Solver

MATLAB implementation for solving the Saint-Venant equations, a set of partial differential equations describing unsteady open-channel water flow. This solver utilizes an explicit finite-difference method, providing a tool for simulating and analyzing the dynamic behavior of water flow over time.

## Usage

1. **Input Parameters:**
   - `L`: Length of the channel or river.
   - `b`: Width of the channel or river.
   - `s0`: Slope.
   - `y0`: Initial water depth.
   - `n`: Manning's roughness coefficient.
   - `time`: Total computation time in minutes.
   - `N`: Number of sections.

2. **Run the Script:**
   - Execute the script in MATLAB.
   - Input the required parameters when prompted.

3. **Results:**
   - The code generates X-Y graphs at different time intervals, illustrating the water depth along the channel.

## Visualization

The script produces visualizations that can help understand the dynamic changes in water depth over the channel length. The graphs are categorized based on special times during the computation.

## Disclaimer

This code is provided as-is, without any warranty. Use it at your own risk.

## Author

Hamed Valaei


## Acknowledgments

- This project was developed as part of the Advanced Hydraulics Programming Project.
