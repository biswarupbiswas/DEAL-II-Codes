# deal.II Relativistic Hydrodynamics (RHD) DG Solver

This project is a high-performance C++ port of a Discontinuous Galerkin (DG) solver for Relativistic Hydrodynamics, originally implemented in C/PETSc. It features a mathematically rigorous port of the **Oscillation-Eliminating (OE) Limiter** based on the research by Peng, Sun & Wu (2024).

## Features
- **Nodal DG Formulation**: High-order accuracy using Legendre-Gauss-Lobatto nodes.
- **Exact RHD Physics**: Robust Newton-solver for primitive variable recovery (`con2prim`).
- **OEDG Limiter**: Full implementation of the multi-index derivative jump damping logic.
- **Parallelized**: Built using deal.II's `WorkStream` to utilize all available CPU cores (Optimized for 12+ cores).
- **RP2D1 Problem**: Pre-configured for the 4-quadrant 2D Riemann Problem.

## Prerequisites
- **deal.II**: Version 9.7.1 (installed in `/home/biswarup-biswas/dealii-install`).
- **CMake**: Version 3.13 or higher.
- **C++ Compiler**: GCC or Clang with C++17 support.
- **Python 3**: With `matplotlib` and `numpy` for visualization.

## Project Structure
- `rhd_dg.cc`: The main C++ source code containing the solver logic and physics kernels.
- `CMakeLists.txt`: Build configuration file.
- `plot_solution.py`: Python script to generate density contour plots from VTK output.

## How to Build
1. Navigate to the project directory:
   ```bash
   cd /home/biswarup-biswas/git_workspace/RHD-ESDG/dealii_rhd
   ```
2. Generate the Makefile using CMake (Force Release mode for maximum speed):
   ```bash
   cmake -DCMAKE_BUILD_TYPE=Release .
   ```
3. Compile using all CPU cores:
   ```bash
   make -j12
   ```

## How to Run
Execute the compiled binary:
```bash
./rhd_dg
```
The solver will output progress to the terminal (e.g., `t=0.05`, `t=0.1`).

## Visualization
Once the simulation is complete, generate a plot of the final state:
```bash
python3 plot_solution.py
```
This will create a high-resolution file named `density_plot.png`.

## Changing Simulation Parameters
All major parameters are located at the top or in the setup functions of `rhd_dg.cc`:

- **Grid Resolution**: In `setup_system()`, change `triangulation.refine_global(4)` to `6` for a $64 \times 64$ grid.
- **Polynomial Degree**: In `main()`, change `RHD_DG(2)` to `RHD_DG(3)` for 4th-order accuracy.
- **Final Time**: In `run()`, modify `while (time < 0.4)`.
- **Initial Conditions**: Modify the `InitialCondition::vector_value` function to change the $\rho, u, p$ values for different Riemann problems.

## References
- Peng, Sun & Wu (2024). *Oscillation-Eliminating Discontinuous Galerkin Methods for Relativistic Hydrodynamics*. arXiv:2310.04807.
