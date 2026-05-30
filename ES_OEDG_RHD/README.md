# deal.II Entropy-Stable OEDG Solver for 2D Relativistic Hydrodynamics

This project implements an **Entropy-Stable Oscillation-Eliminating Discontinuous Galerkin (ES-OEDG)** scheme for the 2D special-relativistic Euler (RHD) equations, built on the [deal.II](https://www.dealii.org/) finite-element library.

The spatial discretization uses the **SBP split-form** (summation-by-parts) volume integral with an entropy-conservative (EC) two-point flux, combined with Rusanov SAT interface corrections, guaranteeing a global discrete entropy inequality. An oscillation-eliminating (OE) modal damping pass and a physical-constraints-preserving (PCP) limiter are applied after each RK stage.

## Features

- **Entropy-stable DGSEM**: SBP split-form volume integral + Rusanov SAT — provably satisfies a discrete entropy inequality.
- **EC two-point flux** (`conflux`): RHD entropy-conservative flux of Bhoriya & Kumar (ZAMP 2020).  
  ⚠ The axis parameters `(dir_x, dir_y)` are **coordinate-axis selectors**, not continuous direction cosines.  Only `(1,0)` (x-axis) and `(0,1)` (y-axis) are valid.
- **OE Limiter**: Multi-index derivative-jump modal damping (Peng, Sun & Wu 2024).
- **PCP Limiter**: Positivity preservation of baryon density and thermodynamic constraint.
- **SSP-RK3** time integration (Shu–Osher 1988).
- **Nodal DG on GLL nodes**: High-order accuracy; diagonal mass matrix (mass-matrix cancellation in ES residual).
- **Parallelized**: Uses deal.II `WorkStream` across all available CPU cores.
- **RP2D1 Problem**: Pre-configured for the 4-quadrant 2D Riemann problem.

## Prerequisites

- **deal.II**: Version 9.7.1 (installed in `/home/biswarup-biswas/dealii-install`).
- **CMake**: Version 3.13 or higher.
- **C++ Compiler**: GCC or Clang with C++17 support.
- **Python 3**: With `matplotlib`, `numpy`, and `pyvista` for visualization. A virtual environment `.venv` is provided.

## Project Structure

```
ES_OEDG_RHD/
├── rhd_dg.cc          Main solver: ES-OEDG scheme, flux helpers, limiters
├── CMakeLists.txt     Build configuration
├── plot_solution.py   Density contour plot (matplotlib)
├── plot_pyvista.py    High-quality visualization (PyVista)
└── Reading/
    └── my-files/
        ├── THEORY.tex  Full theory document (LaTeX source)
        └── THEORY.pdf  Compiled theory document
```

### Key functions in `rhd_dg.cc`

| Function | Description |
|---|---|
| `logavg(a,b)` | Numerically stable logarithmic mean |
| `calflux(U, dx, dy, f)` | Physical RHD flux along coordinate axis |
| `conflux(Ul, Ur, dx, dy, f)` | EC two-point flux — Bhoriya & Kumar (2020) |
| `rusflux(Ul, Ur, dx, dy, f)` | Rusanov (LLF) SAT interface flux |
| `compute_es_rhs(...)` | ES-OEDG residual: SBP volume + Rusanov SAT |
| `apply_oe_limiter(...)` | OE modal damping pass |
| `apply_pcp_limiter(...)` | PCP positivity-preservation pass |

## How to Build

```bash
cd /home/biswarup-biswas/git_workspace/DEAL-II-Codes/ES_OEDG_RHD
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## How to Run

```bash
./build/rhd_dg
```

The solver prints progress (`t=0.05`, `t=0.1`, …) and writes VTK output files.

## Visualization

```bash
source .venv/bin/activate
python3 plot_solution.py      # basic matplotlib density plot
python3 plot_pyvista.py       # high-resolution PyVista plot
```

## Changing Simulation Parameters

All major parameters are in `rhd_dg.cc`:

| Parameter | Location | Default |
|---|---|---|
| Grid resolution | `setup_system()` → `refine_global(N)` | `4` ($16\times16$ cells) |
| Polynomial degree | `main()` → `RHD_DG(p)` | `2` (P2, 3rd order) |
| Final time | `run()` → `while (time < T)` | `0.4` |
| Initial conditions | `InitialCondition::vector_value` | RP2D1 |

## Theory

A detailed derivation covering the RHD equations, SBP-DGSEM formulation, EC flux, entropy stability proof, OE/PCP limiters, and SSP-RK3 is in [`Reading/my-files/THEORY.pdf`](Reading/my-files/THEORY.pdf).

## References

1. **Bhoriya, D. and Kumar, H.** (2020). *Entropy-stable schemes for relativistic hydrodynamics equations*. ZAMP **71**(1), Art. 29. DOI:[10.1007/s00033-020-1250-8](https://doi.org/10.1007/s00033-020-1250-8)
2. **Yang, J. and Fu, G.** (2026). *An entropy-stable oscillation-eliminating DGSEM for the Euler equations on curvilinear meshes*. arXiv:2602.16732.
3. **Peng, M., Sun, Z., and Wu, K.** (2024). *OEDG: Oscillation-eliminating discontinuous Galerkin method for hyperbolic conservation laws*. Math. Comp. DOI:10.1090/mcom/3998.
4. **Carpenter, M. H. et al.** (2014). *Entropy stable spectral collocation schemes for the Navier-Stokes equations*. SIAM J. Sci. Comput. **36**(5):B835–B867.
5. **Chen, T. and Shu, C.-W.** (2017). *Entropy stable high order DG methods with suitable quadrature rules*. J. Comput. Phys. **345**:427–461.
6. **Shu, C.-W. and Osher, S.** (1988). *Efficient implementation of essentially non-oscillatory shock-capturing schemes*. J. Comput. Phys. **77**(2):439–471.
