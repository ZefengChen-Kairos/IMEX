# Semi-implicit IMEX MHD solver with constrained transport

Reference implementation of the scheme described in:

> Zefeng Chen, Riccardo Demattè, Walter Boscheri, Stephen Millmore,
> *A Structure- and Pressure-Positivity-Preserving Semi-implicit IMEX Finite
> Volume Scheme for Ideal MHD at All Acoustic Mach and Alfvén Mach Numbers with
> Generic Equation of State*, arXiv:2608.15837 (2026).
> https://arxiv.org/abs/2608.15837

A 2D finite-volume solver for the ideal compressible MHD equations. The system
is split into three sub-systems by characteristic wave scale: an advective part
for hydrodynamic transport, a magnetic part for the velocity-field coupling, and
a pressure part for the pressure-velocity coupling. The nonlinear advective
terms are integrated explicitly; the other two sub-systems are implicit and are
solved matrix-free with GMRES. The time step is therefore restricted by the
material velocity rather than by the fast magnetosonic speed, so the scheme
stays efficient in gas-pressure dominated, magnetic-pressure dominated and
incompressible limits alike.

Pressure positivity is enforced through a nested-Newton solve of the pressure
sub-system, in which the diagonal term is replaced by a piecewise-linear
approximation below a pressure floor. The equation of state enters only through
that diagonal term and the sound speed, so ideal-gas and Redlich-Kwong closures
share the same machinery; the ideal-gas path is unchanged from the linear
formulation.

The divergence constraint on the magnetic field is handled by constrained
transport, with several variants available (staggered, unstaggered via the
vector potential, and an EMF-consistent staggered variant built from the fluxes
the implicit magnetic sub-system actually uses).

## Building

Requires a C++17 compiler and two header-only libraries:

| dependency | Debian/Ubuntu | macOS (Homebrew) |
|---|---|---|
| Eigen 3 | `apt install libeigen3-dev` | `brew install eigen` |
| Boost.Math | `apt install libboost-math-dev` | `brew install boost` |

Then:

```sh
make            # optimised build -> build/imex_mhd
make debug      # -O0 -g with address and UB sanitizers
make clean
```

The build is clean under `-Wall -Wextra`.

The Makefile searches the usual locations for both libraries. If they are
somewhere else:

```sh
make EIGEN_INC=/opt/homebrew/include/eigen3 BOOST_INC=/opt/homebrew/include
```

## Running

```sh
./build/imex_mhd                    # writes to ./output/...
./build/imex_mhd /path/to/results   # or set MHD_OUTPUT_ROOT
```

Output is written to a directory tree keyed by the run configuration, e.g.

```
output/Ideal Gas/PreservingON/MHD_vortex/PressureFormula/NoDivFree/definition/10x10_IMEXdt1_RK2_Order2_CFL080/
```

containing per-variable `.dat` files (`RK_rho.dat`, `RK_p.dat`, `RK_vx.dat`, …),
a divergence-error and time-step history, and `plt_*.vtk` snapshots.

## Configuration

The run is configured at **compile time** by the globals at the top of
`include/config.hpp`. The main ones:

| variable | meaning |
|---|---|
| `test` | which test case from the `Test` enum to run |
| `nxCells`, `nyCells` | grid resolution |
| `formula` | `EnergyFormula` or `PressureFormula` — which implicit subsystem to solve |
| `magfield` | CT variant: `StaggerCT`, `UnStaggerCT`, `PotentialA`, `NoDivFree` |
| `recon` | reconstruction: `NoRecon`, `Minmod`, `VanLeer`, `ENO3`, `WENO5` |
| `caltimestep` | `IMEXdt`, `Alfvendt` or `Fullydt` time-step restriction |
| `eos_type` | `IDEAL_GAS` or `REDLICH_KWONG` |
| `floor_on`, `p_floor`, `theta` | pressure positivity via the nested-Newton solve |
| `CFL`, `RKs` | Courant number and number of IMEX-RK stages |

Changing any of these requires a rebuild (`make` picks up the header dependency
automatically).

## Layout

```
include/
  config.hpp        compile-time switches, physical parameters, global solver state
  types.hpp         Tensor / Grid / Array containers, Eigen matrix-free operator glue
  declarations.hpp  declarations shared by all translation units
src/
  state.cpp                  primitive/conservative conversions, physical fluxes, wave speeds
  flux.cpp                   shock detector and explicit numerical fluxes
  reconstruction.cpp         TVD reconstruction, limiters, Gauss-Legendre quadrature
  initial_conditions.cpp     initial data for every case in the Test enum
  boundary.cpp               boundary conditions for the state and the vector potential
  timestep.cpp               time-step control, conservative update, eigenvalue helpers
  constrained_transport.cpp  staggered / unstaggered CT, vector-potential init, EMF diagnostics
  implicit_operators.cpp     matrix-free subsystem operators, nested Newton, GMRES
  imex_subsystems.cpp        IMEX subsystem assembly and updates, divergence diagnostic
  ct_flux.cpp                EMF-consistent staggered CT from the implicit magnetic fluxes
  output.cpp                 intermediate plotting and VTK output
  main.cpp                   test configuration and the time-marching driver
```

Global solver state lives in `config.hpp` as C++17 `inline` variables, so there
is a single definition shared across translation units without an `extern`
declaration and a separate definition.

## Known issues

- Configuration is compile-time only; a runtime input file would be a better fit
  for parameter studies.
- In `PreSubsystem_Pressure` (`imex_subsystems.cpp`) the arrays `priL`/`priR`
  are filled in but never read. They are marked `[[maybe_unused]]` rather than
  deleted, since the surrounding block looks like an unfinished flux evaluation
  rather than accidental leftovers.
- Several subsystem functions take reconstruction arguments (`w_xL`, `w_xR`,
  `w_yL`, `w_yR`) and a stage index they never use; the names are commented out
  in the definitions to keep the shared signature these functions are called
  through.
- The extensive inline commentary is in Chinese.

## Citing

```bibtex
@article{chen2026imexmhd,
  title   = {A Structure- and Pressure-Positivity-Preserving Semi-implicit IMEX
             Finite Volume Scheme for Ideal MHD at All Acoustic Mach and
             Alfv\'en Mach Numbers with Generic Equation of State},
  author  = {Chen, Zefeng and Dematt\`e, Riccardo and Boscheri, Walter
             and Millmore, Stephen},
  journal = {arXiv preprint arXiv:2608.15837},
  year    = {2026},
  url     = {https://arxiv.org/abs/2608.15837}
}
```
