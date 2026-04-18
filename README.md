# STRESS Solver

A modern implementation of the structural analysis solver for **CIVILSOFT STRESS Plus** `.DAT` files. This project provides a high-performance C implementation with a shared library and a standalone executable.

## Features

- **High-Performance C Engine**: Built with zero external dependencies.
- **Structural Analysis Engine**:
  - Implements Finite Element Method (FEM) for `PLANE GRID` elements.
  - Accounts for both bending ($EI$) and torsional ($GJ$) stiffness.
  - Supports joint releases (hinges) for complex connectivity.
- **Loading Conditions**:
  - Joint concentrated loads.
  - Member uniform and concentrated loads.
- **Clean Architecture**: Decoupled library logic (`libstress.so`) from the console-based UI.

## Project Structure

- `stress.h / .c`: Core structural analysis library in C.
- `main.c`: C command-line wrapper.
- `CMakeLists.txt`: Build configuration for the C project.
- `build.sh`: Script for building the project with strict compiler flags.

## Getting Started

### Prerequisites

- GCC (with support for C11)
- CMake 3.10+

### Building the Project

```bash
# Using the provided build script
./build.sh

# Or manually using CMake
mkdir -p build && cd build
cmake ..
make
```

### Running the Solver

```bash
# Run the solver on a .DAT file
./build/stress path/to/input.dat
```

## Input File Format

The solver parses `.DAT` files used by the legacy CIVILSOFT STRESS Plus software. These files are text-based and define:
- `JOINT COOR`: Node positions (X, Y).
- `MEMB INCI`: Topology (connectivity).
- `MEMBER PROPERTIES`: Sectional data ($A_z, I_y, I_x$).
- `LOADING`: Forces and moments applied to the structure.

## Accuracy

The solver has been verified against original STRESS output. It correctly handles load distribution and structural equilibrium within high numerical precision.

---
*Developed as part of an Advanced Agentic Coding task.*
