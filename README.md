# Molecular Dynamics Code (Fortran 90)

## Overview

This project implements a modular **Molecular Dynamics (MD)** code written in **Fortran 90**, with a focus on clarity, numerical transparency, and explicit separation of physical and technical concerns.

The code follows a classical MD workflow:
- reading and validating input parameters,
- computing forces and energies,
- integrating equations of motion using the Velocity Verlet algorithm,
- writing trajectories and diagnostic output.

The repository also includes a test suite, a Makefile-based build system, and automatically generated documentation using **Doxygen**.

---

## Repository structure

| Path / File        | Description |
|--------------------|-------------|
| `src/`             | Core Fortran 90 source code implementing the MD engine. |
| `tests/`           | Independent test programs for individual modules. |
| `build/`           | Automatically generated build artifacts. |
| `build/mod/`       | Compiled Fortran module files (`.mod`). |
| `build/obj/`       | Object files (`.o`) generated during compilation. |
| `doc/`             | Project documentation generated with Doxygen. |
| `doc/html/`        | HTML documentation (open `index.html`). |
| `Makefile`         | Build, test, and documentation automation rules. |
| `inp.txt`          | Input file defining simulation parameters. |
| `trajectory.xyz`   | Output trajectory in XYZ format. |
| `MD.out`           | Standard output from the MD run. |
| `program.exe`      | Compiled executable. |
| `doxygen.txt`      | Doxygen configuration file. |
| `README.md`        | Project overview and Doxygen main page. |

---

## Source code ("src")

The `src/` directory contains the full implementation of the MD engine, written in a modular Fortran 90 style.

### Main program

- **`main.f90`**  
  Entry point of the program. Controls the global MD workflow:
  input handling, force evaluation, time integration, and output.

### Physical and numerical modules

- **`acceleration_mod.f90`**  
  Computes forces and particle accelerations.

- **`energies.f90`**  
  Evaluates kinetic, potential, and total energies.

- **`verlet_al.f90`**  
  Implements the Velocity Verlet integration algorithm.

### Input/output and robustness

- **`treat_input.f90`**  
  Reads and validates simulation parameters.

- **`treat_output.f90`**  
  Handles output generation (energies, trajectories, logs).

- **`error_treat.f90`**  
  Centralized error handling and sanity checks.

Each module exposes explicit interfaces and is designed to be testable in isolation.

---

## Test suite ("tests")

The `tests/` directory contains independent Fortran programs to validate individual modules:

- `test_acceleration_mod.f90`
- `test_acceleration_mod_2.f90`
- `test_energies.f90`
- `test_verlet_al.f90`
- `test_treat_input.f90`
- `test_treat_output.f90`
- `test_error_treat.f90`

These tests are intended to ensure numerical correctness and stable behavior of the core components.

---

## Build artifacts ("build")

This directory is generated automatically during compilation.

- **`build/mod/`**  
  Compiled Fortran module files (`.mod`).

- **`build/obj/`**  
  Object files (`.o`) corresponding to each source file.

This directory should not be modified manually.

---

## Documentation ("doc")

HTML documentation is generated using **Doxygen**.

- **`doc/html/`**  
  Contains browsable documentation including:
  - source code listings,
  - module and namespace descriptions,
  - dependency graphs,
  - documentation for test files.

Open `doc/html/index.html` to explore the documentation.

---

## Build system

The project uses a `Makefile` (written using Bash-style rules) to automate:

- compilation and linking,
- test execution,
- documentation generation.

| Command      | Description |
|-------------|-------------|
| `make`      | Compiles the source code and links the executable. |
| `make test` | Builds and runs the module-level test suite. |
| `make help` | Tells the user the different commands associated to make. |
| `make clean` | Clean the files after compilation. |
