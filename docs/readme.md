# Advanced Computational Framework for Worldvolume Hybrid Monte Carlo Methods

This project implements an advanced computational framework based on the **Worldvolume Hybrid Monte Carlo (WV-HMC)** method, as outlined in the paper **"Simplified Algorithm for the Worldvolume HMC and the Generalized-Thimble HMC"** (arXiv:2311.10663). This framework addresses the challenges posed by the sign problem in complex systems by utilizing an efficient algorithm that avoids the need to compute the Jacobian during configuration generation. The framework also incorporates improvements for cost-effective and reliable molecular dynamics steps on constrained manifolds.

## Key Features

- **Simplified RATTLE Algorithm**: Efficiently projects configurations onto constrained manifolds.
- **Fixed-Point Methods**: Streamlines iterative solutions, reducing computational overhead.
- **WV-HMC**: Effective in handling sign and ergodicity problems in path integrals.
- **DOP853 ODE Solver**: Utilizes the explicit 7th-order Runge-Kutta method with advanced dense output capabilities, developed by Ernst Hairer, for solving ordinary differential equations, ensuring exceptional precision and stability in molecular dynamics simulations.

## Computational Steps

The program consists of three main computational phases:

1. **Define the Coefficients in W(t)**:

   - Implements updates to the coefficients of the W(t) function.
   - Handled by the `update_bw.f95` module.

2. **Generate Markov Chain Configurations**:

   - Produces Markov chain samples required for the Monte Carlo process.
   - Implemented in `generate_markov_chain.f95`.

3. **Evaluate Expectation Values**:

   - Computes observables and their statistical properties.
   - Conducted using `virial.f95`.

4. **Run Test Program**:

   - Executes a standalone test program for debugging and validation.
   - Implemented in `test.f95`.

## Project Structure

```plaintext
project/
│
├── src/                  # Source codes
│   ├── mt95.f95          # Core math routines
│   ├── parameters.f95    # Parameter management
│   ├── hmc.f95           # HMC implementation
│   ├── dop853_module.f95 # ODE solver (DOP853 algorithm by Ernst Hairer)
│   ├── markovchain.f95   # Markov chain sampling
│   ├── model.f95         # Model-specific definitions
│   ├── virial.f95        # Expectation value calculations
│   ├── update_bw.f95     # Updating W(t) coefficients
│   ├── generate_markov_chain.f95 # Markov chain generation
│   ├── util.f95          # Utility routines
│
├── data/                 # Input and configuration data
│   ├── bw_parameters.dat # Parameters for W(t)
│   ├── initial_x.dat     # Initial state data
│
├── docs/                 # Documentation and references
│   ├── README.md         # Project README file
│   ├── LICENSE           # Project license file
│
├── build/                # Build-related files
│   ├── makefile          # Build and compile the project
│
├── output/               # Outputs and results
│   └── (generated output files will go here)
├── tests/                # Test programs
    ├── test.f95          # Test program for Hamiltonian conservation
    └── test2.f95         # Test program for Action Derivatives
```

## Prerequisites

1. **Fortran Compiler**: Ensure a Fortran 95 compatible compiler is installed (e.g., `gfortran`).
2. **Mathematical Libraries**: Include BLAS and LAPACK for numerical computations.
3. **Make Utility**: Build the project with the provided `makefile`.

## Setup and Compilation

1. Clone or copy the repository.
2. Navigate to the `build/` directory.
3. Execute the command:
   ```bash
   make
   ```
   This compiles the source code and generates executables in the `build/` directory.

## Running the Program

Run the program with the menu-driven interface:

```bash
make mn
```

This presents options to execute:

1. **Define Coefficients**:

   - Updates W(t) coefficients using `define_coefficients`.
   - Executes the logic in `update_bw.f95`.

2. **Generate Markov Chain**:

   - Creates configurations with `generate_markov_chain`.
   - Relies on `generate_markov_chain.f95`.

3. **Evaluate Expectations**:

   - Calculates observables via `evaluate_expectations`.
   - Utilizes `virial.f95`.

4. **Run Test Program**:

   - Executes the standalone test program `test_program` for validation.
   - Uses `test.f95` for Hamiltonian conservation test.

5. **Run Test Program**:

   - Executes the standalone test program `test_program2` for validation.
   - Uses `test2.f95` for Action Derivatives test.

6. **Exit**: Ends the interface.

Clean build files with:

```bash
make clean
```

Remove all generated files:

```bash
make cleanall
```

## Outputs

The `output/` directory contains:

1. Generated Markov chain configurations.
2. Computed expectation values.
3. Jackknife error analysis.

## Customization

1. **Parameters**: Modify `parameters.f95` to adjust model settings.
2. **Model Definitions**: Customize or extend models in `model.f95`.
3. **Input Data**: Update `bw_parameters.dat` and `initial_x.dat` in the `data/` folder.

## References

- **Paper**: [Simplified Algorithm for the Worldvolume HMC and the Generalized-Thimble HMC](https://arxiv.org/abs/2311.10663v4)
- **Authors**: Masafumi Fukuma
  - **Affiliation**: Department of Physics, Kyoto University, Kyoto 606-8502, Japan.
  - **Email**: [fukuma@gauge.scphys.kyoto-u.ac.jp](mailto\:fukuma@gauge.scphys.kyoto-u.ac.jp).
- **Published Version**: To be published in a forthcoming volume of Progress of Theoretical and Experimental Physics, April 2024.
  - **DOI**: [10.1093/ptep/ptac010](https://doi.org/10.1093/ptep/ptac010).

## License

This project is distributed under the [MIT License](LICENSE).

## Contact

For any questions or issues, please contact:

- **Name**: Chien-Yu Chou
- **Affiliation**: The Graduate University for Advanced Studies (SOKENDAI), Japan
- **Email**: [ccy@post.kek.jp](mailto\:ccy@post.kek.jp)

