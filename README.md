# symmer-magic
Required main package: [Symmer](https://github.com/qmatter-labs/symmer/tree/main).

## Features
1. Calculate stabilizer Renyi entropy for a QuantumState using both exact and numerical sampling methods

## Installation
From the Python environment you want to use for simulations (either Conda or venv) run `pip install git+https://github.com/samalterman/symmer-magic/`

Alternatively, download the package locally, navigate to the symmer-magic root directory and run `pip install .`

## Citation
For Symmer use, cite
> Tim Weaving, Alexis Ralli, Peter J. Love, Sauro Succi, and Peter V. Coveney. *Contextual subspace variational quantum eigensolver calculation of the dissociation curve of molecular nitrogen on a superconducting quantum computer.* [npj Quantum Information, 11(1), 25.](https://www.nature.com/articles/s41534-024-00952-4) (2025).

> Alexis Ralli, Tim Weaving, Andrew Tranter, William M. Kirby, Peter J. Love, and Peter V. Coveney. *Unitary partitioning and the contextual subspace variational quantum eigensolver.* [Phys. Rev. Research 5, 013095](https://doi.org/10.1103/PhysRevResearch.5.013095) (2023).

> Tim Weaving, Alexis Ralli, William M. Kirby, Andrew Tranter, Peter J. Love, and Peter V. Coveney. *A Stabilizer Framework for the Contextual Subspace Variational Quantum Eigensolver and the Noncontextual Projection Ansatz.* [J. Chem. Theory Comput. 2023, 19, 3, 808–821](https://doi.org/10.1021/acs.jctc.2c00910) (2023).

> William M. Kirby, Andrew Tranter, and Peter J. Love, *Contextual Subspace Variational Quantum Eigensolver*, [Quantum 5, 456](https://doi.org/10.22331/q-2021-05-14-456) (2021).

For stabilizer Rényi entropy (SRE), cite
> Lorenzo Leone, Salvatore F. E. Oliviero, and Alioscia Hamma, *Stabilizer Rényi Entropy*, [Phys. Rev. Lett. 128, 050402](https://doi.org/10.1103/PhysRevLett.128.050402) (2022).

For the Metropolis sampling method for calculating SRE, cite
> Xhek Turkeshi, Anatoly Dymarsky, and Piotr Sierant, *Pauli spectrum and nonstabilizerness of typical quantum many-body states*, [Phys. Rev. B 111, 054301](https://doi.org/10.1103/PhysRevB.111.054301) (2025).

For the perfect Paulis sampling method for calculating SRE, cite
> Guglielmo Lami and Mario Collura, *Nonstabilizerness via Perfect Pauli Sampling of Matrix Product States*, [Phys. Rev. Lett. 131, 180401](https://doi.org/10.1103/PhysRevLett.131.180401) (2023).

For the use of the [PauliComposer](https://github.com/sebastianvromero/PauliComposer) algorithm for calculating sparse matrix representations of Pauli strings, cite
>Sebastián Vidal Romero and Juan Santoz-Suárez, *PauliComposer: compute tensor products of Pauli matrices efficiently*, [Quantum Information Processing 22, 449](https://doi.org/10.1007/s11128-023-04204-w) (2023).