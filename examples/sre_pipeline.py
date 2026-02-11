"""
Pipeline: Hamiltonian JSON -> Exact Diagonalization -> SRE of Ground State

Reads a molecular Hamiltonian in Pauli string JSON format, finds the lowest k
eigenstates, and computes the Stabilizer Renyi Entropy (SRE) of the ground state.

Usage:
    # Local test (4 cores)
    mpirun -n 4 python sre_pipeline.py ../hamiltonian_18_or_less/H2O_STO-3G_SINGLET_JW.json

    # With SLEPc solver
    mpirun -n 8 python sre_pipeline.py hamiltonian.json --solver slepc

    # Save output to file
    mpirun -n 16 python sre_pipeline.py hamiltonian.json --output results.txt

    # Cluster (see submit_pipeline.sh)
    sbatch submit_pipeline.sh
"""

import argparse
import json
import numpy as np
import time
import sys

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from symmer.operators import PauliwordOp, QuantumState
from symmermagic.utils import stab_renyi_entropy


def load_hamiltonian(filepath):
    """Load Hamiltonian from JSON file.

    Expected format:
        {
            "hamiltonian": {"PAULI_STRING": [real, imag], ...},
            "data": {"n_qubits": ..., "calculated_properties": {...}, ...}
        }

    Returns:
        H_sparse: scipy sparse matrix of the Hamiltonian
        n_qubits: number of qubits
        metadata: dict of molecular data
    """
    with open(filepath) as f:
        data = json.load(f)

    ham_dict = data['hamiltonian']
    metadata = data.get('data', {})

    pauli_coeff_dict = {k: complex(v[0], v[1]) for k, v in ham_dict.items()}
    H_op = PauliwordOp.from_dictionary(pauli_coeff_dict)
    H_sparse = H_op.to_sparse_matrix

    n_qubits = len(next(iter(ham_dict)))
    n_terms = len(ham_dict)

    return H_sparse, n_qubits, n_terms, metadata


def eigvec_to_quantumstate(eigvec, n_qubits, threshold=1e-12):
    """Convert dense eigenvector to QuantumState.

    Filters out amplitudes below threshold so the SRE computation
    only works with the non-negligible part of the state.
    """
    nonzero_idx = np.where(np.abs(eigvec) > threshold)[0]
    state_dict = {
        format(idx, f'0{n_qubits}b'): complex(eigvec[idx])
        for idx in nonzero_idx
    }
    return QuantumState.from_dictionary(state_dict)


def diag_scipy(H_sparse, k):
    """Find k lowest eigenpairs using scipy Lanczos (single-node)."""
    from scipy.sparse.linalg import eigsh
    eigenvalues, eigenvectors = eigsh(H_sparse, k=k, which='SA')
    idx = np.argsort(eigenvalues)
    return eigenvalues[idx], eigenvectors[:, idx]


def diag_slepc(H_scipy, k, comm):
    """Find k lowest eigenpairs using SLEPc (MPI-distributed).

    Each rank inserts its local rows into a PETSc matrix, then SLEPc
    solves the eigenvalue problem in parallel. Eigenvectors are gathered
    to all ranks.
    """
    try:
        from petsc4py import PETSc
        from slepc4py import SLEPc
    except ImportError:
        raise ImportError(
            "SLEPc solver requires petsc4py and slepc4py.\n"
            "Install with: conda install -c conda-forge petsc4py slepc4py"
        )

    dim = H_scipy.shape[0]
    H_csr = H_scipy.tocsr()

    # PETSc compiled with real scalars requires real-valued matrix.
    # Molecular Hamiltonians are Hermitian with real matrix elements
    # (imaginary parts from Y-Paulis cancel in the full sum).
    if np.issubdtype(H_csr.dtype, np.complexfloating):
        max_imag = np.max(np.abs(H_csr.data.imag)) if H_csr.nnz > 0 else 0
        if max_imag > 1e-10:
            raise ValueError(
                f"Hamiltonian has non-negligible imaginary entries (max={max_imag:.2e}). "
                f"PETSc is compiled with real scalars. Rebuild PETSc with "
                f"--with-scalar-type=complex, or check your Hamiltonian."
            )
        H_csr = H_csr.real.tocsr()

    A = PETSc.Mat().create(comm)
    A.setSizes([dim, dim])
    A.setType('aij')
    nnz_per_row = np.diff(H_csr.indptr)
    A.setPreallocationNNZ(int(nnz_per_row.max()) if len(nnz_per_row) else 1)
    A.setUp()

    # Each rank fills only its local rows
    rstart, rend = A.getOwnershipRange()
    for i in range(rstart, rend):
        rs, re = H_csr.indptr[i], H_csr.indptr[i + 1]
        if rs < re:
            cols = H_csr.indices[rs:re].astype(PETSc.IntType)
            vals = H_csr.data[rs:re].astype(PETSc.ScalarType)
            A.setValues([i], cols, vals)

    A.assemblyBegin()
    A.assemblyEnd()

    # SLEPc eigensolver
    eps = SLEPc.EPS().create(comm)
    eps.setOperators(A)
    eps.setProblemType(SLEPc.EPS.ProblemType.HEP)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
    eps.setDimensions(nev=k)
    eps.setFromOptions()  
    eps.solve()

    nconv = eps.getConverged()
    if nconv < k and comm.Get_rank() == 0:
        print(f"  WARNING: SLEPc converged {nconv}/{k} eigenpairs")

    eigenvalues = []
    eigenvectors = []
    vr, vi = A.createVecs()

    for i in range(min(nconv, k)):
        eigenvalues.append(eps.getEigenvalue(i).real)
        eps.getEigenvector(i, vr, vi)
        # Gather full vector to all ranks
        scatter, vec_all = PETSc.Scatter.toAll(vr)
        scatter.scatter(vr, vec_all)
        eigenvectors.append(np.array(vec_all))
        scatter.destroy()
        vec_all.destroy()

    vr.destroy()
    vi.destroy()
    eps.destroy()
    A.destroy()

    eigenvalues = np.array(eigenvalues)
    eigenvectors = np.column_stack(eigenvectors) if eigenvectors else np.empty((dim, 0))
    idx = np.argsort(eigenvalues)
    return eigenvalues[idx], eigenvectors[:, idx]



def main():
    parser = argparse.ArgumentParser(
        description='Hamiltonian JSON -> Ground State SRE Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  mpirun -n 4  python sre_pipeline.py ../hamiltonian_18_or_less/H2O_STO-3G_SINGLET_JW.json
  mpirun -n 16 python sre_pipeline.py big_hamiltonian.json --solver slepc --k 8
  srun python sre_pipeline.py hamiltonian.json --output results.txt
        """)
    parser.add_argument('hamiltonian', help='Path to Hamiltonian JSON file')
    parser.add_argument('--solver', choices=['scipy', 'slepc'], default='scipy',
                        help='Eigensolver: scipy (single-node Lanczos) or slepc '
                             '(MPI-distributed). Default: scipy')
    parser.add_argument('--k', type=int, default=5,
                        help='Number of lowest eigenstates to compute. '
                             'Extra states help verify ground state stability. Default: 5')
    parser.add_argument('--order', type=int, default=2,
                        help='SRE Renyi order. Default: 2')
    parser.add_argument('--filtered', action='store_true',
                        help='Compute filtered SRE (excludes identity)')
    parser.add_argument('--output', default=None,
                        help='Save results to this file')
    parser.add_argument('--threshold', type=float, default=1e-15,
                        help='Amplitude threshold for eigenvector coefficients. '
                             'Coefficients below this are zeroed and the state is '
                             'renormalized. Reduces SRE computation time for dense '
                             'eigenvectors (e.g. from SLEPc). Default: 1e-15')

    args = parser.parse_args()

    # Load Hamiltonian 
    if rank == 0:
        print(f"Loading: {args.hamiltonian}")

    H_sparse, n_qubits, n_terms, metadata = load_hamiltonian(args.hamiltonian)

    if rank == 0:
        print(f"  n_qubits  = {n_qubits}")
        print(f"  n_terms   = {n_terms}")
        print(f"  dim       = {2**n_qubits}")
        print(f"  nnz       = {H_sparse.nnz}")
        if metadata:
            print(f"  basis     = {metadata.get('basis', 'N/A')}")
        print()

    #  Diagonalize 
    if rank == 0:
        print(f"Diagonalizing ({args.solver}, k={args.k}) ...")

    t0 = time.perf_counter()

    if args.solver == 'scipy':
        # Only rank 0 diagonalizes, then broadcasts results
        if rank == 0:
            eigenvalues, eigenvectors = diag_scipy(H_sparse, args.k)
        else:
            eigenvalues = None
            eigenvectors = None
        eigenvalues = comm.bcast(eigenvalues, root=0)
        eigenvectors = comm.bcast(eigenvectors, root=0)

    elif args.solver == 'slepc':
        # All ranks participate in distributed diagonalization
        eigenvalues, eigenvectors = diag_slepc(H_sparse, args.k, comm)

    t_diag = time.perf_counter() - t0

    if rank == 0:
        print(f"  Time: {t_diag:.2f}s\n")
        print(f"  Eigenvalues:")
        for i, E in enumerate(eigenvalues):
            gap = f"  (gap = {E - eigenvalues[0]:.8f})" if i > 0 else ""
            marker = " <-- ground state" if i == 0 else ""
            print(f"    E_{i} = {E:18.10f}{gap}{marker}")

        # Compare with reference energies if available
        ref = metadata.get('calculated_properties', {})
        if ref:
            print(f"\n  Reference energies:")
            for method in ['HF', 'MP2', 'CCSD', 'FCI']:
                if method in ref:
                    print(f"    {method:6s} = {ref[method]['energy']:18.10f}")
            if 'FCI' in ref:
                err = abs(eigenvalues[0] - ref['FCI']['energy'])
                print(f"    |E0 - FCI| = {err:.2e}")

        # Degeneracy warning
        if len(eigenvalues) > 1:
            gap = eigenvalues[1] - eigenvalues[0]
            if gap < 1e-6:
                print(f"\n  WARNING: near-degenerate ground state! "
                      f"gap = {gap:.2e}")
                print(f"  The ground state eigenvector may be unreliable.")
        print()

    # Build ground state QuantumState
    if rank == 0:
        gs_vec = eigenvectors[:, 0]

        # Truncate small coefficients and renormalize
        raw_nonzero = np.count_nonzero(np.abs(gs_vec) > 1e-15)
        mask = np.abs(gs_vec) > args.threshold
        discarded_weight = np.sum(np.abs(gs_vec[~mask])**2)
        gs_vec[~mask] = 0.0
        norm = np.linalg.norm(gs_vec)
        if norm > 0:
            gs_vec /= norm
        print(f"Coefficient thresholding (threshold={args.threshold:.1e}):")
        print(f"  Before: {raw_nonzero} non-zero amplitudes")
        print(f"  After:  {np.count_nonzero(mask)} non-zero amplitudes")
        print(f"  Discarded weight: {discarded_weight:.2e}")

        ground_state = eigvec_to_quantumstate(gs_vec, n_qubits, args.threshold)
        n_nonzero = len(ground_state.to_dictionary)
        print(f"Ground state: {n_nonzero} non-zero amplitudes "
              f"({100 * n_nonzero / 2**n_qubits:.1f}% of {2**n_qubits})")
    else:
        ground_state = None
    ground_state = comm.bcast(ground_state, root=0)

    #  Compute SRE with MPI
    if rank == 0:
        print(f"Computing SRE (order={args.order}, {size} MPI ranks) ...")

    comm.Barrier()
    t0 = time.perf_counter()
    sre = stab_renyi_entropy(
        ground_state,
        order=args.order,
        filtered=args.filtered,
        parallel='mpi'
    )
    t_sre = time.perf_counter() - t0

    #  Output
    if rank == 0:
        print(f"  Time: {t_sre:.2f}s\n")
        print("=" * 50)
        print(f"  Ground state energy:  {eigenvalues[0]:.10f}")
        print(f"  SRE (order {args.order}):        {sre:.10f}")
        if len(eigenvalues) > 1:
            print(f"  Energy gap (E1-E0):   {eigenvalues[1] - eigenvalues[0]:.10f}")
        print(f"  Diag time:            {t_diag:.2f}s")
        print(f"  SRE time:             {t_sre:.2f}s")
        print(f"  MPI ranks:            {size}")
        print("=" * 50)

        if args.output:
            with open(args.output, 'w') as f:
                f.write(f"hamiltonian   = {args.hamiltonian}\n")
                f.write(f"n_qubits      = {n_qubits}\n")
                f.write(f"n_terms       = {n_terms}\n")
                f.write(f"solver        = {args.solver}\n")
                f.write(f"mpi_ranks     = {size}\n")
                f.write(f"ground_energy = {eigenvalues[0]}\n")
                if len(eigenvalues) > 1:
                    f.write(f"energy_gap    = {eigenvalues[1] - eigenvalues[0]}\n")
                f.write(f"sre_order     = {args.order}\n")
                f.write(f"sre_filtered  = {args.filtered}\n")
                f.write(f"sre           = {sre}\n")
                f.write(f"diag_time_s   = {t_diag:.4f}\n")
                f.write(f"sre_time_s    = {t_sre:.4f}\n")
                for i, E in enumerate(eigenvalues):
                    f.write(f"E_{i}           = {E}\n")
                ref = metadata.get('calculated_properties', {})
                for method, props in ref.items():
                    f.write(f"ref_{method}      = {props.get('energy', 'N/A')}\n")
            print(f"\nResults saved to {args.output}")


if __name__ == '__main__':
    main()
