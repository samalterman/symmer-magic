"""Test and benchmark MPI parallelization for SRE computation.

Usage:
    mpirun -n 4 python tests/test_mpi.py

Compares serial, joblib, and MPI results for correctness, then benchmarks timing.
"""

from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from symmer import QuantumState
from symmermagic.utils import stab_renyi_entropy

ATOL = 1e-6
N_PROC_JOBLIB = 4


def test_correctness():
    """Verify MPI gives the same result as serial and joblib."""
    if rank == 0:
        print("=" * 60)
        print("CORRECTNESS TEST")
        print("=" * 60)

    for n_qubits in [4, 6, 8]:
        # Rank 0 creates the state and broadcasts to all ranks
        if rank == 0:
            state = QuantumState.haar_random(n_qubits=n_qubits)
        else:
            state = None
        state = comm.bcast(state, root=0)

        # MPI (all ranks participate)
        sre_mpi = stab_renyi_entropy(state, parallel='mpi')

        # Serial and joblib only on rank 0
        if rank == 0:
            sre_serial = stab_renyi_entropy(state, parallel=False)
            sre_joblib = stab_renyi_entropy(state, parallel='joblib', n_proc=N_PROC_JOBLIB)

            match_serial = abs(sre_mpi - sre_serial) < ATOL
            match_joblib = abs(sre_mpi - sre_joblib) < ATOL
            status = "PASS" if (match_serial and match_joblib) else "FAIL"

            print(f"\n  {n_qubits} qubits:")
            print(f"    serial  = {sre_serial:.8f}")
            print(f"    joblib  = {sre_joblib:.8f}")
            print(f"    mpi({size}r) = {sre_mpi:.8f}")
            print(f"    [{status}]")


def test_backward_compat():
    """Verify parallel=True still works (maps to joblib)."""
    if rank == 0:
        print("\n" + "=" * 60)
        print("BACKWARD COMPATIBILITY TEST")
        print("=" * 60)

        state = QuantumState.haar_random(n_qubits=6)
        sre_true = stab_renyi_entropy(state, parallel=True, n_proc=4)
        sre_joblib = stab_renyi_entropy(state, parallel='joblib', n_proc=4)

        match = abs(sre_true - sre_joblib) < ATOL
        status = "PASS" if match else "FAIL"
        print(f"  parallel=True  = {sre_true:.8f}")
        print(f"  parallel='joblib' = {sre_joblib:.8f}")
        print(f"  [{status}]")
    comm.Barrier()


def benchmark():
    """Time serial vs joblib vs MPI at increasing qubit counts."""
    if rank == 0:
        print("\n" + "=" * 60)
        print(f"BENCHMARK  (MPI ranks: {size}, joblib procs: {N_PROC_JOBLIB})")
        print("=" * 60)
        header = f"{'qubits':>6}  {'serial':>10}  {'joblib':>10}  {'mpi':>10}  {'speedup':>8}"
        print(header)
        print("-" * len(header))

    for n_qubits in [6, 8, 10, 12]:
        if rank == 0:
            state = QuantumState.haar_random(n_qubits=n_qubits)
        else:
            state = None
        state = comm.bcast(state, root=0)

        # MPI (all ranks)
        comm.Barrier()
        t0 = time.perf_counter()
        sre_mpi = stab_renyi_entropy(state, parallel='mpi')
        t_mpi = time.perf_counter() - t0

        if rank == 0:
            # Serial
            t0 = time.perf_counter()
            sre_serial = stab_renyi_entropy(state, parallel=False)
            t_serial = time.perf_counter() - t0

            # Joblib
            t0 = time.perf_counter()
            sre_joblib = stab_renyi_entropy(state, parallel='joblib', n_proc=N_PROC_JOBLIB)
            t_joblib = time.perf_counter() - t0

            speedup = t_serial / t_mpi if t_mpi > 0 else float('inf')
            print(f"{n_qubits:>6}  {t_serial:>9.4f}s  {t_joblib:>9.4f}s  {t_mpi:>9.4f}s  {speedup:>7.1f}x")

        comm.Barrier()


if __name__ == "__main__":
    test_correctness()
    test_backward_compat()
    benchmark()

    if rank == 0:
        print("\nDone.")
