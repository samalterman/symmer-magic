"""
Template script for computing Stabilizer Renyi Entropy (SRE) on an HPC cluster with MPI.

"""

from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# SRE parameters
ORDER = 2               # Renyi entropy order 
FILTERED = False         # True for filtered SRE 
METHOD = "haar_random"   # options: "haar_random", "sparse_random", "from_dict"

# --- Option 1: Haar random state (dense, for testing/benchmarking) ---
N_QUBITS = 14

# --- Option 2: Sparse random state ---
# N_QUBITS = 16
# N_TERMS = 50  # number of non-zero computational basis states

# --- Option 3: Build from a dictionary {bitstring: amplitude} ---
# This is how you'd use a state from VQE, QSCI, state prep, etc.
# Example:
# STATE_DICT = {
#     "0000": 0.5,
#     "0011": 0.5,
#     "1100": 0.5,
#     "1111": 0.5,
# }

# Output file (written by rank 0 only)
OUTPUT_FILE = "sre_result.txt"


from symmer import QuantumState
from symmermagic.utils import stab_renyi_entropy


def build_state():
    """Build the quantum state on rank 0. Returns None on other ranks."""
    if METHOD == "haar_random":
        return QuantumState.haar_random(n_qubits=N_QUBITS)

    elif METHOD == "sparse_random":
        return QuantumState.random(num_qubits=N_QUBITS, num_terms=N_TERMS)

    elif METHOD == "from_dict":
        amplitudes = np.array(list(STATE_DICT.values()), dtype=complex)
        norm = np.linalg.norm(amplitudes)
        normalized = {k: v / norm for k, v in STATE_DICT.items()}
        return QuantumState(normalized)

    else:
        raise ValueError(f"Unknown method: {METHOD}")


def main():
    # Rank 0 builds the state, then broadcasts to all ranks #
    if rank == 0:
        print(f"Building state (method={METHOD}) ...")
        state = build_state()
        print(f"  n_qubits = {state.n_qubits}")
        print(f"  n_terms  = {len(state.to_dictionary)}")
        print(f"  MPI ranks = {size}")
        print(f"  order    = {ORDER}, filtered = {FILTERED}")
        print()
    else:
        state = None
    state = comm.bcast(state, root=0)

    # Compute SRE with MPI parallelization #
    comm.Barrier()
    t0 = time.perf_counter()
    sre = stab_renyi_entropy(state, order=ORDER, filtered=FILTERED, parallel='mpi')
    elapsed = time.perf_counter() - t0

    # Rank 0 prints and saves the result #
    if rank == 0:
        print(f"SRE (order {ORDER}) = {sre:.10f}")
        print(f"Time: {elapsed:.2f}s  ({size} MPI ranks)")

        with open(OUTPUT_FILE, "w") as f:
            f.write(f"method     = {METHOD}\n")
            f.write(f"n_qubits   = {state.n_qubits}\n")
            f.write(f"n_terms    = {len(state.to_dictionary)}\n")
            f.write(f"order      = {ORDER}\n")
            f.write(f"filtered   = {FILTERED}\n")
            f.write(f"sre        = {sre}\n")
            f.write(f"time_sec   = {elapsed:.4f}\n")
            f.write(f"mpi_ranks  = {size}\n")
        print(f"Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
