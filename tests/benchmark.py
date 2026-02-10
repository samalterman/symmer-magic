"""Benchmark: exact vs symplectic SRE computation for random and sparse states.

Compares: exact, exact parallel, symp, symp parallel, and (optionally) symp GPU.
Produces two log-scale timing plots saved to tests/.
"""

from symmer import QuantumState
from symmermagic.utils import stab_renyi_entropy, stab_entropy_symp
import time
import numpy as np
import matplotlib.pyplot as plt

QUBIT_RANGE = range(2, 11)
K_SPARSE = 20
ORDER = 2
N_PROC = 4
ATOL = 1e-6

# check GPU availability
HAS_GPU = False
try:
    import cupy as cp
    cp.cuda.Device()
    HAS_GPU = True
except Exception:
    pass


def make_sparse_state(n_qubits, k):
    """Create a random state with exactly min(k, 2^n) basis states."""
    k_actual = min(k, 2**n_qubits)
    return QuantumState.random(num_qubits=n_qubits, num_terms=k_actual), k_actual


def bench(fn, state):
    """Time a single SRE computation, return (elapsed_seconds, value)."""
    t0 = time.perf_counter()
    val = fn(state)
    elapsed = time.perf_counter() - t0
    return elapsed, float(np.real(val))


def run_benchmarks():
    methods = {
        'exact':           lambda s: stab_renyi_entropy(s, order=ORDER, approach='exact'),
        'exact_parallel':  lambda s: stab_renyi_entropy(s, order=ORDER, approach='exact', parallel=True, n_proc=N_PROC),
        'symp':            lambda s: stab_entropy_symp(s, order=ORDER),
        'symp_parallel':   lambda s: stab_entropy_symp(s, order=ORDER, parallel=True, n_proc=N_PROC),
    }
    if HAS_GPU:
        methods['symp_gpu'] = lambda s: stab_entropy_symp(s, order=ORDER, gpu=True)

    # ---- collect data ----
    results = {
        'random': {m: {'n': [], 't': []} for m in methods},
        'sparse': {m: {'n': [], 't': []} for m in methods},
    }

    for n in QUBIT_RANGE:
        print(f"\nn = {n} qubits")
        print("-" * 70)

        for state_label, state, k in _make_states(n):
            ref_val = None
            for name, fn in methods.items():
                t, v = bench(fn, state)
                results[state_label][name]['n'].append(n)
                results[state_label][name]['t'].append(t)

                if ref_val is None:
                    ref_val = v
                    status = ""
                else:
                    status = "OK" if abs(v - ref_val) < ATOL else "FAIL"

                tag = f"{state_label}(k={k})" if state_label == "sparse" else state_label
                print(f"  {tag:>14s}  {name:>16s}  {t:10.4f}s  SRE={v:.6f}  {status}")

    # ---- plot ----
    _make_plots(methods, results)


def _make_states(n):
    """Yield (label, state, k) tuples for a given qubit count."""
    dense = QuantumState.haar_random(n_qubits=n)
    yield ("random", dense, 2**n)

    sparse, k_actual = make_sparse_state(n, K_SPARSE)
    yield ("sparse", sparse, k_actual)


def _make_plots(methods, results):
    style = {
        'exact':          dict(marker='o', ls='-',  color='#e74c3c', label='Exact'),
        'exact_parallel': dict(marker='s', ls='--', color='#e67e22', label=f'Exact parallel ({N_PROC} cores)'),
        'symp':           dict(marker='^', ls='-',  color='#2980b9', label='Symplectic + FWHT'),
        'symp_parallel':  dict(marker='D', ls='--', color='#27ae60', label=f'Symplectic + FWHT parallel ({N_PROC} cores)'),
        'symp_gpu':       dict(marker='*', ls='-',  color='#8e44ad', label='Symplectic + FWHT (GPU)', markersize=10),
    }

    for state_label, title_extra in [('random', 'Haar Random (Dense)'),
                                     ('sparse', f'Sparse (k={K_SPARSE})')]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for name in methods:
            ns = results[state_label][name]['n']
            ts = results[state_label][name]['t']
            s = {k: v for k, v in style[name].items()}
            ms = s.pop('markersize', 7)
            ax.plot(ns, ts, **s, markersize=ms, linewidth=2)

        ax.set_yscale('log')
        ax.set_xlabel('Number of qubits', fontsize=13)
        ax.set_ylabel('Time (s)', fontsize=13)
        ax.set_title(f'SRE Computation Time — {title_extra}', fontsize=14)
        ax.set_xticks(list(QUBIT_RANGE))
        ax.legend(fontsize=10)
        ax.grid(True, which='both', alpha=0.3)
        fig.tight_layout()

        fname = f"tests/benchmark_{state_label}.png"
        fig.savefig(fname, dpi=150)
        print(f"\nSaved {fname}")
        plt.close(fig)


if __name__ == "__main__":
    gpu_status = "available" if HAS_GPU else "not available (CuPy not found)"
    print(f"GPU: {gpu_status}")
    print(f"Methods: exact, exact_parallel, symp, symp_parallel" + (", symp_gpu" if HAS_GPU else ""))
    run_benchmarks()
