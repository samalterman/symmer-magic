"""Benchmark: hybrid Monte Carlo vs exact symplectic SRE.

Tests on 13-qubit Haar random state where |X_symps| = 8192.
Shows accuracy vs speed tradeoff for different sample counts.
"""

from symmer import QuantumState
from symmermagic.utils import stab_entropy_symp, stab_entropy_symp_hybrid
import time
import numpy as np
import matplotlib.pyplot as plt

N_QUBITS = 13
N_TRIALS = 5  # repeat hybrid to measure variance
SAMPLE_COUNTS = [50, 100, 200, 500, 1000, 2000, 4000, 8192]


def run():
    print(f"Generating {N_QUBITS}-qubit Haar random state...")
    state = QuantumState.haar_random(n_qubits=N_QUBITS)
    d = 2**N_QUBITS
    print(f"d = {d}, |X_symps| = {d} (dense Haar state)\n")

    # exact baseline
    print("Running exact symp...")
    t0 = time.perf_counter()
    v_exact = stab_entropy_symp(state)
    t_exact = time.perf_counter() - t0
    print(f"  exact:  SRE = {v_exact:.6f}  time = {t_exact:.2f}s\n")

    # hybrid at various sample counts
    results = []
    header = f"{'n_samples':>10}  {'mean SRE':>10}  {'std SRE':>10}  {'rel err':>10}  {'mean time':>10}  {'speedup':>8}  {'% X_symps':>10}"
    print(header)
    print("-" * len(header))

    for ns in SAMPLE_COUNTS:
        times = []
        values = []
        for _ in range(N_TRIALS):
            t0 = time.perf_counter()
            v = stab_entropy_symp_hybrid(state, n_samples=ns)
            t = time.perf_counter() - t0
            times.append(t)
            values.append(v)

        mean_v = np.mean(values)
        std_v = np.std(values)
        mean_t = np.mean(times)
        rel_err = abs(mean_v - v_exact) / abs(v_exact)
        speedup = t_exact / mean_t
        pct = 100 * ns / d

        results.append({
            'ns': ns, 'mean_v': mean_v, 'std_v': std_v,
            'rel_err': rel_err, 'mean_t': mean_t, 'speedup': speedup,
            'pct': pct, 'values': values, 'times': times,
        })

        print(f"{ns:>10}  {mean_v:>10.4f}  {std_v:>10.4f}  {rel_err:>10.4%}  {mean_t:>10.4f}s  {speedup:>7.1f}x  {pct:>9.1f}%")

    # ---- plots ----
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ns_arr = [r['ns'] for r in results]

    # 1) time vs n_samples
    ax = axes[0]
    ax.plot(ns_arr, [r['mean_t'] for r in results], 'o-', color='#2980b9', linewidth=2)
    ax.axhline(t_exact, color='#e74c3c', ls='--', linewidth=2, label=f'Exact ({t_exact:.1f}s)')
    ax.set_xlabel('n_samples', fontsize=12)
    ax.set_ylabel('Time (s)', fontsize=12)
    ax.set_title('Speed', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # 2) relative error vs n_samples
    ax = axes[1]
    ax.plot(ns_arr, [r['rel_err'] * 100 for r in results], 's-', color='#e67e22', linewidth=2)
    ax.set_xlabel('n_samples', fontsize=12)
    ax.set_ylabel('Relative error (%)', fontsize=12)
    ax.set_title('Accuracy', fontsize=13)
    ax.grid(True, alpha=0.3)

    # 3) SRE distribution (box plot)
    ax = axes[2]
    positions = range(len(results))
    bp = ax.boxplot([r['values'] for r in results], positions=positions, widths=0.6, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('#2ecc71')
        patch.set_alpha(0.7)
    ax.axhline(v_exact, color='#e74c3c', ls='--', linewidth=2, label=f'Exact = {v_exact:.4f}')
    ax.set_xticks(positions)
    ax.set_xticklabels([str(r['ns']) for r in results], rotation=45)
    ax.set_xlabel('n_samples', fontsize=12)
    ax.set_ylabel('SRE', fontsize=12)
    ax.set_title(f'Distribution ({N_TRIALS} trials)', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f'Hybrid Monte Carlo — {N_QUBITS}-qubit Haar Random State', fontsize=14, y=1.02)
    fig.tight_layout()
    fname = "tests/benchmark_hybrid.png"
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"\nSaved {fname}")
    plt.close(fig)


if __name__ == "__main__":
    run()
