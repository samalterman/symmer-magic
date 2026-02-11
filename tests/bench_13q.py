"""Quick benchmark: 13-qubit Haar random state with symp method."""

from symmer import QuantumState
from symmermagic.utils import stab_entropy_symp
import time

s = QuantumState.haar_random(n_qubits=13)
print("State generated (13 qubits, dense)")

t0 = time.perf_counter()
v1 = stab_entropy_symp(s)
t1 = time.perf_counter()
print(f"symp serial:       SRE={v1:.6f}  time={t1-t0:.2f}s")

t0 = time.perf_counter()
v2 = stab_entropy_symp(s, parallel=True, n_proc=4)
t1 = time.perf_counter()
print(f"symp parallel (4): SRE={v2:.6f}  time={t1-t0:.2f}s")

print(f"match: {abs(v1-v2) < 1e-6}")
