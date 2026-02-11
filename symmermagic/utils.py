from symmer.operators import PauliwordOp, QuantumState
import numpy as np
import scipy as sp
from itertools import chain, product
import random
from joblib import Parallel, delayed, parallel_config
import warnings
import time

def stab_renyi_entropy(state: QuantumState, order: int=2, filtered : bool = False, approach : str = 'exact', parallel = False, n_proc : int = 4, n_samples : int= 1e6,gpu : bool = False, gpu_device : int = None):
    """Calculates the stabilizer Renyi entropy of the state. See arXiv:2106.12567 for details.

    Args:
        state (QuantumState): the state to calculate the SRE of.
        order (int, optional): the order of SRE to calculate. Default is 2.
        filtered (bool, optional): whether to calculate the filtered stabilizer Renyi entropy by excluding the identity. See arXiv:2312.11631. Default is False.
        approach (str, optional): which calculation approach to use. Valid options are {'exact', 'metropolis'}. Default is 'exact'.
        n_samples (int, optional): the number of samples to use (only applies to 'metropolis'). Default is 1e6.
        parallel (bool or str, optional): parallelization mode. False for serial, True or 'joblib' for
            joblib multiprocessing, 'mpi' for MPI distributed computing (requires mpi4py, launch with
            mpirun/srun). Only applies to 'exact' approach. Default is False.
        n_proc (int, optional): number of processes for joblib parallelization. Ignored for MPI
            (number of ranks is set at launch). Default is 4.
        gpu (bool, optional): whether to use GPU acceleration via CuPy. Can be combined with
            parallel='mpi' for multi-node multi-GPU runs. Default is False.
        gpu_device (int, optional): which GPU device to use. For MPI+GPU, auto-assigned from local
            rank if None. Default is None.

    Returns:
        Mq: the calculated stabilizer Renyi entropy
    """
    Mq=0
    approach=approach.lower()

    if approach=='exact':
        Mq=stab_entropy_symp(state=state,order=order,filtered=filtered,parallel=parallel,n_proc=n_proc,gpu=gpu, gpu_device=gpu_device)
    elif approach=='metropolis':
        state_vec=state.to_sparse_matrix
        Mq=stab_entropy_metropolis(state_vec,order=order,filtered=filtered,n_samples=n_samples)
    else:
        raise ValueError('Unrecognised calculation strategy.')
            
    return Mq

def _fwht(a):
    """Fast Walsh-Hadamard Transform using numpy vectorization.
    Computes f[z] = sum_k a[k] * (-1)^{popcount(k & z)} in O(d log d)."""
    result = np.array(a, dtype=complex)
    h = 1
    while h < len(result):
        result = result.reshape(-1, 2 * h)
        left = result[:, :h].copy()
        right = result[:, h:].copy()
        result[:, :h] = left + right
        result[:, h:] = left - right
        result = result.ravel()
        h *= 2
    return result

def _fwht_batched_gpu(A_gpu):
    """Batched Fast Walsh-Hadamard Transform on GPU using CuPy.
    A_gpu has shape (batch, d). Transforms each row independently."""
    import cupy as cp
    batch, d = A_gpu.shape
    h = 1
    while h < d:
        A_gpu = A_gpu.reshape(batch, -1, 2 * h)
        left = A_gpu[:, :, :h].copy()
        right = A_gpu[:, :, h:].copy()
        A_gpu[:, :, :h] = left + right
        A_gpu[:, :, h:] = left - right
        A_gpu = A_gpu.reshape(batch, -1)
        h *= 2
    return A_gpu

def _build_a_vectors(X_symps_batch, c_states, c_states_set, coeff_dict, conj_dict, d):
    """Build the signal vectors a[s1] = c_{s1} * conj(c_{s1 XOR x}) for a batch of x_symps."""
    A = np.zeros((len(X_symps_batch), d), dtype=np.complex128)
    for j, x_symp in enumerate(X_symps_batch):
        for s1 in c_states:
            s2 = s1 ^ x_symp
            if s2 in c_states_set:
                A[j, s1] = coeff_dict[s1] * conj_dict[s2]
    return A

def stab_entropy_symp(state, order : int = 2, filtered : bool = False, parallel = False, n_proc : int = 4, gpu : bool = False, gpu_device : int = None) -> float:
    """Calculates the exact stabilizer Renyi entropy of the given state by being cheeky in the symplectic representation.
    Args:
        state (QuantumState): the state to calculate the stabilizer entropy for
        order (int): the order of the stabilizer entropy to calculate. default is 2
        filtered (bool): whether to calculate the filtered stabilizer entropy instead of the unfilitered stabilizer entropy. See arXiv:2312.11631 for details. default is False.
        parallel (bool or str, optional): False for serial, True or 'joblib' for joblib multiprocessing,
            'mpi' for MPI distributed computing. Default is False.
        n_proc (int, optional): number of joblib processes. Ignored for MPI. Default is 4.
        gpu (bool, optional): whether to use GPU acceleration via CuPy. Runs batched FWHT on GPU.
            Requires CuPy to be installed and a CUDA-capable GPU available.
            Can be combined with parallel='mpi' for multi-node multi-GPU. Default is False.
        gpu_device (int, optional): which GPU device to use. If None, uses the device set by
            CUDA_VISIBLE_DEVICES or defaults to device 0. For MPI+GPU, auto-assigned from
            local rank if None. Default is None.
    Returns:
        Mq (float): the calculated stabilizer entropy """

    # Normalize parallel parameter for backward compatibility
    if parallel is True:
        _mode = 'joblib'
    elif parallel is False or parallel is None:
        _mode = None
    else:
        _mode = str(parallel).lower()
        if _mode not in ('joblib', 'mpi'):
            raise ValueError(f"Unrecognised parallel mode '{parallel}'. Use False, True, 'joblib', or 'mpi'.")

    # For MPI+GPU, auto-assign GPU device from local MPI rank
    if gpu and _mode == 'mpi' and gpu_device is None:
        gpu_device = _get_mpi_local_rank()

    if gpu:
        try:
            import cupy as cp
        except ImportError:
            raise ImportError("GPU mode requires CuPy. Install with: pip install cupy-cuda12x (or the appropriate CUDA version)")
        try:
            if gpu_device is not None:
                cp.cuda.Device(gpu_device).use()
            dev=cp.cuda.Device()
            dev.compute_capability
        except cp.cuda.runtime.CUDARuntimeError:
            raise RuntimeError(
                f"No CUDA GPU available (are you on a login node?). "
                f"Submit to a GPU node or set CUDA_VISIBLE_DEVICES."
            )
    t1=time.perf_counter()
    n_qubits=state.n_qubits
    d=2**n_qubits
    # build integer-keyed coefficient dict for O(1) lookup without string formatting
    state_dict=state.to_dictionary
    coeff_dict={int(key,2): val for key, val in state_dict.items()}
    conj_dict={k: np.conj(v) for k, v in coeff_dict.items()}
    c_states=list(coeff_dict.keys())
    c_states_set=set(c_states)
    t2=time.perf_counter()
#    print(f'Setup time: {round(t2-t1,6)}s')
    # generate all the X symplectic vectors that could give non-zero contributions
    # sorted() ensures deterministic ordering across MPI ranks
    X_symps=sorted({s1^s2 for s1 in c_states for s2 in c_states})
    t3=time.perf_counter()
 #   print(f'X symps generation time: {round(t3-t2,6)}s')

    if _mode == 'mpi':
        zeta=_symp_mpi(X_symps, c_states, c_states_set, coeff_dict, conj_dict, d, order, gpu=gpu)
    elif gpu:
        zeta=_symp_gpu(cp, X_symps, c_states, c_states_set, coeff_dict, conj_dict, d, order)
    elif _mode == 'joblib':
        def _zeta_for_x(x_symp):
            a=np.zeros(d, dtype=complex)
            for s1 in c_states:
                s2=s1^x_symp
                if s2 in c_states_set:
                    a[s1]=coeff_dict[s1]*conj_dict[s2]
            f=_fwht(a)
            return np.sum(np.abs(f)**(2*order))/d
        batches=max(int(len(X_symps) / n_proc), len(X_symps) // (n_proc * 10))
        with parallel_config(backend='loky'):
            zeta_vals=Parallel(n_jobs=n_proc,batch_size=batches)(delayed(_zeta_for_x)(x) for x in X_symps)
        zeta=sum(zeta_vals)
    else:
        def _zeta_for_x(x_symp):
            a=np.zeros(d, dtype=complex)
            for s1 in c_states:
                s2=s1^x_symp
                if s2 in c_states_set:
                    a[s1]=coeff_dict[s1]*conj_dict[s2]
            f=_fwht(a)
            return np.sum(np.abs(f)**(2*order))/d
        zeta=sum(_zeta_for_x(x) for x in X_symps)
    t4=time.perf_counter()
#    print(f'Zeta calculation time: {round(t4-t3,6)}s')
    if filtered:
        zeta=(zeta-1/d)*d/(d-1)
    Mq=-np.log2(zeta)/(order-1)
    return Mq

def _symp_gpu(cp, X_symps, c_states, c_states_set, coeff_dict, conj_dict, d, order):
    """GPU-accelerated zeta computation. Batches X_symps to fit in GPU memory."""
    num_x=len(X_symps)
    # auto-determine batch size from available GPU memory
    free_mem=cp.cuda.Device().mem_info[0]
    bytes_per_row=d * 16  # complex128
    # use at most 40% of free GPU memory for the working batch (need room for copies in FWHT)
    batch_size=max(1, int(free_mem * 0.4 / (bytes_per_row * 3)))
    batch_size=min(batch_size, num_x)

    zeta=0.0
    for i in range(0, num_x, batch_size):
        batch_x=X_symps[i:i+batch_size]
        # build a-vectors on CPU
        A=_build_a_vectors(batch_x, c_states, c_states_set, coeff_dict, conj_dict, d)
        # transfer to GPU, run batched FWHT, compute contribution
        A_gpu=cp.asarray(A)
        F_gpu=_fwht_batched_gpu(A_gpu)
        zeta+=float(cp.sum(cp.abs(F_gpu)**(2*order)).get())/d
        # free GPU memory between batches
        del A_gpu, F_gpu
        cp.get_default_memory_pool().free_all_blocks()

    return zeta

def _get_mpi_local_rank():
    """Detect MPI local rank from environment variables set by common MPI launchers.
    Used to auto-assign GPU devices in multi-node MPI+GPU runs."""
    import os
    for var in ('OMPI_COMM_WORLD_LOCAL_RANK',   # OpenMPI
                'MV2_COMM_WORLD_LOCAL_RANK',    # MVAPICH2
                'MPI_LOCALRANKID',              # Intel MPI
                'SLURM_LOCALID',                # SLURM
                'LOCAL_RANK'):                   # PyTorch-style / generic
        val = os.environ.get(var)
        if val is not None:
            return int(val)
    # Fallback: use global rank (correct when running on a single node)
    from mpi4py import MPI
    return MPI.COMM_WORLD.Get_rank()

def _symp_mpi(X_symps, c_states, c_states_set, coeff_dict, conj_dict, d, order, gpu=False):
    """MPI-distributed zeta computation. Distributes X_symps across MPI ranks,
    computes partial zeta on each rank (CPU or GPU), and reduces via MPI.SUM.

    All ranks must call this function with the same state data.
    Returns the total zeta on all ranks (broadcast from root)."""
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Distribute X_symps evenly across ranks
    n = len(X_symps)
    chunk = n // size
    remainder = n % size
    if rank < remainder:
        start = rank * (chunk + 1)
        end = start + chunk + 1
    else:
        start = remainder * (chunk + 1) + (rank - remainder) * chunk
        end = start + chunk
    my_X_symps = X_symps[start:end]

    # Compute local zeta contribution
    if len(my_X_symps) == 0:
        local_zeta = 0.0
    elif gpu:
        import cupy as cp
        local_zeta = _symp_gpu(cp, my_X_symps, c_states, c_states_set, coeff_dict, conj_dict, d, order)
    else:
        def _zeta_for_x(x_symp):
            a = np.zeros(d, dtype=complex)
            for s1 in c_states:
                s2 = s1 ^ x_symp
                if s2 in c_states_set:
                    a[s1] = coeff_dict[s1] * conj_dict[s2]
            f = _fwht(a)
            return np.sum(np.abs(f)**(2*order)) / d
        local_zeta = sum(_zeta_for_x(x) for x in my_X_symps)

    # Sum partial zeta across all ranks
    total_zeta = comm.reduce(local_zeta, op=MPI.SUM, root=0)
    # Broadcast so all ranks return the same result
    zeta = comm.bcast(total_zeta, root=0)
    return zeta


def stab_entropy_exact(state_vec, order : int = 2, filtered : bool = False, parallel : bool = False, n_proc : int = 4) -> float:
    """Calculates the exact stabilizer Renyi entropy of the given state by sampling all possible Pauli strings.
    Args: 
        state_vec (csr_matrix): the sparse matrix representation of the state to calculate the stabilizer entropy for
        order (int): the order of the stabilizer entropy to calculate. default is 2
        filtered (bool): whether to calculate the filtered stabilizer entropy instead of the unfilitered stabilizer entropy. See arXiv:2312.11631 for details. default is False.
        parallel (bool, optional): whether to use parallelization approaches. Default is False.
        n_proc (int, optional): if using parallelization, the number of processes to use. Default is 4.
    Returns:
        Mq (float): the calculated stabilizer entropy """
    ###deprecated!
    
    n_qubits=float(np.log2(state_vec.shape[0]))
    assert n_qubits.is_integer(), 'state is wrong shape!'
    state_vec_H=state_vec.getH()
    n_qubits=int(n_qubits)
    d=2**n_qubits
    zeta=0
    pstrings=list(map(lambda plist : ''.join(plist),product(['I','X','Y','Z'],repeat=n_qubits)))

    if parallel:    
        def expval(pstring):
            val=(abs((state_vec_H.dot(pstr_to_sparse(pstring).dot(state_vec))) [0,0])**(2*order))/d
            return val
        
        batches = min(int((d**2)/n_proc),3000)
        with parallel_config(backend='loky'):
            zeta_vals=Parallel(n_jobs=n_proc,return_as='generator_unordered',batch_size=batches)(delayed(expval)(pstring) for pstring in pstrings)
        zeta = accumulator_sum(zeta_vals)

    else:
        for pstring in pstrings:
         sparse_mat=pstr_to_sparse(pstring)
         zeta+=(abs((state_vec_H.dot(sparse_mat.dot(state_vec))) [0,0])**(2*order))/d
    if filtered:
        zeta=(zeta-1/d)*d/(d-1)
    Mq=-np.log2(zeta)/(order-1)
    return Mq

def stab_entropy_metropolis(state_vec, order : int = 2, filtered : bool = False, n_samples : int = 1e6) -> float:
    """Approximates the stabilizer entropy of the given state using a Metropolis-Hastings algorithm. See arXiv:2312.11631 for details.
    Args: 
        state_vec (csr_matrix): the sparse matrix representation of the state to calculate the stabilizer entropy for
        order (int): the order of the stabilizer entropy to calculate. default is 2
        filtered (bool): whether to calculate the filtered stabilizer entropy instead of the unfilitered stabilizer entropy. See arXiv:2312.11631 for details. default is False.
        n_samples (int): the number of samples to use. default is 1e6

    Returns:
        Mq (float): the calculated stabilizer entropy 
    """
    n_qubits=float(np.log2(state_vec.shape[0]))
    assert n_qubits.is_integer(), 'state is wrong shape!'
    state_vec_H=state_vec.getH()
    n_qubits=int(n_qubits)
    d=2**n_qubits

    prob_list=[]
    rng = np.random.default_rng()

    # find starting state which has high enough probability
    while len(prob_list)<1:
        symp_vec=rng.integers(2, size=2*n_qubits)
        # make sure we don't start in the identity state because that will over-sample
        if any(symp_vec):
            pauli_this=PauliwordOp(symp_matrix=symp_vec,coeff_vec=[1])
            sparse_this=pauli_this.to_sparse_matrix
            p_this= abs((state_vec_H.dot(sparse_this.dot(state_vec)))[0,0])**2

            # we want to make sure the state we start in isn't TOO impossible
            if p_this/(d-1) > 1e-8:
                prob_list.append(p_this)

    # initialize the probability list with our initial probability
    prob_list=[p_this]

    # Metropolis-Hastings algorithm
    while len(prob_list)< n_samples:
        # generate hopping positions
        [int1,int2]=rng.integers(2*n_qubits,size=2)
        if int1==int2: # we accidentally drew the same operator twice! we'll forget this happened and try again next time ;)
            continue

        # copy the current PauliworOp symplectic matrix
        symp_next=pauli_this.symp_matrix[0]
        # flip two indicies, equivalent to multiplying by (X/Z)_i and (X/Z)_j  
        symp_next[int1]^=1
        symp_next[int2]^=1

        if not any(symp_next): # this would be taking us to the identity! we'll just forget this happened ;)
            continue
        
        # generate the candidate next PauliwordOp, no operator multiplication required!
        pauli_next=PauliwordOp(symp_matrix=symp_next,coeff_vec=[1])

        # calculate hopping probability
        sparse_next=pauli_next.to_sparse_matrix
        p_next = abs((state_vec_H.dot(sparse_next.dot(state_vec)))[0,0])**2
        hop_prob = p_next/p_this
        # randomly generate an acceptance threshold from 0 to 1
        accept_thresh=random.random()

        # check if the hopping probability is greater than the acceptance threshold
        if hop_prob > accept_thresh:
        # if it is, hop!, if not, we'll sit tight
            pauli_this=pauli_next
            p_this = p_next
        #add p_{k+1} to the list
        prob_list.append(p_this)
    # turn Pauli prob list into zeta
    zeta = sum([p**(order-1) for p in prob_list])/n_samples
    if not filtered:
        zeta=((d-1)*zeta+1)/d
    Mq=-np.log2(zeta)/(order-1)
    return float(Mq)

def stab_linear_entropy(state : QuantumState):
    """Calculates the stabilizer linear entropy of the state. See arXiv:2106.12567 for details.

    Returns:
    Mlin: the calculated stabilizer linear entropy
    """
    zeta=0
    n_qubits=state.n_qubits
    state_vec=state.to_sparse_matrix
    state_vec_H=state.getH()
    pstrings=list(map(lambda plist : ''.join(plist),product(['I','X','Y','Z'],repeat=n_qubits)))
    for pstring in pstrings:
        sparse_mat=pstr_to_sparse(pstring)
        zeta+=(abs(state_vec_H.dot(sparse_mat.dot(state_vec)))[0,0])**4
    Mlin=1-zeta/(2**n_qubits)
    return Mlin

def pauli_spectrum(state : QuantumState) -> float:
    """Calculates the Pauli spectrum of the given state by sampling all possible Pauli strings"""
    state_vec=state.to_sparse_matrix
    n_qubits=float(np.log2(state_vec.shape[0]))
    assert n_qubits.is_integer(), 'state is wrong shape!'
    state_vec_H=state_vec.getH()
    n_qubits=int(n_qubits)
    d=2**n_qubits
    probs=[]
    pstrings=list(map(lambda plist : ''.join(plist),product(['I','X','Y','Z'],repeat=n_qubits)))
    for pstring in pstrings:
        sparse_mat=pstr_to_sparse(pstring)
        probs.append((abs((state_vec_H.dot(sparse_mat.dot(state_vec))) [0,0])**2)/d)
    return probs

def accumulator_sum(generator):
    """Helper method for effectively summing generators"""
    total = 0
    for value in generator:
        total += value
    return total

def pstr_to_sparse(entry: str):
    """Helper method for efficiently converting a Pauli string into a sparse matrix using the PauliComposer algorithm. See https://arxiv.org/abs/2301.00560 for details."""
    BINARY = {'I': '0', 'X': '1', 'Y': '1', 'Z': '0'}
    n = len(entry)
    # Compute some helpful powers
    dim = 1<<n

    # Store the entry converting the Pauli labels into uppercase
    entry = entry.upper()
    #paulis = list(set(entry))

    # (-i)**(0+4m)=1, (-i)**(1+4m)=-i, (-i)**(2+4m)=-1, (-i)**(3+4m)=i
    mat_ent = {0: 1, 1: -1j, 2: -1, 3: 1j}

    # Count the number of ny mod 4
    ny = entry.count('Y') & 3
    init_ent = mat_ent[ny]
    init_entry = init_ent
    iscomplex = np.iscomplex(init_ent)

    # Reverse the input and its 'binary' representation
    rev_entry = entry[::-1]
    rev_bin_entry = ''.join([BINARY[ent] for ent in rev_entry])

    # Column of the first-row non-zero entry
    col_val = int(''.join([BINARY[ent] for ent in entry]), 2)

    # Initialize an empty (2**n x 3)-matrix (rows, columns, entries)
    # row = np.arange(self.dim) 
    col = np.empty(dim, dtype=np.int32)
    # FIXME: storing rows and columns as np.complex64 since NumPy arrays
    # must have the same data type for each entry. Consider using
    # pd.DataFrame?

    col[0] = col_val  # first column
    # The AND bit-operator computes more rapidly mods of 2**n. Check that:
    #    x mod 2**n == x & (2**n-1)
    if iscomplex:
        ent = np.full(dim, init_entry, dtype=np.complex64)
    else:
        ent = np.full(dim, init_entry, dtype=np.int8)

    for ind in range(n):
        p = 1<<int(ind)  # left-shift of bits ('1' (1) << 2 = '100' (4))
        p2 = p<<1
        disp = p if rev_bin_entry[ind] == '0' else -p  # displacements
        col[p:p2] = col[0:p] + disp  # compute new columns
        # col[p:p2] = col[0:p] ^ p  # alternative for computing column

        # Store the new entries using old ones
        if rev_entry[ind] in ['I', 'X']:
            ent[p:p2] = ent[0:p]
        else:
            ent[p:p2] = -ent[0:p]
    row = np.arange(dim)
    return sp.sparse.csr_matrix((ent, (row, col)),
                             shape=(dim, dim))