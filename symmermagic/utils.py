from symmer.operators import PauliwordOp, QuantumState
import numpy as np
import scipy as sp
from itertools import chain, product
import random

def stab_renyi_entropy(state: QuantumState, order: int=2, filtered : bool = False, sampling : bool = False, sampling_approach : str = 'Metropolis', n_samples : int= 1e6):
    """Calculates the stabilizer Renyi entropy of the state. See arXiv:2106.12567 for details.
    
    Args:
        order (int, optional): the order of SRE to calculate. Default is 2.
        filtered (bool, optional): whether to calculate the filtered stabilizer Renyi entropy by excluding the identity. See arXiv:2312.11631. Default is False.
        sampling (bool, optional): (currently experimental) whether to use a sampling approach. Default is False. 
        sampling_approach (str, optional): which sampling approach to use. Valid options are 'Metropolis'. Default is 'Metropolis'.
        n_samples (int, optional): if using a sampling approach, the number of samples to use. Default is 1e6.
        return_xi_vec (bool, optional): whether to also return the vector of probabilities. Defaullt is False.

    Returns:
        Mq: the calculated stabilizer Renyi entropy
    """
    zeta=0
    n_qubits=state.n_qubits
    d=2**n_qubits
    state_vec=state.to_sparse_matrix
    state_vec_H=state_vec.getH()

    if n_qubits > 12 and not sampling:
        print("Warning: Direct computation for large states may take an extremely long time!")

    if sampling:
        if sampling_approach == 'Metropolis':
            zeta=stab_entropy_metropolis(state_vec,order=order,filtered=filtered,n_samples=n_samples)
        else:
            raise ValueError('Unrecognised approximation strategy.')
    else:
        pstrings=list(map(lambda plist : ''.join(plist),product(['I','X','Y','Z'],repeat=n_qubits)))
        for pstring in pstrings:
            sparse_mat=PauliComposer(pstring).to_sparse()
            zeta+=(abs((state_vec_H.dot(sparse_mat.dot(state_vec))) [0,0])**(2*order))/d
        if filtered:
            zeta=(zeta-1/d)*d/(d-1)
    Mq=-np.log2(zeta)/(order-1)
    return Mq


def stab_entropy_metropolis(state_vec, order : int = 2, filtered : bool = False, n_samples : int = 1e6) -> float:
    """Calculates the stabilizer entropy of the given state using a Metropolis-Hastings algorithm. See arXiv:2312.11631 for details.
    Args: 
        state_vec (csr_matrix): the sparse matrix representation of the state to calculate the stabilizer entropy for
        order (int): the order of the stabilizer entropy to calculate. default is 2
        filtered (bool): whether to calculate the filtered stabilizer entropy instead of the unfilitered stabilizer entropy. See arXiv:2312.11631 for details. default is False.
        n_samples (int): the number of samples to use. default is 1e6

    Returns:
        zeta (float): the calculated stabilizer entropy 
    """
    n_qubits=float(np.log2(state_vec.shape[0]))
    assert n_qubits.is_integer(), 'state is wrong shape!'
    state_vec_H=state_vec.getH()
    n_qubits=int(n_qubits)
    d=2**n_qubits
    pool_range=list(range(2*n_qubits))

    prob_list=[]
    loop=True
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

    return float(zeta)

def stab_linear_entropy(state : QuantumState):
    """Calculates the stabilizer linear entropy of the state. See arXiv:2106.12567 for details.
    
    Args:
    return_xi_vec (bool, optional): whether to also return the vector of probabilities. Default is False.

    Returns:
    Mlin: the calculated stabilizer linear entropy
    """
    zeta=0
    n_qubits=state.n_qubits
    symp_list=[list(chain.from_iterable(ps)) for ps in product([[0,0],[0,1],[1,0],[1,1]],repeat=n_qubits)]
    for symp in symp_list:
        pauli_word=PauliwordOp(symp_matrix=symp,coeff_vec=[1])
        exval=np.real(state.dagger * pauli_word * state)
        zeta +=exval**(4)
    Mlin=1-zeta/(2**n_qubits)
    return Mlin


"""
PauliComposer class definition from https://github.com/sebastianvromero/PauliComposer/

See: https://arxiv.org/abs/2301.00560
"""

import warnings
import numpy as np
import scipy.sparse as ss
from numbers import Number

PAULI_LABELS = ['I', 'X', 'Y', 'Z']
NUM2LABEL = {ind: PAULI_LABELS[ind] for ind in range(len(PAULI_LABELS))}
BINARY = {'I': '0', 'X': '1', 'Y': '1', 'Z': '0'}
PAULI = {'I': np.eye(2, dtype=np.uint8),
         'X': np.array([[0, 1], [1, 0]], dtype=np.uint8),
         'Y': np.array([[0, -1j], [1j, 0]], dtype=np.complex64),
         'Z': np.array([[1, 0], [0, -1]], dtype=np.int8)}


class PauliComposer:

    def __init__(self, entry: str, weight: Number = None):
        # Compute the number of dimensions for the given entry
        n = len(entry)
        self.n = n

        # Compute some helpful powers
        self.dim = 1<<n

        # Store the entry converting the Pauli labels into uppercase
        self.entry = entry.upper()
        self.paulis = list(set(self.entry))

        # (-i)**(0+4m)=1, (-i)**(1+4m)=-i, (-i)**(2+4m)=-1, (-i)**(3+4m)=i
        mat_ent = {0: 1, 1: -1j, 2: -1, 3: 1j}

        # Count the number of ny mod 4
        self.ny = self.entry.count('Y') & 3
        init_ent = mat_ent[self.ny]
        if weight is not None:
            # first non-zero entry
            init_ent *= weight
        self.init_entry = init_ent
        self.iscomplex = np.iscomplex(init_ent)

        # Reverse the input and its 'binary' representation
        rev_entry = self.entry[::-1]
        rev_bin_entry = ''.join([BINARY[ent] for ent in rev_entry])

        # Column of the first-row non-zero entry
        col_val = int(''.join([BINARY[ent] for ent in self.entry]), 2)

        # Initialize an empty (2**n x 3)-matrix (rows, columns, entries)
        # row = np.arange(self.dim) 
        col = np.empty(self.dim, dtype=np.int32)
        # FIXME: storing rows and columns as np.complex64 since NumPy arrays
        # must have the same data type for each entry. Consider using
        # pd.DataFrame?

        col[0] = col_val  # first column
        # The AND bit-operator computes more rapidly mods of 2**n. Check that:
        #    x mod 2**n == x & (2**n-1)
        if weight is not None:
            if self.iscomplex:
                ent = np.full(self.dim, self.init_entry)
            else:
                ent = np.full(self.dim, float(self.init_entry))
        else:
            if self.iscomplex:
                ent = np.full(self.dim, self.init_entry, dtype=np.complex64)
            else:
                ent = np.full(self.dim, self.init_entry, dtype=np.int8)

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

        self.col = col
        self.mat = ent

    def to_sparse(self):
        self.row = np.arange(self.dim)
        return ss.csr_matrix((self.mat, (self.row, self.col)),
                             shape=(self.dim, self.dim))

    def to_matrix(self):
        return self.to_sparse().toarray()