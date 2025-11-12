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

    # still experimental so don't really trust this
    if sampling:
        if sampling_approach == 'Metropolis':
            zeta=stab_entropy_metropolis(state_vec,order=order,filtered=filtered,n_samples=n_samples)
        else:
            raise ValueError('Unrecognised approximation strategy.')
    else:
        symps=product((False,True),repeat=2*n_qubits) #generate all the possible symplectic matrices
        sparses=map(lambda symp : PauliwordOp(np.array([symp]),coeff_vec=[1]).to_sparse_matrix, symps)
        # now we go through all of the possible symplectic matrices
        for sparse_mat in sparses:
            #sparse_mat=PauliwordOp(np.array([symp]),coeff_vec=[1]).to_sparse_matrix
            zeta+=(abs((state_vec_H.dot(sparse_mat.dot(state_vec))) [0,0])**(2*order))/d

        #sparse_list=[PauliwordOp(symp_matrix=symp,coeff_vec=[1]).to_sparse_matrix for symp in symp_list]
        #for sparse_mat in sparse_mats:
        #    prob=abs((state_vec_H.dot(sparse_mat.dot(state_vec)))[0,0])**2
        #    zeta +=(prob**order)/d
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