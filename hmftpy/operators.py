from quspin.operators import quantum_operator
import numpy as np

TYPE = np.complex128

def inner_hamiltonian(plaquette, interactions, basis, verbose=False, disorder={}):
    """
    Constructs a QuSpin quantum_operator corresponding to the inner-cluster
        Hamiltonian (i.e. the cluster Hamiltonian with OBC)
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        disorder (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Hi: QuSpin quantum_operator
    """
    bonds = []
    L = plaquette['L']
    for n in interactions:
        if n == 'local':
            for c_op in interactions[n]:
                r = np.ones(L)
                if n in disorder:
                    if c_op in disorder[n]:
                        r = disorder[n][c_op]
                couplings = interactions[n][c_op]*r
                bonds += [[c_op, [[couplings[i], i] for i in range(L)]]]
        else:
            for c_op in interactions[n]:
                r = np.ones((L, L)) # resetting
                if n in disorder:
                    if c_op in disorder[n]:
                        r = disorder[n][c_op]
                for i in range(L):
                    coupling = interactions[n][c_op]
                    bonds += [[c_op, [[coupling*r[i, ni], i, ni] for ni in plaquette['inner'][n][i]]]]
    Hi = quantum_operator({'static': bonds}, basis=basis,
                          check_herm=False, check_symm=False,
                          dtype=TYPE)
    return Hi


def periodic_hamiltonian(plaquette, interactions, basis, verbose=False, disorder={}):
    """
    Constructs a QuSpin quantum_operator corresponding to the cluster
        Hamiltonian with PBC
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        disorder (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Hp: QuSpin quantum_operator
    """
    di = 0
    ld = len(disorder)
    bonds = []
    neighbors = ['nearest', 'n_nearest', 'n_n_nearest']
    L = plaquette['L']
    for n in interactions:
        if n == 'local':
            for c_op in interactions[n]:
                r = np.ones(L)
                if n in disorder:
                    if c_op in disorder[n]:
                        r = disorder[n][c_op]
                couplings = interactions[n][c_op]*r
                bonds += [[c_op, [[couplings[i], i] for i in range(L)]]]
        else:
            for c_op in interactions[n]:
                r = np.ones((L, L)) # resetting
                if n in disorder:
                    if c_op in disorder[n]:
                        r = disorder[n][c_op]
                for i in range(L):
                    coupling = interactions[n][c_op]
                    bonds += [[c_op, [[coupling*r[i,ni], i, ni] for ni in plaquette['inner'][n][i]]]]
                    bonds += [[c_op, [[coupling*r[i,ni], i, ni] for ni in plaquette['outer'][n][i]]]]
    Hp = quantum_operator({'static': bonds}, basis=basis,
                          check_herm=False, check_symm=False,
                          dtype=TYPE)
    return Hp


def outer_hamiltonian(plaquette, mean_fields, interactions, basis, verbose=False, disorder={}):
    """
    Constructs a QuSpin quantum_operator corresponding to the cluster-bath
        Hamiltonian
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        mean_fields (dict): dict of mean fields (see example in README.md)
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        disorder (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Ho: QuSpin quantum_operator
    """
    L = plaquette['L']
    bonds = []
    for n in interactions:
        for c_op in interactions[n]:
            if n != 'local':
                r = np.ones((L, L))
                if n in disorder:
                    if c_op in disorder[n]:
                        r = disorder[n][c_op]
                for i in range(L):
                    coupling = interactions[n][c_op]
                    bonds += [[c_op[1], [[coupling*mean_fields[c_op[0]][ni]*r[i,ni], i] for ni in plaquette['outer'][n][i]]]]
    Ho = quantum_operator({'static': bonds}, basis=basis, check_herm=False,
                          check_symm=False, dtype=TYPE)
    return Ho


def mf_ops(plaquette, basis):
    """
    Constructs a dict of QuSpin quantum operators for all components of spin on
        each site of the cluster
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
    Output:
        ops (list): Format is [{'x': x_op, 'y': y_op, 'z': z_op}] with one
            dict for each site in the cluster
    """
    ops = [{} for i in range(plaquette['L'])]
    for i in range(plaquette['L']):
        ops[i]['x'] = quantum_operator({'static': [['x', [[1.0, i]]]]},
                                           basis=basis, check_herm=False,
                                           check_symm=False, dtype=TYPE)
        ops[i]['y'] = quantum_operator({'static': [['y', [[1.0, i]]]]},
                                           basis=basis, check_herm=False,
                                           check_symm=False, dtype=TYPE)
        ops[i]['z']= quantum_operator({'static': [['z', [[1.0, i]]]]},
                                          basis=basis, check_herm=False,
                                          check_symm=False, dtype=TYPE)
    return ops



def n_nearest_magnetization(sublattices, basis, dir=(1,0,0)):
    ops = []
    for subl in sublattices:
        op_lst = [['x', [[dir[0], i] for i in subl]],
                  ['y', [[dir[1], i] for i in subl]],
                  ['z', [[dir[2], i] for i in subl]]]
        ops += [quantum_LinearOperator(op_lst, basis=basis)]
    return ops

def n_n_nearest_magnetization(sublattices, basis, dir=(1,0,0)):
    ops = []
    for subl in sublattices:
        op_lst = [['x', [[dir[0], i] for i in subl]],
                  ['y', [[dir[1], i] for i in subl]],
                  ['z', [[dir[2], i] for i in subl]]]
        ops += [quantum_LinearOperator(op_lst, basis=basis)]
    return ops


def stripe_magnetization(stripes, basis, dir=(1,0,0)):
    ops = []
    for st in stripes:
        op_lst = [['x', [[dir[0], i] for i in st]],
                  ['y', [[dir[1], i] for i in st]],
                  ['z', [[dir[2], i] for i in st]]]
        ops += [quantum_LinearOperator(op_lst, basis=basis)]
    return ops


def helicity(triangles, basis):
    # specific to 12 site for now
    # doing z projection for now
    ops = []
    for tr in triangles:
        op_lst = [['xy', [[1, tr[0], tr[1]], [-1, tr[1], tr[0]], # s0 x s1
                          [1, tr[1], tr[2]], [-1, tr[2], tr[1]], # s1 x s2,
                          [1, tr[2], tr[0]], [-1, tr[0], tr[2]]]]] # s2 x s3
        ops += [quantum_LinearOperator(op_lst, basis=basis)]
    return ops
