from quspin.operators import quantum_operator
import numpy as np

TYPE = np.complex128

def get_r(L, coeffs, n, c_op):
    """
    Constructes the coefficient array for a given interaction type "c_op" connecting "n"-type neighbors.
    """
    if n in coeffs:
        if c_op in coeffs[n]:
            return coeffs[n][c_op]
    elif n == 'local':
        return np.ones(L)
    else:
        return np.ones((L, L))

def inner_hamiltonian(plaquette, interactions, basis, verbose=False, coeffs={'inner': {}, 'outer': {}},
                      checks=False):
    """
    Constructs a QuSpin quantum_operator corresponding to the inner-cluster
        Hamiltonian (i.e. the cluster Hamiltonian with OBC)
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        coeffs (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Hi: QuSpin quantum_operator
    """
    terms = []
    L = plaquette['L']
    bonds = ['x_bonds', 'y_bonds', 'z_bonds', 'n_bonds', 'nn_bonds', 'a_bonds', 'b_bonds']
    for n in interactions:
        if n == 'local':
            for c_op in interactions[n]:
                r = get_r(L, coeffs['inner'], n, c_op)
                couplings = interactions[n][c_op]*r
                terms += [[c_op, [[couplings[i], i] for i in range(L)]]]
        elif n in bonds:
            for c_op in interactions[n]:
                coupling = interactions[n][c_op]
                terms += [[c_op, [[coupling, b[0], b[1]] for b in plaquette['inner'][n]]]]
        else:
            raise Exception('Unknown interaction type {}'.format(n))
    Hi = quantum_operator({'static': terms}, basis=basis,
                          check_herm=checks, check_symm=checks, check_pcon=checks,
                          dtype=TYPE)
    return Hi


def periodic_hamiltonian(plaquette, interactions, basis, verbose=False,
                         coeffs={'inner': {}, 'outer': {}},
                         checks=False):
    """
    Constructs a QuSpin quantum_operator corresponding to the cluster
        Hamiltonian with PBC
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        coeffs (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Hp: QuSpin quantum_operator
    """
    di = 0
    terms = []
    L = plaquette['L']
    bonds = ['x_bonds', 'y_bonds', 'z_bonds', 'n_bonds', 'nn_bonds', 'a_bonds', 'b_bonds']
    for n in interactions:
        if n == 'local':
            for c_op in interactions[n]:
                r = get_r(L, coeffs['inner'], n, c_op)
                couplings = interactions[n][c_op]*r
                terms += [[c_op, [[couplings[i], i] for i in range(L)]]]
        elif n in bonds:
            for c_op in interactions[n]:
                coupling = interactions[n][c_op]
                terms += [[c_op, [[coupling, b[0], b[1]] for b in plaquette['inner'][n]]]]
                terms += [[c_op, [[coupling, b[0], b[1]] for b in plaquette['outer'][n]]]]
        else:
            raise Exception('Unknown interaction type {}'.format(n))
    Hp = quantum_operator({'static': terms}, basis=basis,
                          check_herm=checks, check_symm=checks, check_pcon=checks,
                          dtype=TYPE)
    return Hp


def outer_hamiltonian(plaquette, mean_fields, interactions, basis, verbose=False,
                      coeffs={'inner': {}, 'outer': {}},
                      checks=False):
    """
    Constructs a QuSpin quantum_operator corresponding to the cluster-bath
        Hamiltonian
    Input:
        plaquette (dict): Contains lists of neighbors, from files in the plaquettes folder
        mean_fields (dict): dict of mean fields (see example in README.md)
        interactions (dict): dict of interactions (see example in README.md)
        basis (QuSpin spin_basis_1d object): basis to construct the Hamiltonian in
        coeffs (list of #s): Factors to multiply each coupling for each site by.
    Output:
        Ho: QuSpin quantum_operator
    """
    L = plaquette['L']
    terms = []
    neighbors = ['nearest', 'n_nearest', 'n_n_nearest']
    bonds = ['x_bonds', 'y_bonds', 'z_bonds', 'n_bonds', 'nn_bonds', 'a_bonds', 'b_bonds']
    for n in interactions:
        if n in bonds:
            for c_op in interactions[n]:
                coupling = interactions[n][c_op]
                terms += [[c_op[0], [[coupling*mean_fields[c_op[1]][b[1]], b[0]] for b in plaquette['outer'][n]]]]
                terms += [[c_op[1], [[coupling*mean_fields[c_op[0]][b[0]], b[1]] for b in plaquette['outer'][n]]]]
        elif n == 'local':
            pass
        else:
            raise Exception('Unknown interaction type {}'.format(n))
    Ho = quantum_operator({'static': terms}, basis=basis, check_herm=checks,
                          check_symm=checks, dtype=TYPE)
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


def initialize_mf_hamiltonian(plaquette, interactions, basis,
                                 checks=False):
    L = plaquette['L']
    terms = {'inner': []}
    bonds = ['x_bonds', 'y_bonds', 'z_bonds', 'n_bonds', 'nn_bonds']
    for n in interactions:
        if n in bonds:
            for c_op in interactions[n]:
                coupling = interactions[n][c_op]
                for b in plaquette['outer'][n]:
                    s01 = c_op[0]+str(b[0])+'_'+c_op[1]+str(b[1])
                    s10 = c_op[1]+str(b[1])+'_'+c_op[0]+str(b[0])
                    if s01 in terms:
                        terms[s01] += [[c_op[0], [[coupling, b[0]]]]]
                    else:
                        terms[s01] = [[c_op[0], [[coupling, b[0]]]]]
                    if s10 in terms:
                        terms[s10] += [[c_op[1], [[coupling, b[1]]]]]
                    else:
                        terms[s10] = [[c_op[1], [[coupling, b[1]]]]]
                terms['inner'] += [[c_op, [[coupling, b[0], b[1]] for b in plaquette['inner'][n]]]]
        elif n == 'local':
            for c_op in interactions[n]:
                couplings = interactions[n][c_op]
                terms['inner'] += [[c_op, [[couplings[i], i] for i in range(L)]]]
        else:
            print('Incorrectly formatted interactions!')
    H = quantum_operator(terms, basis=basis, check_herm=checks, check_symm=checks, dtype=np.complex128)
    return H


def mf_params(mean_fields, interactions, plaquette):
    p_dict = {}
    for n in interactions:
        if n in ['x_bonds', 'y_bonds', 'z_bonds', 'n_bonds', 'nn_bonds']:
            for c_op in interactions[n]:
                coupling = 1. + 0.j # this will be multiplied times the coupling in the original hamiltonian
                for b in plaquette['outer'][n]:
                    s01 = c_op[0]+str(b[0])+'_'+c_op[1]+str(b[1])
                    s10 = c_op[1]+str(b[1])+'_'+c_op[0]+str(b[0])
                    p_dict[s01] = coupling*mean_fields[c_op[1]][b[1]]
                    p_dict[s10] = coupling*mean_fields[c_op[0]][b[0]]
    return p_dict
