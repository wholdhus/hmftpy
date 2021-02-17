import numpy as np
from quspin.basis import spin_basis_1d, tensor_basis
from quspin.operators import hamiltonian, quantum_operator
from tqdm import tqdm
from kithcmb import kitaevhoneycomb

VERBOSE = False
MAXIT = None
TOL = 10**-14
HTOL = 10**-11
TYPE = np.complex128





def log(msg):
    """
    Prints msg only when global variable VERBOSE is true.
    """
    if VERBOSE:
        print(msg)


def inner_hamiltonian(plaquette, interactions, basis):
    bonds = []
    for c_op in interactions['local']:
        coupling = interactions['local'][c_op]
        bonds += [[c_op, [[coupling, i] for i in range(plaquette['L'])]]]
    for i in range(plaquette['L']):
        for c_op in interactions['n']:
            coupling = interactions['n'][c_op]
            bonds+= [[c_op, [[coupling, i, ni] for ni in plaquette['inner_n'][i]]]]
        for c_op in interactions['nn']:
            coupling = interactions['nn'][c_op]
            bonds+= [[c_op, [[coupling, i, nni] for nni in plaquette['inner_nn'][i]]]]
        for c_op in interactions['nnn']:
            coupling = interactions['nnn'][c_op]
            bonds+= [[c_op, [[coupling, i, nnni] for nnni in plaquette['inner_nnn'][i]]]]
    Hi = quantum_operator({'static': bonds}, basis=basis,
                          check_herm=False, check_symm=False,
                          dtype=TYPE)
    return Hi


def periodic_hamiltonian(plaquette, interactions, basis):
    bonds = []
    for c_op in interactions['local']:
        coupling = interactions['local'][c_op]
        bonds += [[c_op, [[coupling, i] for i in range(plaquette['L'])]]]
    for i in range(plaquette['L']):
        for c_op in interactions['n']:
            coupling = interactions['n'][c_op]
            bonds+= [[c_op, [[coupling, i, ni] for ni in plaquette['inner_n'][i]]]]
            bonds+= [[c_op, [[coupling, i, ni] for ni in plaquette['outer_n'][i]]]]
        for c_op in interactions['nn']:
            coupling = interactions['nn'][c_op]
            bonds+= [[c_op, [[coupling, i, nni] for nni in plaquette['inner_nn'][i]]]]
            bonds+= [[c_op, [[coupling, i, nni] for nni in plaquette['outer_nn'][i]]]]
        for c_op in interactions['nnn']:
            coupling = interactions['nnn'][c_op]
            bonds+= [[c_op, [[coupling, i, nnni] for nnni in plaquette['inner_nnn'][i]]]]
            bonds+= [[c_op, [[coupling, i, nnni] for nnni in plaquette['outer_nnn'][i]]]]

    Hi = quantum_operator({'static': bonds}, basis=basis,
                          check_herm=False, check_symm=False,
                          dtype=TYPE)
    return Hi


def outer_hamiltonian(plaquette, mean_fields, interactions, basis):
    bonds = []
    for i in range(plaquette['L']):
        for c_op in interactions['n']:
            coupling = interactions['n'][c_op]
            bonds += [[c_op[1], [[coupling*mean_fields[c_op[0]][ni], i] for ni in plaquette['outer_n'][i]]]]

        for c_op in interactions['nn']:
            coupling = interactions['nn'][c_op]
            bonds += [[c_op[1], [[coupling*mean_fields[c_op[0]][nni], i] for nni in plaquette['outer_nn'][i]]]]

        for c_op in interactions['nnn']:
            coupling = interactions['nnn'][c_op]
            bonds += [[c_op[1], [[coupling*mean_fields[c_op[0]][nnni], i] for nnni in plaquette['outer_nnn'][i]]]]

    Ho = quantum_operator({'static': bonds}, basis=basis, check_herm=False,
                          check_symm=False, dtype=TYPE)
    return Ho


def mf_ops(plaquette, basis):
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


def get_mfs(v, mf_ops):
    L = len(mf_ops)
    xmfs = np.zeros(L, dtype=TYPE)
    ymfs = np.zeros(L, dtype=TYPE)
    zmfs = np.zeros(L, dtype=TYPE)
    for i in range(L):
        xmfs[i] = mf_ops[i]['x'].matrix_ele(v, v)
        ymfs[i] = mf_ops[i]['y'].matrix_ele(v, v)
        zmfs[i] = mf_ops[i]['z'].matrix_ele(v, v)
    mfs = {'x': xmfs, 'y': ymfs, 'z': zmfs}
    return mfs


def do_hmft(plaquette, interactions, basis, max_iter=100, mf0=None, full_diag=False, n_states=1, double_ext=False,
            avoid_zero=False, Hi=None, ops=None):
    L = plaquette['L']
    if Hi is None:
        Hi = inner_hamiltonian(plaquette, interactions, basis)
    if ops is None:
        ops = mf_ops(plaquette, basis)
    if mf0 is None:
        mf0 = {'x': TYPE(2*(np.random.rand(L) - 0.5)),
               'y': TYPE(2*(np.random.rand(L) - 0.5)),
               'z': TYPE(2*(np.random.rand(L) - 0.5))}
    H = Hi + outer_hamiltonian(plaquette, mf0, interactions, basis)
    if double_ext:
        H += outer_hamiltonian(plaquette, mf0, interactions, basis)
    log('Hamiltonian complete!')
    if full_diag:
        e, v = H.eigh()
    else:
        e, v = H.eigsh(k=n_states, which='SA', maxiter=MAXIT, tol=TOL)
    log('Ground state energy with initial seed')
    log(e[0])

    energies = [e[0]]
    vs = [v[:, 0]]
    mf = get_mfs(vs[0], ops)
    iter = 0
    converged = False
    while iter < max_iter and not converged:
        iter += 1
        log('{}th iteration'.format(iter))
        H = Hi + outer_hamiltonian(plaquette, mf, interactions, basis)
        if double_ext:
            H += outer_hamiltonian(plaquette, mf, interactions, basis)

        if full_diag:
            e, v = H.eigh()
            es_degen = e[np.abs(e - e[0]) < TOL]
            log('Ground state degeneracy: {}'.format(len(es_degen)))
        else:
            e, v = H.eigsh(k=n_states, which = 'SA') # just finding lowest energies
        log('Energy: {}'.format(e[0]))
        mf = get_mfs(v[:,0], ops)
        energies += [e[0]]
        vs += [v[:, 0]]
        cvg = np.abs(energies[iter] - energies[iter - 1])
        if cvg < HTOL:
            converged = True
        if avoid_zero:
            if np.mean(np.abs(mf['x'])) < HTOL or np.mean(np.abs(mf['y'])) < HTOL or np.mean(np.abs(mf['z'])) < HTOL:
                print('Zero mean-field! Restarting with random seed!')
                mf = {'x': TYPE(2*(np.random.rand(L) - 0.5)),
                       'y': TYPE(2*(np.random.rand(L) - 0.5)),
                       'z': TYPE(2*(np.random.rand(L) - 0.5))}
                converged = False
    if double_ext:
        H2 = Hi + outer_hamiltonian(plaquette, mf, interactions, basis)
        if full_diag:
            e, v = H2.eigh()
        else:
            e, v = H2.eigsh(k=1, which='SA')
            e0 = e[0]
    else:
        # ep = Hi.matrix_ele(vs[-1], vs[-1])
        # e0 = 0.5*(ep + energies[-1]) # (2 H_i + H_o)/2
        e0 = energies[-1]
    return np.real(e0), vs[-1], mf, converged



def neel_op(L, basis): # Neel OP operator in the z direction
    op_lst = [['z', [[(-1)**i, i] for i in range(L)]]]
    return quantum_operater({'static': op_lst}, basis=basis)


def bond_chiral_z(L, basis):
    pinds = []
    minds = []
    for i in range(L):
        pinds += [[1, i, j] for j in range(L)]
        minds += [[-1, i, j] for j in range(L)]
    op_lst = [['xy', pinds], ['yx', minds]]
    return quantum_operator({'static': op_lst}, basis=basis)


if __name__ == '__main__':
    from plaquettes.triangular import plaq12, plaq12_hex, plaq7
    basis = spin_basis_1d(12, pauli=0)
    ops = mf_ops(plaq12, basis)
    interactions = {'local': {},
                    'n': {'xx': 1, 'yy': 1, 'zz': 1},
                    'nn': {'xx': .1, 'yy': .1, 'zz': .1},
                    'nnn': {}
                    }
    print('12 site')
    print('HMFT')
    e, v, mf, cvg = do_hmft(plaq12, interactions, basis, mf0=None)
    print('Energy density:')
    print(e/24)
    print('Converged?')
    print(cvg)
    Hi = inner_hamiltonian(plaq12, interactions, basis)
    e, v = Hi.eigsh(k=1, which='SA')
    print('ED energy with open b.c.')
    print(e[0]/24)
    Hp = periodic_hamiltonian(plaq12, interactions, basis)
    e, v = Hp.eigsh(k=1, which='SA')
    print('ED energy with periodic b.c.')
    print(e[0]/24)

    print('')
    b7 = spin_basis_1d(7, pauli=0)
    print('7 site')
    print('Running HMFT')
    e, v, mf, cvg = do_hmft(plaq7, interactions, b7, mf0=None)
    print('Energy density:')
    print(e/14)
    print('Converged?')
    print(cvg)
    Hp = periodic_hamiltonian(plaq7, interactions, b7)
    e, v = Hp.eigsh(k=1, which='SA')
    print('ED energy with periodic b.c.')
    print(e[0]/14)
    Hi = inner_hamiltonian(plaq7, interactions, b7)
    e, v = Hi.eigsh(k=1, which='SA')
    print('ED energy with open b.c.')
    print(e[0]/14)
