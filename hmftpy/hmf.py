import numpy as np
from quspin.basis import spin_basis_1d, tensor_basis
from quspin.operators import hamiltonian, quantum_operator, quantum_LinearOperator
from operators import inner_hamiltonian, outer_hamiltonian, mf_ops

VERBOSE = False
TYPE = np.complex128

def log(msg):
    """
    Prints msg only when global variable VERBOSE is true.
    """
    if VERBOSE:
        print(msg)

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
            avoid_zero=False, Hi=None, ops=None, lanczos_tol=10**-15, hmft_tol=10**-13,
            v0_start=False):
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
        e, v = H.eigsh(k=n_states, which='SA', tol=lanczos_tol)
    log('Ground state energy with initial seed')
    log(e[0])

    energies = [e[0]]
    vs = [v[:, 0]]
    mf = get_mfs(vs[0], ops)
    iter = 0
    converged = False
    v0 = None
    while iter < max_iter and not converged:
        iter += 1
        log('{}th iteration'.format(iter))
        H = Hi + outer_hamiltonian(plaquette, mf, interactions, basis)
        if double_ext:
            H += outer_hamiltonian(plaquette, mf, interactions, basis)

        if full_diag:
            e, v = H.eigh()
            es_degen = e[np.abs(e - e[0]) < 5*lanczos_tol]
            log('Ground state degeneracy: {}'.format(len(es_degen)))
        else:
            e, v = H.eigsh(k=n_states, which = 'SA', tol=lanczos_tol, v0=v0) # just finding lowest energies
        if v0_start:
            v0 = v[:,0]
        log('Energy: {}'.format(e[0]))
        mf = get_mfs(v[:,0], ops)
        energies += [e[0]]
        vs += [v[:, 0]]
        cvg = np.abs(energies[iter] - energies[iter - 1])
        if cvg < hmft_tol:
            converged = True
        if avoid_zero:
            if np.mean(np.abs(mf['x'])) < hmft_tol or np.mean(np.abs(mf['y'])) < hmft_tol or np.mean(np.abs(mf['z'])) < hmft_tol:
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
            e, v = H2.eigsh(k=1, which='SA', tol=lanczos_tol)
            e0 = e[0]
    else:
        e0 = energies[-1]
    return np.real(e0), vs[-1], mf, converged




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
