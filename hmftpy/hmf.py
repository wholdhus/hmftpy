import numpy as np
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian, quantum_operator, quantum_LinearOperator
from .operators import inner_hamiltonian, outer_hamiltonian, mf_ops
from .operators import initialize_mf_hamiltonian, mf_params

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
        xmfs[i] = mf_ops[i]['x'].expt_value(v)
        ymfs[i] = mf_ops[i]['y'].expt_value(v)
        zmfs[i] = mf_ops[i]['z'].expt_value(v)
    mfs = {'x': xmfs, 'y': ymfs, 'z': zmfs}
    return mfs


def get_useful_mf_inds(plaquette, interactions):
    good_inds = {'x': [], 'y': [], 'z': []}
    for n in interactions:
        if n != 'local':
            for c_op in interactions[n]:
                for neighbors in plaquette['outer'][n]:
                    for u in ['x', 'y', 'z']:
                        if c_op[0] == u:
                            good_inds[u] += neighbors
                        if c_op[1] == u:
                            good_inds[u] += neighbors
    for u in ['x', 'y', 'z']:
        good_inds[u] = list(set(good_inds[u]))
    return good_inds


def do_hmft(plaquette, interactions, basis, max_iter=100, mf0=None,
            Hi=None, ops=None, lanczos_tol=10**-16, hmft_tol=10**-13,
            coeffs={'inner': {}, 'outer': {}},
            mf_cvg=False, v0=None):
    L = plaquette['L']
    if Hi is None:
        Hi = inner_hamiltonian(plaquette, interactions, basis,
                               coeffs=coeffs)
    if ops is None:
        ops = mf_ops(plaquette, basis)
    if mf0 is None:
        mf0 = {'x': TYPE(2*(np.random.rand(L) - 0.5)),
               'y': TYPE(2*(np.random.rand(L) - 0.5)),
               'z': TYPE(2*(np.random.rand(L) - 0.5))}
    H = Hi + outer_hamiltonian(plaquette, mf0, interactions, basis,
                               coeffs=coeffs)
    log('Hamiltonian complete!')
    e, v = H.eigsh(k=1, which='SA', tol=lanczos_tol, v0=v0)
    v0 = v[:,0]
    log('Ground state energy with initial seed')
    log(e[0])

    energies = []
    vs = []
    mf = mf0
    iter = 0
    converged = False
    if mf_cvg:
        mf_inds = get_useful_mf_inds(plaquette, interactions)

    while iter < max_iter and not converged:
        log('{}th iteration'.format(iter))
        Ho = outer_hamiltonian(plaquette, mf, interactions, basis, coeffs=coeffs)
        if str(Ho) == '':
            log('Warning: mean fields are zero, outer Hamiltonian is null')
            H = Hi
        else:
            H = Hi + Ho
        e, v = H.eigsh(k=1, which = 'SA', tol=lanczos_tol, v0=v0) # just finding lowest energies
        v0 = v[:,0]
        log('Energy: {}'.format(e[0]))
        mf = get_mfs(v[:, 0], ops)
        energies += [e[0]]
        vs += [v[:, 0]]
        cvg = 999
        if iter > 0:
            if mf_cvg:
                mf_r = [mf[u][mf_inds[u]] for u in ['x', 'y', 'z']]
                if iter > 1:
                    cvg = np.sum([np.linalg.norm(mf_r[i] - prev_mf_r[i])
                                  for i in range(3)])
                    log('Total distance from previous mean fields')
                    log(cvg)
                prev_mf_r = mf_r
            else: # basing convergence on energy
                cvg = np.abs(energies[iter] - energies[iter-1])
        if cvg < hmft_tol:
            converged = True
        iter += 1
    Ho = outer_hamiltonian(plaquette, mf, interactions, basis, coeffs=coeffs)
    e0 = 0.5*(Ho+Hi+Hi).expt_value(vs[-1]) # H_in + 0.5*H_out
    return np.real(e0), vs[-1], mf, converged


def outer_hamiltonian_params(mean_fields, interactions, plaquette):
    p_dict = {}
    for n in interactions:
        if n in ['n_bonds', 'nn_bonds']:
            for c_op in interactions[n]:
                for b in plaquette['outer'][n]:
                    s01 = c_op[0]+str(b[0])+'_'+c_op[1]+str(b[1])
                    s10 = c_op[1]+str(b[1])+'_'+c_op[0]+str(b[0])
                    p_dict[s01] = mean_fields[c_op[1]][b[1]]
                    p_dict[s10] = mean_fields[c_op[0]][b[0]]
    return p_dict


def do_hmft_2(plaquette, H, interactions, basis, max_iter=100, mf0=None,
              ops=None, lanczos_tol=10**-16, hmft_tol=10**-13,
              mf_cvg=False):
    L = plaquette['L']
    if ops is None:
        ops = mf_ops(plaquette, basis)
    if mf0 is None:
        mf0 = {'x': TYPE(2*(np.random.rand(L) - 0.5)),
               'y': TYPE(2*(np.random.rand(L) - 0.5)),
               'z': TYPE(2*(np.random.rand(L) - 0.5))}
    params = mf_params(mf0, interactions, plaquette)
    e, v = H.eigsh(params, k=1, which='SA', tol=lanczos_tol)
    log('Ground state energy with initial seed')
    log(e[0])

    energies = [e[0]]
    vs = [v[:, 0]]
    mf = get_mfs(vs[0], ops)
    iter = 0
    converged = False
    v0 = None

    if mf_cvg:
        mf_inds = get_useful_mf_inds(plaquette, interactions)

    while iter < max_iter and not converged:
        iter += 1
        log('{}th iteration'.format(iter))
        params = mf_params(mf, interactions, plaquette)
        e, v = H.eigsh(params, k=1, which = 'SA', tol=lanczos_tol, v0=v0) # just finding lowest energies
        v0 = v[:,0]
        log('Energy: {}'.format(e[0]))
        mf = get_mfs(v0, ops)
        energies += [e[0]]
        vs += [v0]
        if mf_cvg:
            mf_r = [mf[u][mf_inds[u]] for u in ['x', 'y', 'z']]
            cvg = 100
            if iter > 0:
                cvg = np.sum(np.linalg.norm(mf_r[i] - prev_mf_r[i]) for i in range(3))
                log('Total distance from previous mean fields')
                log(cvg)
            prev_mf_r = mf_r
        else: # basing convergence on energy
            cvg = np.abs(energies[iter] - energies[iter-1])
        if cvg < hmft_tol:
            converged = True
    p0 = mf_params({'x': np.zeros(L, dtype=np.complex128),
                    'y': np.zeros(L, dtype=np.complex128),
                    'z': np.zeros(L, dtype=np.complex128)},
                    interactions, plaquette)
    params = mf_params(mf, interactions, plaquette)
    Hi = inner_hamiltonian(plaquette, interactions, basis)
    e0 = 0.5*(Hi.expt_value(vs[-1]) + H.expt_value(vs[-1], pars=params))

    return np.real(e0), vs[-1], mf, converged

if __name__ == '__main__':
    from plaquettes.triangular import plaq12
    import time
    basis = spin_basis_1d(12, pauli=0)
    interactions = {'n_bonds': {'xx': 1., 'yy': 1., 'zz': 1.},
                    'nn_bonds': {'xx': 0.1, 'yy': 0.1, 'zz': 0.1}}
    mf0 = {'x': (np.random.rand(12)-.5), 'y': (np.random.rand(12)-.5), 'z': (np.random.rand(12)-.5)}
    mf0p = outer_hamiltonian_params(mf0, interactions, plaq12)
    mfh = initialize_mf_hamiltonian(plaq12, interactions, basis)

    print('New method: ')
    tic = time.time()
    e, v, mf, cvg = do_hmft_2(plaq12, mfh, interactions, basis, mf0=mf0)
    toc = time.time()
    print('That took {} seconds'.format(toc-tic))
    print('Complete!')
    print('Energy: {}'.format(e))
    print('Converged? ')
    print(cvg)

    print('')
    print('Old method: ')
    tic = time.time()
    e, v, mf, cvg = do_hmft(plaq12, interactions, basis, mf0=mf0)
    toc = time.time()
    print('That took {} seconds'.format(toc-tic))
    print('Complete!')
    print('Energy: {}'.format(e))
    print('Converged? ')
    print(cvg)
