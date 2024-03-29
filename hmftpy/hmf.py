import numpy as np
from quspin.basis import spin_basis_1d, tensor_basis
from quspin.operators import hamiltonian, quantum_operator, quantum_LinearOperator
from .operators import inner_hamiltonian, outer_hamiltonian, mf_ops

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
        good_inds[u]= list(set(good_inds[u]))
    return good_inds


def do_hmft(plaquette, interactions, basis, max_iter=100, mf0=None,
            Hi=None, ops=None, lanczos_tol=10**-16, hmft_tol=10**-13,
            v0_start=False, n_states=1, coeffs={'inner': {}, 'outer': {}},
            mf_cvg=False, every_other=False,
            rescale_e=False):
    L = plaquette['L']
    if 'n_bonds' in interactions or 'x_bonds' in interactions:
        rescale_e = True
    if Hi is None:
        Hi = inner_hamiltonian(plaquette, interactions, basis, coeffs=coeffs, every_other=every_other)
    if ops is None:
        ops = mf_ops(plaquette, basis)
    if mf0 is None:
        mf0 = {'x': TYPE(2*(np.random.rand(L) - 0.5)),
               'y': TYPE(2*(np.random.rand(L) - 0.5)),
               'z': TYPE(2*(np.random.rand(L) - 0.5))}
    H = Hi + outer_hamiltonian(plaquette, mf0, interactions, basis, coeffs=coeffs, every_other=every_other)
    log('Hamiltonian complete!')
    if n_states < 0:
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

    if mf_cvg:
        mf_inds = get_useful_mf_inds(plaquette, interactions)

    while iter < max_iter and not converged:
        iter += 1
        log('{}th iteration'.format(iter))
        Ho = outer_hamiltonian(plaquette, mf, interactions, basis, coeffs=coeffs, every_other=every_other)
        if str(Ho) == '':
            log('Warning: mean fields are zero, outer Hamiltonian is null')
            H = Hi
        else:
            H = Hi + Ho
        if n_states < 0:
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
        if mf_cvg:
            mf_r = [mf[u][mf_inds[u]] for u in ['x', 'y', 'z']]
            cvg = 100
            if iter > 1:
                cvg = np.sum(np.linalg.norm(mf_r[i] - prev_mf_r[i]) for i in range(3))
                log('Total distance from previous mean fields')
                log(cvg)
            prev_mf_r = mf_r
        else: # basing convergence on energy
            cvg = np.abs(energies[iter] - energies[iter-1])
        if cvg < hmft_tol:
            converged = True
    e0 = energies[-1]
    if rescale_e:
        e0 = 0.5*(energies[-1] + Hi.expt_value(vs[-1])) # H_in + 0.5*H_out
    return np.real(e0), vs[-1], mf, converged
