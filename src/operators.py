from quspin.operators import quantum_operator

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
