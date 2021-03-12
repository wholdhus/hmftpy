from quspin.operators import quantum_operator


def inner_hamiltonian(plaquette, interactions, basis):
    bonds = []
    neighbors = ['nearest', 'n_nearest', 'n_n_nearest']
    if 'local' in interactions:
        for c_op in interactions['local']:
            coupling = interactions['local'][c_op]
            bonds += [[c_op, [[coupling, i] for i in range(plaquette['L'])]]]
    for i in range(plaquette['L']):
        for n in neighbors:
            if n in interactions:
                for c_op in interactions[n]:
                    coupling = interactions[n][c_op]
                    bonds += [[c_op, [[coupling, i, ni] for ni in plaquette['inner'][n][i]]]]
    Hi = quantum_operator({'static': bonds}, basis=basis,
                          check_herm=False, check_symm=False,
                          dtype=TYPE)
    return Hi


def periodic_hamiltonian(plaquette, interactions, basis):
    bonds = []
    if 'local' in interactions:
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
