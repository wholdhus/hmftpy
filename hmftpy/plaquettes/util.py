import numpy as np
import matplotlib.pyplot as plt

def test_neighbors(plaquette):
    good = True
    L = plaquette['L']
    icats = ['nearest', 'n_nearest']
    ocats = ['nearest', 'n_nearest']

    inner_correct = {c: [True for i in range(L)] for c in icats}
    outer_correct = {c: [True for i in range(L)] for c in ocats}
    for c in icats:
        for i, sites in enumerate(plaquette['inner'][c]):

            for s in sites:
                if i in plaquette['inner'][c][s]:
                    pass # it's fine
                else:
                    inner_correct[c][i] = False
                    good = False
    for c in ocats:
        for i, sites in enumerate(plaquette['outer'][c]):
            for s in sites:
                if i in plaquette['outer'][c][s]:
                    pass # it's fine
                else:
                    outer_correct[c][i] = False
                    good = False
    return good, inner_correct, outer_correct

def test_bonds(plaquette, coordinations):
    good = True
    L = plaquette['L']
    nearest_neighbors = [0 for i in range(L)]
    next_nearest_neighbors = np.zeros(L)
    for i in range(L):
        for b in plaquette['inner']['n_bonds']:
            if i in b:
                nearest_neighbors[i] += 1
        for b in plaquette['outer']['n_bonds']:
            if i in b:
                nearest_neighbors[i] += 0.5 # because these are counted twice
        for b in plaquette['inner']['nn_bonds']:
            if i in b:
                next_nearest_neighbors[i] += 1
        for b in plaquette['outer']['nn_bonds']:
            if i in b:
                next_nearest_neighbors[i] += 0.5
    nearest_test = nearest_neighbors == coordinations[0]
    next_nearest_test = next_nearest_neighbors == coordinations[1]
    return nearest_test, next_nearest_test

def plot_plaq(plaquette):
    rs = plaquette['rs']
    Rs = plaquette['Rs']
    xs = rs[:,0]
    ys = rs[:,1]
    print(xs)
    print(ys)
    for j, R in enumerate(plaquette['Rs']):
        if j == 0:
            color='black'
        else:
            color='red'
        for i in range(plaquette['L']):
            plt.scatter([xs[i]+R[0]], [ys[i]+R[1]], color=color)
            plt.text(xs[i]+R[0], ys[i]+R[1], str(i))
