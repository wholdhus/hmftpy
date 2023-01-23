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
    nearest_neighbors = np.zeros(L)
    next_nearest_neighbors = np.zeros(L)
    for i in range(L):
        for b in plaquette['inner']['n_bonds']:
            if i in b:
                nearest_neighbors[i] += 1
        for b in plaquette['outer']['n_bonds']:
            if i in b:
                nearest_neighbors[i] += 1
        for b in plaquette['inner']['nn_bonds']:
            if i in b:
                next_nearest_neighbors[i] += 1
        for b in plaquette['outer']['nn_bonds']:
            if i in b:
                next_nearest_neighbors[i] += 1
    nearest_test = nearest_neighbors == coordinations[0]
    next_nearest_test = next_nearest_neighbors == coordinations[1]
    test = nearest_test.all() and next_nearest_test.all()
    return test, nearest_neighbors, next_nearest_neighbors

def plot_plaq(plaq):
    R1 = plaq['Rs'][0]
    R2 = plaq['Rs'][1]
    Rs = [np.zeros(2), R1, R2, -R1, -R2, R1-R2, -R1+R2, R1+R2, -R1-R2]
    for R in Rs:
        outline_rs = [plaq['rs'][o]+R for o in plaq['outline']]
        plt.plot([r[0] for r in outline_rs],
                 [r[1] for r in outline_rs], color='gray', zorder=0)
        plt.fill([r[0] for r in outline_rs],
                 [r[1] for r in outline_rs], color='lightgray', zorder=-1)
        for i, r in enumerate(plaq['rs']):
            ri = r + R
            plt.scatter(ri[0], ri[1], color='black', marker='.')
            plt.text(ri[0], ri[1], i, clip_on=True)
