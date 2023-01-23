import numpy as np

az = .5*np.sqrt(3) # vertical displacement for equilateral triangles
a0 = np.zeros(2)
ax = np.array([1., 0.])
a1 = np.array([0.5, 0.5*np.sqrt(3)])
a2 = np.array([-0.5, 0.5*np.sqrt(3)])


"""
Connection diagram for 3 site
    <0>-<1>-<2>
    / \ / \ / \
  <1>- 2 -<0>-<1>
  / \ / \ / \ /
<2>- 0 - 1 -<2>-
  \ / \ / \ / \
  <1>-<2>-<0>-<1>
"""
plaq3 = {'L': 3,
         'inner': {},
         'outer': {}}
plaq3['rs'] = [np.array([0, 1, 0.5]), np.array([0, 0, az])]
plaq3['vs'] = [np.array([1.5, az]), np.array([0, 2*az])]
plaq3['outline'] = [np.array([0, 1, .5, 0]), np.array([0, 0, az, 0])]
plaq3['inner']['n_bonds'] = [[0,1], [1,2], [2,0]]
plaq3['inner']['nn_bonds'] = []
plaq3['outer']['n_bonds'] = [[1,2], [2,0], [1,0], [2,1], [0,1], [2,0]]
plaq3['outer']['nn_bonds'] = ([[0,0] for i in range(6)]
                              + [[1,1] for i in range(6)]
                              + [[2,2] for i in range(6)])
plaq3['rs'] = [a0, ax, a1] # positions of sites within the cluster
plaq3['Rs'] = [ax+a1, a1+a2] # superlattice basis vectors
plaq3['outline'] = [0,1,2,0]

"""
Connection diagram for 7 site
    <4> <0> <1>
      \ / \ /
       5 - 6 -<2>
      / \ / \ /
     2 - 3 - 4 -<0>
      \ / \ / \
       0 - 1 -<5>

"""
plaq7 = {'L': 7,
         'inner': {},
         'outer': {}}
plaq7['inner']['n_bonds'] = [[0,1], [2,3], [3,4], [5,6],
                             [0,3], [1,4], [2,5], [3,6],
                             [0,2], [1,3], [3,5], [4,6]]
plaq7['inner']['nn_bonds'] = [[0,4], [2,6],
                              [0,5], [1,6],
                              [1,2], [4,5]]
plaq7['outer']['n_bonds'] = [[1,5], [4,0], [6,2],
                             [4,2], [5,0], [6,1],
                             [2,1], [5,4], [6,0]]
plaq7['outer']['nn_bonds'] = [[1,0], [3,2], [4,3], [5,1], [6,5],
                              [2,4], [3,0], [4,1], [5,2], [6,3],
                              [0,6], [2,0], [3,1], [5,3], [6,4]]
plaq7['rs'] = [a0, ax,
               a2, a1, ax+a1,
               a1+a2, 2*a1]
plaq7['Rs'] = [2*ax + a1, 2*a1+a2]
plaq7['outline'] = [0,1,4,6,5,2,0]

"""
Connection diagram for 12 site truncated triangle
     <0> <1> <2>
       \ / \ /
   <6>-10 -11 -<3>
     \ / \ / \ /
 <2>- 7 - 8 - 9 -<0>
   \ / \ / \ / \ /
11>-3 - 4 - 5 - 6 -<10>
   / \ / \ / \ / \
 <9>- 0 - 1 - 2 -<7>
     / \ / \ / \
   <6><10> <11><3>
"""
plaq12 = {'L': 12,
          'inner': {},
          'outer': {}}
plaq12['rs'] = [a0, ax, 2*ax,
                a2, a1, a1+ax, a1+2*ax,
                a1+a2, 2*a1, 2*a1+ax,
                2*a1+a2, 3*a1]
plaq12['Rs'] = [2*(a1+ax), 2*(a1+a2)]
plaq12['inner']['n_bonds'] = [
            [0,1], [1,2], [3,4], [4,5], [5,6], [7,8], [8,9], [10,11],
            [0,4], [1,5], [2,6], [3,7], [4,8], [5,9], [7,10], [8,11],
            [0,3], [1,4], [2,5], [4,7], [5,8], [6,9], [8,10], [9,11]]
plaq12['inner']['nn_bonds'] = [
            [0,5], [1,6], [3,8], [4,9], [7,11],
            [0,7], [1,8], [2,9], [4,10], [5,11],
            [1,3], [2,4], [5,7], [6,8], [9,10]]
plaq12['outer']['n_bonds'] = [
            [2,7], [6,10], [9,0], [11,3],
            [6,0], [9,3], [10,1], [11,2],
            [3,2], [7,6], [10,0], [11,1]]
plaq12['outer']['nn_bonds'] = [
            [2,10], [5,0], [6,1], [8,3], [9,4], [10,2], [11,7],
            [3,6], [6,3], [7,0], [8,1], [9,2], [10,4], [11,5],
            [0,11], [3,1], [4,2], [7,5], [8,6], [10,9], [11,0]]
plaq12['outline'] = [0,1,2,6,9,11,10,7,3,0]

"""
Connection diagram for 12 parallelogram that preserves 3SLS and 4SLS
    <8> <0> <1> <2>
      \ / \ / \ /
       9 -10 -11 -<3
        \ / \ / \ /
         6 - 7 - 8 -<0>
        / \ / \ / \
       3 - 4 - 5 -<9>
        \ / \ / \ /
         0 - 1 - 2 -<6>
"""
plaq12A = {'L': 12,
           'inner': {},
           'outer': {}}
plaq12A['rs'] = [a0, ax, 2*ax,
                 a2, a1, a1+ax,
                 a1+a2, 2*a1, 2*a1+ax,
                 a1+2*a2, 2*a1+a2, 3*a1]
plaq12A['Rs'] = [2*(ax+a1), 2*(ax-a2)]
plaq12A['inner']['n_bonds'] = [
            [0,1], [1,2], [3,4], [4,5], [6,7], [7,8], [9,10], [10,11],
            [0,4], [1,5], [3,6], [4,7], [5,8], [6,10], [7,11],
            [0,3], [1,4], [2,5], [4,6], [5,7], [6,9], [7,10], [8,11]]
plaq12A['inner']['nn_bonds'] = [
            [0,5], [3,7], [4,8], [6,11],
            [0,6], [1,7], [2,8], [3,9], [4,10], [5,11],
            [1,3], [2,4], [5,6], [7,9], [8,10]]
plaq12A['outer']['n_bonds'] = [
            [2,6], [5,9], [8,0], [11,3],
            [2,9], [8,3], [9,0], [10,1], [11,2],
            [3,2], [9,8], [10,0], [11,1]]
plaq12A['outer']['nn_bonds'] = [
            [1,9], [2,10], [5,0], [7,3], [8,4], [9,1], [10,2], [11,6],
            [6,0], [7,1], [8,2], [9,3], [10,4], [11,5],
            [0,11], [3,1], [4,2], [6,5], [9,7], [10,8], [11,0]]
plaq12A['outline'] = [0,1,2,5,8,11,10,9,6,3,0]

"""
Connection diagram for 12 parallelogram that preserves 3SLS, breaks 4SLS
    <2> <0> <1> <2>
      \ / \ / \ /
       9 -10 -11 -<9>
        \ / \ / \ /
         6 - 7 - 8 -<6>
        / \ / \ / \
       3 - 4 - 5 -<3>
        \ / \ / \ /
         0 - 1 - 2 -<0>
"""
plaq12B = {'L': 12,
           'inner': {},
           'outer': {}}
plaq12B['rs'] = [a0, ax, 2*ax,
                 a2, a1, a1+ax,
                 a1+a2, 2*a1, 2*a1+ax,
                 a1+2*a2, 2*a1+a2, 3*a1]
plaq12B['Rs'] = [3*ax, 2*(a1+a2)]
plaq12B['inner']['n_bonds'] = [
            [0,1], [1,2], [3,4], [4,5], [6,7], [7,8], [9,10], [10,11],
            [0,4], [1,5], [3,6], [4,7], [5,8], [6,10], [7,11],
            [0,3], [1,4], [2,5], [4,6], [5,7], [6,9], [7,10], [8,11]]
plaq12B['inner']['nn_bonds'] = [
            [0,5], [3,7], [4,8], [6,11],
            [0,6], [1,7], [2,8], [3,9], [4,10], [5,11],
            [1,3], [2,4], [5,6], [7,9], [8,10]]
plaq12B['outer']['n_bonds'] = [
            [2,0], [5,3], [8,6], [11,9],
            [2,3], [8,9], [9,0], [10,1], [11,2],
            [3,8], [9,2], [10,0], [11,1]]
plaq12B['outer']['nn_bonds'] = [
            [1,3], [2,4], [5,6], [7,9], [8,10],
            [9,1], [10,2], [11,0],
            [6,0], [7,1], [8,2],
            [9,3], [10,4], [11,5],
            [0,5], [3,7], [4,8], [6,11],
            [9,1], [10,2], [11,0]]
plaq12B['outline'] = [0,1,2,5,8,11,10,9,6,3,0]

"""
Connection diagram for 24 site cluster (concurrent with 12-site)
       <0> <1> <2>
         \ / \ /
     <6>-22 -23 -<15>
       \ / \ / \ /
   <2>-19 -20 -21 -<12>
     \ / \ / \ / \ /
<23>-15 -16 -17 -18 -<10>
     / \ / \ / \ / \
   21>-12 -13 -14 -<7>
       / \ / \ / \
    <18>-10 -11 -<3>
       \ / \ / \ /
  <14>- 7 - 8 - 9 -<0>
     \ / \ / \ / \ /
 <11>-3 - 4 - 5 - 6 -<22>
     / \ / \ / \ / \
   <9>- 0 - 1 - 2 -<19>
       / \ / \ / \
     <6>-<22><23><15>
"""
plaq24 = {'L': 24,
          'inner': {},
          'outer': {}}
plaq24['outline'] = [0,1,2,6,9,11,14,18,21,23,22,19,15,12,10,7,3,0]
plaq24['rs'] = plaq12['rs'] + [2*(a1+a2), 2*(a1+a2)+ax, 2*(a1+a2+ax),
                               2*a1+3*a2, 3*a1+2*a2, 4*a1+a2, 5*a1,
                               3*(a1+a2), 3*(a1+a2)+ax, 3*(a1+a2)+2*ax,
                               4*a1+3*a2, 5*a1+2*a2]
plaq24['Rs'] = [2*ax+2*a1, 4*(a1+a2)]
plaq24['inner']['n_bonds'] = (plaq12['inner']['n_bonds']
                             + [[12,13], [13,14], [15,16], [16,17], 
                                [17,18], [19,20], [20,21], [22,23],
                                [10,13], [11,14], [12,16], [13,17], 
                                [14,18], [15,19], [16,20], [17,21],
                                [19,22], [20,23],
                                [10,12], [11,13], [12,15], [13,16],
                                [14,17], [16,19], [17,20], [18,21],
                                [20,22], [21,23]])
plaq24['inner']['nn_bonds'] = (plaq12['inner']['nn_bonds']
                               + [[10,14], [12,17], [13,18], [15,20],
                                  [16,21], [19,23],
                                  [7,12], [8,13], [9,14], [10,16], [11,17],
                                  [12,19], [13,20], [14,21], [16,22], [17,23],
                                  [11,12], [13,15], [14,16], [17,19],
                                  [18,20], [21,22]])
plaq24['outer']['n_bonds'] = [
            [2,19], [6,22], [9,0], [11,3], [14,7], [18,10], [21,12], [23,15],
            [6,0], [9,3], [18,12], [21,15], [22,1], [23,2],
            [3,14], [7,18], [15,2], [19,6], [22,0], [23,1]]
plaq24['outer']['nn_bonds'] = [
            [2,22], [5,0], [6,1], [8,3], [9,4], [11,7], [14,10], [17,12],
            [18,13], [20,15], [21,16], [22,2], [23,19],
            [3,18], [6,3], [15,6], [18,15], [19,0], [20,1], [21,2], 
            [22,4], [23,5],
            [0,11], [3,13], [4,14], [7,17], [8,18], [10,21], [12,23], [15,1],
            [16,2], [19,5], [20,6], [22,9], [23,0]]


def paralellogram(Lx, Ly, offset_x=0, offset_y=0):
    N = Lx * Ly
    sites = np.arange(N).reshape((Lx, Ly))
    plaq = {'L': N,
            'inner': {'nearest': [[] for i in range(N)],
                      'n_nearest': [[] for i in range(N)]},
            'outer': {'nearest': [[] for i in range(N)],
                      'n_nearest': [[] for i in range(N)]}}
    for i in range(N):
        print(i)
        if i % Lx != Lx - 1:
            print('{} has inner horizontal neighbors'.format(i))
            plaq['inner']['nearest'][i] += [i + 1] # right
            plaq['inner']['nearest'][i + 1] += [i] # left
        else:
            plaq['outer']['nearest'][i] += [(i - (1+offset_y)*Lx + 1) % N] # right
            plaq['outer']['nearest'][(i - (1+offset_y)*Lx + 1) % N] += [i] # left
        if i + Lx < N:
            print('{} has inner upper left neighbors'.format(i))
            plaq['inner']['nearest'][i] += [i + Lx] # upper left
            plaq['inner']['nearest'][i + Lx] += [i] # lower right
        else:
            plaq['outer']['nearest'][i] += [(i + Lx - offset_x) % Lx] # upper left
            plaq['outer']['nearest'][(i + Lx - offset_x) % Lx] += [i] # lower right
        if i + Lx + 1 < N and i % Lx != Lx - 1:
            print('{} has inner upper right neighbors'.format(i))
            plaq['inner']['nearest'][i] += [i + Lx + 1] # upper right
            plaq['inner']['nearest'][i + Lx + 1] += [i] # lower left
        elif i % Lx != Lx - 1:
            plaq['outer']['nearest'][i] += [(i + Lx + 1 - offset_x) % Lx] # upper right
            plaq['outer']['nearest'][(i + Lx + 1 - offset_x) % Lx] += [i]
        else:
            plaq['outer']['nearest'][i] += [(i + 1 - offset_y * Lx) % N]
            plaq['outer']['nearest'][(i + 1 - offset_y * Lx) % N] += [i]
        if i + 2 * Lx + 1 < N and i % Lx != Lx - 1:
            print('{} has inner vertical nneighbors'.format(i))
            plaq['inner']['n_nearest'][i] += [i + 2 * Lx + 1]
            plaq['inner']['n_nearest'][i + 2 * Lx + 1] += [i]
        elif i % Lx != Lx - 1:
            plaq['outer']['n_nearest'][i] += [(i + 2 * Lx + 1 - offset_x) % N]
            plaq['outer']['n_nearest'][(i + 2 * Lx + 1 - offset_x) % N] += [i]
        else:
            plaq['outer']['n_nearest'][i] += [(i + (1-offset_y)*Lx + 1) % N]
            plaq['outer']['n_nearest'][(i + (1-offset_y)*Lx + 1) % N] += [i]
        if i < N - Lx and i % Lx != 0: # not top row
            print('{} has inner upper left nneighbors'.format(i))
            plaq['inner']['n_nearest'][i] += [i + Lx - 1]
            plaq['inner']['n_nearest'][i + Lx - 1] += [i]
        elif i % Lx != 0:
            plaq['outer']['n_nearest'][i] += [(i + Lx - 1) % N]
            plaq['outer']['n_nearest'][(i + Lx - 1)% N] += [i]
        else:
            plaq['outer']['n_nearest'][i] += [(i + (2+offset_y) * Lx - 1) % N]
            plaq['outer']['n_nearest'][(i + (2+offset_y) * Lx - 1) % N] += [i]
        if i < N - Lx and i % Lx < Lx - 2:
            print('{} has inner upper right nneighbors'.format(i))
            plaq['inner']['n_nearest'][i] += [i + Lx + 2]
            plaq['inner']['n_nearest'][i + Lx + 2] += [i]
        elif i % Lx < Lx - 2:
            plaq['outer']['n_nearest'][i] += [(i + Lx + 2) % N]
            plaq['outer']['n_nearest'][(i + Lx + 2) % N] += [i]
        else:
            plaq['outer']['n_nearest'][i] += [(i + 2 - offset_y * Lx) % N]
            plaq['outer']['n_nearest'][(i + 2 - offset_y * Lx) % N] += [i]
    return plaq

if __name__ == '__main__':
    from util import test_bonds
    print('Checking plaq3')
    print(test_bonds(plaq3, [6,6]))
    print('Checking plaq7')
    print(test_bonds(plaq7, [6,6]))
    print('Checking plaq12')
    print(test_bonds(plaq12, [6,6]))
    print('Checking plaq12A')
    print(test_bonds(plaq12A, [6,6]))
    print('Checking plaq12B')
    print(test_bonds(plaq12B, [6,6]))
    print('Checking plaq24')
    print(test_bonds(plaq24, [6,6]))