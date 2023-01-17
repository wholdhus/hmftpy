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
plaq3['outer']['n_bonds'] = [[0,1], [0,2], [0,1], [0,2],
                             [1,2], [1,0], [1,2], [1,0],
                             [2,0], [2,1], [2,0], [2,1]]
plaq3['outer']['nn_bonds'] = ([[0,0] for i in range(6)]
                              + [[1,1] for i in range(6)]
                              + [[2,2] for i in range(6)])
plaq3['rs'] = [a0, ax, a1] # positions of sites within the cluster
plaq3['Rs'] = [ax+a1, a1+a2] # superlattice basis vectors

"""
Connection diagram for 7 site
       5 - 6
      / \ / \
     2 - 3 - 4
      \ / \ /
       0 - 1

"""
plaq7 = {'L': 7,
         'inner': {},
         'outer': {}}
plaq7['outline'] = [np.array([0, 1, 1.5, 1, 0, -.5, 0]),
                     np.array([0, 0, az, 2*az, 2*az, az, 0])]
plaq7['rs'] = [np.array([0, 1, -.5, .5, 1.5, 0, 1]),
                np.array([0, 0, az, az, az, 2*az, 2*az])]
plaq7['inner']['nearest'] = [[1,3,2],
                    [4,3,0],
                    [0,3,5],
                    [0,1,2,4,5,6],
                    [6,3,1],
                    [2,3,6],
                    [5,3,4]]
plaq7['inner']['n_nearest'] = [[4,5],
                     [6,2],
                     [1,6],
                     [],
                     [5,0],
                     [0,4],
                     [2,1]]
plaq7['inner']['n_n_nearest'] = [[6],
                      [5],
                      [4],
                      [],
                      [2],
                      [1],
                      [0]]
sites = np.arange(7)
plaq7['outer']['nearest'] = [sites[[s2 not in np.append(plaq7['inner']['nearest'][s], s) for s2 in sites]] for s in sites]
plaq7['outer']['n_nearest'] = [sites[[s2 not in np.append(plaq7['inner']['n_nearest'][s], s) for s2 in sites]] for s in sites]
plaq7['outer']['n_n_nearest'] = [sites[[s2 not in np.append(plaq7['inner']['n_n_nearest'][s], s) for s2 in sites]] for s in sites]
plaq7['inner']['n3_nearest'] = [[] for i in range(12)]
plaq7['inner']['n4_nearest'] = [[] for i in range(12)]
plaq7['inner']['n5_nearest'] = [[] for i in range(12)]
plaq7['inner']['n6_nearest'] = [[] for i in range(12)]
plaq7['inner']['n7_nearest'] = [[] for i in range(12)]

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
plaq12['inner']['n_bonds'] = [
            [0,1], [0,4], [0,3], [1,2], [1,5], [1,4], [2,6], [2,5],
            [3,4], [3,7], [4,5], [4,8], [4,7], [5,6], [5,9], [5,8], [6,9],
            [7,8], [7,10], [8,9], [8,11], [8,10], [9,11], [10,11]]
plaq12['inner']['nn_bonds'] = [
            [0,5], [0,7], [1,6], [1,8], [1,3], [2,9], [2,4],
            [3,8], [4,9], [4,10], [5,11], [5,7], [6,8],
            [7,11], [9,10]]
plaq12['outer']['n_bonds'] = [
            [0,6], [0,10], [0,9], [1,10], [1,11], [2,11], [2,3], [2,7],
            [3,9], [3,11], [3,2], [6,7], [6,10], [6,0],
            [7,2], [7,6], [9,0], [9,3], [10,6],
            [10,0], [10,1], [11,3], [11,2], [11,1]]
plaq12['outer']['nn_bonds'] = [
            [0,7], [0,11], [0,11], [0,5],
            [1,8], [1,3], [1,6],
            [2,9], [2,4], [2,10], [2,10],
            [3,6], [3,6], [3,1], [3,8],
            [4,10], [4,2], [4,9],
            [5,11], [5,7], [5,0],
            [6,3], [6,8], [6,1], [6,3],
            [7,0], [7,5], [7,11],
            [8,3], [8,1], [8,6],
            [9,10], [9,4], [9,2],
            [10,2], [10,4], [10,9], [10,2],
            [11,0], [11,7], [11,5], [11,0]]
plaq12['inner']['nearest'] = [
            [1,3,4], # 0
            [0,2,4,5],
            [1,5,6], # 2
            [0,4,7],
            [0,1,3,5,7,8], # 4
            [1,2,4,6,8,9],
            [2,5,9], # 6
            [3,4,8,10],
            [4,5,7,9,10,11], # 8
            [5,6,8,11],
            [7,8,11], # 10
            [8,9,10]
            ]
plaq12['inner']['n_nearest'] = [
            [5,7], # 0
            [3,6,8],
            [4,9], # 2
            [1,8],
            [2,9,10], # 4
            [0,7,11],
            [1,8], # 6
            [0,5,11],
            [1,3,6], # 8
            [2,4,10],
            [4,9], # 10
            [5,7]]
plaq12['inner']['n_n_nearest'] = [
            [2,8], # 0
            [7,9],
            [0,8], # 2
            [5,10],
            [6,11], # 4
            [3,10],
            [4,11], # 6
            [1,9],
            [0,2], # 8
            [1,7],
            [3,5], # 10
            [4,6]]
plaq12['inner']['n3_nearest'] = [
            [6,9,10], # 0
            [10,11],
            [3,7,11], # 2
            [2,9,11],
            [], # 4
            [],
            [0,7,10], # 6
            [2,6],
            [], # 8
            [0,3],
            [0,1,6], # 10
            [1,2,3]]
plaq12['inner']['n4_nearest'] = [
            [11], # 0
            [],
            [10], # 2
            [6],
            [], # 4
            [],
            [3], # 6
            [],
            [], # 8
            [],
            [2], # 10
            [0]]
plaq12['inner']['n5_nearest'] = [[] for i in range(12)]
plaq12['inner']['n6_nearest'] = [[] for i in range(12)]
plaq12['inner']['n7_nearest'] = [[] for i in range(12)]

plaq12['outer']['nearest'] = [
            [6,9,10],
            [10,11],
            [3,7,11],
            [2,9,11],
            [],
            [],
            [0,7,10],
            [2,6],
            [],
            [0,3],
            [0,1,6],
            [1,2,3]
            ]
plaq12['outer']['n_nearest'] = [
            [5,7,11,11], # 0
            [3,6,8],
            [4,9,10,10], # 2
            [1,6,6,8],
            [2,9,10], # 4
            [0,7,11],
            [1,3,3,8], # 6
            [0,5,11],
            [1,3,6], # 8
            [2,4,10],
            [2,2,4,9], # 10
            [0,0,5,7]]
plaq12['outer']['n_n_nearest'] = [
            [2,2,8,8],
            [7,7,9,9],
            [0,0,8,8],
            [5,5,10,10],
            [6,6,11,11],
            [3,3,10,10],
            [4,4,11,11],
            [1,1,9,9],
            [0,0,2,2],
            [1,1,7,7],
            [3,3,5,5],
            [4,4,6,6]]

plaq12['rs'] = [np.array([0, 1, 2, -.5, .5, 1.5, 2.5, 0, 1, 2, .5, 1.5]),
                np.array([0, 0, 0, az, az, az, az, 2*az, 2*az, 2*az, 3*az, 3*az])]
plaq12['vs'] = [np.array([3, 2*az]), np.array([0, 4*az])]
plaq12['outline'] = [np.array([0, 1, 2, 2.5, 2, 1.5, .5, 0, -.5, 0]),
plaq12['rs'] = np.array([a0, a1, 2*a1,
                         a2-a1, a2, a2+a1, a2+2*a1,
                         2*a2-a1, 2*a2, 2*a2+a1,
                         3*a2-a1, 3*a2
                         ])
plaq12['Rs'] = [a0, 2*a1+2*a2, 4*a2-2*a1, 2*a2-4*a1,
                -2*a2-2*a1, -4*a2+2*a1, 4*a1-2*a2]
plaq12['outline'] = [0,1,2,6,9,11,10,7,3,0]

"""
Connection diagram for 12 zigzag

       9 -10 -11
        \ / \ / \
         6 - 7 - 8
        / \ / \ /
       3 - 4 - 5
        \ / \ / \
         0 - 1 - 2
"""
plaq12z = {'L': 12,
           'inner': {},
           'outer': {}}

plaq12z['inner']['nearest'] = [
            [1,3,4], # 0
            [0,2,4,5],
            [1,5], # 2
            [0,4,6],
            [0,1,3,5,6,7], # 4
            [1,2,4,7,8],
            [3,4,7,9,10], # 6
            [4,5,6,8,10,11],
            [5,7,11], # 8
            [6,10],
            [6,7,9,11], # 10
            [7,8,10]]
plaq12z['inner']['n_nearest'] = [
            [5,6], # 0
            [3,7],
            [4,8], # 2
            [1,7,9],
            [2,8,10], # 4
            [0,6,11],
            [0,5,11], # 6
            [1,3,9],
            [2,4,10], # 8
            [3,7],
            [4,8], # 10
            [5,6]]
plaq12z['inner']['n_n_nearest'] = [
            [2,7], # 0
            [6,8],
            [0,7], # 2
            [5,10],
            [9,11], # 4
            [3,10],
            [1,8], # 6
            [0,2],
            [1,6], # 8
            [4,11],
            [3,5], # 10
            [4,9]]
plaq12z['outer']['nearest'] = [
            [2,9,10],
            [10,11],
            [0,3,9,11], # 2
            [2,5,8],
            [], # 4
            [3],
            [8], # 6
            [],
            [3,6,9], # 8
            [0,2,8,11],
            [0,1], # 10
            [1,2,9]]
plaq12z['outer']['n_nearest'] = [
            [5,6,11,11], # 0
            [3,7,9,9],
            [4,8,10,10], # 2
            [1,7,9,],
            [2,8,10], # 4
            [0,6,11],
            [0,5,11], # 6
            [1,3,9],
            [2,4,10], # 8
            [1,1,3,7],
            [2,2,4,8], # 10
            [0,0,5,6]]
plaq12z['outer']['n_n_nearest'] = []
plaq12z['inner']['n3_nearest'] = [
            [8,9,10], # 0
            [10,11],
            [3,6,11], # 2
            [2,8,11],
            [], # 4
            [9],
            [2], # 6
            [],
            [0,3,9], # 8
            [0,5,8],
            [0,1], # 10
            [1,2,3]]
plaq12z['inner']['n4_nearest'] = [
            [11], # 0
            [9],
            [10], # 2
            [], [], [], [], [], [],
            [1],
            [2], # 10
            [0]]
plaq12z['inner']['n5_nearest'] = [[] for i in range(12)]
plaq12z['inner']['n6_nearest'] = [[] for i in range(12)]
plaq12z['inner']['n6_nearest'][2] = [9]
plaq12z['inner']['n6_nearest'][9] = [2]
plaq12z['inner']['n7_nearest'] = [[] for i in range(12)]

plaq12z['rs'] = [np.array([0, 1, 2, -.5, .5, 1.5, 0, 1, 2, -.5, .5, 1.5]),
                 np.array([0, 0, 0, az, az, az, 2*az, 2*az, 2*az, 3*az, 3*az, 3*az])]
plaq12z['vs'] = [np.array([3, 0]), np.array([0, 4*az])]
plaq12z['outline'] = [np.array([0, 2, 1.5, 2, 1.5, -.5, 0, -.5, 0]),
                      np.array([0, 0, az, 2*az, 3*az, 3*az, 2*az, az, 0])]



"""
Connection diagram for 19 site
           16 -17 -18
           / \ / \ / \
         12 -13 -14 -15
         / \ / \ / \ / \
        7 - 8 - 9 -10 -11
         \ / \ / \ / \ /
          3 - 4 - 5 - 6
           \ / \ / \ /
            0 - 1 - 2
"""
plaq19 = {'L': 19,
          'inner': {},
          'outer': {}}
plaq19['inner']['nearest'] = [
                        [1,3,4], # 0
                        [2,5,4,0], # 1
                        [6,5,1], # 2
                        [0,4,8,7], # 3
                        [0,1,5,9,8,3], # 4
                        [2,6,10,9,4,1], # 5
                        [11,10,5,2], # 6
                        [12,8,3], # 7
                        [3,4,9,13,12,7], # 8
                        [4,5,8,10,13,14], # 9
                        [5,6,11,15,14,9], # 10
                        [6,10,15], # 11
                        [7,8,13,16], # 12
                        [8,9,14,17,16,12], # 13
                        [9,10,15,18,17,13], # 14
                        [10,11,14,18], # 15
                        [12,13,17], # 16
                        [13,14,18,16], # 17
                        [17,14,15]] # 18
plaq19['outer']['nearest'] = [
                        [15,11,16], # 0
                        [16,17],
                        [17,18,7], # 2
                        [18,15],
                        [], # 4
                        [],
                        [7,12], # 6
                        [18,2,6],
                        [], # 8
                        [],
                        [], # 10
                        [12,16,0],
                        [11,6], # 12
                        [],
                        [], # 14
                        [0,3],
                        [1,0,11], # 16
                        [1,2],
                        [3,7,2]] # 18
plaq19['inner']['n_nearest'] = [
                        [5,8], # 0
                        [6,9,3],
                        [10,4], # 2
                        [1,9,12],
                        [2,10,13,7], # 4
                        [11,14,8,0],
                        [15,9,1], # 6
                        [4,13],
                        [0,5,14,16], # 8
                        [1,6,15,17,12,3],
                        [2,18,13,4], # 10
                        [5,14],
                        [3,9,17], # 12
                        [7,4,10,18],
                        [5,11,16,8], # 14
                        [6,9,17],
                        [8,14], # 16
                        [12,9,15],
                        [13,10]] # 18
plaq19['outer']['n_nearest'] = [
                        [18,10,12,17], # 0
                        [11,13,18],
                        [16,14,3,12], # 2
                        [2,14,11],
                        [15,16], # 4
                        [17,7],
                        [18,8,16], # 6
                        [11,5,17,15],
                        [6,18], # 8
                        [],
                        [12,0], # 10
                        [7,13,1,3],
                        [0,10,2], # 12
                        [1,11],
                        [3,2], # 14
                        [16,4,7],
                        [2,15,6,4], # 16
                        [7,5,0],
                        [0,8,6,1]] # 18
plaq19['inner']['n_n_nearest'] = [
                        [2,9,7], # 0
                        [10,8],
                        [11,9,0], # 2
                        [5,13],
                        [6,14,12], # 4
                        [15,13,3],
                        [14,4], # 6
                        [0,9,16],
                        [1,10,17], # 8
                        [0,2,7,11,16,18],
                        [1,17,8], # 10
                        [2,9,18],
                        [4,14], # 12
                        [3,5,15],
                        [4,6,12], # 14
                        [5,13],
                        [7,9,18], # 16
                        [8,10],
                        [9,11,16]] # 18
plaq19['outer']['n_n_nearest'] = [
                        [14,6,13], # 0
                        [15,12,14,7],
                        [13,15,8], # 2
                        [6,16,10,17],
                        [18,11,17], # 4
                        [16,18,12],
                        [17,3,13,0], # 6
                        [10,1,14],
                        [11,2,15], # 8
                        [],
                        [7,16,3], # 10
                        [8,17,4],
                        [1,15,5,18], # 12
                        [2,0,6],
                        [0,7,1], # 14
                        [2,8,1,12],
                        [5,3,10], # 16
                        [3,6,4,11],
                        [4,12,5]] # 18


plaq19['n_nearest_sublattice'] = [[0,5,8,11,14,16],
                                  [1,3,6,9,12,15,17],
                                  [2,4,7,10,13,18]]


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
plaq24['inner']['n_bonds'] = [
            [0,1], [0,4], [0,3], [1,2], [1,5], [1,4], [2,6], [2,5],
            [3,4], [3,7], [4,5], [4,8], [4,7], [5,6], [5,9], [5,8], [6,9],
            [7,8], [7,10], [8,9], [8,11], [8,10], [9,11],
            [10,11], [10,13], [10,12], [11,14], [11,13],
            [12,13], [12,16], [12,15], [13,14], [13,17], [13,16],
                [14,18], [14,17],
            [15,16], [15,19], [16,17], [16,20], [16,19],
                [17,18], [17,21], [17,20], [18,21],
            [19,20], [19,22], [20,21], [20,23], [20,22], [21,23]
            [22,23]]
plaq12['inner']['nn_bonds'] = [
            [0,5], [0,7], [1,6], [1,8], [1,3], [2,9], [2,4],
            [3,8], [4,9], [4,10], [5,11], [5,7], [6,8],
            [7,11], [7,12], [9,14], [9,10],
            [10,14], [10,16], [11,17], [11,12],
            [12,17], [12,19], [13,18], [13,20], [13,15], [14,21], [14,16],
            [15,20], [16,21], [16,22], [17,23], [17,19], [18,20],
            [19,23], [21,22]]
plaq24['outer']['n_bonds'] = [
            [0,9], [0,6], [0,22], [1,22], [1,23], [2,19], [2,23], [2,15],
            [3,14], [3,11], [3,9], [6,22], [6,0], [6,19],
            [7,18], [7,14], [9,0], [9,3],
            [10,18], [11,3],
            [12,21], [12,18], [14,7], [14,3],
            [15,2], [15,23], [15,21], [18,10], [18,12], [18,7],
            [19,6], [19,2], [21,12], [21,15],
            [22,1], [22,0], [22,6], [23,15], [23,2], [23,1]
            ]
plaq24['outer']['nn_bonds'] = [
            [0,7], [0,11], [0,5],
            [1,8], [1,3], [1,6],
            [2,9], [2,4], [2,10], [2,10],
            [3,6], [3,6], [3,1], [3,8],
            [4,10], [4,2], [4,9],
            [5,11], [5,7], [5,0],
            [6,3], [6,8], [6,1], [6,3],
            [7,0], [7,5], [7,11],
            [8,3], [8,1], [8,6],
            [9,10], [9,4], [9,2],
            [10,2], [10,4], [10,9], [10,2],
            [11,0], [11,7], [11,5], [11,0]]


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
