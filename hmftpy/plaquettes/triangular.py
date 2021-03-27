import numpy as np

def test_bonds(plaquette):
    good = True
    L = plaquette['L']
    icats = ['nearest', 'n_nearest', 'n_n_nearest',
             'n3_nearest', 'n4_nearest', 'n5_nearest',
             'n6_nearest', 'n7_nearest']
    ocats = ['nearest', 'n_nearest', 'n_n_nearest']

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

az = .5*np.sqrt(3) # vertical displacement for equilateral triangles




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

plaq3['inner']['nearest'] = [
            [1,2], [0,2], [0,1]]
plaq3['inner']['n_nearest'] = [
            [], [], []]
plaq3['inner']['n_n_nearest'] = [
            [], [], []]
plaq3['outer']['nearest'] = [
            [1,1,2,2],
            [0,0,2,2],
            [0,0,1,1]]
plaq3['outer']['n_nearest'] = [
            [0 for i in range(6)],
            [1 for i in range(6)],
            [2 for i in range(6)]]
plaq3['outer']['n_n_nearest'] = [
            [1,1,1,2,2,2],
            [0,0,0,2,2,2],
            [0,0,0,1,1,1]]

plaq3['inner']['x_bonds'] = [[0,1]]
plaq3['inner']['y_bonds'] = [[1,2]]
plaq3['inner']['z_bonds'] = [[2,0]]
plaq3['outer']['x_bonds'] = [[1,2], [2,0]]
plaq3['outer']['y_bonds'] = [[0,1], [2,0]]
plaq3['outer']['z_bonds'] = [[0,1], [1,2]]

"""
Connection diagram for 7 site
       5 - 6
      / \ / \
     2 - 3 - 4
      \ / \ /
       0 - 1

"""


# Inter couplings: both are real spins
# Outer couplings: 1st is mean field param. index, second is real spin

plaq7 = {'L': 7,
         'inner': {},
         'outer': {}}
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

plaq7['x_stripes'] = [[0,1,5,6],
                      [2,3,4]]
plaq7['y_stripes'] = [[1,4,2,5],
                      [0,3,6]]
plaq7['z_stripes'] = [[0,2,4,6],
                      [1,3,5]]
plaq7['triangles'] = [[3,1,0], [5,3,2], [6,4,3]]
plaq7['n_nearest_sublattice'] = [[0,4,5],
                                  [1,2,6],
                                  [3]]

"""
Connection diagram for 9 site 'pacman'

         <1>- 8 -<0>
           \ / \ / \
            5 - 6 - 7
           / \ / \ /
          2 - 3 - 4
           \ / \ /
            0 - 1

"""


plaq9p = {'L': 9,
          'inner': {},
          'outer': {}}
plaq9p['inner']['nearest'] = [
            [1,2,3], # 0
            [0,3,4],
            [0,3,5], # 2
            [0,1,2,4,5,6],
            [1,3,6,7], # 4
            [2,3,6,8],
            [3,4,5,7,8], # 6
            [4,6],
            [5,6]] # 8
plaq9p['inner']['n_nearest'] = [
            [4,5], # 0
            [2,6],
            [1,6], # 2
            [7,8],
            [0,5], # 4
            [0,4],
            [1,2], # 6
            [3,8],
            [3,7]]
plaq9p['inner']['n_n_nearest'] = [
            [6], # 0
            [5,7],
            [4,8], # 2
            [],
            [2,8], # 4
            [1,7],
            [0], # 6
            [1,5],
            [2,4]] # 8
plaq9p['outer']['nearest'] = [
            [6,7,8], # 0
            [5,7,8],
            [4,7,8], # 2
            [],
            [2,8], # 4
            [1,7],
            [0], # 6
            [0,1,2,5],
            [0,1,2,4]] # 8
plaq9p['outer']['n_nearest'] = [
            [4,4,5,5], # 0
            [2,2,6,6],
            [1,1,6,6], # 2
            [7,7,8,8],
            [0,0,5,5], # 4
            [0,0,4,4],
            [1,1,2,2], # 6
            [3,3,8,8],
            [3,3,7,7]] # 8
plaq9p['outer']['n_n_nearest'] = []
plaq9p['inner']['n3_nearest'] = [[] for i in range(9)]
plaq9p['inner']['n3_nearest'][0] = [7,8]
plaq9p['inner']['n3_nearest'][1] = [8]
plaq9p['inner']['n3_nearest'][2] = [7]
plaq9p['inner']['n3_nearest'][7] = [0,2]
plaq9p['inner']['n3_nearest'][8] = [0,1]
plaq9p['inner']['n4_nearest'] = [[] for i in range(9)]
plaq9p['inner']['n5_nearest'] = [[] for i in range(9)]
plaq9p['inner']['n6_nearest'] = [[] for i in range(9)]
plaq9p['inner']['n7_nearest'] = [[] for i in range(9)]
plaq9p['x_stripes'] = [[0,1,5,6,7],[2,3,4,8]]
plaq9p['y_stripes'] = [[1,4,7,2,5,8],[0,3,6]]
plaq9p['z_stripes'] = [[0,2,4,6,8],[1,3,5,7]]
plaq9p['rs'] = [np.array([0, 1, -.5, .5, 1.5, 0, 1, 2, .5]),
                np.array([0, 0, az, az, az, 2*az, 2*az, 2*az, 3*az])]
plaq9p['vs'] = [np.array([3, 0]), np.array([1.5, 3*az])]
plaq9p['outline'] = [np.array([0, 1, 1.5, 2, 1, .5, 0, -.5, 0]),
                     np.array([0, 0, az, 2*az, 2*az, 3*az, 2*az, az, 0])]
plaq9p['triangles'] = [[3,1,0], [5,3,2], [6,4,3], [8,6,5]]
plaq9p['n_nearest_sublattice'] = [[0,4,5],[1,2,6], [3,7,8]]
plaq9p['inner']['x_bonds'] = [[0,1],
                              [2,3], [3,4],
                              [5,6], [6,7]]
plaq9p['inner']['y_bonds'] = [[0,2], [1,3],
                              [3,5], [4,6],
                              [6,8]]
plaq9p['inner']['z_bonds'] = [[3,0], [4,1],
                              [5,2], [6,3], [7,4],
                              [8,5]]
plaq9p['outer']['x_bonds'] = [[1,8], [4,2], [7,5], [8,0]]
plaq9p['outer']['y_bonds'] = [[2,7], [5,1], [7,0], [8,4]]
plaq9p['outer']['z_bonds'] = [[0,6], [1,7], [2,8]]
"""
Connection diagram for alt. 9 site 'diamond'
            8
           / \
          6 - 7
         / \ / \
        3 - 4 - 5
         \ / \ /
          1 - 2
           \ /
            0

"""
plaq9d = {'L': 9,
          'inner': {},
          'outer': {}}
plaq9d['inner']['nearest'] = [
            [1,2], # 0
            [0,2,3,4],
            [0,1,4,5], # 2
            [1,4,6],
            [1,2,3,5,6,7], # 4
            [2,4,7],
            [3,4,7,8], # 6
            [4,5,6,8],
            [6,7]] # 8
plaq9d['inner']['n_nearest'] = [
            [4], # 0
            [5,6],
            [3,7], # 2
            [2,7],
            [0,8], # 4
            [1,6],
            [1,5], # 6
            [2,3],
            [4]]
plaq9d['inner']['n_n_nearest'] = [
            [3,5], # 0
            [7],
            [6], # 2
            [0,5,8],
            [], # 4
            [0,3,8],
            [2], # 6
            [1],
            [3,5]] # 8
plaq9d['outer']['nearest'] = [
            [3,5,6,7], # 0
            [7,8],
            [6,8], # 2
            [0,5,8],
            [], # 4
            [0,3,8],
            [0,2], # 6
            [0,1],
            [1,2,3,5]] # 8
plaq9d['outer']['n_nearest'] = [
            [4,4,8,8,8], # 0
            [5,5,6,6],
            [3,3,7,7], # 2
            [2,2,7,7],
            [0,0,8,8], # 4
            [1,1,6,6],
            [1,1,5,5], # 6
            [2,2,3,3],
            [0,0,4,4]] # 8
plaq9d['outer']['n_n_nearest'] = []
plaq9d['inner']['n3_nearest'] = [[6,7], # 0
                                 [8],
                                 [8], # 2
                                 [], [], [],
                                 [0], # 6
                                 [0],
                                 [1,2] # 8
                                 ]
plaq9d['inner']['n4_nearest'] = [[] for i in range(9)]
plaq9d['inner']['n4_nearest'][0] = [8]
plaq9d['inner']['n4_nearest'][8] = [0]
plaq9d['inner']['n5_nearest'] = [[] for i in range(9)]
plaq9d['inner']['n6_nearest'] = [[] for i in range(9)]
plaq9d['inner']['n7_nearest'] = [[] for i in range(9)]


plaq9d['rs'] = [np.array([0, -.5, .5, -1, 0, 1, -.5, .5, 0]),
                np.array([0, az, az, 2*az, 2*az, 2*az, 3*az, 3*az, 4*az])]
plaq9d['vs'] = [np.array([3,0]), np.array([1.5, 3*az])]
plaq9d['outline'] = [np.array([0, .5, 1, .5, 0, -.5, -1, -.5, 0]),
                     np.array([0, az, 2*az, 3*az, 4*az, 3*az, 2*az, az, 0])]

plaq9d['x_stripes'] = [[0,3,4,5,8], [1,2,6,7]]
plaq9d['y_stripes'] = [[0,2,5,3,6,8], [1,4,7]]
plaq9d['z_stripes'] = [[0,1,3,5,7,8], [2,4,6]]

plaq9d['triangles'] = [[4,2,1], [6,4,3], [7,5,4], [8,7,6]]
plaq9d['n_nearest_sublattice'] = [[0,4,8],
                                  [1,5,6],
                                  [2,3,7]]

plaq9d['inner']['x_bonds'] = [[1,2],
                             [3,4], [4,5],
                             [6,7]]
plaq9d['inner']['y_bonds'] = [[0,1],
                             [1,3], [2,4],
                             [4,6], [5,7],
                             [7,8]]
plaq9d['inner']['z_bonds'] = [[2,0],
                             [4,1], [5,2],
                             [6,3], [7,4],
                             [8,6]]
plaq9d['outer']['x_bonds'] = [[0,6], [2,8], [5,3], [7,0], [8,1]]
plaq9d['outer']['y_bonds'] = [[3,0], [6,2], [8,5]]
plaq9d['outer']['z_bonds'] = [[0,5], [1,7], [3,8]]

"""
Connection diagram for 12 site truncated triangle
       10 -11
       / \ / \
      7 - 8 - 9
     / \ / \ / \
    3 - 4 - 5 - 6
     \ / \ / \ /
      0 - 1 - 2

"""
plaq12 = {'L': 12,
          'inner': {},
          'outer': {}}
plaq12['x_stripes'] = [[0,1,2,7,8,9],
                       [3,4,5,6,10,11]]
plaq12['y_stripes'] = [[2,6,0,4,8,11],
                       [1,5,9,3,7,10]]
plaq12['z_stripes'] = [[0,3,2,5,8,10],
                       [1,4,7,6,9,11]]

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
                     np.array([0, 0, 0, az, 2*az, 3*az, 3*az, 2*az, az, 0])]

plaq12['triangles'] = [[4,1,0],
                       [5,2,1],
                       [7,4,3],
                       [8,5,4],
                       [9,6,5],
                       [10,8,7],
                       [11,9,8]]
plaq12['n_nearest_sublattice'] = [[0,5,7,11],
                                  [1,3,6,8],
                                  [2,4,9,10]]
plaq12['n_n_nearest_sublattice'] = [[0,2,8],
                                    [1,7,9],
                                    [3,5,10],
                                    [4,6,11]]

plaq12['inner']['x_bonds'] = [[0,1], [1,2],
                              [3,4], [4,5], [5,6],
                              [7,8], [8,9],
                              [10, 11]]
plaq12['inner']['y_bonds'] = [[0,3], [1,4], [2,5],
                              [4,7], [5,8], [6,9],
                              [8,10], [9,11]]
plaq12['inner']['z_bonds'] = [[4,0], [5,1], [6,2],
                              [7,3], [8,4], [9,5],
                              [10,7], [11,8]]
plaq12['outer']['x_bonds'] = [[2,7], [6,10], [9,0], [11,3]]
                              # [0,9], [3,11], [7,2], [10,6] # THESE ARE BACKWARD?
plaq12['outer']['y_bonds'] = [[3,2], [7,6], [10,0], [11,1]]
plaq12['outer']['z_bonds'] = [[0,6], [1,10], [2,11], [3,9]]


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
plaq12z['inner']['n7_nearest'] = [[] for i in range(12)]

plaq12z['rs'] = [np.array([0, 1, 2, -.5, .5, 1.5, 0, 1, 2, -.5, .5, 1.5]),
                 np.array([0, 0, 0, az, az, az, 2*az, 2*az, 2*az, 3*az, 3*az, 3*az])]
plaq12z['vs'] = [np.array([3, 0]), np.array([0, 4*az])]
plaq12z['outline'] = [np.array([0, 2, 1.5, 2, 1.5, -.5, 0, -.5, 0]),
                      np.array([0, 0, az, 2*az, 3*az, 3*az, 2*az, az, 0])]
plaq12z['x_stripes'] = [[0,1,2,6,7,8],
                        [3,4,5,9,10,11]]
plaq12z['y_stripes'] = [[2,0,4,7,11,9],
                        [1,5,8,3,6,10]]
plaq12z['z_stripes'] = [[0,3,2,5,7,10],
                        [1,4,6,9,8,11]]
plaq12z['triangles'] = []
plaq12z['n_nearest_sublattice'] = [[0,5,6,11],
                                   [1,3,7,9],
                                   [2,4,8,10]]
plaq12z['n_n_nearest_sublattice'] = []


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
plaq19['x_stripes'] = [[0,1,2,7,8,9,10,11,16,17,18],
                       [3,4,5,6,12,13,14,15]]
plaq19['y_stripes'] = [[2,6,11,0,4,9,14,18,7,12,16],
                       [1,5,10,15,3,8,13,17]]
plaq19['z_stripes'] = [[0,3,7,2,5,9,13,16,11,15,18],
                       [1,4,8,12,2,5,9,13,16,11,15,18]]
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

plaq19['inner']['n3_nearest'] = [
                        [6,10,12,13], # 0
                        [7,11,13,14],
                        [3,8,14,15], # 2
                        [2,10,14,16],
                        [11,15,16,17], # 4
                        [7,12,17,18],
                        [0,8,13,18], # 6
                        [1,5,14,17],
                        [2,6,15,18], # 8
                        [],
                        [0,3,12,16], # 10
                        [1,4,13,17],
                        [0,5,10,18], # 12
                        [0,1,6,11],
                        [1,2,3,7], # 14
                        [2,4,8,16],
                        [3,4,10,15], # 16
                        [4,5,7,11],
                        [5,6,8,12]] # 18
plaq19['inner']['n4_nearest'] = [
                        [14], # 0
                        [12,15],
                        [13], # 2
                        [6,17],
                        [18], # 4
                        [16],
                        [3,17], # 6
                        [10],
                        [11], # 8
                        [],
                        [7], # 10
                        [8],
                        [1,15], # 12
                        [2],
                        [0], # 14
                        [1,12],
                        [5], # 16
                        [3,6],
                        [4]] # 18
plaq19['inner']['n5_nearest'] = [
                        [11,16], # 0
                        [17],
                        [7,18], # 2
                        [15],
                        [], # 4
                        [],
                        [12], # 6
                        [2,18],
                        [], # 8
                        [],
                        [], # 10
                        [0,16],
                        [6], # 12
                        [],
                        [], # 14
                        [3],
                        [0,11], # 16
                        [1],
                        [2,7]] # 18
plaq19['inner']['n6_nearest'] = [
                        [15,17], # 0
                        [16,18],
                        [12,17], # 2
                        [11,18],
                        [], # 4
                        [],
                        [7,16], # 6
                        [6,15],
                        [], # 8
                        [],
                        [], # 10
                        [3,12],
                        [2,11], # 12
                        [],
                        [], # 14
                        [0,7],
                        [1,6], # 16
                        [0,2],
                        [1,3]] # 18
plaq19['inner']['n7_nearest'] = [
                        [18], # 0
                        [],
                        [16], # 2
                        [],
                        [], # 4
                        [],
                        [], # 6
                        [11],
                        [], # 8
                        [],
                        [], # 10
                        [7],
                        [], # 12
                        [],
                        [], # 14
                        [],
                        [2], # 16
                        [],
                        [0]] # 18

plaq19['n_nearest_sublattice'] = [[0,5,8,11,14,16],
                                  [1,3,6,9,12,15,17],
                                  [2,4,7,10,13,18]]


"""
Connection diagram for 21 site


             20
             / \
           16 -17 -18 -19
           / \ / \ / \ /
         12 -13 -14 -15
         / \ / \ / \ / \
        7 - 8 - 9 -10 -11
         \ / \ / \ / \ /
          3 - 4 - 5 - 6
           \ / \ / \ /
            0 - 1 - 2

"""
plaq21 = {'L': 21,
          'inner': {},
          'outer': {}}
plaq21['xs'] = np.array([])
plaq21['ys'] = np.array([])
plaq21['inner']['nearest'] = [
            [1,3,4], # 0
            [0,2,4,5], # 1
            [1,5,6], # 2
            [0,4,7,8], # 3
            [0,1,3,5,8,9], # 4
            [1,2,4,6,9,10], # 5
            [2,5,10,11], # 6
            [3,8,12], # 7
            [3,4,7,9,12,13], # 8
            [4,5,8,10,13,14], # 9
            [5,6,9,11,14,15], # 10
            [6,10,15], # 11
            [7,8,13,16], # 12
            [8,9,12,14,16,17], # 13
            [9,10,13,15,17,18], # 14
            [10,11,14,18,19], # 15
            [12,13,17,20], # 16
            [13,14,16,18,20], # 17
            [14,15,17,19], # 18
            [15,18],
            [16,17]] # 20
plaq21['outer']['nearest'] = []
plaq21['inner']['n_nearest'] = []
plaq21['outer']['n_nearest'] = []
plaq21['inner']['n_n_nearest'] = []
plaq21['outer']['n_n_nearest'] = []


"""
27 site

             24 -25 -26
             / \ / \ / \
           20 -21 -22 -23
           / \ / \ / \ / \
         15 -16 -17 -18 -19
         / \ / \ / \ / \ / \
        9 -10 -11 -12 -13 -14
         \ / \ / \ / \ / \ /
          4 - 5 - 6 - 7 - 8
           \ / \ / \ / \ /
            0 - 1 - 2 - 3

"""
plaq21 = {'L': 27,
          'inner': {},
          'outer': {}}
plaq21['xs'] = np.array([])
plaq21['ys'] = np.array([])
plaq21['inner']['nearest'] = []
plaq21['outer']['nearest'] = []
plaq21['inner']['n_nearest'] = []
plaq21['outer']['n_nearest'] = []
plaq21['inner']['n_n_nearest'] = []
plaq21['outer']['n_n_nearest'] = []
