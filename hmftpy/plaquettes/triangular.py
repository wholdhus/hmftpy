import numpy as np

def test_bonds(plaquette):
    good = True
    L = plaquette['L']
    cats = ['nearest', 'n_nearest', 'n_n_nearest']
    inner_correct = {c: [True for i in range(L)] for c in cats}
    outer_correct = {c: [True for i in range(L)] for c in cats}
    for c in cats:
        for i, sites in enumerate(plaquette['inner'][c]):
            for s in sites:
                if i in plaquette[c][s]:
                    pass # it's fine
                else:
                    inner_correct[c][i] = False
                    good = False
        for i, sites in enumerate(plaquette['outer'][c]):
            for s in sites:
                if i in plaquette[c][s]:
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
Connection diagram for 9 site 'pacman'
            <7>-<5>
            / \ / \
      <3>-<4>-<2>-<3>-<4>
        \ / \ / \ / \ /
    <0>-<1>- 8 -<0>-<1>-<8>
      \ / \ / \ / \ / \ /
  <6>-<7>- 5 - 6 - 7 -<5>-<6>
    \ / \ / \ / \ / \ / \ /
<3>-<4>- 2 - 3 - 4 -<2>-<3>-<4>
  \ / \ / \ / \ / \ / \ / \ /
  <1>-<8>- 0 - 1 -<8>-<0>-<1>
      / \ / \ / \ / \ / \
    <5>-<6>-<7>-<5>-<6>-<7>
          \ / \ / \ /
          <4>-<2>-<3>
"""


plaq9p = {'L': 9,
          'inner': {},
          'outer': {}}
plaq9p['rs'] = [np.array([0, 1, -.5, .5, 1.5, 0, 1, 2, .5]),
                np.array([0, 0, az, az, az, 2*az, 2*az, 2*az, 3*az])]
plaq9p['vs'] = [np.array([3, 0]), np.array([1.5, 3*az])]
plaq9p['outline'] = [np.array([0, 1, 1.5, 2, 1, .5, 0, -.5, 0]),
                     np.array([0, 0, az, 2*az, 2*az, 3*az, 2*az, az, 0])]
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
plaq9p['inner']['n_n_nearest'] = []
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
plaq9p['stripes'] = [[0,1,5,6,7],[2,3,4,8]]
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
plaq9d['rs'] = [np.array([0, -.5, .5, -1, 0, 1, -.5, .5, 0]),
                np.array([0, az, az, 2*az, 2*az, 2*az, 3*az, 3*az, 4*az])]
plaq9d['vs'] = [np.array([3,0]), np.array([1.5, 3*az])]
plaq9d['outline'] = [np.array([0, .5, 1, .5, 0, -.5, -1, -.5, 0]),
                     np.array([0, az, 2*az, 3*az, 4*az, 3*az, 2*az, az, 0])]
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
plaq9d['inner']['n_n_nearest'] = []
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
plaq9d['stripes'] = [[3,6,8,0,2,5],[1,4,7]]
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
            <4>-<5>
            / \ / \
      <9>-<0>-<1>-<2>-<7>
        \ / \ / \ / \ /
    <5>-<6>- a - b -<3>-<4>
      \ / \ / \ / \ / \ /
  <1>-<2>- 7 - 8 - 9 -<0>-<1>
    \ / \ / \ / \ / \ / \ /
<a>-<b>- 3 - 4 - 5 - 6 -<a>-<b>
  \ / \ / \ / \ / \ / \ / \ /
  <8>-<9>- 0 - 1 - 2 -<7>-<8>
      / \ / \ / \ / \ / \
    <5>-<6>-<a>-<b>-<3>-<4>
          \ / \ / \ /
          <7>-<8>-<9>


"""
plaq12 = {'L': 12,
          'inner': {},
          'outer': {}}
plaq12['rs'] = [np.array([0, 1, 2, -.5, .5, 1.5, 2.5, 0, 1, 2, .5, 1.5]),
                np.array([0, 0, 0, az, az, az, az, 2*az, 2*az, 2*az, 3*az, 3*az])]
plaq12['vs'] = [np.array([3, 2*az]), np.array([0, 4*az])]
plaq12['outline'] = [np.array([0, 1, 2, 2.5, 2, 1.5, .5, 0, -.5, 0]),
                     np.array([0, 0, 0, az, 2*az, 3*az, 3*az, 2*az, az, 0])]
plaq12['inner']['nearest'] = [
            [1,3,4],
            [0,2,4,5],
            [1,5,6],
            [0,4,7],
            [0,1,3,5,7,8],
            [1,2,4,6,8,9],
            [2,5,9],
            [3,4,8,10],
            [4,5,7,9,10,11],
            [5,6,8,11],
            [7,8,11],
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
            [2,8],
            [7,9],
            [0,8],
            [5,10],
            [6,11],
            [3,10],
            [4,11],
            [1,9],
            [0,2],
            [1,7],
            [3,5],
            [4,6]]
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
plaq12['stripes'] = [[0,3,2,5,8,10],
                     [1,4,7,6,9,11]]
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
    <2>-<0>-<1>-<2>
    / \ / \ / \ / \
  <b>- 9 - a - b -<9>-<a>
    \ / \ / \ / \ / \ / \
    <8>- 6 - 7 - 8 -<6>-<7>
  \ / \ / \ / \ / \ / \
  <5>- 3 - 4 - 5 -<3>-<4>
    \ / \ / \ / \ / \ / \ /
<1>-<2>- 0 - 1 - 2 -<0>-<1>
      \ / \ / \ / \ / \
      <9>-<a>-<b>-<9>-<a>


"""
plaq12z = {'L': 12,
           'inner': {},
           'outer': {}}
plaq12z['rs'] = [np.array([0, 1, 2, -.5, .5, 1.5, 0, 1, 2, -.5, .5, 1.5]),
                 np.array([0, 0, 0, az, az, az, 2*az, 2*az, 2*az, 3*az, 3*az, 3*az])]
plaq12z['vs'] = [np.array([3, 0]), np.array([0, 4*az])]
plaq12z['outline'] = [np.array([0, 2, 1.5, 2, 1.5, -.5, 0, -.5, 0]),
                      np.array([0, 0, az, 2*az, 3*az, 3*az, 2*az, az, 0])]
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
plaq12z['inner']['n_n_nearest'] = []
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
plaq12z['stripes'] = [[0,3,2,5,7,10],
                      [1,4,6,9,8,11]]
plaq12z['triangles'] = []
plaq12z['n_nearest_sublattice'] = [[0,5,6,11],
                                   [1,3,7,9],
                                   [2,4,8,10]]
plaq12z['n_n_nearest_sublattice'] = []

"""
Connection diagram for 21 site

             20
             / \
           17 -18 -19
           / \ / \ / \
         12 -13 -14 -15 -16
         / \ / \ / \ / \ /
        7 - 8 - 9 -10 -11
         \ / \ / \ / \ /
          3 - 4 - 5 - 6
           \ / \ / \ /
            0 - 1 - 2

"""
plaq21 = {'L': 21,
          'inner': {},
          'outer': {}}
plaq21['xs'] = np.array([0, 1, 2, 
                         -.5, .5, 1.5, 2.5, 
                         -1, 0, 1, 2, 3,
                         -.5, .5, 1.5, 2.5, 3.5,
                         0, 1, 2,
                         .5])
plaq21['ys'] = np.array([0, 0, 0, 
                         az, az, az, az, 
                         2*az, 2*az, 2*az, 2*az, 2*az,
                         3*az, 3*az, 3*az, 3*az, 3*az,
                         4*az, 4*az, 4*az,
                         5*az])
plaq21['inner']['nearest'] = [
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
            [6,10,15,16], # 11
            [7,8,13,17], # 12
            [8,9,12,14,17,18], # 13
            [9,10,13,15,18,19], # 14
            [10,11,14,16,19], # 15
            [11,15], # 16
            [12,13,18,20], # 17
            [13,14,17,19,20], # 18
            [14,15,18],
            [17,18]] # 20
plaq21['outer']['nearest'] = []
plaq21['inner']['n_nearest'] = []
plaq21['outer']['n_nearest'] = []
plaq21['inner']['n_n_nearest'] = []
plaq21['outer']['n_n_nearest'] = []
