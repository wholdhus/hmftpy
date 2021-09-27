"""
Connection diagram for 4-site square (symmetric tiling)

  <1>-<0>-<1>-<0>
   |   |   |   |
  <2>- 3 - 2 -<3>
   |   |   |   |
  <1>- 0 - 1 -<0>
   |   |   |   |
  <2>-<3>-<2>-<3>
"""
plaq4 = {'L': 4,
          'inner': {},
          'outer': {}}
plaq4['inner']['nearest'] = [
            [1,3], # 0
            [0,2],
            [1,3], # 2
            [0,2]]
plaq4['inner']['n_nearest'] = [
            [2], # 0
            [3],
            [0], # 2
            [1]]
plaq4['outer']['nearest'] = [
            [1,3], # 0
            [0,2],
            [1,3], # 2
            [0,2]] 
plaq4['outer']['n_nearest'] = [
            [2, 2, 2],
            [3, 3, 3],
            [0, 0, 0],
            [1, 1, 1]]
plaq4['inner']['x_bonds'] = [[0,1], [3,2]]
plaq4['inner']['y_bonds'] = [[0,3], [1,2]]
plaq4['inner']['a_bonds'] = [[0,2]]
plaq4['inner']['b_bonds'] = [[1,3]]
plaq4['outer']['x_bonds'] = [[0,1], [1,0], [3,2], [2,3]]
plaq4['outer']['y_bonds'] = [[0,3], [3,0], [1,2], [2,1]]
plaq4['outer']['a_bonds'] = [[0,2], [0,2], [0,2],
                             [2,0], [2,0], [2,0]]
plaq4['outer']['b_bonds'] = [[1,3], [1,3], [1,3],
                             [3,1], [3,1], [3,1]]

plaq4['inner']['n_bonds'] = [[0,1], [3,2], [0,3], [1,2],
                            ]
plaq4['inner']['nn_bonds'] = [[0,2], [1,3],
                             ]
plaq4['outer']['n_bonds'] = [[0,1], [0,3],
                             [1,0], [1,2],
                             [2,1], [2,3],
                             [3,0], [3,2]]
plaq4['outer']['nn_bonds'] = [[0,2], [0,2], [0,2], 
                              [1,3], [1,3], [1,3],
                              [2,0], [2,0], [2,0],
                              [3,1], [3,1], [3,1]]