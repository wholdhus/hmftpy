"""
Connection diagram for 4-site square (symmetric tiling)

      <0>-<1>
       |   |
  <3>- 2 - 3 -<2>
   |   |   |   |
  <1>- 0 - 1 -<0>
       |   |
      <2>-<3>
"""
plaq4 = {'L': 4,
          'inner': {},
          'outer': {}}
plaq4['inner']['nearest'] = [[1,2], # 0
                             [0,3],
                             [0,3], # 2
                             [1,2]]
plaq4['inner']['n_nearest'] = [[3], # 0
                               [2],
                               [1], # 2
                               [0]]
plaq4['outer']['nearest'] = [[1,2], # 0
                             [0,3],
                             [0,3], # 2
                             [1,2]] # 8
plaq4['outer']['n_nearest'] = [[3, 3, 3],
                               [2, 2, 2],
                               [1, 1, 1],
                               [0, 0, 0]]
plaq4['inner']['x_bonds'] = [[0,1], [3,2]]
plaq4['inner']['y_bonds'] = [[1,3], [2,0]]
plaq4['outer']['x_bonds'] = [[0,1], [3,2]]
plaq4['outer']['y_bonds'] = [[1,3], [2,0]]
# plaq4['inner']['x_bonds'] = [[0,1], [2,3]]
# plaq4['inner']['y_bonds'] = [[0,2], [1,3]]
# plaq4['outer']['x_bonds'] = [[1,0], [3,2]]
# plaq4['outer']['y_bonds'] = [[2,0], [3,1]]


"""
Connection diagram for 9-site square (symmetric tiling)
      <0>-<1>-<2>
       |   |   |
  <8>- 6 - 7 - 8 -<6>
   |   |   |   |   |
  <5>- 3 - 4 - 5 -<3>
   |   |   |   |   |
  <2>- 0 - 1 - 2 -<0>
       |   |   |
      <6>-<7>-<8>
"""
plaq9 = {'L': 9,
          'inner': {},
          'outer': {}}
plaq9['inner']['nearest'] = [[]]
plaq9['inner']['x_bonds'] = [[0,1], [1,2], 
                             [3,4], [4,5],
                             [6,7], [7,8]]
plaq9['inner']['y_bonds'] = [[0,3], [1,4], [2,5],
                             [3,6], [4,7], [5,8]]
plaq9['outer']['x_bonds'] = [[2,0], [5,3], [8,6]]
plaq9['outer']['y_bonds'] = [[6,0], [7,1], [8,2]]


"""
Connection diagram for 16-site square (symmetric tiling)
      <0>-<1>-<2>-<3>
       |   |   |   |
  15>-12 -13 -14 -15 -12>
       |   |   |   |   |
  11>- 8 - 9 -10 -11 -<8>
   |   |   |   |   |   |
  <7>- 4 - 5 - 6 - 7 -<4>
   |   |   |   |   |   |
  <0>- 0 - 1 - 2 - 3 -<0>
       |   |   |   |
      12>-13>-14>-15>
"""
plaq16 = {'L': 16,
          'inner': {},
          'outer': {}}
plaq16['inner']['x_bonds'] = [[0,1], [1,2], [2,3],
                              [4,5], [5,6], [6,7],
                              [8,9], [9,10], [10,11],
                              [12,13], [13,14], [14,15]]
plaq16['inner']['y_bonds'] = [[0,4], [1,5], [2,6], [3,7],
                              [4,8], [5,9], [6,10], [7,11],
                              [8,12], [9,13], [10,14], [11,15]]
plaq16['outer']['x_bonds'] = [[3,0], [7,4], [11,8], [15,12]]
plaq16['outer']['y_bonds'] = [[12,0], [13,1], [14,2], [15,3]]