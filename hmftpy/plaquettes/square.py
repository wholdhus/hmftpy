"""
Connection diagram for 4-site square (symmetric tiling)

      <0>-<1>-<0>
       |   |   |
  <2>- 3 - 2 -<3>
   |   |   |   |
  <1>- 0 - 1 -<0>
       |   |   |
      <3>-<2>-<3>
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
            [2],
            [3],
            [0],
            [1]]
plaq4['inner']['n_bonds'] = [[0,1], [0,3], [2,1], [2,3],
                            ]
plaq4['inner']['nn_bonds'] = [[0,2], [1,3],
                             ]
plaq4['outer']['n_bonds'] = [[0,1], [0,3],
                             [1,0], [1,2],
                             [2,1], [2,3],
                             [3,0], [3,2]]
plaq4['outer']['nn_bonds'] = [[0,2], [1,3], [2,0], [3,1]]


"""
Connection diagram for 4-site square (symmetric tiling)

      <0>-<1>
       |   |
  <1>- 2 - 3 -<0>
   |   |   |   |
  <3>- 0 - 1 -<2>
       |   |
      <2>-<3>
"""
plaq4_os = {'L': 4,
          'inner': {},
          'outer': {}}
plaq4_os['inner']['nearest'] = [
            [1,2], # 0
            [0,3],
            [0,3], # 2
            [1,2]]
plaq4_os['inner']['n_nearest'] = [
            [3], # 0
            [2],
            [1], # 2
            [0]]
plaq4_os['outer']['nearest'] = [
            [2,3], # 0
            [2,3],
            [0,1], # 2
            [1,0]]
plaq4_os['outer']['n_nearest'] = [
            [1],
            [0],
            [3],
            [2]]