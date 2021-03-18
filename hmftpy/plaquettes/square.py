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
plaq4['inner']['nearest'] = [
            [1,2], # 0
            [0,3],
            [0,3], # 2
            [1,2]]
plaq4['inner']['n_nearest'] = [
            [3], # 0
            [2],
            [1], # 2
            [0]]
plaq4['outer']['nearest'] = [
            [1,2], # 0
            [0,3],
            [0,3], # 2
            [1,2]] # 8
plaq4['outer']['n_nearest'] = [
            [3],
            [2],
            [1],
            [0]]
