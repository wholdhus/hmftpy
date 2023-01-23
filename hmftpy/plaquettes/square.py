import numpy as np

a0 = np.zeros(2)
ax = np.array([0,1])
ay = np.array([1,0])


"""
Connection diagram for 4-site square (symmetric tiling)

  <1>-<0>-<1>-<0>
   |   |   |   |
  <3>- 2 - 3 -<2>
   |   |   |   |
  <1>- 0 - 1 -<0>
   |   |   |   |
  <3>-<2>-<3>-<2>
"""
plaq4 = {'L': 4,
         'inner': {},
         'outer': {}}
plaq4['rs'] = [a0, ax, ay, ax+ay]
plaq4['Rs'] = [2*ax, 2*ay]
plaq4['inner']['n_bonds'] = [[0,1], [2,3], [0,2], [1,3]]
plaq4['inner']['nn_bonds'] = [[0,3], [1,2]]
plaq4['outer']['n_bonds'] = [[1,0], [3,2], [2,0], [3,1]]
plaq4['outer']['nn_bonds'] = [[1,2], [3,0], [2,1],
                              [0,3], [2,1], [3,0]]


"""
Connection diagram for 9-site square (symmetric tiling)
  <2>-<0>-<1>-<2>-<0>
   |   |   |   |   |
  <8>- 6 - 7 - 8 -<6>
   |   |   |   |   |
  <5>- 3 - 4 - 5 -<3>
   |   |   |   |   |
  <2>- 0 - 1 - 2 -<0>
   |   |   |   |   |
  <8>-<6>-<7>-<8>-<6>
"""
plaq9 = {'L': 9,
          'inner': {},
          'outer': {}}
plaq9['rs'] = [a0, ax, 2*ax,
               ay, ax+ay, 2*ax+ay,
               2*ay, ax+2*ay, 2*ax+2*ay]
plaq9['rs'] = [3*ax, 3*ay]
plaq9['inner']['n_bonds'] = [[0,1], [1,2], [3,4], [4,5], [6,7], [7,8],
                             [0,3], [1,4], [2,5], [3,6], [4,7], [5,8]]
plaq9['inner']['nn_bonds'] = [[0,4], [1,5], [3,7], [4,8],
                              [1,3], [2,4], [4,6], [5,7]]
plaq9['outer']['n_bonds'] = [[2,0], [5,3], [8,6], [6,0], [7,1], [8,2]]
plaq9['outer']['nn_bonds'] = [[2,3], [5,6], [6,1], [7,2], [8,0],
                              [0,5], [3,8], [6,2], [7,0], [8,1]]


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
  <3>- 0 - 1 - 2 - 3 -<0>
       |   |   |   |
      12>-13>-14>-15>
"""
plaq16 = {'L': 16,
          'inner': {},
          'outer': {}}
plaq16['rs'] = [a0, ax, 2*ax, 3*ax,
                ay, ax+ay, 2*ax+ay, 3*ax+ay,
                2*ay, ax+2*ay, 2*(ax+ay), 3*ax+2*ay,
                3*ay, ax+3*ay, 2*ax+3*ay, 3*(ax+ay)]
plaq16['Rs'] = [4*ax, 4*ay]
plaq16['inner']['n_bonds'] = [[0,1], [1,2], [2,3], [4,5], [5,6], [6,7],
                              [8,9], [9,10], [10,11], [12,13], [13,14], [14,15],
                              [0,4], [1,5], [2,6], [3,7], [4,8], [5,9], [6,10], [7,11],
                              [8,12], [9,13], [10,14], [11,15]]
plaq16['inner']['nn_bonds'] = [[0,5], [1,6], [2,7], [4,9], [5,10], [6,11],
                               [8,13], [9,14], [10,15],
                               [1,4], [2,5], [3,6], [5,8], [6,9], [7,10],
                               [9,12], [10,13], [11,14]]
plaq16['outer']['n_bonds'] = [[3,0], [7,4], [11,8], [15,12],
                              [12,0], [13,1], [14,2], [15,3]]
plaq16['outer']['nn_bonds'] = [[3,4], [7,8], [11,12], [12,1], [13,2], [14,3], [15,0],
                               [0,7], [4,11], [8,15], [12,3], [13,0], [14,1], [15,2]]

if __name__ == '__main__':
    from util import test_bonds
    print('Checking plaq4')
    print(test_bonds(plaq4, [4,4]))
    print('Checking plaq9')
    print(test_bonds(plaq9, [4,4]))
    print('Checking plaq16')
    print(test_bonds(plaq16, [4,4]))