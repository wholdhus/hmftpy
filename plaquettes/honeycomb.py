"""
Creating lists of the bonds of various types
in the Kitaev model. Each entry is a list of the indices of the
two spins connected by the bond.
"""

"""
Connection diagram for "small" (6 site)

         <0>
          |
  <5>     3     <1>
    \   /   \   /
      2       4
      |       |
      1       5
    /   \   /   \
  <4>     0     <2>
          |
         <3>
"""


# Inter couplings: both are real spins
# Outer couplings: 1st is mean field param. index, second is real spin
plaq6 = {'L': 6,
              'inter_x': [[5,0], [2,3]], # , [0,5], [3,2]],
              'inter_y': [[0,1], [3,4]], #, [1,0], [4,3]],
              'inter_z': [[1,2], [4,5]], # , [2,1], [5,4]],
              'outer_x': [[1,4], [4,1]],
              'outer_y': [[2,5], [5,2]],
              'outer_z': [[0,3], [3,0]]
              }
"""
Connection diagram for "Large" (18 site)

                     <0>
                      |
              <12>    16    <10>
                \   /   \   /
         <11>     17      15
          |       |       |
  <14>    3       6       14
    \   /   \   /   \   /   \
      2       4       7       <2>
      |       |       |
      1       5       8       <1>
    /   \   /   \   /   \   /
  <13>    0       9       13
          |       |       |
         <16>     10      12
                /   \   /   \
              <15>    11    <17>
                      |
                     <3>
sps_even = [[0,2,4],[4,7,9],[7,17,15],[9,11,13]]
sps_odd = [[1,3,5], [5,6,8],[10,8,12],[6,16,14]]
"""
plaq18 = {'L': 18,
              'inter_x': [[5,0], [2,3],
                          [4,6], [8,9],
                          [11, 12],
                          [7,14], [17,16]
                          ],
              'inter_y': [[0,1], [3,4],
                          [6,7], [9,5],
                          [8,13],[11,10],
                          [16,15]
                          ],
              'inter_z': [[1,2], [4,5],
                          [7,8],
                          [10,9], [13,12],
                          [6,17], [15,14]
                          ],
              'outer_x': [[1,13], [13,1],
                          [10,15], [15,10]],
              'outer_y': [[2,14], [14,2],
                          [12,17],[17,12]],
              'outer_z': [[3,11], [11,3],
                          [0,16], [16,0]]
              }
"""
Connection diagram for "Largest" (24 site)

              |       |
       <10>   4       6    <16>
         \  /   \   /  \   /
          3       5      7
          |       |      |
  <12>    2      20      8    <0>
    \   /   \   /   \   /  \  /
      1       19      21    9
      |       |       |     |
      0       18      22    10
    /   \   /   \   /   \  /  \
  <9>     17     23      11    <3>
          |      |       |
          16     14      12
         /   \  /  \   /   \
       <7>    15      13     <1>
              |       |

sps_even = [[0,2,18], [18,20,22], [22,8,10], [16,18,14], [14,22,12], [2,4,20], [20,6,8]]
sps_odd = [[0,2,18], [18,20,22], [22,8,10], [16,18,14], [14,22,12], [2,4,20], [20,6,8]]

"""

plaq24 = {'L': 24,
          'inter_x': [[1,2], [3,4], # outer ring
                      [5,6], [11,10], # outer
                      [13,12], [15,14], # outer
                      [19,20], [23,22], # inner hex
                      [21,8], [17,18]
                      ],
          'inter_y': [[4,5], [6,7],
                      [8,9], [14,13],
                      [16,15],[0,17],
                      [20,21], [18, 23],
                      [2, 19], [22, 11]
                      ],
          'inter_z': [[0,1], [2,3],
                      [8,7], [10,9],
                      [12,11], [16,17],
                      [18, 19], [22, 21],
                      [20, 5], [14, 23]
                      ],
          'outer_x': [[9,0], [0,9],
                      [7,16], [16,7]],
          'outer_y': [[3,10], [10,3],
                      [12,1], [1,12]],
          'outer_z': [[4,15], [15,4],
                      [6,13], [13,6]]
          }
