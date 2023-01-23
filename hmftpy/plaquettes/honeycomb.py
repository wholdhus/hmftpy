"""
6-site hexagonal cluster
       <1>     <0>
         \     /
          4 - 5
         /     \
   <3>- 2       3 -<2>
         \     /
          0 - 1
         /     \
       <5>     <4>
"""
plaq6 = {'L': 6, 'inner': {}, 'outer': {}}
plaq6['inner']['n_bonds'] = [[0,1], [4,5],
                             [1,3], [2,4],
                             [0,2], [3,5]]
plaq6['inner']['nn_bonds'] = [[0,3], [2,5],
                              [0,4], [1,5],
                              [1,2], [3,4]]
plaq6['outer']['n_bonds'] = [[3,2], [5,0], [4,1]]
plaq6['outer']['nn_bonds'] = [[1,2], [3,4], [4,0], [5,1],
                              [2,1], [3,0], [4,3], [5,2],
                              [0,3], [2,5], [4,0], [5,1]]

if __name__ == '__main__':
    from util import test_bonds
    print('Checking plaq6')
    print(test_bonds(plaq6, [3,6]))