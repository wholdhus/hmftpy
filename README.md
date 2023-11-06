# hmftpy: Hierarchical mean-field theory implemented in Python
Implements hierarchical mean-field theory, a cluster mean-field theory that systematically preserves
and breaks symmetries of the Hamiltonian to understand the phase diagrams of strongly-correlated models.

This package relies on `QuSpin`, an exact-diagonalization package. To install `QuSpin`, make sure you have the
`conda` package manager installed and then use command
    conda install -c weinbe58 quspin

## Usage
For example usage, see the Jupyter notebook `demo.ipynb`.

## Input/output dictionary formats
Functions in this package use the following formats for input and output dictionaries:

    cluster = {'L': # of sites in cluster,
                 'inner': {
                    'nearest': [[intra-cluster nearest neighbors of site i]
                                for i in cluster],
                    'n_nearest': [[intra-cluster next-nearest neighbors of site i]
                                  for i in cluster],
                    'n_n_nearest': [[intra-cluster next-next-nearest neighbors of site i]
                                    for i in cluster]}
                 'outer': {
                    'nearest': [[inter-cluster nearest neighbors of site i]
                                for i in cluster],
                    'n_nearest': [[inter-cluster next-nearest neighbors of site i]
                                  for i in cluster],
                    'n_n_nearest': [[inter-cluster next-next-nearest neighbors of site i]
                                    for i in cluster]}
                 }

    interactions = {'local': {'z': -2},
                    'nearest': {'xx': 1, 'yy': 1},
                    'n_nearest': {'xy': 1, 'yx': -1},
                    'n_n_nearest': {'yx': -1, 'yx', -1}}

    mean_fields = {'x': [List of <sigma_i^x> for i in cluster],
                   'y': [List of <sigma_i^y> for i in cluster],
                   'z': [List of <sigma_i^z> for i in cluster]}

    coeffs = {'inner': {local': {'z': 1-D NUMPY ARRAY OF LENGTH L},
                        'nearest': {'xx': 2-D NUMPY ARRAY OF LENGTH L, 'yy': ...},
                        'n_nearest': {'xy': ..., 'yx': ...},
                        'n_n_nearest': {'yx': ..., 'yx', ...}}
              'outer': {local': {'z': ...},
                        'nearest': {'xx': ..., 'yy': ...},
                        'n_nearest': {'xy': ..., 'yx': ...},
                        'n_n_nearest': {'yx': ..., 'yx', ...}}}
