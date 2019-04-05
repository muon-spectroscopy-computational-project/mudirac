# Useful script to generate tables of generalised Laguerre polynomials with SciPy

import os
import numpy as np
from scipy.special import genlaguerre

lag_params = '{'

min_n = 1
max_n = 3
x = np.linspace(0, 5)

lvals = []
for n in range(min_n, max_n):
    for alpha in range(0, n):
        lag_params += '{{ {0}, {1} }},\n'.format(n, alpha)
        lp = genlaguerre(n, alpha)
        lvals.append(lp(x))

lag_params += '}'
lag_values = '{'

lvals = np.array(lvals)

for i in range(0, len(x)):
    lag_values += '{{ {0}, {1} }},\n'.format(x[i], ', '.join(map(str,
                                                                 lvals[:, i])))
lag_values += '}'

path = os.path.split(__file__)[0]
with open(os.path.join(path, 'genlag.h'), 'w') as f:
    f.write('int laguerreParams[{0}][2] = '.format(
        lvals.shape[0]) + lag_params + ';\n')
    f.write('double laguerreValues[{0}][{1}] = '.format(
        lvals.shape[1], lvals.shape[0]+1) + lag_values + ';')
