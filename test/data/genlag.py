# Useful script to generate tables of generalised Laguerre polynomials with SciPy

import os
import numpy as np
from scipy.special import genlaguerre

header = ''

x = np.linspace(0, 5)
lvals = []
for n in range(1, 3):
    for alpha in range(0, n):
        header += '# {0} {1}\n'.format(n, alpha)
        lp = genlaguerre(n, alpha)
        lvals.append(lp(x))

lvals = np.array(lvals)

fout = open('{0}.dat'.format(os.path.splitext(__file__)[0]), 'w')
fout.write(header)
fout.write('## END HEADER ## \n')
for i in range(0, len(x)):
    fout.write('{0} {1}\n'.format(x[i], ' '.join(map(str, lvals[:,i]))))
fout.close()


