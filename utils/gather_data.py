#!/usr/bin/env python

"""
A utility script meant to download and parse atomic mass and spin data from
official sources (AME 2016 and NUBASE) and create header files to store
this data for the main program.
"""

import os
import re
import urllib.request

dir = os.path.split(__file__)[0]


def get_file(url):
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    return data.decode('utf8')

# Load isotopic abundance data
iso_abund = {}
with open(os.path.join(dir, 'abundant.dat')) as f:
    for l in f.readlines():
        if l[0] == '#':
            continue
        Z, A = map(int, l.split())
        iso_abund[Z] = A

# Load nuclear radii
iso_radii = {}
with open(os.path.join(dir, 'nuclear_radii.dat')) as f:
    for l in f.readlines():
        if l[0] == '#':
            continue
        Z, A = map(int, l.split()[:2])
        iso_radii[(Z, A)] = float(l.split()[-1])

ame_url = 'https://www-nds.iaea.org/amdc/ame2016/mass16.txt'
nubase_url = 'https://www-nds.iaea.org/amdc/ame2016/nubase2016.txt'

# Start with ame
try:
    ame_path = os.path.join(dir, 'ame2016.dat')
    ame_txt = open(ame_path).read()
except IOError:
    ame_txt = get_file(ame_url)
    with open(ame_path, 'w') as f:
        f.write(ame_txt)
ame_lines = ame_txt.split('\n')
# Discard header
for i, l in enumerate(ame_lines):
    if l[:4] == '1N-Z':
        ame_lines = ame_lines[i+2:]
        break

atomic_data = {}

for l in ame_lines:
    lspl = l.split()
    if len(lspl) == 0:
        break
    # Skip the initial 0 that groups by A
    if l[0] == '0':
        lspl = lspl[1:]
    # Skip the decay type
    try:
        _ = float(lspl[5])
    except ValueError:
        del(lspl[5])
    Z = int(lspl[2])
    A = int(lspl[3])
    if ((int(lspl[1]) != A-Z) or int(lspl[0]) != A-2*Z):
        raise RuntimeError('Invalid line in ame2016 data')
    el = lspl[4]
    m = float(lspl[-3])+float(lspl[-2].replace('#', ''))*1e-6
    if abs(m-A) > 0.5:
        raise RuntimeError('Invalid mass in ame2016 data')

    # Skip the neutron, we don't need it
    if el == 'n':
        continue

    atomic_data[el] = atomic_data.get(el, {'Z': Z, 'isos': {}})
    atomic_data[el]['isos'][A] = {'m': m, 'spin': None}

# Now nubase
try:
    nubase_path = os.path.join(dir, 'nubase2016.dat')
    nubase_txt = open(nubase_path).read()
except IOError:
    nubase_txt = get_file(nubase_url)
    with open(nubase_path, 'w') as f:
        f.write(nubase_txt)

elems = sorted(atomic_data.keys(), key=lambda x: len(x), reverse=True)

el_re = re.compile('([0-9]+)({0})'.format('|'.join(elems)))
spin_re = re.compile(r'([0-9]/)*([0-9])(\+|-)')
nubase_lines = nubase_txt.split('\n')[1:]  # First line is neutron
for l in nubase_lines:
    if len(l) == 0:
        break
    lspl = l.split()
    A, el = el_re.match(lspl[2]).groups()
    A = int(A)
    if A != int(lspl[0]):
        raise RuntimeError('Invalid line in NUBASE data')
    if not (A in atomic_data[el]['isos']):
        print(
            'Isotope {0}{1} in NUBASE data has no ame2016 information; skipping.'.format(A, el))
    try:
        spin = spin_re.findall(l)[0]
        spin = (float(spin[1]) if spin[0] == '' else float(
            spin[0][:-1])/float(spin[1]))*(1 if spin[2] == '+' else -1)
        atomic_data[el]['isos'][A]['spin'] = spin
    except IndexError:
        continue

# Compile that all as a C++ initialiser string
cpp_string = 'const map<string, element> atomic_data = {'
cpp_Z_string = 'const map<int, string> atomic_Z_lookup = {'
for el, data in atomic_data.items():
    Z = data['Z']
    cpp_string += '{{ "{0}", {{ {1}, {2}, {{\n'.format(el, Z, iso_abund[Z])
    for A, iso in data['isos'].items():
        cpp_string += '{{ {0}, {{ {m}, {spin}, {R} }} }},\n'.format(A, m=iso['m'],
                                                             spin=(iso['spin'] if iso['spin']
                                                                   is not None else 'NAN'), 
                                                             R=iso_radii.get((Z, A), 'NAN'))
    cpp_string += '} } },\n'

    cpp_Z_string += '{{ {0}, "{1}" }},'.format(data['Z'], el)
cpp_string += '};\n'
cpp_Z_string += '};'

with open(os.path.join(dir, 'elements.in.hpp')) as f:
    ftxt = f.read()
    ftxt = ftxt.replace('//{HERE GOES THE ACTUAL DATA}//',
                        cpp_string + cpp_Z_string)
    fout = open(os.path.join(dir, '../lib/elements.hpp'), 'w')
    fout.write(ftxt)
    fout.close()
