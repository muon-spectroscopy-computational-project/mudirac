import os
import sys
import json
import argparse
import tempfile
import warnings
import subprocess as sp

fdir = os.path.split(os.path.abspath(__file__))[0]

class MuDiracError(Exception):
    pass

# Convert a quantum number N to a spectral convention letter
def n2shell(n):
    return chr(ord('K')+n-1)

def write_mudirac_input(element, isotope, path='.', options={}):

    elconf = options.get('electronic_config')
    max_shell = options.get('max_shell', 6)
    max_calc_shell = options.get('max_calc_shell', max_shell)

    max_line = n2shell(max_shell) + str(2*max_shell-1)

    # Build the input file as a dictionary
    infile = {
        'element': element,
        'isotope': isotope,
        'nuclear_model': options.get('nuclear_model', 'FERMI2'),
        'uehling_correction': 'T' if options.get('uehling_correction', True) else 'F',
        'electronic_config': elconf if elconf else element,
        'xr_lines': 'K1:{0}-K1:{0}'.format(max_line),
        'ideal_atom_minshell': n2shell(max_calc_shell+1)
    }

    # Save the input file in the path
    with open(os.path.join(path, 'mudirac.in'), 'w') as f:
        for k, v in infile.items():
            f.write('{0}: {1}\n'.format(k, v))

def run_mudirac(command='mudirac', path='.'):
    
    proc = sp.Popen([command, 'mudirac.in'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, 
                    cwd=path)
    proc.communicate()

    outfile = os.path.join(path, 'mudirac.xr.out')

    # Read the lines
    if not os.path.isfile(outfile):
        raise MuDiracError('Calculation failed')

    log = open(os.path.join(path, 'mudirac.log')).read()
    err = open(os.path.join(path, 'mudirac.err')).read()

    results = []
    with open(outfile) as f:
        lines = f.readlines()
        for l in lines[2:]:
            xline, DE, W = l.strip().split()
            s1, s2 = xline.split('-')
            results.append({
                'from': s1,
                'to': s2,
                'DE': float(DE),
                'W': float(W)
            })

    return results, log, err
            

parser = argparse.ArgumentParser(description='Compile a text file database of MuDirac X-Ray lines')
parser.add_argument('--command', type=str, default='mudirac',
                    help='MuDirac running command')
parser.add_argument('--elemlist', type=str, default=os.path.join(fdir, 'all_elements.dat'),
                    help='The path to a file containing a list of [Symbol, Isotope] pairs to use as list of elements.')
parser.add_argument('--nucmodel', type=str, default='FERMI2', choices=['POINT', 'SPHERE', 'FERMI2'], 
                    help='The nuclear charge density model to use for the atom')
parser.add_argument('--no-uehling', action='store_true', 
                    help='Turn off Uehling correction (speeds up considerably, not recommended for heavy atoms)')
parser.add_argument('--elec-config', type=str, default=None,
                    help='Electronic config to use for all atoms')
parser.add_argument('--max-shell', type=int, default=6,
                    help='Quantum number n of the highest shell to consider for transitions')
parser.add_argument('--max-calc-shell', type=int, default=6,
                    help='Quantum number n of the highest shell to compute numerically (idealised '
                    'hydrogen-like solutions will be used for all higher shells)')
parser.add_argument('--min-DE', type=float, default=1e-4,
                    help='Minimum energy difference (in eV) between two states to consider the line worth storing')
parser.add_argument('--outfile', '-o', type=str, default=None,
                    help='Name of the file to save the database to')
parser.add_argument('--save-logs', action='store_true', 
                    help='Save all .log and .err files for each calculation run')

args = parser.parse_args()

# Compile options
options = {
    'nuclear_model': args.nucmodel,
    'uehling_correction': not args.no_uehling,
    'electronic_config': args.elec_config,
    'max_shell': args.max_shell,
    'max_calc_shell': args.max_calc_shell
}
command = args.command

# Grab the element-isotope pairs
to_calc = []
with open(args.elemlist) as f:
    for l in f.readlines():
        el, A = l.strip().split()[:2]
        to_calc.append((el, A))

N = len(to_calc)
print('{0} isotopes queued')
print('Starting MuDirac calculations...')
db_entries = {}
for i, (el, iso) in enumerate(to_calc):

    sys.stdout.write('\rCalculating: {0}/{1}                 '.format(i+1, len(to_calc)))

    with tempfile.TemporaryDirectory() as tmpdir:

        write_mudirac_input(el, iso, tmpdir, options)
        try:
            results, log, err = run_mudirac(command, tmpdir)
        except MuDiracError:
            warnings.warn('Calculation failed for {0}{1}'.format(iso, el))

        if args.save_logs:
            with open('{0}-{1}.log'.format(el, iso), 'w') as f:
                f.write(log)
            with open('{0}-{1}.err'.format(el, iso), 'w') as f:
                f.write(err)

        # Filter out the lines
        results = [r for r in results if r['DE'] > args.min_DE]

        db_entries[el] = db_entries.get(el, {})
        db_entries[el][iso] = results

db_string = json.dumps(db_entries, indent=2)

if args.outfile:
    with open(args.outfile, 'w') as f:
        f.write(db_string)
else:
    print(db_string)