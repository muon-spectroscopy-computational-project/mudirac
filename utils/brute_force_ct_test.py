import os
import argparse
import tempfile
import subprocess
import itertools
import sys
import numpy as np
from pathlib import Path

EXPERIMENTAL_DATA = {
    "element": "C",
    "isotope": 12,
    "transitions": [
        ("K1-L2", 75248.0, 15000.0),
        ("K1-L3", 89212.0, 15000.0),
        ("K1-N2", 94025.0, 15000.0),
    ]
}

NUM_SAMPLE_POINTS = 3

def  generate_energy_samples(energy,error, num_points):
    return np.linspace(energy- error, energy + error, num_points)

def write_main_config(path, element, isotope, transitions):
    xr_lines = ", ".join([t[0] for t in transitions])

    config = f"""#Mudirac main configuration for {element} {isotope}
element: {element}
isotope: {isotope}
nuclear_model: FERMI2
uehling_correction: TRUE
reduced_mass: TRUE
electronic_config: {element}
optimise_fermi_parameters: TRUE
xr_lines: {xr_lines}
verbosity: 1
output: 1
2pF_coords: CT
min_2pF_algorithm: lm
"""
    config_path = os.path.join(path, "mudirac.in")
    with open(config_path, "w") as f:
        f.write(config)
    return config_path

def write_experimentl_config(path, transitions,energies):
    xr_lines = ", ".join(transitions)
    xr_energies = ", ".join([f"{energies[t][0]:.6f}" for t in transitions])
    xr_errors = ", ".join([f"{energies[t][1]:.6f}" for t in transitions])
    config = f""" Experimental mesure for optimpistaion brute force
xr_lines: {xr_lines}
xr_energies: {xr_energies}
xr_errors: {xr_errors}
"""
    config_path = os.path.join(path, "experimental.in")
    with open(config_path, "w") as f:
        f.write(config)
    return config_path

def run_mudirac(command, work_dir):
    try:
        proc =  subprocess.run(
        [command, "mudirac.in", "experimental.in"],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=600
        )

        #Parse the ouput
        fermi_file = os.path.join(work_dir, "mudiracfermi_parameters.out")
        if not os.path.exists(fermi_file):
            print("Fermi parameters file not found.")
            return None
        
        with open(fermi_file, "r") as f:
            lines = f.readlines()
        
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) >= 5:
                return {
                    "fermi_c": float(parts[0]),
                    "fermi_t": float(parts[1]),
                    "rms_radius": float(parts[2]),
                    "theta": float(parts[3]),
                    "mse": float(parts[4])
                }
        return None
    except subprocess.TimeoutExpired:
        print("Mudirac run timed out.")
        return None
    except Exception as e:
        print(f"An error occurred while running Mudirac: {e}")
        return None
    

def run_brute_force(mudirac_cmd, exp_data, num_points, output_file=None):
    transitions = exp_data['transitions']
    element = exp_data['element']
    isotope = exp_data['isotope']

    energy_samples = {}

    for name , central, error in transitions:
        energy_samples[name] = generate_energy_samples(central, error, num_points)
    print(f"Generated {num_points} sample points for each transition.")
    print(f" Samples: {energy_samples[name]}")

    #create combinations
    transition_names = [t[0] for t in transitions]
    transition_errors = {t[0]: t[2] for t in transitions}
    all_combinations =list(itertools.product(*[energy_samples[name] for name in transition_names]))

    total_runs = len(all_combinations)
    print(f"Total combinations to test: {total_runs}")
    print("="*50)

    results = []

    for i, energy_combo in enumerate(all_combinations):
        
        energies = {
        name: (energy_combo, transition_errors[name])
        for name, energy in zip(transition_names, energy_combo)
        }

        print

        with tempfile.TemporaryDirectory() as temp_dir:
            write_main_config(temp_dir, element, isotope, transitions)
            write_experimentl_config(temp_dir, transition_names, energy_samples)

            result =  run_mudirac(mudirac_cmd, temp_dir)
            if result:
                result["input_energies"] = dict(energies)
                result["energy_deltas"] = {
                    name: energy_combo[j] - transitions[j][1]
                    for j, name in enumerate(transition_names)
                }
                results.append(result)
                print(f"   -> c = {result['fermi_c']:.6f}, t = {result['fermi_t']:.6f}, rms = {result['rms_radius']:.6f}, mse = {result['mse']:.6f}")
            else:
                print("   -> Mudirac run failed or produced no output.")

    print("\n" + "="*50)
    print("SUMMARY OF RESULTS")
    print("="*50)
    if results:
        c_values = [r['fermi_c'] for r in results]
        t_values = [r['fermi_t'] for r in results]
        rms_values = [r['rms_radius'] for r in results]
        mse_values = [r['mse'] for r in results]

        print(f"|\n Successful runs: {len(results)} / {total_runs}")
        print(f"\n Fermu c parameter:")
        print(f" Range : {min(c_values):.6f} - {max(c_values):.6f} fm")
        print(f" Mean  : {np.mean(c_values):.6f} fm")
        print(f" Std   : {np.std(c_values):.6f} fm")    

        print(f"\n Fermi t parameter:")
        print(f" Range : {min(t_values):.6f} - {max(t_values):.6f} fm")
        print(f" Mean  : {np.mean(t_values):.6f} fm")
        print(f" Std   : {np.std(t_values):.6f} fm")    

        print(f"\n RMS radius parameter:")
        print(f" Range : {min(rms_values):.6f} - {max(rms_values):.6f} fm")
        print(f" Mean  : {np.mean(rms_values):.6f} fm")
        print(f" Std   : {np.std(rms_values):.6f} fm")    

        print(f"\n MSE parameter:")
        print(f" Range : {min(mse_values):.6f} - {max(mse_values):.6f}")
        print(f" Mean  : {np.mean(mse_values):.6f}")
        print(f" Std   : {np.std(mse_values):.6f}")    

    if output_file:
        save_results(results, transitions, output_file)
        print(f"\nResults saved to {output_file}")

    return results


def save_results(results, transitions, output_file):
        transition_names = [t[0] for t in transitions]
        with open(output_file, "w") as f:
            #Header
            header = ["run"]
            for name in transition_names:
                header.append(f"E_{name}")
                header.append(f"delta_E_{name}")
            header.extend(["fermi_c", "fermi_t", "rms_radius", "mse"])
            f.write(",".join(header) + "\n")

            #data rows
            for i, res in enumerate(results):
                row = [str(i+1)]
                for name in transition_names:
                    row.append(f"{res['input_energies'][name][0]:.6f}")
                    row.append(f"{res['energy_deltas'][name]:.6f}")
                row.extend([
                    f"{res['fermi_c']:.6f}",
                    f"{res['fermi_t']:.6f}",
                    f"{res['rms_radius']:.6f}",
                    f"{res['mse']:.6f}"
                ])
                f.write(",".join(row) + "\n")   
    
def main():
    parser = argparse.ArgumentParser(description="Brute-force optimization of Fermi parameters using Mudirac.")
    parser.add_argument("--mudirac", "-m", required=True, default="mudirac", help="Path to the Mudirac executable.")
    parser.add_argument("--points", "-n", type=int, default=NUM_SAMPLE_POINTS, help="Number of sample points per transition.")
    parser.add_argument("--output", "-o", default="brute_force_mudirac_results.csv", help="Output CSV file to save results.")
    args = parser.parse_args()

    print("="*50)
    print("Mudirac Brute-Force Fermi Parameter Optimization test")
    print("="*50)
    print(f"\n Element:", EXPERIMENTAL_DATA['element']+"-"+str(EXPERIMENTAL_DATA['isotope']))
    print(f"Transitionss: {len(EXPERIMENTAL_DATA['transitions'])}")
    print(f" Sample points per transition: {args.points}")
    print(f" Total combinations to test: {args.points ** len(EXPERIMENTAL_DATA['transitions'])}")
    print("="*50 + "\n")    

    results = run_brute_force(
        mudirac_cmd=args.mudirac,
        exp_data=EXPERIMENTAL_DATA,
        num_points=args.points,
        output_file=args.output
    )
    return 0 if results else 1

if __name__ == "__main__":
    sys.exit(main())

            