/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * mudirac.cpp
 *
 * Main binary of the program
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "mudirac.hpp"

int main(int argc, char *argv[]) {
  string seed = "mudirac";
  MuDiracInputFile config;

  chrono::high_resolution_clock::time_point t0, t1;
  t0 = chrono::high_resolution_clock::now();

  if (argc < 2) {
    cout << "Input file missing\n";
    cout << "Please use the program as `mudirac <input_file>`\n";
    cout << "Quitting...\n";
    return -1;
  }

  seed = splitString(argv[1], ".")[0];
  try {
    config.parseFile(argv[1]);
  } catch (runtime_error e) {
    cout << "Invalid configuration file:\n";
    cout << e.what() << "\n";
    return -1;
  }

  int output_verbosity = config.getIntValue("output");

  // Set up logging
  AixLog::Severity log_verbosity;
  switch (config.getIntValue("verbosity")) {
    case 1:
      log_verbosity = AixLog::Severity::info;
      break;
    case 2:
      log_verbosity = AixLog::Severity::debug;
      break;
    case 3:
      log_verbosity = AixLog::Severity::trace;
      break;
    default:
      log_verbosity = AixLog::Severity::info;
      break;
  }

  // Sets up logging to the .log file according to user specified verbosity,
  // and writes errors out to the .err file
  AixLog::Log::init({ make_shared<AixLog::SinkFile>(log_verbosity, AixLog::Type::normal, seed + ".log"),
                      make_shared<AixLog::SinkFile>(AixLog::Severity::warning, AixLog::Type::special, seed + ".err") });

  LOG(INFO) << "MuDirac, a muonic atomic solver\n";
  LOG(INFO) << "by Simone Sturniolo\n";
  LOG(INFO) << "Released under the MIT License (2019)\n";
  LOG(INFO) << " \n";
  LOG(INFO) << "Please cite the following paper:\n";
  LOG(INFO) << "Sturniolo, S, Hillier, A., \"Mudirac: A Dirac equation solver for elemental analysis with muonic X‐rays.\"\n";
  LOG(INFO) << "X‐Ray Spectrom. 2020; 1– 17. https://doi.org/10.1002/xrs.3212\n";
  LOG(INFO) << " \n";

  // Are we running any debug tasks?
  string debugtask = config.getStringValue("devel_debug_task");

  if (debugtask == "EdEscan") {
    LOG(INFO) << "Running debug task: E=>dE scan\n";
    runEdEscan(config);
    return 0;
  }


  // Here we construct the atom and pass the transition data
  DiracAtom da = config.makeAtom();
  vector<TransLineSpec> transqnums;

  // Print out potential at high levels of verbosity
  if (output_verbosity >= 2 && (da.getPotentialFlags() && da.HAS_ELECTRONIC)) {
    writeEConfPotential(da.getPotentialElectronic(), seed + ".epot.dat");
  }

  // This contains both state names, the actual dirac states, and the transition matrix
  vector<TransitionData> transitions;

  // Check if we would like to optimise our fermi parameters using least squares
  // The idea is that we will calculate all transition energies and rates for many different pairs
  // of the fermi parameters (c,t)
  // We can then use this to perform least squares optimisation and finally obtain the rms nuclear radius
  if (config.getBoolValue("optimise_fermi_parameters")) {
    string coord_system_2pF;
    string min_2pF_algo;
    config.validateOptimisation(argc, coord_system_2pF, min_2pF_algo);
    double opt_time = 0;

    // try to read the experimental results file
    ExperimentalResultFile measurements;

    // switch to continue with optimisation if read is successful
    bool xr_read = false;
    try {
      measurements.parseFile(argv[2]);

    } catch (runtime_error e) {
      cout << "Invalid experimental measurements file:\n";
      cout << e.what() << "\n";
      return -1;
    }

    // read the measured transition lines
    measurements.validate(); 
    vector<string> xr_lines_measured = measurements.getStringValues("xr_lines");
    vector<double> xr_energies = measurements.getDoubleValues("xr_energy");
    vector<double> xr_errors = measurements.getDoubleValues("xr_error");
    LOG(INFO) << "Successfully read xray measurements input file \n";

    transqnums = measurements.parseXRLines();
    // data structure for storing best parameters.
    OptimisationData best_fermi_parameters;
    double MSE =0;

    // set all optimization values in Dirac Atom
    da.setExpOptData(coord_system_2pF, transqnums, xr_lines_measured, xr_energies, xr_errors, best_fermi_parameters, opt_time);

    // perform the optimisation routine
    optFermi2(da, min_2pF_algo, best_fermi_parameters, opt_time);

    // output file containing best fermi parameters and the associated MSE
    writeFermiParameters(da, best_fermi_parameters, opt_time,  seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
  
  }

  // Default mudirac behaviour

  //reset the transition quantum numbers to those from the config file
  // Here we read in the user specific transition lines and store a vector
  // containing the quantum numbers for each state in each transition
  
  transqnums = config.parseXRLines();
  da.transqnums = transqnums;
  transitions = da.getAllTransitions();

  // Sort transitions by energy if requested
  if (config.getBoolValue("sort_by_energy")) {
    sort(transitions.begin(), transitions.end(),
    [](TransitionData t1, TransitionData t2) {
      return (t1.ds2.E - t1.ds1.E) > (t2.ds2.E - t2.ds1.E);
    });
  }

  // Now create output files
  if (output_verbosity >= 1) {
    // Save a file for all lines
    ofstream out(seed + ".xr.out");
    out << "# Z = " << da.getZ() << ", A = " << da.getA() << " amu, m = " << da.getm() << " au\n";
    out << "Line\tDeltaE (eV)\tW_12 (s^-1)\n";
    out << fixed;

    if (config.getIntValue("xr_print_precision") > -1) {
      out << setprecision(config.getIntValue("xr_print_precision"));
    } else {
      out << setprecision(15); //Setting the maximum precision
    }

    for (int i = 0; i < transitions.size(); ++i) {
      double dE = (transitions[i].ds2.E - transitions[i].ds1.E);
      double tRate = transitions[i].tmat.totalRate();
      if (dE <= 0 || tRate <= 0)
        continue; // Transition is invisible
      out << transitions[i].name << '\t' << dE / Physical::eV;
      out << "\t\t" << tRate * Physical::s << '\n';
    }

    if (config.getBoolValue("write_spec")) {
      // Write a spectrum
      writeSimSpec(transitions, config.getDoubleValue("spec_step"), config.getDoubleValue("spec_linewidth"), config.getDoubleValue("spec_expdec"),
                   seed + ".spec.dat");
    }

    out.close();
  }

  if (output_verbosity >= 2) {
    vector<string> saved_states;
    // Save each individual state
    for (int i = 0; i < transitions.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
        DiracState ds = (j == 0 ? transitions[i].ds1 : transitions[i].ds2);
        string sname = ds.name();
        string fname = seed + "." + sname + ".out";

        if (vectorContains(saved_states, sname)) {
          continue;
        }

        LOG(DEBUG) << "Printing out state file for state " << sname << "\n";

        writeDiracState(ds, fname, config.getIntValue("state_print_precision"));

        saved_states.push_back(sname);
      }
      string fname = seed + "." + transitions[i].name + ".tmat.out";
      writeTransitionMatrix(transitions[i].tmat, fname);
    }
  }

  t1 = chrono::high_resolution_clock::now();

  LOG(INFO) << "Calculation completed in " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1.0e3 << " seconds\n";
}


