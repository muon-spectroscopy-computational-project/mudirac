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

  // Here we read in the user specific transition lines and return a vector
  // containing the quantum numbers for each state in each transition
  vector<TransLineSpec> transqnums = config.parseXRLines();

  // Here we construct the atom
  DiracAtom da = config.makeAtom();

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

    // check the experimental results input file is provided
    if (argc < 3) {
      cout << "Experimental results input file missing\n";
      cout << "When optimise_fermi_parameters is True, an additional experimental results file is expected\n";
      cout << "Please use the program as `mudirac <input_file> <experimental_results_input_file>`\n";
      cout << "If experimental results cannot be provided, optimise_fermi_parameters should be set to False \n";
      cout << "Quitting...\n";
      return -1;
    }

    // check the nuclear model is suitable for optimisation
    if (config.getStringValue("nuclear_model") != "FERMI2") {
      cout << "nuclear model parameters can only be optimised for the 2 parameter Fermi model\n";
      cout << "Please add the line `nuclear_model: FERMI2` to your first input file\n";
      cout << "Quitting...\n";
      return -1;
    }

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
    xr_read = measurements.validate();

    if (xr_read) {
      
      LOG(INFO) << "Successfully read xray measurements input file \n";
      vector<string> xr_lines_measured = measurements.getStringValues("xr_lines");
      vector<double> xr_energies = measurements.getDoubleValues("xr_energy");
      vector<double> xr_errors = measurements.getDoubleValues("xr_error");
      
      // data structure for storing best parameters.
      OptimisationData best_fermi_parameters;
      double MSE =0;

      string coords;
      // get the 2pF optimsation coordinate system (dev)
      if (da.getA() < 5) {
        LOG(INFO) << "using ct coordinate system as polar parameterisation no longer holds for A < 5 \n";
        coords = "ct";
      }
      else {
        coords = config.getStringValue("2pF_coords");
      }
      
      const string coord_system_2pF = coords;
    
      if (!((coord_system_2pF == "ct")||(coord_system_2pF == "polar"))){
        cout << "Invalid 2pF coordinate system choice for minimsation\n";
        cout << "please use \"ct\" or \"polar\" (default is \"polar\") \n";
        cout << "You used: \""<<coord_system_2pF<<"\" \n";
        cout << "Quitting...\n";
        LOG(ERROR) << "Invalid 2pF coordinate system choice for minimsation\n";
        return -1;

      }
      
      // implement the choice of minimisation algorithm
      const string min_2pF_algo = config.getStringValue("min_2pF_algorithm");
      double opt_time = 0;

      // set all optimization values in Dirac Atom
      da.setExpOptData(coord_system_2pF, transqnums, xr_lines_measured, xr_energies, xr_errors, best_fermi_parameters, opt_time);

      if (min_2pF_algo == "bfgs"){
        optimizeFermiParameters(da, best_fermi_parameters, opt_time);
      }
      else if (min_2pF_algo =="trust"){
        opt_2pF_model opt_obj(da);
        optimizeFermiParameters(opt_obj, da, best_fermi_parameters, opt_time);
      }
      else if (min_2pF_algo=="global"){
        globalOptimizeFermiParameters(da, best_fermi_parameters, opt_time);
      }
      else {
        cout << "Invalid 2pF optimisation algorithm choice for minimsation\n";
        cout << "please use \"bfgs\" or \"trust\" (default is \"bfgs\") \n";
        cout << "You used: \""<<min_2pF_algo<<"\" \n";
        cout << "Quitting...\n";
        LOG(ERROR) << "Invalid 2pF optimisation algorithm choice for minimsation: \""<<min_2pF_algo<<"\"\n";
        return -1;
      }
  
      // output file containing best fermi parameters and the associated MSE
      writeFermiParameters(da, best_fermi_parameters, opt_time,  seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
    }

  }

  // Default mudirac behaviour
  // Wrapped the calculation of the states, their energies and the transition probabilities into here,
  // so that we can easily loop over it for least squares optimisation
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


