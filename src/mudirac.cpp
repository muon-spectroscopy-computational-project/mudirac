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
  vector<TransLineSpec> transqnums = parseXRLines(config);

  // We have to decide on how to obtain the experimental energies
  // One suggestion I have is to define a new input file which will look something like:
  //
  // transition energy energy uncertainty
  // i.e
  // K1-L2 300000 1
  // K1-L3 500000 10

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
  if (config.getBoolValue("optimise_fermi_parameters") && config.getStringValue("nuclear_model") == "FERMI2") {

    // variables to store experimental values passed by reference to readXrayMeasurements function
    // experimental measurements input file object
    ExperimentalResultFile measurements;

    // switch to continue with optimisation if read is successful
    bool xr_measurement_read_success = false;

    // vectors to store the experimental data
    vector<string> xr_lines_measured;
    vector<double> xr_energies;
    vector<double> xr_errors;

    // read in experimental results file and populate the associated vectors
    try {
      readXrayMeasurements(measurements, xr_measurement_read_success,  xr_lines_measured,  xr_energies, xr_errors, argv);
    } catch (runtime_error e) {
      cout << "Invalid experimental measurements file:\n";
      cout << e.what() << "\n";
      return -1;
    }


    if (xr_measurement_read_success) {
      LOG(INFO) << "Successfully read xray measurements input file \n";
      // data structure for tracking best parameters.
      OptimisationData optimal_fermi_parameter;
      optimal_fermi_parameter.mse = 1.0;
      vector<OptimisationData> valid_fermi_parameters;

      // find the valid and optimal fermi parameters
      optimiseFermiParameters(xr_lines_measured, xr_energies, xr_errors, transqnums, valid_fermi_parameters, optimal_fermi_parameter, transitions, config);
      // transitions updated by reference for output later

      // output file containing all valid fermi parameters and the associated MSE
      writeFermiParameters(da, valid_fermi_parameters, seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
    }

  } else {
    // Default mudirac behaviour
    // Wrapped the calculation of the states, their energies and the transition probabilities into here,
    // so that we can easily loop over it for least squares optimisation
    transitions = getAllTransitions(transqnums, da);
  }

  // Sort transitions by energy if requested
  if (config.getBoolValue("sort_byE")) {
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

// Calculate the states and transition probabilities
vector<TransitionData> getAllTransitions(vector<TransLineSpec> transqnums, DiracAtom da) {

  vector<TransitionData> transitions;
  vector<string> failconv_states; // Store states whose convergence has failed already, so we don't bother any more
  for (int i = 0; i < transqnums.size(); ++i) {
    int n1, l1, n2, l2;
    bool s1, s2;
    bool success = true;
    TransitionData tdata;

    n1 = transqnums[i].n1;
    l1 = transqnums[i].l1;
    s1 = transqnums[i].s1;
    n2 = transqnums[i].n2;
    l2 = transqnums[i].l2;
    s2 = transqnums[i].s2;

    tdata.sname1 = printIupacState(n1, l1, s1);
    tdata.sname2 = printIupacState(n2, l2, s2);
    tdata.name = tdata.sname1 + "-" + tdata.sname2;

    // Have these been tried before?
    if (vectorContains(failconv_states, tdata.sname1)) {
      LOG(INFO) << "Skipping line " << tdata.name << " because " << tdata.sname1 << " failed to converge before\n";
      continue;
    }
    if (vectorContains(failconv_states, tdata.sname2)) {
      LOG(INFO) << "Skipping line " << tdata.name << " because " << tdata.sname2 << " failed to converge before\n";
      continue;
    }

    LOG(INFO) << "Computing transition " << tdata.name << "\n";

    try {
      LOG(INFO) << "Computing state " << tdata.sname1 << "\n";

      // Here we actually solve for the first state
      tdata.ds1 = da.getState(n1, l1, s1);
      LOG(INFO) << "Computing state " << tdata.sname2 << "\n";

      // Here we actually solve for the final state
      tdata.ds2 = da.getState(n2, l2, s2);
    } catch (AtomErrorCode aerr) {
      LOG(ERROR) << SPECIAL << "Transition energy calculation for line " << tdata.name << " failed with AtomErrorCode " << aerr << "\n";
      success = false;
    } catch (const exception &e) {
      LOG(ERROR) << SPECIAL << "Unknown error: " << e.what() << "\n";
      success = false;
    }
    if (!success) {
      LOG(INFO) << "Convergence of one state failed for line " << tdata.name << ", skipping\n";
      if (!tdata.ds1.converged) {
        failconv_states.push_back(tdata.sname1);
      } else {
        failconv_states.push_back(tdata.sname2);
      }
      continue;
    }

    // Compute transition probability
    tdata.tmat = da.getTransitionProbabilities(n2, l2, s2, n1, l1, s1);

    LOG(INFO) << "Transition energy = " << (tdata.ds2.E - tdata.ds1.E) / (Physical::eV * 1000) << " kEv\n";

    transitions.push_back(tdata);
  }

  return transitions;


}
/** This takes in the input config and returns a vector of transitiondata
 */
vector<TransLineSpec> parseXRLines(MuDiracInputFile config) {

  // First we unravel the user specified string
  vector<string> xr_lines = config.getStringValues("xr_lines");

  // Convert the user specified strings into quantum numbers for the start
  // and end states in the transition
  vector<TransLineSpec> transqnums;

  for (int i = 0; i < xr_lines.size(); ++i) {
    vector<string> ranges = splitString(xr_lines[i], "-");
    vector<int> n1range, n2range, l1range, l2range;
    vector<bool> s1range, s2range;

    LOG(TRACE) << "Parsing XR line specification " << xr_lines[i] << "\n";

    if (ranges.size() != 2) {
      LOG(ERROR) << SPECIAL << "Line " << xr_lines[i] << " can not be interpreted properly\n";
      throw invalid_argument("Invalid spectral line in input file");
    }

    vector<int> nr, lr;
    vector<bool> sr;

    parseIupacRange(ranges[0], nr, lr, sr);
    n1range.insert(n1range.end(), nr.begin(), nr.end());
    l1range.insert(l1range.end(), lr.begin(), lr.end());
    s1range.insert(s1range.end(), sr.begin(), sr.end());

    parseIupacRange(ranges[1], nr, lr, sr);
    n2range.insert(n2range.end(), nr.begin(), nr.end());
    l2range.insert(l2range.end(), lr.begin(), lr.end());
    s2range.insert(s2range.end(), sr.begin(), sr.end());

    for (int j = 0; j < n1range.size(); ++j) {
      for (int k = 0; k < n2range.size(); ++k) {
        TransLineSpec tnums;
        tnums.n1 = n1range[j];
        tnums.l1 = l1range[j];
        tnums.s1 = s1range[j];

        tnums.n2 = n2range[k];
        tnums.l2 = l2range[k];
        tnums.s2 = s2range[k];

        if (tnums.n2 < tnums.n1 || abs(tnums.l2 - tnums.l1) != 1) {
          continue;
        }

        transqnums.push_back(tnums);

        LOG(TRACE) << "Identified transition: " << tnums.n1 << ", " << tnums.l1 << ", " << tnums.s1 << "\t";
        LOG(TRACE) << tnums.n2 << ", " << tnums.l2 << ", " << tnums.s2 << "\n";
      }
    }
  }

  return transqnums;
}


void optimiseFermiParameters(const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, const vector<TransLineSpec> &transqnums,  vector<OptimisationData> &valid_fermi_parameters, OptimisationData &optimal_fermi_parameter,vector<TransitionData> &optimal_transitions, MuDiracInputFile &config ) {
  double opt_fermi_c, opt_fermi_t;
  //double uniform_radius = 1.25 * cbrt((double) config.getIntValue("isotope"));
  // define RMS_0
  // double rms_radius_0 = sqrt(3.0/5.0) * uniform_radius;

  // set variables for rms radius loop limits and increments
  double rms_radius_min = config.getDoubleValue("rms_radius_min");
  double rms_radius_max = config.getDoubleValue("rms_radius_max");
  double rms_iterations_factor = (double) pow(10, config.getIntValue("rms_radius_decimals"));
  double rms_radius_increment = 1.0/rms_iterations_factor;
  double rms_radius = rms_radius_min;
  double rms_radius_iterations = (rms_radius_max - rms_radius_min)*rms_iterations_factor;
  // set variables for theta loop iterations
  double theta;
  int theta_iterations = config.getIntValue("theta_iterations");
  double theta_min = M_PI * config.getDoubleValue("theta_min");
  double theta_max =  M_PI * config.getDoubleValue("theta_max");
  double theta_increment = (theta_max-theta_min)/((double)(theta_iterations));

  // estimate iterations required for the logs
  int total_2pF_iterations = theta_iterations * (1 + (int)((rms_radius_max - rms_radius)*rms_iterations_factor));

  // ouput optimisation domain and iterations to log
  LOG(INFO) << "Starting scan for optimal fermi parameters \n";
  LOG(INFO) << "search domain: rms_radius: [" <<rms_radius << ", "<< rms_radius_max << "]\n";
  LOG(INFO) << "search domain: theta: [" <<theta_min << ", "<< theta_max << "]\n";
  LOG(INFO) << "total iterations: "<< total_2pF_iterations <<"\n";

  // optimisation loops

  for (int i=0; i < theta_iterations; ++i) {
    theta = theta_min + (i * theta_increment);

    //while (rms_radius < rms_radius_max){
    for (int j=0; j < rms_radius_iterations; ++j) {
      rms_radius = rms_radius_min + ((double) j )* rms_radius_increment;
      // get fermi parameters from rms_radius, theta
      tie(opt_fermi_c, opt_fermi_t) = fermiParameters(rms_radius, theta);

      // set new iteration of fermi parameters in config and get transitions
      config.defineDoubleNode("fermi_t", InputNode<double>(opt_fermi_t));
      config.defineDoubleNode("fermi_c", InputNode<double>(opt_fermi_c));
      LOG(DEBUG) << "creating atom with fermi parameters: " << opt_fermi_c << ", " << opt_fermi_t;
      LOG(DEBUG) << " RMS radius: " << rms_radius << " theta: "<< theta << "\n";
      DiracAtom opt_da = config.makeAtom();
      vector<TransitionData> transitions_iteration = getAllTransitions(transqnums, opt_da);

      // declare validity variables for each iteration of fermi parameters
      double mean_square_error = 0;
      double total_square_error = 0;
      bool fermi_parameters_are_valid = true;

      // store MSE for each transition in a vector, RMS, theta, MSE
      OptimisationData fermi_parameter_iteration;

      // loop calculating MSE for each iteration of fermi parameters
      LOG(DEBUG) << "MSE loop \n";
      for (int k = 0; k < transitions_iteration.size(); ++k) {

        // calculate transition energy and rate
        double dE = (transitions_iteration[k].ds2.E - transitions_iteration[k].ds1.E);
        double tRate = transitions_iteration[k].tmat.totalRate();

        // square error for each transitions calculated
        double square_error = 0;

        if (dE <= 0 || tRate <= 0)
          continue; // Transition is invisible

        // check transition allign with experimental transitions
        if (transitions_iteration[k].name == xr_lines_measured[k]) {
          // convert to eV
          double transition_energy = dE / Physical::eV;

          // calculate the square error of each transition
          double square_deviation = (transition_energy-xr_energies[k])*(transition_energy-xr_energies[k]);
          double valid_uncertainty = (xr_errors[k])*(xr_errors[k]);
          square_error = square_deviation/valid_uncertainty;

          // output square error to LOG
          LOG(DEBUG) << transitions_iteration[k].name << " SE: "<< square_error << "\n";
        }

        // break MSE loop if c, t are invalid for any of the transitions_iteration
        if (square_error >= 1) {
          LOG(DEBUG) << "c,t doesnt fit\n";
          fermi_parameters_are_valid = false;
          break;
        }
        // sum the square error for all transitions
        total_square_error += square_error;
      }

      // store the fermi parameters if valid (square_error < 1 for all transitions)
      if (fermi_parameters_are_valid) {
        LOG(DEBUG) << "fermi parameters valid\n";

        //store the parameters in the datastructure
        fermi_parameter_iteration.rms_radius = rms_radius;
        fermi_parameter_iteration.theta = theta;
        mean_square_error = total_square_error / transitions_iteration.size();
        fermi_parameter_iteration.mse = mean_square_error;
        tie(fermi_parameter_iteration.fermi_c, fermi_parameter_iteration.fermi_t) = fermiParameters(rms_radius, theta);

        // add the valid parameter datastructure to a list containing all valid parameters
        valid_fermi_parameters.push_back(fermi_parameter_iteration);

        // check if current iteration is the best iteration so far
        if (mean_square_error < optimal_fermi_parameter.mse) {

          // store parameters and transitions of the current best iteration
          optimal_fermi_parameter = fermi_parameter_iteration;
          optimal_transitions = transitions_iteration;
        }
      }
    }
  }

  // if there are no valid fermi parameters
  if (optimal_fermi_parameter.mse >= 1.0) {
    LOG(DEBUG) << "no valid fermi parameters found\n";
    exit(-1);
  } else {
    // end of loop want optimum c,t with MSE < 1
    // rms radius scanned from lowest to highest so front/back of valid parameters vector contain min/max valid rms radius values
    double max_valid_rms_radius = valid_fermi_parameters.back().rms_radius;
    double min_valid_rms_radius = valid_fermi_parameters.front().rms_radius;
    double rms_radius_uncertainty = (max_valid_rms_radius - min_valid_rms_radius)/2;

    // output valid rms radius range, and optimal fermi parameters to LOG
    LOG(INFO) << "valid rms radius range [" << min_valid_rms_radius << ", " << max_valid_rms_radius << "]\n";
    LOG(INFO) << "found optimal fermi polar parameters " << optimal_fermi_parameter.rms_radius;
    LOG(INFO) << ", "<< optimal_fermi_parameter.theta << "\n";
    LOG(INFO) << "found optimal fermi parameters: " << optimal_fermi_parameter.fermi_c;
    LOG(INFO) << ", " << optimal_fermi_parameter.fermi_t << "\n";
  }
}