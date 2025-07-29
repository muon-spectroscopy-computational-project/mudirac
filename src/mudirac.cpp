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

  chrono::high_resolution_clock::time_point t0, t1;
  t0 = chrono::high_resolution_clock::now();
  string seed = "mudirac";
  MuDiracInputFile config;
  config.validate(argc, argv, seed);
  int output_verbosity = config.getIntValue("output");
  
  setupLogging(config, seed);

  // Are we running any debug tasks?
  debugTasks(config);

  // Here we construct the atom
  DiracAtom da = config.makeAtom();

  // Print out potential at high levels of verbosity
  writeEConfV(da, seed, output_verbosity);

  // data structure to contain transition quantum numbers
  vector<TransLineSpec> transqnums;

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

    // initialise experimental results file object
    ExperimentalResultFile measurements;

    // switch to continue with optimisation if read is successful
    bool xr_read = false;

    // read the measured transition lines and check they are valid
    measurements.validate(argv[2]); 
    vector<string> xr_lines_measured = measurements.getStringValues("xr_lines");
    vector<double> xr_energies = measurements.getDoubleValues("xr_energy");
    vector<double> xr_errors = measurements.getDoubleValues("xr_error");
    LOG(INFO) << "Successfully read xray measurements input file \n";

    // get transition quantum numbers for measured transitions
    transqnums = measurements.parseXRLines();
    // data structure for storing best parameters.
    OptimisationData best_fermi_parameters;
    double MSE =0;

    // set all optimization values in Dirac Atom
    da.setExpOptData(coord_system_2pF, transqnums, xr_lines_measured, xr_energies, xr_errors, best_fermi_parameters, opt_time);

    // perform the 2pF nuclear charge model optimisation routine
    optFermi2(da, min_2pF_algo, best_fermi_parameters, opt_time);

    // output file containing best fermi parameters and the associated MSE
    writeFermiParameters(da, best_fermi_parameters, opt_time,  seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
  
  }

  // Default mudirac behaviour

  // set the new transitions quantum numbers in the dirac atom
  da.transqnums = config.parseXRLines();

  // get the final transition data calculated by mudirac
  transitions = da.getAllTransitions();

  writeOutputFiles(seed, config, da, transitions);
  
  t1 = chrono::high_resolution_clock::now();

  LOG(INFO) << "Calculation completed in " << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1.0e3 << " seconds\n";
}

void setupLogging(MuDiracInputFile & config, const string & seed){
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

  printInitLogMessage();
}


void debugTasks(MuDiracInputFile & config){
  string debugtask = config.getStringValue("devel_debug_task");
  if (debugtask == "EdEscan") {
    LOG(INFO) << "Running debug task: E=>dE scan\n";
    runEdEscan(config);
    exit(0);
  }
}

void writeEConfV(DiracAtom & da, const string & seed, const int & out_verbosity){
  if (out_verbosity >= 2 && (da.getPotentialFlags() && da.HAS_ELECTRONIC)) {
    writeEConfPotential(da.getPotentialElectronic(), seed + ".epot.dat");
  }
}