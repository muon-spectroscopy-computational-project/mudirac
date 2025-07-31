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
  setupLogging(config.getIntValue("verbosity"), seed);
  debugTasks(config);   // runs debug tasks if configured
  DiracAtom da = config.makeAtom();   // construct the atom

  // Print out electronic potential if high levels of verbosity are chosen 
  writeEConfV(da, seed, config.getIntValue("output"));

  // contains both state names, actual dirac states, transition matrix
  vector<TransitionData> transitions;     

  // performs least squares optimisation of 2pF parameters and find rms radius
  if (config.getBoolValue("optimise_fermi_parameters")) {
    runFermiModelOptimisation(config, argc, argv, da, seed);
  }

  // Default muDirac behaviour:
  da.transqnums = config.parseXRLines();    // set qnums for config transitions
  transitions = da.getAllTransitions();     // calculate transition energies
  writeOutputFiles(seed, config, da, transitions);    // output results
  t1 = chrono::high_resolution_clock::now();
  double t_tot = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();
  LOG(INFO) << "Calculation completed in " << t_tot / 1.0e3 << " seconds\n";
}

void setupLogging(const int verbosity, const string & seed){
  AixLog::Severity log_verbosity;
  switch (verbosity) {
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

