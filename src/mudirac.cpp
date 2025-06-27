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
    bool xr_measurement_read_success = false;
    try {
      measurements.parseFile(argv[2]);

    } catch (runtime_error e) {
      cout << "Invalid experimental measurements file:\n";
      cout << e.what() << "\n";
      return -1;
    }

    // read the measured transition lines
    vector<string> xr_lines_measured = measurements.getStringValues("xr_lines");
    LOG(DEBUG) << "Reading experimental Xray measurments for transitions: ";
    for (auto transition: xr_lines_measured) {
      LOG(DEBUG) << transition << ", ";
    }
    LOG(DEBUG) << "\n";

    // read the measured transition energies
    vector<double> xr_energies = measurements.getDoubleValues("xr_energy");
    LOG(DEBUG) << "Reading experimental Xray energies: ";
    for (auto transition_energy: xr_energies) {
      LOG(DEBUG) << transition_energy << ", ";
    }
    LOG(DEBUG) << "\n";

    // read the measured transition errors
    vector<double> xr_errors = measurements.getDoubleValues("xr_error");
    LOG(DEBUG) << "Reading experimental Xray energy errors: ";
    for (auto transition_energy_error: xr_errors) {
      LOG(DEBUG) << transition_energy_error << ", ";
    }
    LOG(DEBUG) << "\n";

    // checking that the file has contents and not the default values
    LOG(DEBUG) << "Validating experimental results input \n";
    if (xr_lines_measured[0] == "") {
      cout << "Experimental results input file is empty\n";
      cout << "Please check the filename of the experimental results input file \n";
      cout << "Quitting...\n";
      return -1;
    }

    // check that the data provided is complete: all transitions measured have energies and errors
    if (xr_lines_measured.size() == xr_energies.size() && xr_energies.size() == xr_errors.size()) {
      xr_measurement_read_success = true;
    } else {
      cout << "Invalid experimental measurements file: Missing input values\n";
      cout << "please check energies and errors are listed for each xray transition line \n";
      cout << "Quitting...\n";
      return -1;

    }

    if (xr_measurement_read_success) {
      LOG(INFO) << "Successfully read xray measurements input file \n";
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

      if (min_2pF_algo == "bfgs"){
        optimizeFermiParameters(config, coord_system_2pF, da, transqnums, xr_lines_measured, xr_energies, xr_errors, best_fermi_parameters, opt_time);
      }
      else if (min_2pF_algo =="trust"){
        opt_2pF_model opt_obj(config, coord_system_2pF, transqnums, xr_lines_measured, xr_energies, xr_errors);
        optimizeFermiParameters(opt_obj, coord_system_2pF, config, da, best_fermi_parameters, opt_time);
      }
      else if (min_2pF_algo=="global"){
        globalOptimizeFermiParameters(config, coord_system_2pF, da, transqnums, xr_lines_measured, xr_energies, xr_errors, best_fermi_parameters, opt_time);
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
      writeFermiParameters(da, best_fermi_parameters,iteration_counter_2pF, opt_time,  seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
    }

  }

  // Default mudirac behaviour
  // Wrapped the calculation of the states, their energies and the transition probabilities into here,
  // so that we can easily loop over it for least squares optimisation
  transitions = da.getAllTransitions(transqnums);

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


void init2pFModelParams(DiracAtom & da, const string coord_system, column_vector & init_params){

  double fermi_c_init, fermi_t_init;
  // get c and t after it has either been defined by the user or default values by the marinova tables
  LOG(DEBUG) << " getting initial fermi parameters \n";
  tie(fermi_c_init, fermi_t_init) =  da.getFermi2();

  // set starting optimisation paramters based on the coordinate system we optimize over.
  if (coord_system == "ct"){
    LOG(DEBUG) << " optimising fermi paramaters using c, t coordinate system \n";
    LOG(DEBUG) << " initial fermi parameters in fm: ("<< fermi_c_init <<", " << fermi_t_init <<  ") \n";
    init_params = {fermi_c_init, fermi_t_init};
  }
  else if (coord_system == "polar"){
    // get rms radius estimate from marinova tables
    LOG(DEBUG) << " optimising fermi paramaters using polar coordinate system \n";
    LOG(DEBUG) << " getting initial polar fermi parameters \n";
    double rms_radius_estimate = rmsRadius(fermi_c_init, fermi_t_init);
    double theta_estimate = atan(fermi_t_init/fermi_c_init);
    LOG(DEBUG) << " initial polar fermi parameters: ("<< rms_radius_estimate <<" fm, " << theta_estimate << ") \n";

    init_params = {rms_radius_estimate, theta_estimate};
  }
}


void configureNuclearModel(const column_vector& m, const string coord_system, MuDiracInputFile &config, DiracAtom & da, OptimisationData &fermi_parameters) {
  double fermi_c, fermi_t, rms_radius, theta;
  if (coord_system == "polar"){
    LOG(DEBUG) << " configuring dirac atom using polar coordinate system \n";
    rms_radius = m(0);
    theta = m(1);
    // calculate new c and t
    tie(fermi_c, fermi_t) = fermiParameters(rms_radius, theta);
  }
  else if(coord_system=="ct"){
    LOG(DEBUG) << " configuring dirac atom using c, t coordinate system \n";
    fermi_c = m(0);
    fermi_t = m(1);
    rms_radius=rmsRadius(fermi_c, fermi_t);
    theta=0;
  }

  // populate the fermi parameters structure
  fermi_parameters.rms_radius = rms_radius;
  fermi_parameters.theta = theta;
  fermi_parameters.fermi_c = fermi_c;
  fermi_parameters.fermi_t = fermi_t;

  // set new iteration of fermi parameters in config and get transitions
  config.defineDoubleNode("fermi_t", InputNode<double>(fermi_t));
  config.defineDoubleNode("fermi_c", InputNode<double>(fermi_c));
  LOG(DEBUG) << "creating atom with fermi parameters: " << fermi_c << " fm, " << fermi_t << " fm\n";
  LOG(DEBUG) << " RMS radius: " << rms_radius << " fm, theta: "<< theta << "\n";

  // make the new Dirac atom with the new configuration
  da = config.makeAtom();
}


double calculateMSE(const column_vector& m, const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors) {
  
  // increment iteration counter
  iteration_counter_2pF++;
  // create a new dirac atom
  DiracAtom da;
  OptimisationData iteration_parameters;

  // configure the new dirac atom with the new fermi parameters in the given coordinate system
  configureNuclearModel(m, coord_system, config, da, iteration_parameters);

  // calculate the energies of the transitions with the new fermi parameters.
  vector<TransitionData> transitions_iteration = da.getAllTransitions(transqnums);

  LOG(DEBUG) << "MSE loop \n";
  double MSE = 0;
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
      MSE += square_error;
    }

  }
  MSE = MSE/transitions_iteration.size();
  // output MSE to LOG
  LOG(DEBUG) << "MSE: "<< MSE << "\n";
  return MSE;
}

void globalOptimizeFermiParameters(MuDiracInputFile &config, const string coord_system, DiracAtom & da, const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters, double & opt_time) {
  // initialise starting parameters for optimisation based on the coordinate system
  column_vector init_params;
  init2pFModelParams(da, coord_system, init_params);

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  column_vector bound_lower, bound_upper;
  if (coord_system =="polar"){
    // set bounds for rms radius and theta close to Marinova table values
    double rms_lower = 0.95*init_params(0);
    double rms_upper = 1.05*init_params(0);
    double theta_lower = 0;
    double theta_upper = M_PI*0.5;
    bound_lower = {rms_lower, theta_lower};
    bound_upper = {rms_upper, theta_upper};
  }
  else if (coord_system =="ct"){
    double c_lower = 0.5*init_params(0);
    double c_upper = 1.5*init_params(0);
    double t_lower = 0.01;
    double t_upper = 4;
    bound_lower = {c_lower, t_lower};
    bound_upper = {c_upper, t_upper};
  }

  // lambda function for MSE required for dlib find min global
  auto MSE_lambda = [&](double rms_radius, double theta){
    column_vector x = {rms_radius, theta};

    return calculateMSE(x,coord_system, config, transqnums, xr_lines_measured, xr_energies, xr_errors);
  };
  
  // global optimisation where result.x is opt rms radius and coords, result.y is MSE
  // max run time 10 minutes can be changed
  auto result = dlib::find_min_global(
    MSE_lambda,
    bound_lower,
    bound_upper,
    std::chrono::milliseconds(60*1000*10));

  opt_t1 = chrono::high_resolution_clock::now();
  auto final_params = result.x;
  double MSE = result.y;

  // update fermi_paramters structure
  fermi_parameters.mse = MSE;

  if (coord_system =="polar"){
    fermi_parameters.rms_radius = final_params(0);
    fermi_parameters.theta = final_params(1);
    tie(fermi_parameters.fermi_c, fermi_parameters.fermi_t) = fermiParameters(fermi_parameters.rms_radius, fermi_parameters.theta);
  }
  else if (coord_system=="ct"){
    fermi_parameters.fermi_c =  final_params(0); 
    fermi_parameters.fermi_t = final_params(1);
    fermi_parameters.rms_radius = rmsRadius(final_params(0), final_params(1));
  }

  LOG(INFO) << "minimised with MSE: "<< MSE << " and fermi " << coord_system << " parameters: "<< final_params <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << iteration_counter_2pF <<" iterations from MuDirac objective function \n";

  // repeat the final configuration of the nuclear model
  configureNuclearModel(final_params, coord_system, config, da, fermi_parameters);
}


void optimizeFermiParameters(MuDiracInputFile &config, const string coord_system, DiracAtom & da, const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters, double & opt_time) {
  LOG(INFO) << "Starting minimisation for fermi model \n";

  // initialise starting parameters for optimisation based on the coordinate system
  column_vector init_params;
  init2pFModelParams(da, coord_system, init_params);

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  double MSE;
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  MSE = dlib::find_min(
          dlib::bfgs_search_strategy(),
          dlib::objective_delta_stop_strategy(1e-2).be_verbose(),  // gradient change < 0.01
          std::bind(&calculateMSE, std::placeholders::_1, coord_system, config, transqnums, xr_lines_measured, xr_energies, xr_errors),
          std::bind(&MSE_2pF_derivative, std::placeholders::_1, coord_system, config, transqnums, xr_lines_measured, xr_energies, xr_errors),
          init_params,
          -1);
  
  opt_t1 = chrono::high_resolution_clock::now();
  LOG(INFO) << "minimised with MSE: "<< MSE << " and fermi " << coord_system << " parameters: "<< init_params <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << iteration_counter_2pF <<" iterations from MuDirac objective function \n";

  // repeat the final configuration of the nuclear model
  configureNuclearModel(init_params, coord_system, config, da, fermi_parameters);
  fermi_parameters.mse = MSE;
}


void optimizeFermiParameters(opt_2pF_model &opt_obj, const string coord_system, MuDiracInputFile & config, DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time){
  LOG(INFO) << "Starting minimisation for fermi model using trust region method\n";

  // initialise starting parameters for optimisation based on the coordinate system
  column_vector init_params;
  init2pFModelParams(da, coord_system, init_params);

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  double MSE;

  // start time of minimisation
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  MSE = dlib::find_min_trust_region(dlib::objective_delta_stop_strategy(1e-2),
                                    opt_obj,
                                    init_params,
                                    0.1);

  opt_t1 = chrono::high_resolution_clock::now();
  
  LOG(INFO) << "minimised with MSE: "<< MSE << " and "<< coord_system << "fermi parameters: "<< init_params <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << iteration_counter_2pF <<" iterations from MuDirac objective function \n";
  // repeat the final configuration of the nuclear model
  configureNuclearModel(init_params, coord_system, config, da, fermi_parameters);
  fermi_parameters.mse = MSE;
}


const column_vector MSE_2pF_derivative( const column_vector &m, const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors ) {

  // compute gradient
  // bind calculate mse as function of just polar fermi parameters column vector
  auto MSE_2pF = std::bind(&calculateMSE, std::placeholders::_1, coord_system, config, transqnums, xr_lines_measured, xr_energies, xr_errors);
  LOG(DEBUG) << "computing derivative at " << coord_system <<" fermi parameters (" << m(0) << ", " << m(1) <<") \n";

  // get the derivative by central differences
  auto res_func = dlib::derivative(MSE_2pF, 1e-7);
  //LOG(DEBUG) << " computing derivative, config fermi parameters: (" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  // get the values of the derivative at m as a column vector?
  auto res = res_func(m);
  //LOG(DEBUG) << " derivative computed, config fermi parameters: (" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  return res;
}



dlib::matrix<double> MSE_2pF_hessian(const column_vector &m, const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors){
  dlib::matrix<double> res(2,2);
  
  // choose derivative step size
  double d_r = 1e-7;
  double d_theta = 1e-7;
  column_vector delta_r = {d_r, 0};
  column_vector delta_theta = {0, d_theta};
  LOG(DEBUG) << "computing hessian at" << coord_system << "fermi parameters (" << m(0) << ", " << m(1) <<") \n";
  //LOG(DEBUG) << "computing hessian at config fermi parameters(" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  auto MSE_derivative = std::bind(&MSE_2pF_derivative, std::placeholders::_1, coord_system, config, transqnums, xr_lines_measured, xr_energies, xr_errors); 

  auto hess_r_component = (MSE_derivative(m + delta_r) - MSE_derivative(m - delta_r));
  auto hess_theta_component = (MSE_derivative(m + delta_theta) - MSE_derivative(m - delta_theta));
  //LOG(DEBUG) << "computed hessian, config fermi parameters(" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  // hessian df/dc^2 component
  res(0,0) = hess_r_component(0,0)/(2*d_r);

  // hessian df/dt^2 component
  res(1,1) = hess_theta_component(0,1)/(2*d_theta);

  // hessian df/dcdt components
  res(0,1) = hess_theta_component(0,0)/(2*d_theta);
  res(1,0) = hess_r_component(0,1)/(2*d_r);

  // return the finite differenced matrix
  return res;

}
