/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * optimisation.cpp
 *
 * Functions for optimising the dirac atom nuclear model parameters.
 *
 * @author Milan Kumar
 * @version 1.0 30/06/2025
 */


#include "optimisation.hpp"

void globalOptimizeFermiParameters(DiracAtom & da, double & opt_time) {
  // initialise starting parameters for optimisation based on the coordinate system
  array<double,2> init_params = da.getFermi2(da.coord_system);

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  double c1_lower_coef, c1_upper_coef, c2_lower_coef, c2_upper_coef;
  if (da.coord_system =="polar") {
    // set bounds for rms radius and theta close to Marinova table values
    c1_lower_coef = 0.95;
    c1_upper_coef = 1.05;
    c2_lower_coef = 0;
    c2_upper_coef = M_PI*0.5;
  } else if (da.coord_system =="ct") {
    // set bounds for c and t within reasonable range of original values
    c1_lower_coef = 0.5;
    c1_upper_coef = 1.5;
    c2_lower_coef = 0.01;
    c2_upper_coef = 4;
  }

  // lambda function for MSE required for dlib find min global
  auto MSE_lambda = [&](double rms_radius, double theta) {
    return da.calculateMSE(rms_radius, theta);
  };

  // global optimisation where result.x is opt rms radius and coords, result.y is MSE
  // max run time 10 minutes can be changed
  auto result = dlib::find_min_global(
                  MSE_lambda,
  {c1_lower_coef*init_params[0], c2_lower_coef*init_params[1]}, //lower bounds
  {c1_upper_coef*init_params[0], c2_upper_coef*init_params[1]}, //upper bounds
  std::chrono::milliseconds(60*1000*10));

  opt_t1 = chrono::high_resolution_clock::now();
  auto final_params = result.x;
  double MSE = result.y;

  // update fermi_paramters structure
  da.fermi2.mse = MSE;
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  finaliseFermi2(da, da.coord_system, final_params, opt_time, MSE);
}


void bfgsOptimizeFermiParameters(DiracAtom & da, double & opt_time) {
  LOG(INFO) << "Starting minimisation for fermi model \n";

  // initialise starting parameters for optimisation based on the coordinate system
  array<double, 2> fermi_coords = da.getFermi2(da.coord_system);
  // change to dlib column vector type
  column_vector init_params = {fermi_coords.at(0), fermi_coords.at(1)};

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  double MSE;
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  auto MSE_lambda = [&](column_vector x) {
    return da.calculateMSE(x(0), x(1));
  };

  auto der_lambda = [&](column_vector x) {
    return MSE_2pF_derivative(x, da);
  };

  MSE = dlib::find_min(
          dlib::bfgs_search_strategy(),
          dlib::objective_delta_stop_strategy(1e-2).be_verbose(),  // gradient change < 0.01
          MSE_lambda,
          der_lambda,
          init_params,
          -1);

  opt_t1 = chrono::high_resolution_clock::now();
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  finaliseFermi2(da, da.coord_system, init_params, opt_time, MSE);
}


void trustOptimizeFermiParameters(const opt_2pF_model &opt_obj, DiracAtom & da, double & opt_time) {
  LOG(INFO) << "Starting minimisation for fermi model using trust region method\n";

  // initialise starting parameters for optimisation based on the coordinate system
  array<double, 2> fermi_coords = da.getFermi2(da.coord_system);
  column_vector init_params = {fermi_coords.at(0), fermi_coords.at(1)};

  // start time of minimisation
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  double MSE = dlib::find_min_trust_region(dlib::objective_delta_stop_strategy(1e-2),
                                    opt_obj,
                                    init_params,
                                    0.1);
  da.iteration_counter_2pF = opt_obj.opt_iterations;
  opt_t1 = chrono::high_resolution_clock::now();
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  finaliseFermi2(da, da.coord_system, init_params, opt_time, MSE);
}

void finaliseFermi2(DiracAtom & da, const string coord_sys, column_vector final_fermi_params, double opt_time, double MSE) {

  // output final optimisation values to LOG
  LOG(INFO) << "minimised with MSE: "<< MSE << " and "<< coord_sys << "fermi parameters: "<< final_fermi_params <<" \n";
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << da.iteration_counter_2pF <<" iterations from MuDirac objective function \n";

  // ensure final fermi parameters are set
  da.setFermi2(final_fermi_params(0), final_fermi_params(1), coord_sys);

  //set all the optimisation parameter values
  array<double, 2> polar_coords = da.getFermi2("polar");
  da.fermi2.mse = MSE;
}


column_vector MSE_2pF_derivative(const column_vector &m, DiracAtom & da) {

  // compute gradient
  // bind calculate mse as function of just polar fermi parameters column vector
  auto MSE_lambda = [&](column_vector x) {
    return da.calculateMSE(x(0), x(1));
  };
  LOG(DEBUG) << "computing derivative at " << da.coord_system <<" fermi parameters (" << m(0) << ", " << m(1) <<") \n";

  // get the derivative by central differences
  auto res_func = dlib::derivative(MSE_lambda, 1e-7);
  //LOG(DEBUG) << " computing derivative, config fermi parameters: (" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  // get the values of the derivative at m as a column vector?
  auto res = res_func(m);
  //LOG(DEBUG) << " derivative computed, config fermi parameters: (" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  return res;
}


dlib::matrix<double> MSE_2pF_hessian(const column_vector & m, DiracAtom & da) {
  dlib::matrix<double> res(2,2);

  // choose derivative step size
  double d_r = 1e-7;
  double d_theta = 1e-7;
  column_vector delta_r = {d_r, 0};
  column_vector delta_theta = {0, d_theta};
  LOG(DEBUG) << "computing hessian at " << da.coord_system << " fermi parameters (" << m(0) << ", " << m(1) <<") \n";
  //LOG(DEBUG) << "computing hessian at config fermi parameters(" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  auto der_lambda = [&](column_vector x) {
    return MSE_2pF_derivative(x, da);
  };

  auto hess_r_component = (der_lambda(m + delta_r) - der_lambda(m - delta_r));
  auto hess_theta_component = (der_lambda(m + delta_theta) - der_lambda(m - delta_theta));
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

void optFermi2(DiracAtom & da, const string algo, double & opt_time) {
  if (algo == "bfgs") {
    bfgsOptimizeFermiParameters(da, opt_time);
  } else if (algo == "trust") {
    opt_2pF_model opt_obj(da);
    trustOptimizeFermiParameters(opt_obj, da, opt_time);
  } else if (algo=="global") {
    globalOptimizeFermiParameters(da, opt_time);
  } else {
    cout << "Invalid 2pF optimisation algorithm choice for minimsation\n";
    cout << "please use \"bfgs\" or \"trust\" (default is \"bfgs\") \n";
    cout << "You used: \""<<algo<<"\" \n";
    cout << "Quitting...\n";
    LOG(ERROR) << "Invalid 2pF optimisation algorithm choice for minimsation: \""<<algo<<"\"\n";
  }
}


void runFermiModelOptimisation(MuDiracInputFile & config, const int & argc, char * argv[], DiracAtom &da, const string & seed) {
  string coord_system_2pF;
  string min_2pF_algo;
  config.validateOptimisation(argc, coord_system_2pF, min_2pF_algo);

  // initialise experimental results file object, read xr input file and validate
  ExperimentalResultFile measurements;
  measurements.validate(argv[2]);
  LOG(INFO) << "Successfully read xray measurements input file \n";

  // get transition quantum numbers for measured transitions
  da.transqnums = measurements.parseXRLines();

  // set all optimization values in Dirac Atom
  da.setExpOptData(coord_system_2pF, measurements.getStringValues("xr_lines"), measurements.getDoubleValues("xr_energy"), measurements.getDoubleValues("xr_error"));

  double opt_time = 0;

  // perform the 2pF nuclear charge model optimisation routine
  optFermi2(da, min_2pF_algo, opt_time);

  // output file containing best fermi parameters and the associated MSE
  writeFermiParameters(da, opt_time,  seed + "fermi_parameters.out", config.getIntValue("rms_radius_decimals"));
}

