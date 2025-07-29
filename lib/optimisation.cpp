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

void globalOptimizeFermiParameters(DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time) {
  // initialise starting parameters for optimisation based on the coordinate system
  const string coord_sys = da.coord_system;
  vector<double> init_params = da.getFermi2(coord_sys);


  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

  column_vector bound_lower, bound_upper;
  if (coord_sys =="polar"){
    // set bounds for rms radius and theta close to Marinova table values
    double rms_lower = 0.95*init_params.at(0);
    double rms_upper = 1.05*init_params.at(0);
    double theta_lower = 0;
    double theta_upper = M_PI*0.5;
    bound_lower = {rms_lower, theta_lower};
    bound_upper = {rms_upper, theta_upper};
  }
  else if (coord_sys =="ct"){
    double c_lower = 0.5*init_params.at(0);
    double c_upper = 1.5*init_params.at(0);
    double t_lower = 0.01;
    double t_upper = 4;
    bound_lower = {c_lower, t_lower};
    bound_upper = {c_upper, t_upper};
  }

  // lambda function for MSE required for dlib find min global
  auto MSE_lambda = [&](double rms_radius, double theta){
    return da.calculateMSE(rms_radius, theta);
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

  if (coord_sys =="polar"){
    fermi_parameters.rms_radius = final_params(0);
    fermi_parameters.theta = final_params(1);
    tie(fermi_parameters.fermi_c, fermi_parameters.fermi_t) = fermiParameters(fermi_parameters.rms_radius, fermi_parameters.theta);
  }
  else if (coord_sys=="ct"){
    fermi_parameters.fermi_c =  final_params(0); 
    fermi_parameters.fermi_t = final_params(1);
    fermi_parameters.rms_radius = rmsRadius(final_params(0), final_params(1));
  }

  LOG(INFO) << "minimised with MSE: "<< MSE << " and fermi " << coord_sys << " parameters: "<< final_params(0)<<", " << final_params(1) <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << da.iteration_counter_2pF <<" iterations from MuDirac objective function \n";

  // repeat the final configuration of the nuclear model
  //configureNuclearModel(final_params, coord_sys, da);
  da.setFermi2(final_params(0), final_params(1), coord_sys);
}


void optimizeFermiParameters(DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time) {
  LOG(INFO) << "Starting minimisation for fermi model \n";

  // initialise starting parameters for optimisation based on the coordinate system
  const string coord_sys = da.coord_system;
  vector<double> fermi_coords = da.getFermi2(coord_sys);
  // change to dlib column vector type
  column_vector init_params = {fermi_coords.at(0), fermi_coords.at(1)};

  // dlib functions for minimisation only finds minimum, no bayesian uncertainty  analysis.
  double MSE;
  LOG(INFO) << "Minimising the MSE over the fermi parameters using the bgfs search strategy \n";
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();

   auto MSE_lambda = [&](column_vector x){
    return da.calculateMSE(x(0), x(1));
  };

  auto der_lambda = [&](column_vector x){
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
  LOG(INFO) << "minimised with MSE: "<< MSE << " and fermi " << coord_sys << " parameters: "<< init_params <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << da.iteration_counter_2pF <<" iterations from MuDirac objective function \n";

  // repeat the final configuration of the nuclear model
  // init_params have been updated by the minimisation to the optimising params
  da.setFermi2(init_params(0), init_params(1), coord_sys);

  vector<double> ct_coords = da.getFermi2("ct");
  vector<double> polar_coords = da.getFermi2("polar");
  fermi_parameters.fermi_c = ct_coords.at(0);
  fermi_parameters.fermi_t = ct_coords.at(1);
  fermi_parameters.rms_radius = polar_coords.at(0);
  fermi_parameters.theta = polar_coords.at(1);
  fermi_parameters.mse = MSE;
}


void optimizeFermiParameters(const opt_2pF_model &opt_obj, DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time){
  LOG(INFO) << "Starting minimisation for fermi model using trust region method\n";

  // initialise starting parameters for optimisation based on the coordinate system
  const string coord_sys = da.coord_system;
  vector<double> fermi_coords = da.getFermi2(coord_sys);
  column_vector init_params = {fermi_coords.at(0), fermi_coords.at(1)};

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
  
  LOG(INFO) << "minimised with MSE: "<< MSE << " and "<< coord_sys << "fermi parameters: "<< init_params <<" \n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";
  LOG(INFO) << "minimised using " << da.iteration_counter_2pF <<" iterations from MuDirac objective function \n";
  // repeat the final configuration of the nuclear model
  // init_params have been updated by the minimisation to the optimising params
  da.setFermi2(init_params(0), init_params(1), coord_sys);
  fermi_parameters.mse = MSE;
  vector<double> ct_coords = da.getFermi2("ct");
  vector<double> polar_coords = da.getFermi2("polar");
  fermi_parameters.fermi_c = ct_coords.at(0);
  fermi_parameters.fermi_t = ct_coords.at(1);
  fermi_parameters.rms_radius = polar_coords.at(0);
  fermi_parameters.theta = polar_coords.at(1);
}


column_vector MSE_2pF_derivative(const column_vector &m, DiracAtom & da) {

  // compute gradient
  // bind calculate mse as function of just polar fermi parameters column vector
  auto MSE_lambda = [&](column_vector x){
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



dlib::matrix<double> MSE_2pF_hessian(const column_vector & m, DiracAtom & da){
  dlib::matrix<double> res(2,2);
  
  // choose derivative step size
  double d_r = 1e-7;
  double d_theta = 1e-7;
  column_vector delta_r = {d_r, 0};
  column_vector delta_theta = {0, d_theta};
  LOG(DEBUG) << "computing hessian at " << da.coord_system << " fermi parameters (" << m(0) << ", " << m(1) <<") \n";
  //LOG(DEBUG) << "computing hessian at config fermi parameters(" << config.getDoubleValue["fermi_c"] << ", " << config.getDoubleValue["fermi_t"] <<") \n";
  auto der_lambda = [&](column_vector x){
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

void optFermi2(DiracAtom & da, const string algo, OptimisationData & best_fermi_parameters, double & opt_time){
  if (algo == "bfgs"){
      optimizeFermiParameters(da, best_fermi_parameters, opt_time);
      }
  else if (algo == "trust"){
      opt_2pF_model opt_obj(da);
      optimizeFermiParameters(opt_obj, da, best_fermi_parameters, opt_time);
  }
  else if (algo=="global"){
      globalOptimizeFermiParameters(da, best_fermi_parameters, opt_time);
  }
  else {
      cout << "Invalid 2pF optimisation algorithm choice for minimsation\n";
      cout << "please use \"bfgs\" or \"trust\" (default is \"bfgs\") \n";
      cout << "You used: \""<<algo<<"\" \n";
      cout << "Quitting...\n";
      LOG(ERROR) << "Invalid 2pF optimisation algorithm choice for minimsation: \""<<algo<<"\"\n";
  }
}