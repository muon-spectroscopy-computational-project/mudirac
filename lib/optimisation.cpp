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

#define INFO CERES_MINILOG_INFO
#define WARNING CERES_MINILOG_WARNING
#define ERROR CERES_MINILOG_ERROR
#define FATAL CERES_MINILOG_FATAL
#include "ceres/ceres.h"
#include "ceres/covariance.h"
#undef INFO
#undef WARNING
#undef ERROR
#undef FATAL
#undef LOG

#include "optimisation.hpp"



void finaliseFermi2(DiracAtom & da, Fermi2CoordinateSystem coord_sys, double final_c1, double final_c2, double opt_time, double MSE, double sigma_c1, double sigma_c2, double cov_c1_c2) {

  // ensure final fermi parameters are set
  da.setFermi2Femto(final_c1, final_c2, coord_sys);
  //set all the optimisation parameter values
  da.fermi2.mse = MSE;
  if(coord_sys == CT){
    da.fermi2.sigma_c = sigma_c1;
    da.fermi2.sigma_t = sigma_c2;
    da.fermi2.cov_c_t = cov_c1_c2;

    double rms_radius = da.fermi2.rms_radius;
    if (da.fermi2.rms_radius > 0){
      double pifactor = M_PI / (4.0 * log(3.0));
      double drms_dc = (3.0/5.0) * da.fermi2.c /rms_radius;
      double drms_dt = (7.0/5.0) * pow(pifactor, 2) * da.fermi2.t / rms_radius;
      double var_rms_radius = pow(drms_dc, 2) * pow(sigma_c1, 2) + pow(drms_dt, 2) * pow(sigma_c2, 2)
       + 2 * drms_dc * drms_dt * cov_c1_c2;
      da.fermi2.sigma_rms_radius = sqrt(max(var_rms_radius, 0.0)); //max with 0 to avoid numerical issues with negative variance
    }

  } else if (coord_sys == POLAR) {
    da.fermi2.sigma_rms_radius = sigma_c1;
    
    double k = (7.0 * M_PI * M_PI )/ (48.0 * log(3.0) * log(3.0));
    double D = sqrt(pow(cos(da.fermi2.theta), 2) + k * pow(sin(da.fermi2.theta), 2));
    double sqrt_5_3 = sqrt(5.0/3.0);
    double sin_theta = sin(da.fermi2.theta);
    double cos_theta = cos(da.fermi2.theta);
    double R = sqrt_5_3 * da.fermi2.rms_radius;

    double dc_drms = sqrt_5_3 * cos_theta / D;
    double dt_drms = sqrt_5_3 * sin_theta / D;
    double dc_dtheta = -R * k * sin_theta / pow(D,3);
    double dt_dtheta = R * cos_theta / pow(D,3);

    double var_rms = pow(sigma_c1,2);
    double var_theta = pow(sigma_c2,2);
    double cov_rms_theta = cov_c1_c2;

    double var_c = pow(dc_drms, 2) * var_rms + pow(dc_dtheta, 2) * var_theta + 2 * dc_drms * dc_dtheta * cov_rms_theta;
    double var_t = pow(dt_drms, 2) * var_rms + pow(dt_dtheta, 2) * var_theta + 2 * dt_drms * dt_dtheta * cov_rms_theta;
    double cov_c_t = dc_drms * dt_drms * var_rms + dc_dtheta * dt_dtheta * var_theta +
      (dc_drms * dt_dtheta + dt_drms * dc_dtheta) * cov_rms_theta;

    da.fermi2.sigma_c = sqrt(max(var_c, 0.0));
    da.fermi2.sigma_t = sqrt(max(var_t, 0.0));
    da.fermi2.cov_c_t = cov_c_t;
  } 
  // output final optimisation values to LOG
  array<double, 2> final_ct_params = da.getFermi2Femto(CT);
  array<double, 2> final_polar_params = da.getFermi2Femto(POLAR);
  LOG(INFO) << "minimised with MSE: "<< MSE << "\n";
  LOG(INFO) << "(c , t) fermi parameters: "<< final_ct_params[0] << "+/-" << da.fermi2.sigma_c << ", " << final_ct_params[1] << "+/-" << da.fermi2.sigma_t <<" \n";
  LOG(INFO) << "cov(c, t) = " << da.fermi2.cov_c_t << "\n";
  LOG(INFO) << "(rms_radius, theta) fermi parameters: "<< final_polar_params[0] << "+/-" << da.fermi2.sigma_rms_radius << ", " << final_polar_params[1] <<" \n";
  LOG(INFO) << "2pF optimisation completed in " << opt_time << " seconds\n";


}


void optFermi2(DiracAtom & da, const string algo, double & opt_time) {
  if ((algo=="lm")|| (algo=="ls")) {
    ceresOptimizeFermiParameters(da, opt_time, algo);
  } else {
    cout << "Invalid 2pF optimisation algorithm choice for minimsation\n";
    cout << "please use \"lm\", \" or ls\" (default is \"lm\") \n";
    cout << "You used: \""<<algo<<"\" \n";
    cout << "Quitting...\n";
    LOG(ERROR) << "Invalid 2pF optimisation algorithm choice for minimsation: \""<<algo<<"\"\n";
  }
}


void runFermiModelOptimisation(MuDiracInputFile & config, const int & argc, char * argv[], DiracAtom &da, const string & seed) {

  string min_2pF_algo;
  Fermi2CoordinateSystem coord_system_2pF;
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

void ceresOptimizeFermiParameters(DiracAtom & da, double & opt_time, const string & algo) {

  LOG(INFO) << "optimizing fermi parameters using the ceres software \n" ;
  // Get initial guess
  array<double, 2> fermi_coords = da.getFermi2Femto(da.coord_system);
  double  c1 = fermi_coords[0], c2 = fermi_coords[1];
  double c2_upper_bound;
  if (da.coord_system == CT) {
    c2_upper_bound = 3;
  } else {
    c2_upper_bound = M_PI / 4.0;
  }

  // define the cost function
  int num_residuals = da.xr_energies.size();
  ceres::Problem problem;
  auto * cost_function =
    new ceres::NumericDiffCostFunction<CostFunctor, ceres::CENTRAL, ceres::DYNAMIC, 1, 1>(new CostFunctor(da),ceres::TAKE_OWNERSHIP,num_residuals);
  problem.AddResidualBlock(cost_function, nullptr, &c1, &c2);

  // set minimisation algorithm
  ceres::MinimizerType minimizer;
  if (algo == "lm") {   // levenberg marquardt
    minimizer = ceres::TRUST_REGION;

    // set bounds for problem
    problem.SetParameterLowerBound(&c1,0,0);
    problem.SetParameterLowerBound(&c2,0,0);
    problem.SetParameterUpperBound(&c1,0, 7);
    problem.SetParameterUpperBound(&c2,0,c2_upper_bound);
    LOG(INFO) << "optimizing fermi parameters using the levenberg marquardt algorithm" ;
  } else if (algo == "ls") {  // bfgs
    minimizer = ceres::LINE_SEARCH;
    LOG(INFO) << "optimizing fermi parameters using the bfgs algorithm" ;
  }

  // start time of minimisation
  chrono::high_resolution_clock::time_point opt_t0, opt_t1;
  opt_t0 = chrono::high_resolution_clock::now();


  // set options and solve minimisation
  ceres::Solver::Options options;
  options.minimizer_type = minimizer;
  options.gradient_tolerance =1e-5;
  options.parameter_tolerance = 1e-5;
  options.min_relative_decrease = 1e-2;
  options.linear_solver_type = ceres::DENSE_QR;

  options.minimizer_progress_to_stdout = true;

  options.line_search_interpolation_type = ceres::QUADRATIC;
  options.max_num_line_search_step_size_iterations = 10;
  options.line_search_sufficient_curvature_decrease = 0.90;
  options.max_line_search_step_contraction = 0.1;
  options.min_line_search_step_contraction = 0.5;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  opt_t1 = chrono::high_resolution_clock::now();
  std::cout << summary.FullReport() << "\n";
  opt_time = chrono::duration_cast<chrono::milliseconds>(opt_t1 - opt_t0).count() / 1.0e3;


  double sigma_c1 = 0.0, sigma_c2 = 0.0, cov_c1_c2 = 0.0;
  if (summary.termination_type == ceres::CONVERGENCE) {
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::DENSE_SVD;
    cov_options.null_space_rank = -1;
    ceres::Covariance covariance(cov_options);

    std::vector<std::pair<const double*, const double*>> covariance_blocks;
    covariance_blocks.push_back(std::make_pair(&c1, &c1));
    covariance_blocks.push_back(std::make_pair(&c2, &c2));
    covariance_blocks.push_back(std::make_pair(&c1, &c2));
    
    if (covariance.Compute(covariance_blocks, &problem)) {
      double cov_c1c1 = 0.0, cov_c2c2 = 0.0, cov_c1c2_val = 0.0;
      covariance.GetCovarianceBlock(&c1, &c1, &cov_c1c1);
      covariance.GetCovarianceBlock(&c2, &c2, &cov_c2c2);
      covariance.GetCovarianceBlock(&c1, &c2, &cov_c1c2_val);
      
      sigma_c1 = sqrt(cov_c1c1);
      sigma_c2 = sqrt(cov_c2c2);
      cov_c1_c2 = cov_c1c2_val;

      if (cov_c1c1 <= 0 || cov_c2c2 <= 0) {
        LOG(WARNING) << "Covariance diagonal non-positive; uncertainties unreliable\n";
      }

    } else {
      LOG(WARNING) << "Failed to compute covariance for 2pF parameters\n";
    }
  } else {
    LOG(WARNING) << "2pF optimisation did not converge successfully, final parameters may not be optimal \n";
  }

  finaliseFermi2(da, da.coord_system, c1, c2, opt_time, summary.final_cost*2.0/double(num_residuals), sigma_c1, sigma_c2, cov_c1_c2);

}
