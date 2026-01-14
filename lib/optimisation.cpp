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


void finaliseFermi2(DiracAtom & da, Fermi2CoordinateSystem coord_sys, double final_c1, double final_c2, double opt_time, double MSE) {

  // ensure final fermi parameters are set
  da.setFermi2Femto(final_c1, final_c2, coord_sys);
  //set all the optimisation parameter values
  da.fermi2.mse = MSE;
  // output final optimisation values to LOG
  array<double, 2> final_ct_params = da.getFermi2Femto(CT);
  array<double, 2> final_polar_params = da.getFermi2Femto(POLAR);
  LOG(INFO) << "minimised with MSE: "<< MSE << "\n";
  LOG(INFO) << "(c , t) fermi parameters: "<< final_ct_params[0] << ", " << final_ct_params[1] <<" \n";
  LOG(INFO) << "(rms_radius, theta) fermi parameters: "<< final_polar_params[0] << ", " << final_polar_params[1] <<" \n";
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
  finaliseFermi2(da, da.coord_system, c1, c2, opt_time, summary.final_cost*2.0/double(num_residuals));

}
