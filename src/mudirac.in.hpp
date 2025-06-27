/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * mudirac.hpp
 *
 * Main binary of the program (header)
 *
 * @author Simone Sturniolo
 * @version ${PROJECT_VERSION} 20/03/2020
 */

#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <dlib/optimization.h>
#include "../lib/config.hpp"
#include "../lib/experiment.hpp"
#include "../lib/atom.hpp"
#include "../lib/output.hpp"
#include "../lib/utils.hpp"
#include "../lib/elements.hpp"
#include "../lib/constants.hpp"
#include "../lib/debugtasks.hpp"

#include "../vendor/aixlog/aixlog.hpp"

#define PROJECT_VERSION "${PROJECT_VERSION}"

using namespace std;


typedef dlib::matrix<double,0,1> column_vector;

double calculateMSE(const column_vector& m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);
void optimizeFermiParameters(MuDiracInputFile &config, DiracAtom & da,const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters);


const column_vector MSE_2pF_derivative( const column_vector &m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);

dlib::matrix<double> MSE_2pF_hessian(const column_vector &m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);


class opt_2pF_model
{
  /*!
    This object is a "function model" which can be used with the
    find_min_trust_region() routine.  
  !*/

  public:
    typedef ::column_vector column_vector;
    typedef dlib::matrix<double> general_matrix;
    MuDiracInputFile config;
    vector<TransLineSpec> transqnums;
    vector<string> xr_lines_measured;
    vector<double> xr_energies;
    vector<double> xr_errors;

    // constructor
    opt_2pF_model(MuDiracInputFile cfg, const vector<TransLineSpec> tqn, const vector<string> xr_lines, const vector<double> xr_e, const vector<double> xr_er){
      config = cfg;
      transqnums = tqn;
      xr_lines_measured = xr_lines;
      xr_energies = xr_e;
      xr_errors = xr_er;
    }

    double operator() (
      const column_vector& x
    ) const {return calculateMSE(x, config, transqnums, xr_lines_measured, xr_energies, xr_errors);}

    void get_derivative_and_hessian (
      const column_vector& x,
      column_vector& der,
      general_matrix & hess
    ) const
    {
      der = MSE_2pF_derivative(x, config, transqnums, xr_lines_measured, xr_energies, xr_errors);
      hess = MSE_2pF_hessian(x, config, transqnums, xr_lines_measured, xr_energies, xr_errors);
    }
};