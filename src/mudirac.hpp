/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * mudirac.hpp
 *
 * Main binary of the program (header)
 *
 * @author Simone Sturniolo
 * @version 1.0.1 20/03/2020
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
#include <dlib/global_optimization.h>
#include "../lib/config.hpp"
#include "../lib/experiment.hpp"
#include "../lib/atom.hpp"
#include "../lib/output.hpp"
#include "../lib/utils.hpp"
#include "../lib/elements.hpp"
#include "../lib/constants.hpp"
#include "../lib/debugtasks.hpp"

#include "../vendor/aixlog/aixlog.hpp"

#define PROJECT_VERSION "1.0.1"

using namespace std;

struct TransLineSpec {
  int n1, n2;
  int l1, l2;
  bool s1, s2;
};

vector<TransitionData> getAllTransitions(vector<TransLineSpec> transqnums, DiracAtom da);
vector<TransLineSpec> parseXRLines(MuDiracInputFile config);

typedef dlib::matrix<double,0,1> column_vector;

int iteration_counter_2pF = 0;


/**
 *
 *
 * @brief  A function which finds the default c, t rms radius and theta of an atom as initial values for optimisation
 * @note   This function determines initial values 2pF domain coordinates as a vector with either the ct or polar coordinates
 * based on the coordinate system. The inital values are extracted from the dirac atom, which are either supplied by the user or from the default
 * parameters taken from marinova tables.
 *
 *
 * @param da:     previous dirac atom to be reconfigured
 *  * @param coord_system: "ct" or "polar" coordinate system used to configure the nuclear model

 * @retval initial fermi parameters: in coordinates defined by the coordinate system
 *
 */
void init2pFModelParams(DiracAtom & da, const string coord_system, column_vector & init_params);

/**
 *
 *
 * @brief  A function whichs updates 2pF parameters in a config object and creates a new dirac atom
 * @note   A function which takes polar fermi parameters and updates a config object by reference to generate a dirac atom
 * with a new 2pF model. The convetional and polar fermi parameters are also passed by reference into a structure to contain
 * the optimal fermi parameters in a minimisation process.
 *
 * @param m:    polar fermi parameters (rms_radius, theta)
 * @param coord_system: "ct" or "polar" coordinate system used to configure the nuclear model
 * @param config:     config object for MuDirac
 * @param da:     previous dirac atom to be reconfigured
 * @param fermi_parameters:      data structure to contain the polar and conventional fermi parameters and the MSE when calcualted as part
 * of optimisation
 * @retval None
 *
 */
void configureNuclearModel(const column_vector& m,const string coord_system, MuDiracInputFile &config, DiracAtom & da, OptimisationData &fermi_parameters);


/**
 *
 *
 * @brief  A function which calculates the Mean Square Error between a MuDirac simulation with a 2pF model and experimental measurements.
 * @note   A function which calculates the non linear Mean Square Error between MuDirac simualted energies and energies experimentally measured.
 * This function configures a dirac atom with new fermi parameters and calculates the square error between mudirac and experimental values for each muonic xray transition.
 * The Mean Square error is taken over all measured transitions provided.
 *
 * @param m:    polar fermi parameters (rms_radius, theta)
 * @param coord_system: "ct" or "polar" coordinate system used to configure the nuclear model
 * @param config:     config object for MuDirac
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 *
 * @retval MSE: Mean Square error between transition energies calculated by MuDirac and measured experimentally.
 *
 */
double calculateMSE(const column_vector& m,const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);


/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib bgfs optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. The starting rms radius for each atom is taken from the
 * default values for the Dirac atom. The starting theta value is set at 0.3 radians.
 *
 *
 * @param config:     config object for MuDirac updated in every iteration
 * @param coord_system: "ct" or "polar" coordinate system the optimisation routine optimises over to configure the nuclear model of each iteration.
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void optimizeFermiParameters(MuDiracInputFile &config,const string coord_system, DiracAtom & da,const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters, double & opt_time);



/**
 *
 *
 * @brief  A function which finds the derivative with respect to the fermi parameters of the mean square error of MuDirac energies compared to experimental energies for selecte transitions. 
 * @note   This function uses the dlib derivative method to numerically calculate the derivate of the calculateMSE function. The derivative value is returned as a dlib column vector and can be used in minimisation routines.
 *
 * @param m: fermi parameters in ct or polar coordinates
 * @param coord_system: "ct" or "polar" determines the coordinate system the derivative function will use.
 * @param config:     config object for MuDirac updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
*
 *
 * @retval derivative: dlib::column_vector (length 2)
 *
 */
const column_vector MSE_2pF_derivative(const column_vector &m, const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);


/**
 *
 *
 * @brief  A function which finds the hessian with respect to the fermi parameters of the mean square error of MuDirac energies compared to experimental energies for selecte transitions. 
 * @note   This function uses a central differences finite difference method on the MSE_2pF_derivative to find the hessian with respect to the fermi parameters.
 * The Hessian is left un symmetrised as this yielded better optimisation results. 
 *
 * @param m: fermi parameters in ct or polar coordinates
 * @param coord_system: "ct" or "polar" determines the coordinate system the hessian function will use.
 * @param config:     config object for the initial dirac atom
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
*
 *
 * @retval hessian matrix: dlib::matrix (2x2)
 *
 */
dlib::matrix<double> MSE_2pF_hessian(const column_vector &m, const string coord_system, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);

/**
 * 
 * @brief functor class which allows the means square error of the energies simulated by mudirac compared with experimental energies to be optimised with respect to the fermi parameters.
 * @note This class uses a numerically defined derivative and hessian. This functor exists to be used in dlib::find_min_trust_function.
 */
class opt_2pF_model
{
  /*!
    This object is a "function model" which can be used with the
    find_min_trust_region() routine.  
  !*/

  public:
    // common parameters for the objective, derivative and hessian functions
    typedef ::column_vector column_vector;
    typedef dlib::matrix<double> general_matrix;
    MuDiracInputFile config;
    vector<TransLineSpec> transqnums;
    vector<string> xr_lines_measured;
    vector<double> xr_energies;
    vector<double> xr_errors;
    string coord_sys;

    // constructor
    opt_2pF_model(MuDiracInputFile cfg,const string coord_system, const vector<TransLineSpec> tqn, const vector<string> xr_lines, const vector<double> xr_e, const vector<double> xr_er){
      config = cfg;
      transqnums = tqn;
      xr_lines_measured = xr_lines;
      xr_energies = xr_e;
      xr_errors = xr_er;
      coord_sys = coord_system;
    }

    double operator() (
      const column_vector& x
    )  const {return calculateMSE(x, coord_sys, config, transqnums, xr_lines_measured, xr_energies, xr_errors);}

    // function for the dlib minisation routine to get the derivative and hessian
    void get_derivative_and_hessian (
      const column_vector& x,
      column_vector& der,
      general_matrix & hess
    ) const 
    {
      der = MSE_2pF_derivative(x,coord_sys, config,  transqnums, xr_lines_measured, xr_energies, xr_errors);
      hess = MSE_2pF_hessian(x, coord_sys, config,  transqnums, xr_lines_measured, xr_energies, xr_errors);
    }
};

/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib trust region optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. In polar coordinates, the starting rms radius for each atom is taken from the
 * default values for the Dirac atom. The starting theta value is set at 0.3 radians.
 *@param opt_obj: object defined for dlib optmisation routines which has the experimental values as attributes already defined.
 The operator function of this functor is the objective function for minimisation and the functor also has functions for the derivative and hessian.
 * @param coord_system:   either "ct" or "polar" sets the coordinate system ised by the minisation algorithm.
 * @param config:     config object for MuDirac updated in every iteration
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void optimizeFermiParameters(opt_2pF_model &opt_obj, const string coord_system, MuDiracInputFile & config, DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time);


/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib global optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. The domain is bounded close the rms radius values near those
 * given by default. The optimization is limited to 10 minutes.
 *
 *
 * @param config:     config object for MuDirac updated in every iteration
 * @param coord_system: "ct" or "polar" coordinate system the optimisation routine optimises over to configure the nuclear model of each iteration.
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void globalOptimizeFermiParameters(MuDiracInputFile &config, const string coord_system, DiracAtom & da, const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters, double & opt_time);
