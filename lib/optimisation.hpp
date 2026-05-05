/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * optimisation.cpp
 *
 * Functions for optimising the dirac atom nuclear model parameters - header file
 *
 * @author Milan Kumar
 * @version 1.0 30/06/2025
 */
#include "ceres/ceres.h"
#undef LOG
#include "atom.hpp"
#include "config.hpp"
#include "experiment.hpp"
#include "output.hpp"



/**
 * @brief a function which outputs optimisation results to LOG and sets finalised optimisation data values in a dirac atom.
 * @note This function outputs the fermi parameters in the selected coordinate system to log along with the optimisation time
 * and the MSE. The fermi parameters in the other coordinate system are them also stored in the fermi2 datastructure of the Dirac atom.
 * @param da: Dirac atom after optimisation
 * @param coord_sys: the selected coordinate system of the optimisation
 * @param final_c1: the optimised first fermi parameters
 * @param final_c2: the optimised second fermi parameters
 * @param opt_time: the time taken to perform just the optimisation algorithm
 * @param MSE: the mean square error of new xray energies calculated using the optimised parameters
 *  with respect to the experimental energies.
 */
void finaliseFermi2(DiracAtom & da, Fermi2CoordinateSystem coord_sys, double final_c1, double final_c2, double opt_time, double MSE);



/**
 * @brief: Function which wraps the other 2pF optimisation functions and validates the algorithm choice
 *
 * @param da: dirac atom used for the minimisation calculations
 * @param algo: the selected minimisation algorithm
 * @param best_fermi_parameters: stores the values of c, t, rms radius, theta and MSE of the optimal fermi parameters.
 * @param opt_time: time taken to complete the optimisation.
 */
void optFermi2(DiracAtom & da, const string algo, double & opt_time);


/**
 * @brief  Highest level function which runs 2pF optimisation
 * @note   validates both mudirac input files and parses them to get the
 * experimental transition lines and energies. This function then uses other
 * lower level functions to run the optimisation with algorithm and coordinate
 * choices defined in the config file. Finally the outputs are written to file
 * with a name based on the seed.
 *
 * @param  config: object containing parameters which allow optimisation
 * @param argc: number of arguments passed with mudirac for validation
 * @param argv: arguments passed to mudirac including config filenames
 * @param da: the pre configured dirac atom
 * @param seed: seed name for the output file.
 * @retval None
 */
void runFermiModelOptimisation(MuDiracInputFile & config, const int & argc, char * argv[], DiracAtom &da, const string & seed);

// Ceres optimisation problem things
struct CostFunctor {
  mutable DiracAtom ma;

  CostFunctor(DiracAtom &da_in) : ma(da_in) {}


  bool operator() (const double* const c1, const double* const c2, double* residual) const {
    ma.setFermi2Femto(c1[0], c2[0], ma.coord_system);
    vector<TransitionData> transitions_iteration = ma.getAllTransitions();

    if ((int)transitions_iteration.size() < (int)ma.xr_lines_measured.size())
      return false; // Not enough transitions to compare with experimental data

    for (int k = 0; k < transitions_iteration.size(); ++k) {
      // calculate transition energy and rate
      double dE = (transitions_iteration[k].ds2.E - transitions_iteration[k].ds1.E);
      double transition_energy = dE / Physical::eV;
      residual[k] = (transition_energy - ma.xr_energies[k])/ma.xr_errors[k];
    }
    return true;
  }
};

/**
 * @brief uses the ceres optimisation library to perform the optimisation of the 2pF parameters
 * @note This funciton uses the ceres solver to minimize the Mean Square Error (MSE) between simulated energies
 * calculated using 2pf parameters and experimental xray transition energies
 * @param da: dirac atom with 2pf model to be optimised
 * @param opt_time: time taken to complete the optimisation.
 * @param algo: the selected minimisation algorithm
 */
void ceresOptimizeFermiParameters(DiracAtom & da, double & opt_time,const  string & algo);