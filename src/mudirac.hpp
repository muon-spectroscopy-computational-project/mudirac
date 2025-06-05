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

/**
 *
 *
 * @brief  A function whichs updates 2pF parameters in a config object and creates a new dirac atom
 * @note   A function which takes polar fermi parameters and updates a config object by reference to generate a dirac atom
 * with a new 2pF model. The convetional and polar fermi parameters are also passed by reference into a structure to contain
 * the optimal fermi parameters in a minimisation process.
 *
 * @param m:    polar fermi parameters (rms_radius, theta)
 * @param config:     config object for MuDirac
 * @param da:     previous dirac atom to be reconfigured
 * @param fermi_parameters:      data structure to contain the polar and conventional fermi parameters and the MSE when calcualted as part
 * of optimisation
 * @retval None
 *
 */
void configureNuclearModel(const column_vector& m, MuDiracInputFile &config, DiracAtom & da, OptimisationData &fermi_parameters);


/**
 *
 *
 * @brief  A function which calculates the Mean Square Error between a MuDirac simulation with a 2pF model and experimental measurements.
 * @note   A function which calculates the non linear Mean Square Error between MuDirac simualted energies and energies experimentally measured.
 * This function configures a dirac atom with new fermi parameters and calculates the square error between mudirac and experimental values for each muonic xray transition.
 * The Mean Square error is taken over all measured transitions provided.
 *
 * @param m:    polar fermi parameters (rms_radius, theta)
 * @param config:     config object for MuDirac
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 *
 * @retval MSE: Mean Square error between transition energies calculated by MuDirac and measured experimentally.
 *
 */
double calculateMSE(const column_vector& m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);


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
 * * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void optimizeFermiParameters(MuDiracInputFile &config, DiracAtom & da,const vector<TransLineSpec> &transqnums, const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, OptimisationData &fermi_parameters);