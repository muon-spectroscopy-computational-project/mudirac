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

struct TransLineSpec {
  int n1, n2;
  int l1, l2;
  bool s1, s2;
};

vector<TransitionData> getAllTransitions(vector<TransLineSpec> transqnums, DiracAtom da);
vector<TransLineSpec> parseXRLines(MuDiracInputFile config);

/**
 * 
 * @brief  A function to fit fermi parameters of a 2pF model to experimental data
 * @note   A function which scans over polar coordinates RMS_radius and theta (which map to fermi half density radius and skin thickness),
 * calculating muonic transition energies and comparing them to the experimental transition energies. fermi parameters which generate energies
 * with mean square error < 1 when compared to the experimental energies are listed and the fermi parameters with the lowest MSE is stored as
 * the optimal fermi parameters. All parameters are passed by reference.
 * 
 * @param xr_lines_measured:      labels of xray lines measured provided in the experimental data input file.
 * @param xr_energies:      energies of the measured xrays.
 * @param xr_errors:      energy errors of the measured xrays.
 * @param transqnums:     quantum numbers of the transitions provided in the first mudirac input file.
 * @param valid_fermi_parameters:     a list to contain structures containing valid fermi parameters and the MSE
 * @param optimal_fermi_parameter:      a structure to contain the fermi parameters with the lowest MSE
 * @param optimal_transitions:      transition data calculated using the optimal fermi parameters
 * @param config:     contains the general config of the atom input by the user. fermi parameters are updated with each iteration.
 */
void optimiseFermiParameters(const vector<string> &xr_lines_measured, const vector<double> &xr_energies, const vector<double> &xr_errors, const vector<TransLineSpec> &transqnums, vector<OptimisationData> &valid_fermi_parameters, OptimisationData &optimal_fermi_parameter,vector<TransitionData> &optimal_transitions, MuDiracInputFile &config);
typedef dlib::matrix<double,0,1> column_vector;
double testFunction(const column_vector& m);
double minimise_MSE(const column_vector& m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);