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
void configureNuclearModel(const column_vector& m, MuDiracInputFile &config, DiracAtom & da, OptimisationData &fermi_parameters);
double minimise_MSE(const column_vector& m, MuDiracInputFile config, const vector<TransLineSpec> transqnums, const vector<string> xr_lines_measured, const vector<double> xr_energies, const vector<double> xr_errors);
