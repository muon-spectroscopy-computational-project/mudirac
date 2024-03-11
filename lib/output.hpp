/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * output.hpp
 *
 * Functions to write output files - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "atom.hpp"
#include "potential.hpp"
#include "constants.hpp"

using namespace std;

#ifndef MUDIRAC_OUTPUT
#define MUDIRAC_OUTPUT

void writeDiracState(DiracState ds, string fname);
void writeTransitionMatrix(TransitionMatrix tmat, string fname);
void writeEConfPotential(EConfPotential epot, string fname);
void writeSimSpec(vector<TransitionData> transitions, double dE, double lw, double expd, string fname);

// Debug tasks
void writeEdEscan(vector<double> Es, vector<double> dEs, vector<int> nodes, string fname="EdEscan.dat");

#endif
