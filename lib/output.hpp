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
#include "config.hpp"

using namespace std;

#ifndef MUDIRAC_OUTPUT
#define MUDIRAC_OUTPUT

void writeFermiParameters(DiracAtom &da, const double opt_time, string fname, int output_precision);
void writeDiracState(DiracState ds, string fname, int output_precision=-1);
void writeTransitionMatrix(TransitionMatrix tmat, string fname);
void writeEConfPotential(EConfPotential epot, string fname);
void writeSimSpec(vector<TransitionData> transitions, double dE, double lw, double expd, string fname);

// outputs refactored from main loop

/**
 * @brief  Prints initial log message
 * @retval None
 */
void printInitLogMessage();

/**
 * @brief  Writes final output files of default MuDirac calculation
 * @note   Writes output files for the transitions, energies and intensities
 *  with options to sort by energy and output full transition matrices defined
 * in the config.
 * @param seed: seed for the output file
 * @param  config:  config object used to make a dirac atom and run an EdEscan
 * @param da: dirac atom containing data about the atom
 * @param transitions: calculated transitions which need to be output
 * @retval None
 */
void writeOutputFiles(string seed, MuDiracInputFile & config, DiracAtom & da, vector<TransitionData> & transitions);
// Debug tasks
void writeEdEscan(vector<double> Es, vector<double> dEs, vector<int> nodes, string fname="EdEscan.dat");

#endif
