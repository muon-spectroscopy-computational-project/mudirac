/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * output.hpp
 * 
 * Functions to write output files - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <fstream>
#include <string>
#include "atom.hpp"
#include "constants.hpp"

using namespace std;

#ifndef MUDIRAC_OUTPUT
#define MUDIRAC_OUTPUT

void writeDiracState(DiracState ds, string fname);

#endif
