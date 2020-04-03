/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * debugtasks.hpp
 * 
 * Methods or classes used to perform debug specific tasks - header file
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <tuple>
#include <vector>
#include "atom.hpp"
#include "output.hpp"
#include "config.hpp"


#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_DEBUGTASKS
#define MUDIRAC_DEBUGTASKS

void runEdEscan(MuDiracInputFile infile);

#endif