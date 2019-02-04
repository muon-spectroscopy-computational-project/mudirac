/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * hydrogenic.hpp
 * 
 * Computes known hydrogenic atom solutions for reference - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <vector>
#include "constants.hpp"

using namespace std;

double hydrogenicSchroEnergy(double Z=1.0, double mu=1.0, int n=1);