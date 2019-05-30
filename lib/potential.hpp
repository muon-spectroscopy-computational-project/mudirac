/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * potential.hpp
 * 
 * Classes for various types of potentials - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <cmath>
#include <vector>
#include <stdexcept>
#include "integrate.hpp"

using namespace std;

#ifndef MUDIRAC_POTENTIAL
#define MUDIRAC_POTENTIAL

class Potential
{
public:
  virtual double V(double r){};
};

class CoulombPotential : Potential
{
public:
  CoulombPotential(double Z = 1.0, double R = -1);
  double V(double r) override;

protected:
  double R, R3, VR, Z;
};

#endif