/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * atom.hpp
 * 
 * Atom classes, designed to be a high level interface to solution of quantum mechanical equations - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <cmath>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include "utils.hpp"
#include "constants.hpp"
#include "hydrogenic.hpp"
#include "integrate.hpp"
#include "boundary.hpp"
#include "elements.hpp"
#include "potential.hpp"

using namespace std;

// Exception classes
class AtomConvergenceException : exception
{
public:
  enum ACEType
  {
    NAN_ENERGY,
    NODES_WRONG,
    MAXIT_REACHED
  };

  AtomConvergenceException(ACEType t = NAN_ENERGY);
  ~AtomConvergenceException(void) {}
  ACEType getType() const throw() { return type; };
  const char *what() const throw() { return msg.c_str(); };

private:
  ACEType type;
  string msg;
};

enum NuclearRadiusModel
{
  POINT,
  SPHERE
};

// Main classes
class State
{
public:
  bool init = false;
  int nodes = 0;
  double E = 0;
  pair<int, int> grid_indices;
  vector<double> grid;
  vector<double> loggrid;
  vector<double> V;

  State();
};

class SchroState : public State
{
public:
  vector<double> R;
  int l;

  SchroState(const SchroState &s);
};

class DiracState : public State
{
public:
  vector<double> Q;
  vector<double> P;
  int k;
  int nodesQ = 0;

  DiracState(int N = 0);
  DiracState(double x0, double x1, int N);
  DiracState(double rc, double dx, int i0, int i1);
  DiracState(const DiracState &s);
};

class Atom
{
public:
  // Tolerances and other details
  double Etol = 1e-10, Esearch = 1.2;
  int maxit = 100;

  // Radius models
  double static sphereNuclearModel(double A);

protected:
  // Fundamental properties
  double Z, A;  // Nuclear charge and mass
  double m, mu; // Mass of the orbiting particle (e.g. muon, electron) and effective mass of the system
  double R;     // Nuclear radius
  // Grid
  double rc = 1.0;   // Central radius
  double dx = 0.005; // Step
  CoulombPotential V;

  vector<double> recalcPotential(vector<double> r); // Recalculate V for a certain grid

public:
  Atom(double Z = 1, double m = 1, double A = -1, NuclearRadiusModel radius_model = POINT,
       double fc = 1.0, double dx = 0.005);

  // Basic getters
  double getZ() { return Z; };
  double getA() { return A; };
  double getm() { return m; };
  double getmu() { return mu; };
  double getR() { return R; };
  double getV(double r) { return V.V(r); };
};

class DiracAtom : public Atom
{
private:
  double restE; // Rest energy
  // Eigenstates
  map<tuple<int, int, bool>, DiracState> states;

public:
  double out_eps = 1e-5;
  double in_eps = 1e-5;
  int min_n = 1000;

  enum GridLimitsFailcode
  {
    OK,
    UNBOUND,
    SMALL_GAMMA,
  };

  DiracAtom(double Z = 1, double m = 1, double A = -1, NuclearRadiusModel radius_model = POINT,
            double fc = 1.0, double dx = 0.005);

  double getRestE() { return restE; };

  pair<int, int> gridLimits(double E, int k, GridLimitsFailcode &failcode);

  void resetStates();
  void calcState(int n, int l, bool s, bool force = false);
  void calcAllStates(int max_n, bool force = false);
  // Convergence
  DiracState convergeState(double E0, int k = -1);
  DiracState getState(int n, int l, bool s);
};