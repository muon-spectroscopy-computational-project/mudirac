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

using namespace std;

// Exception classes
class AtomConvergenceException : exception
{
public:
  static const int NAN_ENERGY = 0;
  static const int NODES_WRONG = 1;
  static const int MAXIT_REACHED = 2;

  AtomConvergenceException(int t = NAN_ENERGY);
  ~AtomConvergenceException(void) {}
  int getType() const throw() { return type; };
  const char *what() const throw() { return msg.c_str(); };

private:
  int type;
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
  vector<double> grid;
  vector<double> loggrid;

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
  double Z, A;         // Nuclear charge and mass
  double m, mu;        // Mass of the orbiting particle (e.g. muon, electron) and effective mass of the system
  double R;            // Nuclear radius
  vector<double> V;    // Radial potential
  vector<double> bkgQ; // Background charge
  // Grid
  int N = 1000;                // Number of points
  double r0 = 1e-2, r1 = 5e1;  // Extremes
  vector<vector<double>> grid; // Grid (x, r)
  double dx;                   // Step

  void recalcPotential(); // Recalculate V

public:
  Atom(double Z_in = 1, double m_in = 1, double A_in = -1, NuclearRadiusModel radius_model = POINT);

  // Basic getters
  double getZ() { return Z; };
  double getA() { return A; };
  double getm() { return m; };
  double getmu() { return mu; };
  double getR() { return R; };
  double getN() { return N; };

  // Grid and potential
  void setGrid(double r0_in = 1e-2, double r1_in = 5e1, int N_in = 1000);
  void setGridRelative(double r0_in = 1e-2, double r1_in = 5e1, int N_in = 1000);
  void setBackgroundCharge(vector<double> bkgQ_in);
  vector<double> getGrid(bool log = false);
  vector<double> getPotential();
};

class DiracAtom : public Atom
{
private:
  // Eigenstates
  map<tuple<int, int, bool>, DiracState> states;

  // Convergence
  DiracState convergeState(double E0, int k = -1);

public:
  DiracAtom(double Z_in = 1, double m_in = 1, double A_in = -1, NuclearRadiusModel radius_model = POINT);

  void resetStates();
  void calcState(int n, int l, bool s, bool force = false);
  void calcAllStates(int max_n, bool force = false);
  DiracState getState(int n, int l, bool s);
};