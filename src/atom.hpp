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

#include <math.h>
#include <vector>
#include <tuple>
#include <map>
#include "utils.hpp"
#include "constants.hpp"
#include "hydrogenic.hpp"
#include "integrate.hpp"
#include "boundary.hpp"

class State
{
public:
  bool init = false;
  int nodes = 0;
  double E = 0;

  State();
};

class SchroState : State
{
public:
  vector<double> R;
  int l;

  SchroState(const SchroState &s);
};

class DiracState : State
{
public:
  vector<double> Q;
  vector<double> P;
  int k;

  DiracState();
  DiracState(const DiracState &s);
};

class Atom
{
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
  Atom(double Z_in = 1, double m_in = 1, double A_in = -1, double R_in = -1);

  // Basic getters
  double getZ() { return Z; };
  double getA() { return A; };
  double getm() { return m; };
  double getmu() { return mu; };
  double getR() { return R; };
  double getN() { return N; };

  // Grid and potential
  void setGrid(double r0_in = 1e-2, double r1_in = 5e1, int N_in = 1000);
  void setBackgroundCharge(vector<double> bkgQ_in);
  vector<double> getGrid(bool log = false);
  vector<double> getPotential();
};

class DiracAtom : public Atom
{
private:
  // Eigenstates
  map<tuple<int, int, bool>, DiracState *> states;

public:
  DiracAtom(double Z_in = 1, double m_in = 1, double A_in = -1, double R_in = -1);

  void resetStates();
  void calcState(int n, int l, bool s, bool force = false);
  void calcAllStates(int max_n, bool force = false);
  DiracState getState(int n, int l, bool s);
};