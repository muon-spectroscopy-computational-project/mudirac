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
#include "state.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "hydrogenic.hpp"
#include "integrate.hpp"
#include "boundary.hpp"
#include "elements.hpp"
#include "econfigs.hpp"
#include "potential.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_ATOM
#define MUDIRAC_ATOM
// Error codes
enum AtomErrorCode
{
  UNBOUND_STATE = 1, // State is unbound, E - mc^2 > 0
  RMAX_SMALL,        // Upper bound of grid is smaller than turning point
  RMIN_LARGE,        // Lower bound of grid is bigger than turning point
  SMALL_GAMMA,       // Gamma parameter too small (won't happen for known elements)
  NODES_HIGH,        // State has converged but does not contain the right amount of nodes (too many)
  NODES_LOW,         // State has converged but does not contain the right amount of nodes (too few)
};

enum NuclearRadiusModel
{
  POINT,
  SPHERE
};

// Main classes
class TransitionMatrix
{
public:
  int k1, k2;
  vector<double> m1;
  vector<double> m2;
  vector<vector<double>> T;

  TransitionMatrix(int k1, int k2);

  double totalRate();
};

class Atom
{
public:
  // Tolerances and other details
  double Etol = 1e-7, Edamp = 0.5;
  double max_dE_ratio = 1e-1;
  double nodetol = 1e-6;
  int maxit_E = 100;
  int maxit_nodes = 100;
  int maxit_state = 100;

  // Radius models
  double static sphereNuclearModel(int Z, int A);

protected:
  // Fundamental properties
  int Z, A;        // Nuclear charge and mass number
  double m, M, mu; // Mass of the orbiting particle (e.g. muon, electron), of the nucleus, and effective mass of the system
  double R;        // Nuclear radius
  // Grid
  double rc = 1.0;   // Central radius
  double dx = 0.005; // Step
  // Potential
  CoulombSpherePotential V;
  // Additional potential terms
  bool use_uehling = false;
  UehlingSpherePotential V_uehling;
  bool use_econf = false;
  EConfPotential V_econf;

public:
  Atom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
       double fc = 1.0, double dx = 0.005);

  // Basic getters
  double getZ() { return Z; };
  double getA() { return A; };
  double getm() { return m; };
  double getmu() { return mu; };
  double getM() { return M; };
  double getR() { return R; };
  double getV(double r);
  vector<double> getV(vector<double> r);
  double getrc() { return rc; };
  double getdx() { return dx; };

  // Potential getters
  static const uint HAS_UEHLING = 1;
  static const uint HAS_ELECTRONIC = 2;
  uint getPotentialFlags() { return HAS_UEHLING * use_uehling + HAS_ELECTRONIC * use_econf; };
  CoulombSpherePotential getPotentialCoulomb() { return V; };
  UehlingSpherePotential getPotentialUehling() { return V_uehling; };
  EConfPotential getPotentialElectronic() { return V_econf; };

  void setgrid(double rc, double dx);

  // Additional potential terms get/setters
  bool getUehling() { return use_uehling; };
  void setUehling(bool s, int usteps = 1000, double cut_low = 0, double cut_high = INFINITY);
  // Electronic background
  void setElectBkgConfig(bool s, ElectronicConfiguration econf, double rho_eps = 1e-5, double max_r0 = -1, double min_r1 = -1);

  // Clear computed states
  virtual void reset(){};
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

  DiracAtom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
            double fc = 1.0, double dx = 0.005);

  double getRestE() { return restE; };

  void reset() override;

  void calcState(int n, int l, bool s, bool force = false);
  void calcAllStates(int max_n, bool force = false);

  // Convergence
  pair<double, double> energyLimits(int nodes = 0, int k = -1);
  pair<int, int> gridLimits(double E, int k);
  DiracState initState(double E, int k = -1);
  void integrateState(DiracState &state, TurningPoint &tp);
  void integrateState(DiracState &state, TurningPoint &tp, double &dE);
  void convergeNodes(DiracState &state, TurningPoint &tp, int targ_nodes, double &minE, double &maxE);
  void convergeE(DiracState &state, TurningPoint &tp, double &minE, double &maxE);
  DiracState convergeState(int n = 1, int k = -1);
  DiracState getState(int n, int l, bool s);
  TransitionMatrix getTransitionProbabilities(int n1, int l1, bool s1, int n2, int l2, bool s2, bool approx_j0 = false);
};

// A class used mainly for debugging purposes, works as DiracAtom but uses the analytical hydrogen-like solution
class DiracIdealAtom : public DiracAtom
{
public:
  DiracIdealAtom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
                 double fc = 1.0, double dx = 0.005);

private:
  DiracState convergeState(int n = 1, int k = -1);
};

#endif