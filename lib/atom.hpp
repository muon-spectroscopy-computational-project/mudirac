/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * atom.hpp
 *
 * Atom classes, designed to be a high level interface to solution of quantum
 * mechanical equations - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "../vendor/aixlog/aixlog.hpp"
#include "boundary.hpp"
#include "constants.hpp"
#include "econfigs.hpp"
#include "elements.hpp"
#include "hydrogenic.hpp"
#include "integrate.hpp"
#include "potential.hpp"
#include "state.hpp"
#include "utils.hpp"
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

#ifndef MUDIRAC_ATOM
#define MUDIRAC_ATOM
// Error codes
enum AtomErrorCode {
  UNBOUND_STATE = 1, // State is unbound, E - mc^2 > 0
  RMAX_SMALL,        // Upper bound of grid is smaller than turning point
  RMIN_LARGE,        // Lower bound of grid is bigger than turning point
  SMALL_GAMMA,       // Gamma parameter too small (won't happen for known elements)
  NODES_HIGH,        // State has converged but does not contain the right amount of nodes (too many)
  NODES_LOW,         // State has converged but does not contain the right amount of nodes (too few)
};

enum NuclearRadiusModel {
  POINT,
  SPHERE,
  FERMI2
};

// Main classes
class TransitionMatrix {
 public:
  int k1, k2;
  vector<double> m1;
  vector<double> m2;
  vector<vector<double>> T;

  TransitionMatrix() {
    return;
  };
  TransitionMatrix(int k1, int k2);

  double totalRate();
};

struct TransLineSpec {
  int n1, n2;
  int l1, l2;
  bool s1, s2;
};

struct TransitionData {
  string name;
  string sname1;
  string sname2;
  DiracState ds1;
  DiracState ds2;
  TransitionMatrix tmat;
};


/**
 * @brief Data structure to store a set of conventional and polar fermi parameters and related mean square error
 *
 */
struct OptimisationData {
  double rms_radius;
  double theta;
  double mse;
  double fermi_c;
  double fermi_t;
};


class Atom {
 public:
  // Tolerances and other details
  double Etol = 1e-7, Edamp = 0.5;
  double max_dE_ratio = 1e-1;
  double nodetol = 1e-6;
  int maxit_E = 100;
  int maxit_nodes = 100;
  int maxit_state = 100;

  bool write_debug =
    false; // If true, files will be printed to store details of the process

  // Radius models
  double static sphereNuclearModel(int Z, int A);

 protected:

  // Fundamental properties
  int Z, A;                  // Nuclear charge and mass number
  double m, M, mu;           // Mass of the orbiting particle (e.g. muon, electron), of the nucleus, and effective mass of the system
  double R;                  // Nuclear radius
  NuclearRadiusModel rmodel; // Nuclear model
  bool reduced_mass;        // turn on reduced mass

  //Grid
  double rc = 1.0;   // Central radius
  double dx = 0.005; // Step

  //Potential
  CoulombSpherePotential *V_coulomb;

  // Additional potential terms
  bool use_uehling = false;
  UehlingSpherePotential V_uehling;
  bool use_econf = false;
  EConfPotential V_econf;

 public:
   // fermi 2pF parameters
  double fermi_c=0, fermi_t=0;

  Atom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
       double radius = -1, double fc = 1.0,double dx = 0.005, bool reduced_mass = true);

  // Basic getters
  double getZ() {
    return Z;
  };
  double getA() {
    return A;
  };
  double getm() {
    return m;
  };
  double getmu() {
    return mu;
  };
  double getM() {
    return M;
  };
  double getR() {
    return R;
  };
  double getV(double r);
  vector<double> getV(vector<double> r);
  double getrc() {
    return rc;
  };
  double getdx() {
    return dx;
  };

  // Potential getters
  static const uint HAS_UEHLING = 1;
  static const uint HAS_ELECTRONIC = 2;
  uint getPotentialFlags() {
    return HAS_UEHLING * use_uehling + HAS_ELECTRONIC * use_econf;
  };
  CoulombSpherePotential *getPotentialCoulomb() {
    return V_coulomb;
  };
  UehlingSpherePotential getPotentialUehling() {
    return V_uehling;
  };
  EConfPotential getPotentialElectronic() {
    return V_econf;
  };

  void setgrid(double rc, double dx);

  // Additional potential terms get/setters
  void setFermi2(double thickness = Physical::fermi2_t, double fermi2_potential = -1);

  /**
   * @brief sets the 2 parameter fermi model parameters of the atoms nuclear model.
   * @note sets the 2pF parameters in ct or polar coordinates. units are in fm.
   * @param coord_1: half density radius c or rms radius 
   * @param coord_2: skin thickness t or theta
   * @param coord_sys: coordinate system "ct" or "polar"
   * @retval None 
   */
  void setFermi2(const double coord_1, const double coord_2, const string coord_sys);
  
  /**
   * @brief gets the 2 parameter fermi model parameters of the atoms nuclear model.
   */
  vector<double> getFermi2(const string coord_sys);

  bool getUehling() {
    return use_uehling;
  };
  void setUehling(bool s, int usteps = 1000, double cut_low = 0,
                  double cut_high = INFINITY);
  // Electronic background
  void setElectBkgConfig(bool s, ElectronicConfiguration econf,
                         double rho_eps = 1e-5, double max_r0 = -1,
                         double min_r1 = -1);

  // Clear computed states
  virtual void reset() {};
};

class DiracAtom : public Atom {
 private:
  double restE; // Rest energy
  // Eigenstates
  map<tuple<int, int, bool>, DiracState> states;
  int idshell = -1;

 public:
  double out_eps = 1e-5;
  double in_eps = 1e-5;
  int min_n = 1000;

  // 2pF optimisation attributes
  int  iteration_counter_2pF =0;
  string coord_system;
  vector<TransLineSpec> transqnums;
  vector<string> xr_lines_measured;
  vector<double> xr_energies;
  vector<double> xr_errors;
  OptimisationData fermi_parameters;
  double  opt_time;

  DiracAtom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
            double radius = -1, double fc = 1.0, double dx = 0.005, int ideal_minshell = -1, bool reduced_mass = true);

  double getRestE() {
    return restE;
  };

  void reset() override;

  void calcState(int n, int l, bool s, bool force = false);
  void calcAllStates(int max_n, bool force = false);

  // Convergence
  pair<double, double> energyLimits(int nodes = 0, int k = -1);
  pair<int, int> gridLimits(double E, int k);
  DiracState initState(double E, int k = -1);
  void integrateState(DiracState &state, TurningPoint &tp);
  void integrateState(DiracState &state, TurningPoint &tp, double &dE);
  void convergeNodes(DiracState &state, TurningPoint &tp, int targ_nodes,
                     double &minE, double &maxE);
  void convergeE(DiracState &state, TurningPoint &tp, double &minE,
                 double &maxE);
  DiracState convergeState(int n = 1, int k = -1);
  DiracState getState(int n, int l, bool s);
  TransitionMatrix getTransitionProbabilities(int n1, int l1, bool s1, int n2,
      int l2, bool s2, bool approx_j0 = false);

  // Full energy calculation

  vector<TransitionData> getAllTransitions();

  // optimisation

  void setExpOptData(string coord_sys, const vector<string> xr_lines, const vector<double> xr_e, const vector<double> xr_err, OptimisationData fermi_ps){
    iteration_counter_2pF =0;
    coord_system = coord_sys;
    xr_lines_measured = xr_lines;
    xr_energies = xr_e;
    xr_errors = xr_err;
    fermi_parameters = fermi_ps;
    opt_time = 0.0;
  }

  double calculateMSE(double coord_1, double coord_2);
};

// A class used mainly for debugging purposes, works as DiracAtom but uses only
// the analytical hydrogen-like solution
class DiracIdealAtom : public DiracAtom {
 public:
  DiracIdealAtom(int Z = 1, double m = 1, int A = -1, NuclearRadiusModel radius_model = POINT,
                 double radius = -1, double fc = 1.0, double dx = 0.005);
};


#endif