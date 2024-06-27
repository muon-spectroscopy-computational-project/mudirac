/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * state.hpp
 *
 * State classes, designed to store the solution of quantum mechanical equations - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include<vector>
#include "integrate.hpp"

using namespace std;

#ifndef MUDIRAC_STATE
#define MUDIRAC_STATE

class State {
 public:
  bool converged = false; // Used to flag states that are "good" to use
  int nodes = 0;
  double E = 0;
  pair<int, int> grid_indices;
  vector<double> grid;
  vector<double> loggrid;
  vector<double> V;

  State();

  virtual double norm() = 0;
  virtual int size() = 0;
  virtual void resize(int i0, int i1) {
    return;
  };
};

class SchroState : public State {
 public:
  vector<double> R;
  int l;

  SchroState(const SchroState &s);

  int getn();
  double norm() override;
  int size() override {
    return R.size();
  };
  void resize(int i0, int i1) override;
};

class DiracState : public State {
 public:
  vector<double> Q;
  vector<double> P;
  int k;
  int nodesQ = 0;
  double m;

  DiracState(int N = 0);
  DiracState(double x0, double x1, int N);
  DiracState(double rc, double dx, int i0, int i1);
  DiracState(const DiracState &s);

  int getn();
  int getl();
  bool gets();
  double restE();
  double bindingE();
  double norm() override;
  int size() override {
    return Q.size();
  };
  void resize(int i0, int i1) override;
  string name();

  void continuify(TurningPoint tp);
  void findNodes(double tol = 1e-6);
  void normalize();
};
#endif