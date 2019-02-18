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
#include "utils.hpp"
#include "constants.hpp"

class Atom
{
  private:
    // Fundamental properties
    double Z, A;      // Nuclear charge and mass
    double m, mu;     // Mass of the orbiting particle (e.g. muon, electron) and effective mass of the system
    double R;         // Nuclear radius
    vector<double> V; // Radial potential
    // Grid
    int N = 1000;                // Number of points
    double r0 = 1e-2, r1 = 5e1;  // Extremes
    vector<vector<double>> grid; // Grid (x, r)

    void recalcPotential(); // Recalculate V

  public:
    Atom(double Z_in = 1, double m_in = 1, double A_in = -1, double R_in = -1);

    void setGrid(double r0_in = 1e-2, double r1_in = 5e1, int N_in = 1000);
    vector<double> getGrid(bool log = false);
    vector<double> getPotential();
};