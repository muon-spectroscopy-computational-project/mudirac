/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * output.cpp
 *
 * Functions to write output files
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "output.hpp"

/**
  * @brief  Write a DiracState object to a text file
  * @note   Write down a full DiracState object in an ASCII file,
  * saving potential and P and Q component values
  *
  * @param  ds:         DiracState to save
  * @param  fname:      Filename
  * @retval None
 */
void writeDiracState(DiracState ds, string fname) {
  ofstream out(fname);

  // Start with writing a header
  out << "#####################################################\n";
  out << "# DiracState with n = " << ds.getn() << ", l = " << ds.getl() << ", s = " << ds.gets() << "\n";
  out << "# E = " << ds.bindingE() / Physical::eV << " + mc^2 = " << ds.E / Physical::eV << " eV\n";
  out << "# nodes = " << ds.nodes << ", " << ds.nodesQ << "\n";
  out << "#####################################################\n";

  // Now write down the content
  for (int i = 0; i < ds.grid.size(); ++i) {
    out << ds.grid[i] << '\t' << ds.V[i] / Physical::eV << '\t' << ds.P[i] << '\t' << ds.Q[i] << '\n';
  }

  out.close();
}

/**
 * @brief  Write a TransitionMatrix object to a text file
 * @note   Write down a full TransitionMatrix object in an ASCII file,
 * putting down each transition probability by component
 *
 * @param  tmat:        TransitionMatrix to save
 * @param  fname:       Filename
 * @retval None
 */
void writeTransitionMatrix(TransitionMatrix tmat, string fname) {
  ofstream out(fname);

  // Header
  out << "#####################################################\n";
  out << "# TransitionMatrix from state with k = " << tmat.k1 << " to state with k = " << tmat.k2 << "\n";
  out << "# Total rate = " << tmat.totalRate() * Physical::s << " s^-1\n";
  out << "#####################################################\n";

  for (int i = -1; i < (int)tmat.m1.size(); ++i) {
    for (int j = -1; j < (int)tmat.m2.size(); ++j) {
      if (j == -1 && i == -1) {
        out << setw(5) << ' ';
      } else if (j == -1) {
        out << setw(5) << tmat.m1[i];
      } else if (i == -1) {
        out << setw(15) << tmat.m2[j];
      } else {
        out << setw(15) << tmat.T[i][j];
      }
    }
    out << '\n';
  }

  out.close();
}

/**
 * @brief  Write an EConfPotential object to a text file
 * @note   Write down a full Electronic Configuration
 * Potential in ASCII format, as grid, charge density,
 * and potential.
 *
 * @param  epot:        EConfPotential to print
 * @param  fname:       Filename
 * @retval None
 */
void writeEConfPotential(EConfPotential epot, string fname) {
  ofstream out(fname);

  int i0 = epot.getGridLimits().first;
  vector<double> r = epot.getGrid()[1];

  out << "# Q = " << epot.getQ() << "\n";

  for (int i = 0; i < r.size(); ++i) {
    out << r[i] << '\t' << epot.getrho()[i] << '\t' << epot.Vgrid(i0 + i) << '\n';
  }

  out.close();
}

/**
 * @brief  Write a simulated spectrum as an ASCII file
 * @note   Write a simulated spectrum from the given transition
 * energies and probabilities. The spectrum will be composed of a series of Gaussian
 * lines centred at the transition energies, with height proportional to the probabilities.
 * An exponential decay can be added to account for the potential loss of sensitivity
 * at high energy of the detector.
 *
 * @param  transitions:         A vector of TransitionData structs, containing all the necessary information
 * @param  dE:                  Step of the points in the spectrum
 * @param  lw:                  Line width for the Gaussians
 * @param  expd:                Exponential decay factor
 * @param  fname:               Name of the file to save
 * @retval None
 */
void writeSimSpec(vector<TransitionData> transitions, double dE, double lw, double expd, string fname) {
  ofstream out(fname);

  int N = transitions.size();
  double minE, maxE; // Looking for the boundaries
  vector<double> energies, intensities;

  minE = INFINITY;
  maxE = 0;

  for (int i = 0; i < N; ++i) {
    double E = (transitions[i].ds2.E - transitions[i].ds1.E) / Physical::eV;
    minE = min(E, minE);
    maxE = max(E, maxE);

    energies.push_back(E);
    intensities.push_back(transitions[i].tmat.totalRate() * transitions[i].tmat.m1.size());
  }

  // Now adjust the boundaries
  minE = max(0.0, minE - 5 * lw);
  maxE += 5 * lw;

  for (double E = minE; E <= maxE; E += dE) {
    out << E << '\t';
    double S = 0.0;
    for (int i = 0; i < N; ++i) {
      S += exp(-0.5 * pow((energies[i] - E) / lw, 2)) * intensities[i];
    }
    if (expd > 0) {
      S *= exp(-E / expd);
    }
    out << S << '\n';
  }

  out.close();
}

// Debug tasks

void writeEdEscan(vector<double> Es, vector<double> dEs, vector<int> nodes, string fname) {
  ofstream out(fname);
  int N = Es.size();

  for (int i = 0; i < N; ++i) {
    out << Es[i]*Physical::eV << "\t" << dEs[i]*Physical::eV << "\t" << nodes[i] << "\n";
  }

  out.close();
}