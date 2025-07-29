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
  * @brief  Write valid 2 parameter Fermi Distribution parameters to a file
  * @note   Write a tsv of valid fermi parameters which fit the input experimental results in an ASCII file.
  * The file has columns of fermi_c, fermi_t, rms_radius, theta, mean_sq_error, where fermi_c and fermi_t are
  * the parameters for the Fermi distribution, rms_radius is the root mean square radius of that distribution,
  * theta is a polar coordinate which maps the rms_radius to the fermi parameters, and mean_sq_error is the
  * mean square error between the experimental xray transition energies and the transition energies simulated
  * by mudirac using the fermi parameters in the row.
  *
  * @param  da:         Dirac Atom containing required configuration data (Z, A, m)
  * @param  fermi_parameters:         vector containing valid fermi parameter data to be output
  * @param  fname:      Filename
  * @param  output_precision: output decimal places
  *
  * @retval None
 */
void writeFermiParameters(DiracAtom &da, const OptimisationData &fermi_parameters, const double opt_time, string fname, int output_precision) {

  // output file containing all valid fermi parameters and the associated MSE
  ofstream out(fname);
  out << "# Z = " << da.getZ() << ", A = " << da.getA() << " amu, m = " << da.getm() << " au\n";
  out << "fermi_c\tfermi_t\trms_radius\ttheta\tmean_sq_error\titerations\ttime\n";
  out << fixed;
  out << setprecision(output_precision);

  // output fermi_c, fermi_c, rms radius, theta, MSE, n iterations, optimisation time
  out << fermi_parameters.fermi_c << '\t' << fermi_parameters.fermi_t << '\t';
  out << fermi_parameters.rms_radius << '\t' << fermi_parameters.theta  << '\t';
  out << fermi_parameters.mse <<'\t' << da.iteration_counter_2pF << '\t' << opt_time << '\n';


  out.close();
}

/**
  * @brief  Write a DiracState object to a text file
  * @note   Write down a full DiracState object in an ASCII file,
  * saving potential and P and Q component values
  *
  * @param  ds:         DiracState to save
  * @param  fname:      Filename
  * @retval None
 */
void writeDiracState(DiracState ds, string fname, int output_precision) {
  ofstream out(fname);

  // Start with writing a header
  out << "#####################################################\n";
  out << "# DiracState with n = " << ds.getn() << ", l = " << ds.getl() << ", s = " << ds.gets() << "\n";
  out << "# E = " << ds.bindingE() / Physical::eV << " + mc^2 = " << ds.E / Physical::eV << " eV\n";
  out << "# nodes = " << ds.nodes << ", " << ds.nodesQ << "\n";
  out << "#####################################################\n";

  if (output_precision > -1) {
    out << fixed;
    out << setprecision(output_precision);
  }

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

// main output functions 

void printInitLogMessage(){
  LOG(INFO) << "MuDirac, a muonic atomic solver\n";
  LOG(INFO) << "by Simone Sturniolo\n";
  LOG(INFO) << "Released under the MIT License (2019)\n";
  LOG(INFO) << " \n";
  LOG(INFO) << "Please cite the following paper:\n";
  LOG(INFO) << "Sturniolo, S, Hillier, A., \"Mudirac: A Dirac equation solver for elemental analysis with muonic X‐rays.\"\n";
  LOG(INFO) << "X‐Ray Spectrom. 2020; 1– 17. https://doi.org/10.1002/xrs.3212\n";
  LOG(INFO) << " \n";
}

void writeOutputFiles(string seed, MuDiracInputFile & config, DiracAtom & da, vector<TransitionData> & transitions){
  
  // Sort transitions by energy if requested
  if (config.getBoolValue("sort_by_energy")) {
    LOG(DEBUG) << "sorting transitions by energy \n";
    sort(transitions.begin(), transitions.end(),
    [](TransitionData t1, TransitionData t2) {
      return (t1.ds2.E - t1.ds1.E) > (t2.ds2.E - t2.ds1.E);
    });
  }

  // Now create output files
  int output_verbosity = config.getIntValue("output");
  if (output_verbosity >= 1) {
    // Save a file for all lines
    LOG(DEBUG) << "writing xray transitions output file \n";
    ofstream out(seed + ".xr.out");
    out << "# Z = " << da.getZ() << ", A = " << da.getA() << " amu, m = " << da.getm() << " au\n";
    out << "Line\tDeltaE (eV)\tW_12 (s^-1)\n";
    out << fixed;

    if (config.getIntValue("xr_print_precision") > -1) {
      out << setprecision(config.getIntValue("xr_print_precision"));
    } else {
      out << setprecision(15); //Setting the maximum precision
    }

    for (int i = 0; i < transitions.size(); ++i) {
      double dE = (transitions[i].ds2.E - transitions[i].ds1.E);
      double tRate = transitions[i].tmat.totalRate();
      if (dE <= 0 || tRate <= 0)
        continue; // Transition is invisible
      out << transitions[i].name << '\t' << dE / Physical::eV;
      out << "\t\t" << tRate * Physical::s << '\n';
    }

    // could be a config method in output.hpp?
    if (config.getBoolValue("write_spec")) {
      // Write a spectrum
      writeSimSpec(transitions, config.getDoubleValue("spec_step"), config.getDoubleValue("spec_linewidth"), config.getDoubleValue("spec_expdec"),
                   seed + ".spec.dat");
    }

    out.close();
  }

  if (output_verbosity >= 2) {
    vector<string> saved_states;
    // Save each individual state
    for (int i = 0; i < transitions.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
        DiracState ds = (j == 0 ? transitions[i].ds1 : transitions[i].ds2);
        string sname = ds.name();
        string fname = seed + "." + sname + ".out";

        if (vectorContains(saved_states, sname)) {
          continue;
        }

        LOG(DEBUG) << "Printing out state file for state " << sname << "\n";

        writeDiracState(ds, fname, config.getIntValue("state_print_precision"));

        saved_states.push_back(sname);
      }
      string fname = seed + "." + transitions[i].name + ".tmat.out";
      writeTransitionMatrix(transitions[i].tmat, fname);
    }
  }
}