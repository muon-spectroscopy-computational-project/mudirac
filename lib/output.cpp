/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * output.cpp
 * 
 * Functions to write output files
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
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
void writeDiracState(DiracState ds, string fname)
{
    ofstream out(fname);

    // Start with writing a header
    out << "#####################################################\n";
    out << "# DiracState with n = " << ds.getn() << ", l = " << ds.getl() << ", s = " << ds.gets() << "\n";
    out << "# E = " << ds.E / Physical::eV << " eV\n";
    out << "# nodes = " << ds.nodes << ", " << ds.nodesQ << "\n";
    out << "#####################################################\n";

    // Now write down the content
    for (int i = 0; i < ds.grid.size(); ++i)
    {
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
void writeTransitionMatrix(TransitionMatrix tmat, string fname)
{
    ofstream out(fname);

    // Header
    out << "#####################################################\n";
    out << "# TransitionMatrix from state with k = " << tmat.k1 << " to state with k = " << tmat.k2 << "\n";
    out << "# Total rate = " << tmat.totalRate() * Physical::s << " s^-1\n";
    out << "#####################################################\n";

    for (int i = -1; i < (int)tmat.m1.size(); ++i)
    {
        for (int j = -1; j < (int)tmat.m2.size(); ++j)
        {
            if (j == -1 && i == -1)
            {
                out << setw(5) << ' ';
            }
            else if (j == -1)
            {
                out << setw(5) << tmat.m1[i];
            }
            else if (i == -1)
            {
                out << setw(15) << tmat.m2[j];
            }
            else
            {
                out << setw(15) << tmat.T[i][j];
            }
        }
        out << '\n';
    }
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
void writeEConfPotential(EConfPotential epot, string fname)
{
    ofstream out(fname);

    int i0 = epot.getGridLimits().first;
    vector<double> r = epot.getGrid()[1];

    out << "# Q = " << epot.getQ() << "\n";

    for (int i = 0; i < r.size(); ++i)
    {
        out << r[i] << '\t' << epot.getrho()[i] << '\t' << epot.Vgrid(i0 + i) << '\n';
    }
}