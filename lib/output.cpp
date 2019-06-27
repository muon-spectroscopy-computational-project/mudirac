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
        out << ds.grid[i] << '\t' << ds.V[i]/Physical::eV << '\t' << ds.P[i] << '\t' << ds.Q[i] << '\n';
    }

    out.close();
}