/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * debugtasks.cpp
 * 
 * Methods or classes used to perform debug specific tasks
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "debugtasks.hpp"

void runEdEscan(MuDiracInputFile infile)
{
    DiracAtom da = infile.makeAtom();

    int k = infile.getIntValue("devel_EdEscan_k");
    double minE = infile.getDoubleValue("devel_EdEscan_minE");
    double maxE = infile.getDoubleValue("devel_EdEscan_maxE");

    pair<double, double> limE = da.energyLimits(0, k);

    limE.first = max(limE.first-da.getRestE(), minE)/Physical::eV;
    limE.second = min(limE.second-da.getRestE(), maxE)/Physical::eV;

    // What's the range?
    int nE = infile.getIntValue("devel_EdEscan_steps");
    bool logE = infile.getBoolValue("devel_EdEscan_log");
    double stepE = logE ? pow(limE.second - limE.first, 1.0 / (nE - 1)) : (limE.second - limE.first) / (nE - 1);

    vector<double> Erange(nE);

    LOG(TRACE) << "EdE scan limits: " << limE.first << " <= E <= " << limE.second << "\n";

    for (int i = 0; i < nE; ++i)
    {
        Erange[i] = logE ? limE.second - (pow(stepE, nE - i - 1)) : limE.first + stepE * i;
    }

    DiracState ds;
    TurningPoint tp;
    double dE;

    ds.k = k;

    vector<double> dEs(nE);
    vector<int> nodes(nE);
    
    for (int i = 0; i < nE; ++i) {
        ds = da.initState(Erange[i]+da.getRestE(), k);
        da.integrateState(ds, tp, dE);
        
        ds.continuify(tp);
        ds.normalize();
        ds.findNodes();

        dEs[i] = dE;
        nodes[i] = ds.nodes;
    }

    string fname = "EdEscan_" + to_string(k) + ".dat";

    writeEdEscan(Erange, dEs, nodes, fname);
}