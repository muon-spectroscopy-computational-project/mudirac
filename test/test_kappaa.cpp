#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "datapath.h"
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main(int argc, char *argv[])
{

    // Open the existing data and use it to compare with computed quantities
    vector<int> dataZ;
    vector<double> dataA, dataKa, calcKa, calcUehlKa;
    string line;
    ifstream datafile(string(CURRENT_DATAPATH) + "/data/kappa_alpha_exp.dat");
    float Z, A, Ka;

    AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace, AixLog::Type::normal);

    while (getline(datafile, line))
    {
        istringstream iss(line);
        if (line[0] == '#')
            continue;
        iss >> Z >> A >> Ka;
        dataZ.push_back(Z);
        dataA.push_back(A);
        dataKa.push_back(Ka);
    }

    // Now, simulations
    for (int i = 0; i < dataZ.size(); ++i)
    {
        DiracAtom da = DiracAtom(dataZ[i], Physical::m_mu, dataA[i], SPHERE);
        DiracState p2 = da.getState(2, 1, false);
        DiracState s1 = da.getState(1, 0, false);

        Ka = (p2.E - s1.E) / (1e3*Physical::eV);        

        calcKa.push_back(Ka);

        // Now with Uehling potential
        da.setUehling(true);
        p2 = da.getState(2, 1, false);
        s1 = da.getState(1, 0, false);

        Ka = (p2.E - s1.E) / (1e3*Physical::eV);        

        calcUehlKa.push_back(Ka);
    }

    // Print out
    if (argc > 1)
    {
        ofstream outfile(argv[1]);
        for (int i = 0; i < dataZ.size(); ++i)
        {
            outfile << dataZ[i] << '\t' << dataA[i] << '\t' << dataKa[i] << '\t' << calcKa[i] << '\t' << (calcKa[i] - dataKa[i]) / dataKa[i] << '\t';
            outfile << calcUehlKa[i] << '\t' << (calcUehlKa[i] - dataKa[i]) / dataKa[i] << '\n';
        }
    }
}