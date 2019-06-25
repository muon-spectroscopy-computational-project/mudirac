/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 *
 * wavefunction.cpp
 *
 * Functions useful for operations on wavefunctions
 *
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "wavefunction.hpp"

double braKetLog(vector<double> psiL, vector<double> psiR, vector<double> r, double dx)
{
    double tot;
    int N = psiL.size();
    vector<double> ker(N);

    if (psiR.size() != N)
    {
        throw invalid_argument("Bra and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = psiL[i] * psiR[i] * r[i];
    }

    return trapzInt(dx, ker);
}

double braKetLog(vector<vector<double>> psiL, vector<vector<double>> psiR, vector<double> r, double dx)
{
    double tot;
    int N = psiL[0].size();
    vector<double> ker(N);

    if (psiL.size() != 2 || psiR.size() != 2 || psiR[0].size() != N)
    {
        throw invalid_argument("Bra and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = (psiL[0][i] * psiR[0][i] + psiL[1][i] * psiR[1][i]) * r[i];
    }

    return trapzInt(dx, ker);
}