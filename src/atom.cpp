/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * atom.cpp
 * 
 * Atom classes, designed to be a high level interface to solution of quantum mechanical equations
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "atom.hpp"

/**
 * @brief  Initialise an Atom class instance
 * @note   Creates an Atom object defined by the given properties
 * 
 * @param  Z_in: Atomic number (nuclear charge, can be fractional)
 * @param  m_in: Mass of the orbiting particle (e.g. electron)
 * @param  A_in: Atomic mass (amus, ignored if -1)
 * @param  R_in: Atomic radius (treated as point-like if <= 0)
 * @retval 
 */
Atom::Atom(double Z_in, double m_in, double A_in, double R_in)
{
    // Set the properties
    Z = Z_in;
    A = A_in;
    m = m_in;
    R = R_in;

    if (A > 0)
    {
        mu = effectiveMass(m, A * Physical::amu);
    }
    else
    {
        mu = m;
    }

    // Compute the grid
    setGrid();
}

void Atom::recalcPotential()
{
    vector<double> r = grid[1];
    vector<double> Vbkg;
    double V0 = R <= 0 ? 0 : -1.5 * Z / R;
    double R3 = pow(R, 3.0);

    V = vector<double>(N, 0);

    // Additional potential contribution from the background charge
    shootNumerov(Vbkg, vector<double>(N, 0), bkgQ, dx);

    for (int i = 0; i < N; ++i)
    {
        if (r[i] <= R)
        {
            V[i] = Z * pow(r[i], 2) / (2 * R3) + V0;
        }
        else
        {
            V[i] = -Z / r[i];
        }
    }
}

void Atom::setGrid(double r0_in, double r1_in, int N_in)
{
    r0 = r0_in;
    r1 = r1_in;
    N = N_in;

    grid = logGrid(r0, r1, N);
    dx = grid[0][1] - grid[0][0];
    recalcPotential();
}

void Atom::setBackgroundCharge(vector<double> bkgQ_in)
{
    bkgQ = vector<double>(bkgQ_in);
    recalcPotential();
}

vector<double> Atom::getGrid(bool log)
{
    return vector<double>(grid[!log]);
}

vector<double> Atom::getPotential()
{
    return vector<double>(V);
}