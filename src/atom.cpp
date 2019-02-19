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
#include <iostream>

State::State()
{
    E = 0;
    nodes = 0;
}

SchroState::SchroState(const SchroState &s)
{
    nodes = s.nodes;
    E = s.E;
    l = s.l;
    R = vector<double>(s.R);
}

DiracState::DiracState()
{
}

DiracState::DiracState(const DiracState &s)
{
    nodes = s.nodes;
    E = s.E;
    k = s.k;
    Q = vector<double>(s.Q);
    P = vector<double>(s.P);
}

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
    double V0 = R <= 0 ? 0 : -1.5 * Z / R;
    double R3 = pow(R, 3.0);

    V = vector<double>(N, 0);

    // Start with the background charge contribution
    shootPotentialLog(V, bkgQ, dx);

    for (int i = 0; i < N; ++i)
    {
        if (r[i] <= R)
        {
            V[i] += Z * pow(r[i], 2) / (2 * R3) + V0;
        }
        else
        {
            V[i] += -Z / r[i];
        }
    }

    // Now send it to zero at the end
    for (int i = 0; i < N; ++i)
    {
        V[i] -= V[N - 1];
    }
}

void Atom::setGrid(double r0_in, double r1_in, int N_in)
{
    r0 = r0_in;
    r1 = r1_in;
    N = N_in;

    grid = logGrid(r0, r1, N);
    dx = grid[0][1] - grid[0][0];
    bkgQ = vector<double>(N, 0); // Reset background charge
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

DiracAtom::DiracAtom(double Z_in, double m_in, double A_in, double R_in) : Atom(Z_in, m_in, A_in, R_in)
{
}

void DiracAtom::calcState(int n, int l, bool s, bool force)
{
    int k = s ? l : -l - 1;
    int maxit = 100;
    double E, err;
    DiracState *state = new DiracState();
    TurningPoint tp;

    // First, check if it's already calculated
    if (!force && !(states[{n, l, s}] == NULL))
    {
        return;
    }

    // Then start with a guess for the energy
    E = hydrogenicDiracEnergy(Z, mu, n, k);

    for (int it = 0; it < maxit; ++it)
    {
        // Start by applying boundary conditions
        boundaryDiracCoulomb(state->Q, state->P, grid[1], E, k, mu, Z);
        tp = shootDiracLog(state->Q, state->P, grid[1], V, E, k, mu, dx);
    }

    cout << E << '\n';
}