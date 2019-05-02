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

AtomConvergenceException::AtomConvergenceException(int t)
{

    type = t;

    switch (type)
    {
    case NAN_ENERGY:
        msg = "Atomic energy could not converge";
        break;
    case NODES_WRONG:
        msg = "Nodal theorems conditions not satisfied";
        break;
    case MAXIT_REACHED:
        msg = "Maximum numbers of iterations reached";
        break;
    default:
        msg = "Unknown atomic convergence error";
        break;
    }
}

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

/**
 * @brief  Initialise a DiracState instance
 * @note   Creates a DiracState with given grid size
 * 
 * @param  N: Desired grid size for Q and P
 * @retval 
 */
DiracState::DiracState(int N)
{
    grid = vector<double>(N, 0);
    loggrid = vector<double>(N, 0);
    Q = vector<double>(N, 0);
    P = vector<double>(N, 0);
}

/**
 * @brief  Initialise a DiracState instance with a grid
 * @note   Creates a DiracState with given grid size and initialises the
 * spatial grid itself
 * 
 * @param  x0:  Grid starting point
 * @param  x1:  Grid ending point
 * @param  N:   Number of steps
 * @retval 
 */
DiracState::DiracState(double x0, double x1, int N) : DiracState(N)
{
    vector<vector<double>> grids = logGrid(x0, x1, N);
    loggrid = grids[0];
    grid = grids[1];
}

DiracState::DiracState(const DiracState &s)
{
    init = true;
    nodes = s.nodes;
    E = s.E;
    k = s.k;
    grid = vector<double>(s.grid);
    loggrid = vector<double>(s.loggrid);
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
    setGridRelative();
}

/**
 * @brief  Recalculate the electrostatic potential
 * @note   Recalculate the electrostatic potential for an atom. Done automatically 
 * after changes in grid or background charge
 * 
 * @retval None
 */
void Atom::recalcPotential()
{
    vector<double> r = grid[1];
    double V0 = R <= r0 ? 0 : -1.5 * Z / R;
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
}

/**
 * @brief  Set logarithmic integration grid
 * @note   Set parameters of desired logarithmic integration grid for this Atom
 * 
 * @param  r0_in: Inferior boundary of the grid
 * @param  r1_in: Superior boundary of the grid
 * @param  N_in:  Number of points
 * @retval None
 */
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

/**
 * @brief Set logarithmic integration grid relative to a physical scale
 * @note  Set parameters of desired logarithmic integration grid for this Atom,
 * using a meaningful physical scale which corresponds to the expected radius of
 * the Schroedinger 1s orbital. The parameters are turned to distances following:
 * 
 * r = f/(Z*mu)
 * 
 * with mu effective mass of the atom. This scales down the dimension to what is
 * reasonable for the atom in question, allowing one to use homogeneous parameters
 * for multiple atoms.
 * 
 * @param   f0:       Inferior boundary of the grid
 * @param   f1:       Superior boundary of the grid
 * @param   N:        Number of points
 * @retval  None
 */
void Atom::setGridRelative(double f0, double f1, int N)
{
    setGrid(f0 / (Z * mu), f1 / (Z * mu), N);
}

/**
 * @brief  Set the background charge
 * @note   Set the background charge for this Atom. The charge must be already 
 * expressed as amount of charge per spherical shell (so for example a constant 
 * spatial charge q should be here 4pi*r^2*q)
 * 
 * @param  bkgQ_in: Vector of new background charge
 * @retval None
 */
void Atom::setBackgroundCharge(vector<double> bkgQ_in)
{
    bkgQ = vector<double>(bkgQ_in);
    recalcPotential();
}

/**
 * @brief  Get the Atom object's grid
 * @note   Get either r or log(r/r0) for this Atom's grid
 * 
 * @param  log: If true, return log(r/r0) instead of r (default = false)
 * @retval 
 */
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

/**
 * @brief  Converge a state with given k and initial energy guess E0
 * @note   Converge iteratively a Dirac orbital for this atom from a given k and
 * energy starting guess. Will fail if convergence can't be achieved or if nodal theorems
 * are violated. If successful, will return the state found
 * 
 * @param  E0: Initial energy guess
 * @param  k:  Quantum number k
 * @retval     Found DiracState
 */
DiracState DiracAtom::convergeState(double E0, int k)
{
    int Qn, Pn;
    double E, dE, err, norm;
    vector<double> y(N), zetai(N), zetae(N);
    DiracState state = DiracState(N);
    TurningPoint tp;

    E = E0;

    for (int it = 0; it < maxit; ++it)
    {
        // Start by applying boundary conditions
        boundaryDiracCoulomb(state.Q, state.P, grid[1], E, k, mu, Z, R > r0);
        // cout << state.P[0] << '\t' << state.P[N - 1] << '\n';
        // cout << state.Q[0] << '\t' << state.Q[N - 1] << '\n';
        // Integrate here
        tp = shootDiracLog(state.Q, state.P, grid[1], V, E, k, mu, dx);
        err = tp.Qi / tp.Pi - tp.Qe / tp.Pe;
        // Compute the derivative of the error in dE
        for (int i = 0; i < N; ++i)
        {
            y[i] = state.Q[i] / state.P[i];
        }

        // First the forward version
        y[tp.i] = tp.Qi / tp.Pi;
        shootDiracErrorDELog(zetai, y, grid[1], V, tp.i, E, k, mu, dx);
        // Then the backwards one
        y[tp.i] = tp.Qe / tp.Pe;
        boundaryDiracErrorDECoulomb(zetae, E, k, mu);
        shootDiracErrorDELog(zetae, y, grid[1], V, tp.i, E, k, mu, dx, 'b');

        dE = err / (zetai[tp.i] - zetae[tp.i]);
        // cout << tp.Qe << "\n";
        // cout << tp.Pe << "\n";
        // cout << zetai[tp.i] << " " << zetae[tp.i] << "\n";
        // cout << (it + 1) << '\t' << E - mu * pow(Physical::c, 2) << '\t' << dE << '\n';
        if (!std::isnan(dE) && abs(dE) < Etol)
        {
            break;
        }
        E = E - dE;
        if (std::isnan(E))
        {
            // Something bad happened
            throw AtomConvergenceException(AtomConvergenceException::NAN_ENERGY);
        }
    }

    // Make states continuous
    for (int i = tp.i; i < N; ++i)
    {
        state.P[i] *= tp.Pi / tp.Pe;
        state.Q[i] *= tp.Pi / tp.Pe;
    }
    // Now normalise
    for (int i = 0; i < N; ++i)
    {
        y[i] = (pow(state.P[i], 2) + pow(state.Q[i] * Physical::alpha, 2)) * grid[1][i];
    }
    norm = sqrt(trapzInt(grid[0], y));
    for (int i = 0; i < N; ++i)
    {
        state.P[i] /= norm;
        state.Q[i] /= norm;
    }

    // Count nodes
    Pn = countNodes(state.P);
    Qn = countNodes(state.Q);

    if (Qn - Pn != (R > r0))
    {
        throw AtomConvergenceException(AtomConvergenceException::NODES_WRONG);
    }

    state.nodes = Pn;
    state.nodesQ = Qn;
    state.E = E;
    state.init = true;
    state.k = k;

    return state;
}

/**
 * @brief  Search for an orbital with given set of quantum numbers
 * @note   Search for a Dirac orbital for this Atom with a given set of
 * quantum numbers. Will start with a guess equal to the one for the hydrogenic solution
 * and search from there. If any other states are found by accident during the search they
 * are stored for future use.
 * 
 * @param  n: Principal quantum number
 * @param  l: Orbital quantum number
 * @param  s: Spin quantum number (true = 1/2, false = -1/2)
 * @param  force: If true, force recalculation of the orbital even if already present
 * @retval None
 */
void DiracAtom::calcState(int n, int l, bool s, bool force)
{
    int k = ((s && l > 0) ? l : -l - 1);
    int dnode;
    double E0;
    bool found = false;
    DiracState state;
    TurningPoint tp;

    // First, check if it's already calculated
    if (!force && states[make_tuple(n, l, s)].init)
    {
        return;
    }

    // Then start with a guess for the energy
    E0 = hydrogenicDiracEnergy(Z, mu, n, k);

    for (int it = 0; it < maxit; ++it)
    {
        state = convergeState(E0, k);
        // Is the nodes condition respected?
        dnode = state.nodes - (n - l - 1);
        if (dnode == 0)
        {
            break;
        }
        else
        {
            // Still save the state for future use
            states[make_tuple(state.nodes + l + 1, l, s)] = state;
            E0 = dnode > 0 ? E0 / Esearch : E0 * Esearch;
        }
    }

    states[make_tuple(n, l, s)] = state;
}

/**
 * @brief  Return an orbital with given set of quantum numbers
 * @note   Search for a Dirac orbital for this Atom with a given set of
 * quantum numbers. If the state has already been calculated and stored, return it.
 * Otherwise, calculate it, then return it.  
 * 
 * @param  n: Principal quantum number
 * @param  l: Orbital quantum number
 * @param  s: Spin quantum number (true = 1/2 / false = -1/2)
 * @retval Requested orbital
 */
DiracState DiracAtom::getState(int n, int l, bool s)
{
    calcState(n, l, s);
    DiracState st = states[make_tuple(n, l, s)];

    if (!st.init)
    {
        throw AtomConvergenceException(AtomConvergenceException::MAXIT_REACHED);
    }

    return DiracState(st);
}