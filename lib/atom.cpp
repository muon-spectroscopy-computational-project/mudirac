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

AtomConvergenceException::AtomConvergenceException(ACEType t)
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
    V = vector<double>(N, 0);
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
DiracState::DiracState(double rc, double dx, int i0, int i1) : DiracState(i1 - i0 + 1)
{
    vector<vector<double>> grids = logGrid(rc, dx, i0, i1);
    loggrid = grids[0];
    grid = grids[1];

    grid_indices.first = i0;
    grid_indices.second = i1;
}

DiracState::DiracState(const DiracState &s)
{
    init = true;
    nodes = s.nodes;
    E = s.E;
    k = s.k;
    grid_indices = pair<int, int>(s.grid_indices);
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
 * @param  fc:   Central point of the grid (corresponding to i = 0), as a fraction
 *               of 1/(Z*mu), the 1s orbital radius for this atom (default = 1)
 * @param  dx:   Logarithmic step of the grid (default = 0.005)
 * @retval 
 */
Atom::Atom(double Z, double m, double A, NuclearRadiusModel radius_model,
           double fc, double dx)
{
    // Set the properties
    this->Z = Z;
    this->A = A;
    this->m = m;

    // Sanity checks
    if (Z <= 0)
    {
        throw invalid_argument("Z must be positive");
    }
    if (m <= 0)
    {
        throw invalid_argument("Mass can not be negative");
    }
    if (fc <= 0 || dx <= 0)
    {
        throw invalid_argument("Invalid grid parameters passed to Atom");
    }

    if (A > 0)
    {
        mu = effectiveMass(m, A * Physical::amu);
    }
    else
    {
        mu = m;
    }

    // Define radius
    if (A == -1)
    {
        R = -1;
    }
    else
    {
        switch (radius_model)
        {
        case POINT:
            R = -1;
            break;
        case SPHERE:
            R = sphereNuclearModel(A);
            break;
        default:
            R = -1;
            break;
        }
    }

    // Grid
    rc = fc / (Z * mu);
    this->dx = dx;

    // Potential
    V = CoulombPotential(Z, R);
}

/**
 * @brief  Recalculate the electrostatic potential
 * @note   Recalculate the electrostatic potential for an atom. Done automatically
 * after changes in grid or background charge
 *
 * @param r:        Grid to compute the potential on
 * @retval          Computed potential
 */
vector<double> Atom::recalcPotential(vector<double> r)
{
    int N = r.size();
    vector<double> Vout(N, 0);

    for (int i = 0; i < N; ++i)
    {
        Vout[i] = V.V(r[i]);
    }

    return Vout;
}

// Nuclear radius models

/**
 * @brief  Spherical nuclear radius model
 * @note   Describe the nucleus' finite size as if it was a sphere,
 * assuming a model by which R = 1.2 A^(1/3) fm.
 * 
 * @param  A:   Atomic mass
 * @retval      Nuclear radius
 */
double Atom::sphereNuclearModel(double A)
{
    return 1.2 * Physical::fm * pow(A, 1.0 / 3.0);
}

DiracAtom::DiracAtom(double Z, double m, double A, NuclearRadiusModel radius_model, double fc, double dx) : Atom(Z, m, A, radius_model, fc, dx)
{
}

pair<int, int> DiracAtom::gridLimits(double E, int k, GridLimitsFailcode &failcode)
{
    double B;
    double K = pow(mu * Physical::c, 2) - pow(E / Physical::c, 2);
    double gamma = pow(k, 2) - pow(Z * Physical::alpha, 2);
    double r_out, r_in;
    int i_out, i_in;

    failcode = GridLimitsFailcode::OK;

    if (K < 0)
    {
        failcode = GridLimitsFailcode::UNBOUND;
        return {0, 0};
    }
    if (gamma < 0)
    {
        // Unlikely while we're in the periodic table...
        failcode = GridLimitsFailcode::SMALL_GAMMA;
        return {0, 0};
    }

    K = sqrt(K);
    gamma = sqrt(gamma);
    B = E - mu * pow(Physical::c, 2);

    // Upper limit
    if (out_eps > 1 || out_eps < 0)
    {
        throw runtime_error("Invalid value for out_eps in DiracAtom; must be 0 < out_eps < 1");
    }
    r_out = Z / abs(B) - log(out_eps) / K;

    // Lower limit
    if (in_eps > 1 || in_eps < 0)
    {
        throw runtime_error("Invalid value for in_eps in DiracAtom; must be 0 < in_eps < 1");
    }
    r_in = pow(in_eps, 1.0 / gamma) / M_E * gamma / K;

    // Now get these as integer numbers of steps.
    i_out = ceil(log(r_out / rc) / dx);
    i_in = floor(log(r_in / rc) / dx);

    return {i_in, i_out};
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
    int N, Qn, Pn;
    double E, dE, err, norm;
    pair<int, int> glimits;
    vector<double> y, zetai, zetae;
    DiracState state;
    TurningPoint tp;
    GridLimitsFailcode fcode;

    E = E0;

    for (int it = 0; it < maxit; ++it)
    {
        glimits = gridLimits(E, k, fcode);
        switch (fcode)
        {
        case GridLimitsFailcode::OK:
            break;
        case GridLimitsFailcode::UNBOUND:
            // LOG(INFO) << "Convergence failed - State with B = " << E-mu*pow(Physical::c, 2) << " is unbound";            
            continue;
            break;
        case GridLimitsFailcode::SMALL_GAMMA:
            break;
        default:
            throw runtime_error("Unknown error occurred in grid estimation in convergeState");
            break;
        }
        state = DiracState(rc, dx, glimits.first, glimits.second);
        N = state.grid.size();
        y = vector<double>(N, 0);
        zetai = vector<double>(N, 0);
        zetae = vector<double>(N, 0);
        // Start by applying boundary conditions
        boundaryDiracCoulomb(state.Q, state.P, state.grid, E, k, mu, Z, R > state.grid[0]);
        // Potential
        state.V = recalcPotential(state.grid);
        // Integrate here
        tp = shootDiracLog(state.Q, state.P, state.grid, state.V, E, k, mu, dx);
        err = tp.Qi / tp.Pi - tp.Qe / tp.Pe;
        // Compute the derivative of the error in dE
        for (int i = 0; i < N; ++i)
        {
            y[i] = state.Q[i] / state.P[i];
        }

        // First the forward version
        y[tp.i] = tp.Qi / tp.Pi;
        shootDiracErrorDELog(zetai, y, state.grid, state.V, tp.i, E, k, mu, dx);
        // Then the backwards one
        y[tp.i] = tp.Qe / tp.Pe;
        boundaryDiracErrorDECoulomb(zetae, E, k, mu);
        shootDiracErrorDELog(zetae, y, state.grid, state.V, tp.i, E, k, mu, dx, 'b');

        dE = err / (zetai[tp.i] - zetae[tp.i]);
        if (!std::isnan(dE) && abs(dE) < Etol)
        {
            break;
        }
        E = E - dE;
        if (std::isnan(E))
        {
            // Something bad happened
            throw AtomConvergenceException(AtomConvergenceException::ACEType::NAN_ENERGY);
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
        y[i] = (pow(state.P[i], 2) + pow(state.Q[i] * Physical::alpha, 2)) * state.grid[i];
    }
    norm = sqrt(trapzInt(state.loggrid, y));
    for (int i = 0; i < N; ++i)
    {
        state.P[i] /= norm;
        state.Q[i] /= norm;
    }

    // Count nodes
    Pn = countNodes(state.P);
    Qn = countNodes(state.Q);

    if (Qn - Pn != (R > state.grid[0]))
    {
        throw AtomConvergenceException(AtomConvergenceException::ACEType::NODES_WRONG);
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
 * @param  n:       Principal quantum number
 * @param  l:       Orbital quantum number
 * @param  s:       Spin quantum number (true = 1/2, false = -1/2)
 * @param  force:   If true, force recalculation of the orbital even if already present
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
 * @brief  Calculate all states up to a given n
 * @note   Calculate all states up to a given quantum number n,
 * including all orbital and spin quantum numbers.
 * 
 * @param  max_n:   Maximum value of principal quantum number
 * @param  force:   If true, force recalculation of the orbital even if already present
 * @retval None
 */
void DiracAtom::calcAllStates(int max_n, bool force)
{
    for (int n = 1; n <= max_n; ++n)
    {
        for (int l = 0; l < n; ++l)
        {
            for (int s = 0; s < 2; ++s)
            {
                calcState(n, l, bool(s), force);
            }
        }
    }
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
        throw AtomConvergenceException(AtomConvergenceException::ACEType::MAXIT_REACHED);
    }

    return DiracState(st);
}