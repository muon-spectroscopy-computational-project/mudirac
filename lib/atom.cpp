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
 * @brief   Compute the norm of this eigenstate
 * @note    Compute the norm of this eigenstate,
 * defined as:
 * 
 *              ( int_0^infty R^2 dr )^(1/2)
 * 
 * @retval      Norm 
 */
double SchroState::norm()
{
    int N = R.size();
    vector<double> rho(N);

    for (int i = 0; i < N; ++i)
    {
        rho[i] = pow(R[i], 2) * grid[i];
    }

    return sqrt(trapzInt(loggrid, rho));
}

/**
 * @brief  Get principal quantum number n
 * 
 * @retval      n
 */
int SchroState::getn()
{
    int n;
    qnumNodes2Principal(nodes, l, n);
    return n;
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
    converged = s.converged;
    nodes = s.nodes;
    nodesQ = s.nodesQ;
    E = s.E;
    k = s.k;
    grid_indices = pair<int, int>(s.grid_indices);
    grid = vector<double>(s.grid);
    loggrid = vector<double>(s.loggrid);
    Q = vector<double>(s.Q);
    P = vector<double>(s.P);
    V = vector<double>(s.V);
}

/**
 * @brief   Compute the norm of this eigenstate
 * @note    Compute the norm of this eigenstate,
 * defined as:
 * 
 *             ( int_0^infty (P^2+Q^2) dr )^(1/2)
 * 
 * @retval      Norm 
 */
double DiracState::norm()
{
    int N = P.size();
    vector<double> rho(N);

    for (int i = 0; i < N; ++i)
    {
        rho[i] = (pow(P[i], 2) + pow(Q[i], 2)) * grid[i];
    }

    return sqrt(trapzInt(loggrid, rho));
}

/**
 * @brief  Make the wavefunction continuous at the turning point
 * @note   Make the wavefunction continuous with the information
 * from a TurningPoint object from integration
 * 
 * @param  tp:              Turning point
 * @retval None
 */
void DiracState::continuify(TurningPoint tp)
{
    int N = P.size();
    double f = tp.Pi / tp.Pe;

    for (int i = tp.i; i < N; ++i)
    {
        P[i] *= f;
        Q[i] *= f;
    }
}

/**
 * @brief  Find and count the wavefunction's nodes
 * @note   Find and count the wavefunction's nodes
 * 
 * @param  tol:         Relative tolerance used to define when the function counts as non-zero for node counting
 * @retval None
 */
void DiracState::findNodes(double tol)
{
    nodes = countNodes(P, tol);
    nodesQ = countNodes(Q, tol);
}

/**
 * @brief  Normalize the wavefunction
 * @note   Normalize the wavefunction
 * 
 * @retval None
 */
void DiracState::normalize()
{
    int N = P.size();
    double n = norm();

    for (int i = 0; i < N; ++i)
    {
        P[i] /= n;
        Q[i] /= n;
    }
}

/**
 * @brief  Get principal quantum number n
 * 
 * @retval      n
 */
int DiracState::getn()
{
    int l = getl();
    int n;
    qnumNodes2Principal(nodes, l, n);
    return n;
}

/**
 * @brief  Get orbital quantum number l
 * 
 * @retval      l
 */
int DiracState::getl()
{
    int l;
    bool s;
    qnumDirac2Schro(k, l, s);
    return l;
}

/**
 * @brief  Get spin quantum number s (true for +1/2)
 * 
 * @retval      s
 */

bool DiracState::gets()
{
    int l;
    bool s;
    qnumDirac2Schro(k, l, s);
    return s;
}

/**
 * @brief  Initialise a TransitionMatrix class instance
 * @note   Initialise a TransitionMatrix class instance. 
 * Uses the given quantum numbers to compute the number of
 * m states to allocate.
 * 
 * @param  k1: Dirac quantum number for state 1
 * @param  k2: Dirac quantum number for state 2
 * @retval 
 */
TransitionMatrix::TransitionMatrix(int k1, int k2)
{
    // Remember:
    // j = abs(k)-1/2
    // 2j+1 = 2*abs(k)

    int mn1 = 2 * abs(k1); // Number of m states for j
    int mn2 = 2 * abs(k2);

    m1 = vector<double>(mn1);
    m2 = vector<double>(mn2);
    T = vector<vector<double>>(mn1, vector<double>(mn2, 0));

    for (int i = 0; i < mn1; ++i)
    {
        m1[i] = -(mn1 - 1.0) / 2.0 + i;
    }
    for (int i = 0; i < mn2; ++i)
    {
        m2[i] = -(mn2 - 1.0) / 2.0 + i;
    }
}

/**
 * @brief  Total transition rate between two shells
 * @note   Total transition rate between two shells, 
 * assuming the particle is in the first one and the
 * second is empty. In this case the rate is summed 
 * over the destination states, and averaged over
 * the origin states.
 * 
 * @retval The total transition rate
 */
double TransitionMatrix::totalRate()
{
    double tot = 0;

    for (int i = 0; i < m1.size(); ++i) {
        for (int j = 0; j < m2.size(); ++j) {
            tot += T[i][j];
        }
    }

    tot /= m1.size();
    return tot;
}

/**
 * @brief  Initialise an Atom class instance
 * @note   Creates an Atom object defined by the given properties
 * 
 * @param  Z_in: Atomic number (nuclear charge, can be fractional)
 * @param  m_in: Mass of the orbiting particle (e.g. electron)
 * @param  A_in: Atomic mass (amus, ignored if -1)
 * @param  R_in: Nuclear radius (treated as point-like if <= 0)
 * @param  fc:   Central point of the grid (corresponding to i = 0), as a fraction
 *               of 1/(Z*mu), the 1s orbital radius for this atom, or of
 *               the nuclear radius, depending on which one is bigger (default = 1)
 * @param  dx:   Logarithmic step of the grid (default = 0.005)
 * @retval 
 */
Atom::Atom(int Z, double m, int A, NuclearRadiusModel radius_model,
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
        M = getIsotopeMass(Z, A);
        mu = effectiveMass(m, M * Physical::amu);
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
            R = sphereNuclearModel(Z, A);
            break;
        default:
            R = -1;
            break;
        }
    }

    // Grid
    rc = fc * max(1 / (Z * mu), R);
    this->dx = dx;

    // Potential
    V = CoulombSpherePotential(Z, R);

    // Logging
    LOG(INFO) << "Created atom with Z = " << Z << ", A = " << A << "\n";
    LOG(INFO) << "Particle mass = " << m << " au, effective mass = " << mu << " au\n";
    LOG(INFO) << "Atomic radius = " << R << " au\n";
    LOG(INFO) << "Logarithmic grid parameters rc = " << rc << " au, dx = " << dx << "\n";
}

/**
  * @brief  Set parameters for the Uehling potential term
  * @note   Set up the Uehling potential term, activating/deactivating
  * it and setting the number of steps used for its integration. 
  * Calling this function resets all computed states.
  * 
  * @param  s:          Whether the Uehling potential should be on/off
  * @param  usteps:     Number of integration steps used for it (default = 1000)
  * @retval None
 */
void Atom::setUehling(bool s, int usteps, double cut_low, double cut_high)
{
    use_uehling = s;
    if (s)
    {
        V_uehling = UehlingSpherePotential(Z, R, usteps);
        V_uehling.set_exp_cutoffs(cut_low, cut_high);
    }
    reset();
}

/**
  * @brief  Set grid parameters
  * @note   Set parameters defining the logarithmic grid, rc and dx.
  * The grid will be made of points at
  * 
  * r = rc*exp(i*dx)
  * 
  * where i can be any integer number. Calling this function resets
  * all computed states.
  * 
  * @param  rc:         Central radius of the grid
  * @param  dx:         Logarithmic step of the grid
  * @retval None
 */
void Atom::setgrid(double rc, double dx)
{
    this->rc = rc;
    this->dx = dx;

    reset();
}

/**
 * @brief  Recalculate the electrostatic potential
 * @note   Recalculate the electrostatic potential for an atom. Done automatically
 * after changes in grid or background charge
 *
 * @param r:        Grid to compute the potential on
 * @retval          Computed potential
 */
vector<double> Atom::getV(vector<double> r)
{
    int N = r.size();
    vector<double> Vout(N, 0);

    for (int i = 0; i < N; ++i)
    {
        Vout[i] = V.V(r[i]);
        if (use_uehling)
        {
            Vout[i] += V_uehling.V(r[i]);
        }
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
double Atom::sphereNuclearModel(int Z, int A)
{
    try
    {
        return Physical::fm * getIsotopeRadius(Z, A);
    }
    catch (invalid_argument e)
    {
        LOG(TRACE) << "Isotope not found; falling back on default model for nuclear radius";
        return 1.2 * Physical::fm * pow(A, 1.0 / 3.0);
    }
}

DiracAtom::DiracAtom(int Z, double m, int A, NuclearRadiusModel radius_model, double fc, double dx) : Atom(Z, m, A, radius_model, fc, dx)
{
    restE = mu * pow(Physical::c, 2);
    LOG(DEBUG) << "Rest energy = " << restE / Physical::eV << " eV\n";
}

void DiracAtom::reset()
{
    states.clear();
}

/**
 * @brief  Bounds for the energy of a state with given k and nodes
 * @note   Lower and upper bound for the energy of a state with given
 * k and nodes, taking into account states that have been already
 * found with the same k
 * 
 * @param  nodes:   Target number of nodes in P
 * @param  k:       Quantum number k
 * @retval          Lower and upper bounds for the energy
 */
pair<double, double> DiracAtom::energyLimits(int nodes, int k)
{
    int n, l;
    bool s;
    double minE, maxE;
    map<tuple<int, int, bool>, DiracState>::iterator it;

    qnumDirac2Schro(k, l, s);
    qnumNodes2Principal(nodes, l, n);

    // Not below the lowest potential energy, and also use the hydrogen-like energy as a reasonable lower bound scale
    minE = max(V.V(0) + restE, hydrogenicDiracEnergy(Z, mu, 1, -1, true) * 2 + restE);
    // Required for the state to be bound
    maxE = restE;

    for (it = states.begin(); it != states.end(); it++)
    {
        int itn, itl;
        bool its;

        if (!it->second.converged)
            continue;

        itn = get<0>(it->first);
        itl = get<1>(it->first);
        its = get<2>(it->first);

        if (itl == l && its == s)
        {
            if (itn <= n)
            {
                minE = max(minE, it->second.E);
            }
            else
            {
                maxE = min(maxE, it->second.E);
            }
        }
    }

    return {minE, maxE};
}

/**
  * @brief  Converge a state to fall within an attraction basin with the required number of nodes
  * @note   Perform a preliminary, rough bisection search to find an energy that produces a
  * wavefunction with the desired number of nodes for the given state. This is then used as a 
  * starting point for full convergence.
  * 
  * @param  &state:     DiracState to integrate
  * @param  &tp:        TurningPoint object to store turning point info
  * @param  targ_nodes: Target number of nodes
  * @param  &minE:      Minimum energy (boundary will be updated through search)
  * @param  &maxE:      Maximum energy (boundary will be updated through search)
  * @retval None
  */
void DiracAtom::convergeNodes(DiracState &state, TurningPoint &tp, int targ_nodes, double &minE, double &maxE)
{
    int k;
    int nl = -1, nr = -1;
    double El, Er, oldEl = maxE + 1, oldEr = maxE + 1;
    pair<int, int> glim;

    k = state.k;
    El = minE + (maxE - minE) / 3.0;
    Er = maxE - (maxE - minE) / 3.0;

    LOG(DEBUG) << "Running convergeNodes to search energy with solution with " << targ_nodes << " nodes\n";

    for (int it = 0; it < maxit_nodes; ++it)
    {
        LOG(DEBUG) << "Iteration " << (it + 1) << ", El = " << El - restE << "+mc2, nl = " << nl << ", Er = "
                   << Er - restE << "+mc2, nr = " << nr << "\n";
        if (El != oldEl)
        {
            state = initState(El, k);
            integrateState(state, tp);
            state.continuify(tp);
            state.normalize();
            state.findNodes(nodetol);
            nl = state.nodes;
            if (nl == targ_nodes)
            {
                LOG(TRACE) << "State with " << targ_nodes << " nodes found at E = " << El - restE << "+mc2\n";
                return;
            }
        }

        if (Er != oldEr)
        {
            state = initState(Er, k);
            integrateState(state, tp);
            state.continuify(tp);
            state.normalize();
            state.findNodes(nodetol);
            nr = state.nodes;
            if (nr == targ_nodes)
            {
                LOG(TRACE) << "State with " << targ_nodes << " nodes found at E = " << Er - restE << "+mc2\n";
                return;
            }
        }

        LOG(TRACE) << "Nodes count: nl = " << nl << ", nr = " << nr << "\n";

        // Otherwise, what are their signs?
        int dl = (nl - targ_nodes);
        int dr = (nr - targ_nodes);

        if (dl > 0 && dr > 0)
        {
            // Both are too high
            oldEr = Er;
            Er = El;
            El = (minE + El) / 2.0;
            maxE = Er;
        }
        else if (dl < 0 && dr < 0)
        {
            // Both are too low
            oldEl = El;
            El = Er;
            Er = (maxE + Er) / 2.0;
            minE = El;
        }
        else if (dl < 0 && dr > 0)
        {
            // It's in between!
            oldEl = El;
            minE = El;
            El = (El + Er) / 2.0;
        }
        else
        {
            // Doesn't make sense
            throw runtime_error("convergeNodes failed - higher number of nodes for lower energy");
        }
    }

    throw runtime_error("convergeNodes failed to find a suitable state - maximum iterations hit");
}

void DiracAtom::convergeE(DiracState &state, TurningPoint &tp, double &minE, double &maxE)
{
    int k;
    double E, dE;
    double Edamp_eff = abs(Edamp);
    pair<int, int> glim;

    k = state.k;
    E = state.E;

    LOG(DEBUG) << "Running convergeE to search energy from starting value of " << E - restE << " + mc2\n";
    LOG(DEBUG) << "Energy limits: " << minE - restE << " + mc2 < E < " << maxE - restE << " + mc2\n";

    for (int it = 0; it < maxit_E; ++it)
    {
        LOG(TRACE) << "Iteration " << (it + 1) << ", E = " << E - restE << " + mc2\n";

        state = initState(E, k);
        integrateState(state, tp, dE);

        LOG(TRACE) << "Integration complete, computed error dE = " << dE << "\n";

        if (std::isnan(dE))
        {
            throw runtime_error("Invalid dE value returned by integrateState");
        }

        if (abs(dE) < Etol)
        {
            E = E - dE;
            LOG(TRACE) << "Convergence complete after " << (it + 1) << " iterations\n";
            state.E = E;
            state.continuify(tp);
            state.normalize();
            state.findNodes(nodetol);
            return;
        }
        // Apply maximum step ratio
        if (abs(dE / E) > max_dE_ratio)
        {
            dE = abs(E) * max_dE_ratio * (dE > 0 ? 1 : -1);
            LOG(TRACE) << "Step exceeds maximum allowed dE/E ratio, resized to " << dE << "\n";
        }
        E = E - dE * Edamp_eff;
        if (E > maxE)
        {
            // Something has gone wrong. Try to go back to a more reasonable search
            E = (maxE + E + dE * Edamp) / 2.0;
            LOG(TRACE) << "New energy exceeds maxE, resized to " << E - restE << " + mc2, reduced damping\n";
            Edamp_eff /= 2;
        }
        else if (E < minE)
        {
            // As above
            E = (minE + E + dE * Edamp) / 2.0;
            LOG(TRACE) << "New energy below minE, resized to " << E - restE << " + mc2, reduced damping\n";
            Edamp_eff /= 2;
        }
    }

    throw runtime_error("Convergence failed in given numer of iterations");
}

/**
 * @brief  Compute grid limits for given E and k
 * @note   Compute ideal indices to use as grid limits for
 * given E and k. This is done by making use of the fact
 * that the function decays exponentially past the 
 * turning point (for outer limit) and knowledge of its 
 * power law behaviour near the centre (for the inner one).
 * 
 * @param  E:           Energy of the state we need the grid for
 * @param  k:           Quantum number k of the state
 * @param  &failcode:   Error code in case of failure
 * @retval              Pair of limit indices for the grid {inner, outer}
 */
pair<int, int> DiracAtom::gridLimits(double E, int k)
{
    double B;
    double K = pow(mu * Physical::c, 2) - pow(E / Physical::c, 2);
    double gamma = pow(k, 2) - pow(Z * Physical::alpha, 2);
    double r_out, r_in, r_tp;
    int i_out, i_in;

    if (K < 0)
    {
        throw AtomErrorCode::UNBOUND_STATE;
        return {0, 0};
    }
    if (gamma < 0)
    {
        throw AtomErrorCode::SMALL_GAMMA;
        return {0, 0};
    }

    K = sqrt(K);
    gamma = sqrt(gamma);
    B = E - restE;

    r_tp = Z / abs(B); // Coulombic turning point radius

    LOG(TRACE) << "Computing optimal grid size for state with E = " << E - restE << " + mc2, k = " << k << "\n";
    LOG(TRACE) << "K = " << K << ", gamma = " << gamma << ", r_tp = " << r_tp << "\n";

    // Upper limit
    if (out_eps > 1 || out_eps < 0)
    {
        throw runtime_error("Invalid value for out_eps in DiracAtom; must be 0 < out_eps < 1");
    }
    r_out = r_tp - log(out_eps) / K;

    LOG(TRACE) << "Outer grid radius = " << r_out << "\n";

    // Lower limit
    if (in_eps > 1 || in_eps < 0)
    {
        throw runtime_error("Invalid value for in_eps in DiracAtom; must be 0 < in_eps < 1");
    }
    r_in = pow(in_eps, 1.0 / gamma) / M_E * gamma / K;

    LOG(TRACE) << "Inner grid radius = " << r_in << "\n";

    if (r_in > r_tp)
    {
        LOG(ERROR) << "Inner grid radius " << r_in << " is smaller than turning point radius " << r_tp << "; please decrease in_eps\n";
        throw runtime_error("Inner grid radius is too small for given atom and state; please decrease in_eps");
    }

    // Now get these as integer numbers of steps.
    i_out = ceil(log(r_out / rc) / dx);
    i_in = floor(log(r_in / rc) / dx);

    return {i_in, i_out};
}

/**
 * @brief  Initialise a DiracState based on E and k
 * @note   Initialise a DiracState based on energy E 
 * and quantum number k
 * 
 * @param  E:           Energy
 * @param  k:           Quantum number k
 * @retval              Initialised state
 */
DiracState DiracAtom::initState(double E, int k)
{
    DiracState state;
    pair<int, int> glimits;

    LOG(TRACE) << "Initialising state with E = " << E - restE << "+mc2, k = " << k << "\n";

    glimits = gridLimits(E, k);
    state = DiracState(rc, dx, glimits.first, glimits.second);
    state.k = k;
    state.E = E;
    state.V = getV(state.grid);

    return state;
}

/**
 * @brief  Integrate a DiracState of given E, k and V
 * @note   Perform a single integration of a DiracState,
 * given its E, k and V (which must be set in the DiracState
 * object itself).
 * 
 * @param  &state:  DiracState to integrate
 * @param  &tp:     TurningPoint object to store turning point info
 * @retval          
 */
void DiracAtom::integrateState(DiracState &state, TurningPoint &tp)
{
    int N;

    N = state.grid.size();
    if (N == 0)
    {
        throw runtime_error("Can not integrate state with zero-sized grid");
    }
    LOG(TRACE) << "Integrating state with grid of size " << N << "\n";
    // Start by applying boundary conditions
    boundaryDiracCoulomb(state.Q, state.P, state.grid, state.E, state.k, mu, Z, R > state.grid[0] ? R : -1);
    tp = shootDiracLog(state.Q, state.P, state.grid, state.V, state.E, state.k, mu, dx);
    LOG(TRACE) << "Integration complete, turning point found at " << tp.i << "\n";

    return;
}

/**
 * @brief  Integrate a DiracState of given E, k and V
 * @note   Perform a single integration of a DiracState,
 * given its E, k and V (which must be set in the DiracState
 * object itself). Computes also a suggested correction for the energy
 * 
 * @param  &state:  DiracState to integrate
 * @param  &tp:     TurningPoint object to store turning point info
 * @param  &dE:     Energy correction
 * @retval          
 */
void DiracAtom::integrateState(DiracState &state, TurningPoint &tp, double &dE)
{
    int N;
    double err;
    vector<double> y, zetai, zetae;

    integrateState(state, tp);

    N = state.grid.size();
    y = vector<double>(N, 0);
    zetai = vector<double>(N, 0);
    zetae = vector<double>(N, 0);

    err = tp.Qi / tp.Pi - tp.Qe / tp.Pe;
    // Compute the derivative of the error in dE
    for (int i = 0; i < N; ++i)
    {
        y[i] = state.Q[i] / state.P[i];
    }
    // First the forward version
    y[tp.i] = tp.Qi / tp.Pi;
    shootDiracErrorDELog(zetai, y, state.grid, state.V, tp.i, state.E, state.k, mu, dx);
    // Then the backwards one
    y[tp.i] = tp.Qe / tp.Pe;
    boundaryDiracErrorDECoulomb(zetae, state.E, state.k, mu);
    shootDiracErrorDELog(zetae, y, state.grid, state.V, tp.i, state.E, state.k, mu, dx, 'b');

    dE = err / (zetai[tp.i] - zetae[tp.i]);

    return;
}

DiracState DiracAtom::convergeState(int n, int k)
{
    int l;
    bool s;
    int targ_nodes;
    double minE, maxE;
    pair<double, double> Elim;
    DiracState state;
    TurningPoint tp;

    // Compute the required number of nodes
    qnumDirac2Schro(k, l, s);
    qnumPrincipal2Nodes(n, l, targ_nodes);

    // Find energy limits
    Elim = energyLimits(targ_nodes, k);
    minE = Elim.first;
    maxE = Elim.second;

    LOG(DEBUG) << "Converging state with n = " << n << ", k = " << k << "\n";
    LOG(DEBUG) << "Energy limits: " << minE - restE << " + mc2 < E < " << maxE - restE << " + mc2\n";

    for (int it = 0; it < maxit_state; ++it)
    {
        LOG(TRACE) << "Iteration " << (it + 1) << ", minE = " << minE - restE << "+mc2, maxE = " << maxE - restE << "+mc2\n";
        state.k = k;
        // Find appropriate basin
        convergeNodes(state, tp, targ_nodes, minE, maxE);
        // Now converge energy
        double hydroE = hydrogenicDiracEnergy(Z, mu, n, k);
        if (it == 0 && hydroE > minE && hydroE < maxE)
        {
            // We only try this the first time; if it fails, it ain't working any better later...
            LOG(TRACE) << "Using hydrogenic Dirac energy " << hydroE - restE << " + mc2 as starting guess\n";
            state.E = hydroE; // Speeds up things a lot when we got a broad interval
        }
        convergeE(state, tp, minE, maxE);

        // Check node condition
        if (state.nodes != targ_nodes)
        {
            LOG(TRACE) << "Converged state contains " << state.nodes << " nodes instead of " << targ_nodes << "\n";
            LOG(TRACE) << "Converged state has E = " << state.E - restE << "+mc2\n";
            if (state.nodes > targ_nodes)
            {
                maxE = min(maxE, state.E);
            }
            else
            {
                minE = max(minE, state.E);
            }
            // Store it for the future
            state.normalize();
            state.converged = true;
            states[make_tuple(n, l, s)] = state;
            state = DiracState();
        }
        else
        {
            state.normalize();
            state.converged = true;

            LOG(TRACE) << "Convergence achieved at E = " << state.E - restE << " + mc2\n";

            // And return
            return state;
        }
    }

    throw runtime_error("Failed to converge with convergeState");
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
    int k;
    int dnode;
    double E0;
    bool found = false;
    DiracState state;
    TurningPoint tp;

    qnumSchro2Dirac(l, s, k);

    // First, check if it's already calculated
    if (!force && states[make_tuple(n, l, s)].converged)
    {
        LOG(DEBUG) << "State with n = " << n << ", k = " << k << " already calculated";
        return;
    }

    try
    {
        state = convergeState(n, k);
    }
    catch (runtime_error re)
    {
        LOG(ERROR) << "Convergence failed with error: " << re.what() << "\n";
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

    if (!st.converged)
    {
        throw runtime_error("State is not converged");
    }

    return DiracState(st);
}

/**
 * @brief  Return a matrix of transition probabilities between two states
 * @note   Return a matrix containing all spontaneous transition probabilities 
 * between two states of given quantum numbers n, l and s; the matrix corresponds
 * to all possible transition between m states. The transition is expected to be 
 * from 1 to 2: in other words, if the energy of 2 is higher than that of 1, 
 * all probabilities will be returned as zero.
 * 
 * @param  n1: Principal quantum number of state 1
 * @param  l1: Orbital quantum number of state 1
 * @param  s1: Spin quantum number (true = 1/2 / false = -1/2) of state 1
 * @param  n2: Principal quantum number of state 2
 * @param  l2: Orbital quantum number of state 1
 * @param  s2: Spin quantum number (true = 1/2 / false = -1/2) of state 2
 * @param  approx_j0: If true, approximate the Bessel function j0(K*r) as 1 
 * @retval Transition matrix
 */
TransitionMatrix DiracAtom::getTransitionProbabilities(int n1, int l1, bool s1, int n2, int l2, bool s2, bool approx_j0)
{
    int k1, k2;

    qnumSchro2Dirac(l1, s1, k1);
    qnumSchro2Dirac(l2, s2, k2);

    TransitionMatrix tmat(k1, k2);

    // Get the relevant states
    DiracState psi1 = getState(n1, l1, s1);
    DiracState psi2 = getState(n2, l2, s2);

    float DE = psi1.E - psi2.E;
    float K = DE / Physical::c;

    if (DE < 0 || abs(l2 - l1) != 1)
    {
        // Invalid: state 2 has a higher energy, or transition forbidden
        return tmat;
    }

    // Now the integrals
    int i0 = max(psi1.grid_indices.first, psi2.grid_indices.first);
    int i1 = min(psi1.grid_indices.second, psi2.grid_indices.second);
    int delta1 = max(i0 - psi1.grid_indices.first, 0);
    int delta2 = max(i0 - psi2.grid_indices.first, 0);

    vector<double> intgrid = logGrid(rc, dx, i0, i1)[1];
    vector<double> kerP1Q2(i1 - i0 + 1), kerP2Q1(i1 - i0 + 1);

    for (int i = 0; i < intgrid.size(); ++i)
    {
        double j0 = (approx_j0? 1.0 : sinc(K * intgrid[i]));
        kerP1Q2[i] = psi1.P[i + delta1] * psi2.Q[i + delta2] * j0 * intgrid[i];
        kerP2Q1[i] = psi1.Q[i + delta1] * psi2.P[i + delta2] * j0 * intgrid[i];
    }

    double J12 = trapzInt(dx, kerP1Q2);
    double J21 = trapzInt(dx, kerP2Q1);

    int sgk1 = (k1 < 0 ? -1 : 1);
    int sgk2 = (k2 < 0 ? -1 : 1);

    // Now on to the full matrix elements
    for (int im1 = 0; im1 < tmat.m1.size(); ++im1)
    {
        for (int im2 = 0; im2 < tmat.m2.size(); ++im2)
        {
            double m1 = tmat.m1[im1];
            double m2 = tmat.m2[im2];

            double u1 = cgCoeff(k1, m1, true);
            double u2 = cgCoeff(k1, m1, false);
            double u3 = cgCoeff(-k1, m1, true);
            double u4 = cgCoeff(-k1, m1, false);
            double v1 = cgCoeff(k2, m2, true);
            double v2 = cgCoeff(k2, m2, false);
            double v3 = cgCoeff(-k2, m2, true);
            double v4 = cgCoeff(-k2, m2, false);

            double MSx = ((u1 * v4 * (m1 == m2 + 1) + u2 * v3 * (m1 + 1 == m2)) * (l1 == l2 - sgk2) * J12 +
                          (u3 * v2 * (m1 == m2 + 1) + u4 * v1 * (m1 + 1 == m2)) * (l1 - sgk1 == l2) * J21);
            double MSy = ((u1 * v4 * (m1 == m2 + 1) - u2 * v3 * (m1 + 1 == m2)) * (l1 == l2 - sgk2) * J12 +
                          (u3 * v2 * (m1 == m2 + 1) - u4 * v1 * (m1 + 1 == m2)) * (l1 - sgk1 == l2) * J21);
            double MSz = 2 * (u1 * v3 * (l1 == l2 - sgk2) * J12 + u3 * v1 * (l1 - sgk1 == l2) * J21) * (m1 == m2);

            tmat.T[im1][im2] = 4.0 / 3.0 * K * (pow(MSx, 2) + pow(MSy, 2) + pow(MSz, 2));
        }
    }

    /*
    CURRENT OBSERVED PROBLEMS:

    Based on comparison of results with Table III in Wilbur Payne's thesis

    * j = l+1/2 and j = l-1/2 states seem swapped
    * j = l+1/2 state (which now is returned by s = false instead of true) gives a rate approximately twice the correct one
     
     
    */

    return tmat;
}