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
    E = s.E;
    k = s.k;
    grid_indices = pair<int, int>(s.grid_indices);
    grid = vector<double>(s.grid);
    loggrid = vector<double>(s.loggrid);
    Q = vector<double>(s.Q);
    P = vector<double>(s.P);
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
 * @retval None
 */
void DiracState::findNodes()
{
    nodes = countNodes(P);
    nodesQ = countNodes(Q);
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
    V = CoulombSpherePotential(Z, R);
}

void Atom::setUehling(bool s, int usteps)
{
    use_uehling = s;
    if (s)
    {
        V_uehling = UehlingSpherePotential(Z, R, usteps);
    }
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
double Atom::sphereNuclearModel(double A)
{
    return 1.2 * Physical::fm * pow(A, 1.0 / 3.0);
}

DiracAtom::DiracAtom(double Z, double m, double A, NuclearRadiusModel radius_model, double fc, double dx) : Atom(Z, m, A, radius_model, fc, dx)
{
    restE = mu * pow(Physical::c, 2);
    LOG(TRACE) << "Dirac Atom created with rest energy = " << restE << "\n";
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

    // Not below the lowest potential energy, and not an antiparticle
    minE = max(V.V(0) + restE, -restE);
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
    int nl, nr;
    double El, Er, oldEl = maxE + 1, oldEr = maxE + 1;
    pair<int, int> glim;

    k = state.k;
    El = minE + (maxE - minE) / 3.0;
    Er = maxE - (maxE - minE) / 3.0;

    LOG(TRACE) << "Running convergeNodes to search energy with solution with " << targ_nodes << " nodes\n";

    for (int it = 0; it < maxit; ++it)
    {
        LOG(TRACE) << "Iteration " << (it + 1) << ", El = " << El - restE << "+mc2, Er = " << Er - restE << "+mc2\n";
        if (El != oldEl)
        {
            state = initState(El, k);
            integrateState(state, tp);
            state.continuify(tp);
            state.findNodes();
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
            state.findNodes();
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

    LOG(TRACE) << "Running convergeE to search energy from starting value of " << E - restE << " + mc2\n";
    LOG(TRACE) << "Energy limits: " << minE - restE << " + mc2 < E < " << maxE - restE << " + mc2\n";

    for (int it = 0; it < maxit; ++it)
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
            state.findNodes();
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

    LOG(TRACE) << "Computing optimal grid size for state with E = " << E -restE << " + mc2, k = " << k << "\n";
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
    boundaryDiracCoulomb(state.Q, state.P, state.grid, state.E, state.k, mu, Z, R > state.grid[0]);
    LOG(TRACE) << "Boundary conditions applied\n";
    tp = shootDiracLog(state.Q, state.P, state.grid, state.V, state.E, state.k, mu, dx);
    LOG(TRACE) << "Integration complete, turning point at " << tp.i << "\n";

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

    LOG(TRACE) << "Converging state with n = " << n << ", k = " << k << "\n";
    LOG(TRACE) << "Energy limits: " << minE-restE << " + mc2 < E < " << maxE-restE << " + mc2\n";

    for (int it = 0; it < maxit; ++it)
    {
        state.k = k;
        // Find appropriate basin
        convergeNodes(state, tp, targ_nodes, minE, maxE);
        // Now converge energy
        double hydroE = hydrogenicDiracEnergy(Z, mu, n, k);
        if (hydroE > minE && hydroE < maxE)
        {
            LOG(TRACE) << "Using hydrogenic Dirac energy " << hydroE - restE << " + mc2 as starting guess\n";
            state.E = hydroE; // Speeds up things a lot when we got a broad interval
        }
        convergeE(state, tp, minE, maxE);

        // Check node condition
        if (state.nodes != targ_nodes)
        {
            LOG(TRACE) << "Converged state contains " << state.nodes << " nodes instead of " << targ_nodes << "\n";
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

    // First, check if it's already calculated
    if (!force && states[make_tuple(n, l, s)].converged)
    {
        return;
    }

    qnumSchro2Dirac(l, s, k);
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