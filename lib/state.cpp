/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * state.cpp
 * 
 * State classes, designed to store the solution of quantum mechanical equations
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "state.hpp"

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