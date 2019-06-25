/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 *
 * utils.cpp
 *
 * Utility functions of various sorts
 *
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "utils.hpp"

/**
 * @brief  Compute the effective mass of a two-body system
 *
 * @param  m1:  First mass
 * @param  m2:  Second mass
 * @retval      Effective mass
 */
double effectiveMass(double m1, double m2) { return (m1 * m2) / (m1 + m2); }

/**
 * @brief  Factorial of n
 *
 * @param  n:   Argument
 * @retval      n!
 */
int factorial(int n)
{
    int fn = 1;

    if (n < 0)
    {
        throw "Factorial of a negative number";
    }

    for (int i = 1; i <= n; ++i)
    {
        fn *= i;
    }

    return fn;
}

/**
 * @brief Sinc function, or sin(x)/x
 * 
 * @param x         Argument
 * @retval          sin(x)/x 
 */
double sinc(double x)
{
    return x <= 0 ? 1 : sin(x) / x;
}

/**
  * @brief  Linear interpolation
  * @note   Lineral interpolation between two values with a
  * parameter t
  * 
  * @param  a:      First value
  * @param  b:      Second value
  * @param  t:      Interpolation parameter (t=0 -> a, t=1 -> b)
  * @retval         Interpolated value
 */
double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}

/**
 * @brief  Generate a linear grid
 * @note   Generate a linearly spaced grid between two given points.
 *
 * @param  x0: Starting point.
 * @param  x1: End point.
 * @param  n: Number of points (default = 100).
 * @retval x: Vector containing the grid.
 */
vector<double> linGrid(double x0, double x1, int n)
{
    double dx = (x1 - x0) / (n - 1.0);
    vector<double> x(n);

    for (int i = 0; i < n; ++i)
    {
        x[i] = i * dx + x0;
    }

    return x;
}

/**
 * @brief  Generate a logarithmic grid
 * @note   Generate a logarithmically spaced grid between two given points.
 *
 * @param  x0: Starting point.
 * @param  x1: End point.
 * @param  n: Number of points (default = 100).
 * @retval grid: Vectors containing the grids, [x, exp(x)].
 */
vector<vector<double>> logGrid(double x0, double x1, int n)
{
    double dlog = log(x1 / x0) / (n - 1.0);
    vector<double> x(n), ex(n);

    for (int i = 0; i < n; ++i)
    {
        x[i] = i * dlog;
        ex[i] = x0 * exp(x[i]);
    }

    return vector<vector<double>>({x, ex});
}

/**
 * @brief  Generate a logarithmic grid
 * @note   Generate a logarithmically spaced grid from a central point
 * and a step, ranging between the given indices.
 *
 * @param  xc:  Central point
 * @param  dx:  Logarithmic step
 * @param  i0:  Starting index
 * @param  i1:  End index
 * @retval grid: Vectors containing the grids, [x, exp(x)].
 */
vector<vector<double>> logGrid(double xc, double dx, int i0, int i1)
{
    int n = i1 - i0 + 1;
    if (n <= 0)
    {
        throw invalid_argument("i1 must be greater or equal than i0 in logGrid");
    }
    vector<double> x(n), ex(n);

    for (int i = i0; i <= i1; ++i)
    {
        x[i - i0] = i * dx;
        ex[i - i0] = xc * exp(x[i - i0]);
    }

    return vector<vector<double>>({x, ex});
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at x
 *
 * @param  x:       Argument of the polynomial
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
double genLaguerrePoly(double x, int n, double alpha)
{
    // We do this recursively. Not the most efficient way, but since we only need
    // small numbers, it should be sufficient.

    double L = 1.0;

    if (n < 0)
    {
        throw "Parameter k of generalised Laguerre polynomial must be >= 0.";
    }

    switch (n)
    {
    case 0:
        return L;
        break;
    case 1:
        return (1.0 + alpha - x);
        break;
    default:
        return ((2 * n - 1 + alpha - x) * genLaguerrePoly(x, n - 1, alpha) -
                (n - 1 + alpha) * genLaguerrePoly(x, n - 2, alpha)) /
               n;
        break;
    }
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at multiple x
 *
 * @param  x:       Arguments of the polynomial
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
vector<double> genLaguerrePoly(vector<double> x, int n, int alpha)
{
    // Just iteratively go over all x
    vector<double> L(x.size());

    for (int i = 0; i < x.size(); ++i)
    {
        L[i] = genLaguerrePoly(x[i], n, alpha);
    }

    return L;
}

/**
 * @brief  Count the nodes in a given function
 * @note   Count the nodes (times the value crosses zero) in a
 * function. Includes only an interval comprised between two extremes that are
 * above a given tolerance.
 *
 * @param  v:   Vector containing the function
 * @param  tol: Tolerance
 * @retval      Number of nodes
 */
int countNodes(vector<double> v, double tol)
{
    int nc = 0;
    int i0, i1;

    for (i0 = 0; i0 < v.size(); ++i0)
    {
        if (abs(v[i0]) > tol)
        {
            break;
        }
    }
    for (i1 = v.size() - 1; i1 > i0; --i1)
    {
        if (abs(v[i1]) > tol)
        {
            break;
        }
    }
    for (int i = i0 + 1; i <= i1; ++i)
    {
        nc += ((v[i] * v[i - 1]) < 0);
    }

    return nc;
}

/**
 * @brief  Convert quantum numbers l and s to k
 * @note   Convert quantum numbers l and s as used in
 * the radial Schroedinger equation to their equivalent k
 * in the Dirac one
 * 
 * @param  l:       Orbital quantum number
 * @param  s:       Spin quantum number
 * @param  &k:      Dirac quantum number
 * @retval None
 */
void qnumSchro2Dirac(int l, bool s, int &k)
{
    k = ((!s && l > 0) ? l : -l - 1);
}

/**
 * @brief  Convert quantum number k to l and s
 * @note   Convert quantum number k as used in
 * the radial Dirac equation to its equivalent l
 * and s in the Schroedinger one.
 * 
 * @param  k:      Dirac quantum number
 * @param  &l:       Orbital quantum number
 * @param  &s:       Spin quantum number
 * @retval None
 */
void qnumDirac2Schro(int k, int &l, bool &s)
{
    s = (k < 0);
    l = s ? -(k + 1) : k;
}

/**
 * @brief  Convert from nodes and l to n
 * @note   Convert from number of nodes in the wavefunction and 
 * orbital quantum number l to principal quantum number
 * 
 * @param  nodes:   Number of nodes
 * @param  l:       Orbital quantum number
 * @param  &n:      Principal quantum number
 * @retval None
 */
void qnumNodes2Principal(int nodes, int l, int &n)
{
    n = nodes + l + 1;
}

/**
 * @brief  Convert from n and l to nodes
 * @note   Convert from principal quantum number
 * and orbital quantum number l to number of nodes in the wavefunction
 * 
 * @param  n:      Principal quantum number
 * @param  l:       Orbital quantum number
 * @param  &nodes:   Number of nodes
 * @retval None
 */
void qnumPrincipal2Nodes(int n, int l, int &nodes)
{
    nodes = n - l - 1;
}

/**
 * @brief  Split a string with a given separator
 * @note   Split a string s using each occurrence of a separator sep
 *
 * @param  s:       String to split
 * @param  sep:     Separator
 * @param  merge:   If true, merge multiple consecutive occurrences of sep
 * @retval          Vector of string fragments
 */
vector<string> splitString(string s, string sep, bool merge)
{
    int pos = -1;
    vector<string> fragments;

    while (s.size() > 0)
    {
        pos = s.find_first_of(sep);
        while (pos == 0 && merge)
        {
            s = s.substr(1);
            pos = s.find_first_of(sep);
        }
        if (s.size() == 0)
            break;
        fragments.push_back(s.substr(0, pos));
        s = s.substr(pos + 1);
        if (pos == string::npos)
            break;
    }

    return fragments;
}

/**
 * @brief   Strip a string of the given characters on both ends
 * @note    Strip a string of any instances of the given characters
 * on both ends. If no characters are passed, spaces are used.
 *
 * @param   s:          String to strip
 * @param   strip:      Characters to strip
 * @retval              Stripped string
 */
string stripString(string s, string strip)
{
    int pos;
    pos = s.find_first_not_of(strip);
    if (pos == -1)
    {
        // All to strip
        return "";
    }
    s = s.substr(pos);
    pos = s.find_last_not_of(strip);
    s = s.substr(0, pos + 1);

    return s;
}
