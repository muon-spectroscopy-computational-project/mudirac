/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * boundary.cpp
 * 
 * Boundary conditions
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "boundary.hpp"

/**
 * @brief  Impose boundary conditions to a Dirac wavefunction based on a Coulomb potential
 * @note   Impose boundary conditions to the Q and P components of a Dirac radial wavefunction based on a Coulomb potential.
 * The conditions take the form detailed in Richard R Silbar and T Goldman, "Solving the radial dirac equations: a
 * numerical odyssey" - European Journal of Physics, 32(1):217–233, Dec 2010, and I. P. Grant, "B-spline methods for
 * radial Dirac equations" - J. Phys. B: At. Mol. Opt. Phys., 42:055002, Jan 2009, and found in the Development Guide.
 * This method is used before passing Q and P to a routine such as shootDiracLog for integration.
 * 
 * @param  &Q: Vector for Q. Must have size greater than 4.
 * @param  &P: Vector for P. Must have size equal to Q.
 * @param  r:  Radial grid
 * @param  E:  Energy (binding + mc^2)
 * @param  k:  Quantum number (default = -1)
 * @param  m:  Mass of the particle (default = 1)
 * @param  Z:  Nuclear charge (default = 1)
 * @param  R:  Nuclear radius; considered point-like if <= 0 (default = -1)
 * @retval None
 */
void boundaryDiracCoulomb(vector<double> &Q, vector<double> &P, vector<double> r, double E, int k, double m, double Z, double R)
{
    int l; // Angular momentum number
    bool s;
    int N = Q.size();
    double K;

    qnumDirac2Schro(k, l, s);

    if (P.size() != N || r.size() != N || N < 4)
    {
        throw invalid_argument("Invalid array size passed to boundaryDiracCoulomb");
    }

    K = pow(m * Physical::c, 2) - pow(E * Physical::alpha, 2);
    if (K < 0)
    {
        throw invalid_argument("Can't compute boundary conditions for non-bound state");
    }
    K = sqrt(K);

    LOG(TRACE) << "Computing Coulomb boundary conditions for Dirac wavefunction, k = " << k << ", K = " << K << "\n";

    // r = 0 limit
    // Depends on whether we consider the nucleus of finite size or point-like
    if (R <= 0)
    {
        double gamma = k * k - pow(Physical::alpha * Z, 2.0);
        if (gamma < 0)
        {
            throw invalid_argument("Can't compute boundary conditions for state with point-like nucleus and negative gamma");
        }
        gamma = sqrt(gamma);
        // Point like
        // We set P, the major component, as having a prefactor of 1, and scale Q accordingly
        for (int i = 0; i < 2; ++i)
        {
            P[i] = pow(r[i], gamma);
            Q[i] = Z * Physical::alpha / (gamma - k) * P[i];
        }
    }
    else
    {
        // Finite size
        for (int i = 0; i < 2; ++i)
        {
            if (k < 0)
            {
                P[i] = pow(r[i], -k);
                Q[i] = pow(r[i], -k + 1) * 1.5 * Z * Physical::alpha / (R * (-2 * k + 1));
            }
            else
            {
                P[i] = pow(r[i], k + 2);
                Q[i] = -pow(r[i], k + 1) * Physical::c * (2 * k + 1) * R / (1.5 * Z);
            }
        }
    }

    LOG(TRACE) << "Boundary conditions at r => 0 (" << r[0] << "), P = [" << P[0] << "," << P[1] << "], Q = [" << Q[0] << "," << Q[1] << "]\n";

    // r = inf limit
    // Same as above
    for (int i = 1; i < 3; ++i)
    {
        P[N - i] = exp(-K * r[N - i]);
        Q[N - i] = -K / (m * Physical::c + E * Physical::alpha) * P[N - i];

        if (P[N - i] == 0)
        {
            // We went below numerical precision!
            LOG(ERROR) << SPECIAL << "Boundary conditions give zero at r = " << r[N - i] << "; a smaller grid is necessary\n";
            throw runtime_error("Boundary conditions give zero at the outside edge - you may need a smaller grid");
        }
    }

    LOG(TRACE) << "Boundary conditions at r => inf (" << r[N - 1] << "), P = [" << P[N - 2] << "," << P[N - 1] << "], Q = [" << Q[N - 2] << "," << Q[N - 1] << "]\n";
}

/**
 * @brief  Impose boundary conditions to d/dE (Q/P) for a Dirac wavefunction based on a Coulomb potential
 * @note   Impose boundary conditions on zeta = d/dE (Q/P) for a Dirac wavefunction based on a Coulomb potential. 
 * Zeta is used to find the optimal energy correction at each step when converging a solution. This function is 
 * used before passing zeta to a routine such as shootDiracErrorDELog for integration.
 * 
 * @param  &zeta: Vector for zeta. Must have size greater than 4.
 * @param  E:  Energy (binding + mc^2)
 * @param  k:  Quantum number (default = -1)
 * @param  m:  Mass of the particle (default = 1)
 * @retval None
 */
void boundaryDiracErrorDECoulomb(vector<double> &zeta, double E, int k, double m)
{
    int l; // Angular momentum number
    bool s;
    int N = zeta.size();
    double K, gp;

    qnumDirac2Schro(k, l, s);

    if (N < 4)
    {
        throw "Invalid array size passed to boundaryDiracErrorDECoulomb";
    }

    // In the r=0 limit, it's fine to have it be 0
    zeta[0] = 0;
    zeta[1] = 0;

    // r = inf limit
    K = pow(m * Physical::c, 2) - pow(E * Physical::alpha, 2);
    gp = (m * Physical::c + E * Physical::alpha);
    if (K < 0)
    {
        throw "Can't compute boundary conditions for non-bound state";
    }
    K = sqrt(K);
    zeta[N - 1] = E / (K * gp) * pow(Physical::alpha, 2) + K / pow(gp, 2) * Physical::alpha;
    zeta[N - 2] = zeta[N - 1];
}