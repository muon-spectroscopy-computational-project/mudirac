/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * hydrogenic.cpp
 * 
 * Computes known hydrogenic atom solutions for reference
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "hydrogenic.hpp"

/**
 * @brief Eigenenergy of a given hydrogenic Schroedinger state
 * @note Returns the known eigenenergy for a hydrogenic atom computed with
 * the Schroedinger equation and for given charge, mass and quantum number n.
 * 
 * E_n = -mu*Z^2/(2*n^2)
 * 
 * @param Z:    Nuclear charge (default 1)
 * @param mu:   Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param n:    Principal quantum number (default 1)
 * @return:     Energy
 */
double hydrogenicSchroEnergy(double Z, double mu, int n) {
    return -mu*Z*Z/(2*n*n);
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Schroedinger state
 * @note  Retyrns the known eigenenfunction for a hydrogenic atom computed at distance r with
 * the Schroedinger equation and for given charge, mass and quantum numbers n and l.
 * 
 *  R_nl(r) = [(2*Z*mu/n)^3*(n-l-1)!/(2*n*(n+l)!)]^0.5*exp(-Z*r*mu/n)*(2*Z*r*mu/n)^l*L_{n-l-1}^{2*l+1}(2*Z*r*mu/n)
 * 
 * @param  r:   Distance at which to compute the wavefunction
 * @param  Z:   Nuclear charge (default 1)
 * @param  mu:  Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param  n:   Principal quantum number (default 1)
 * @param  l:   Orbital quantum number (default 0)
 * @retval      Radial wavefunction value
 */
double hydrogenicSchroWavefunction(double r, double Z, double mu, int n, int l) {
    
    double arg = Z*mu/n;
    double R;

    if (n < 1 || l < 0 || l >= n) {
        throw "Invalid quantum numbers for hydrogenic wavefunction.";
    }

    R = sqrt(pow(2*arg, 3)*factorial(n-l-1)/(2.0*n*factorial(n+l)))*exp(-arg*r)*pow(2*arg*r, l)*genLaguerrePoly(2*arg*r, n-l-1, 2*l+1);

    return R;
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Schroedinger state
 * @note  Retyrns the known eigenenfunction for a hydrogenic atom computed at multiple distances r with
 * the Schroedinger equation and for given charge, mass and quantum numbers n and l.
 * 
 *  R_nl(r) = [(2*Z*mu/n)^3*(n-l-1)!/(2*n*(n+l)!)]^0.5*exp(-Z*r*mu/n)*(2*Z*r*mu/n)^l*L_{n-l-1}^{2*l+1}(2*Z*r*mu/n)
 * 
 * @param  r:   Distances at which to compute the wavefunction
 * @param  Z:   Nuclear charge (default 1)
 * @param  mu:  Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param  n:   Principal quantum number (default 1)
 * @param  l:   Orbital quantum number (default 0)
 * @retval      Radial wavefunction value
 */
vector<double> hydrogenicSchroWavefunction(vector<double> r, double Z, double mu, int n, int l) {

    vector<double> R(r.size());

    for (int i = 0; i < r.size(); ++i) {
        R[i] = hydrogenicSchroWavefunction(r[i], Z, mu, n, l);
    }
}