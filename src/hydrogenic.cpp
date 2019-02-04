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
 * @brief  Eigenfunction of a given hydrogenic Schroedinger state
 * @note  he known eigenenfunction for a hydrogenic atom computed with
 * the Schroedinger equation and for given charge, mass and quantum numbers n and l.
 * 
 * 
 * 
 * @param  r: 
 * @param  Z: 
 * @param  mu: 
 * @param  n: 
 * @param  l: 
 * @retval 
 */
vector<double> hydrogenicSchroWavefunction(vector<double> r, double Z, double mu, int n, int l) {

}