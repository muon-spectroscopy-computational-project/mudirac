/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * hydrogenic.cpp
 *
 * Computes known hydrogenic atom solutions for reference
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "hydrogenic.hpp"

/**
 * @brief Energy of a given hydrogenic Schroedinger state
 * @note Returns the known eigenenergy for a hydrogenic atom computed with
 * the Schroedinger equation and for given nuclear charge, mass and quantum number n.
 *
 * E_n = -mu*Z^2/(2*n^2)
 *
 * @param Z:    Nuclear charge (default 1)
 * @param mu:   Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param n:    Principal quantum number (default 1)
 * @return:     Energy
 */
double hydrogenicSchroEnergy(double Z, double mu, int n) {
  if (n < 0) {
    throw "Invalid principal quantum number for hydrogenic energy.";
  }
  return -mu * Z * Z / (2 * n * n);
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Schroedinger state
 * @note  Retyrns the known eigenenfunction for a hydrogenic atom computed at distance r with
 * the Schroedinger equation and for given nuclear charge, mass and quantum numbers n and l.
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

  double arg = Z * mu / n;
  double R;

  if (n < 1 || l < 0 || l >= n) {
    throw "Invalid quantum numbers for hydrogenic wavefunction.";
  }

  R = sqrt(pow(2 * arg, 3) * factorial(n - l - 1) / (2.0 * n * factorial(n + l))) * r * exp(-arg * r) * pow(2 * arg * r, l) *
      genLaguerrePoly(2 * arg * r, n - l - 1, 2 * l + 1);

  return R;
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Schroedinger state
 * @note  Returns the known eigenenfunction for a hydrogenic atom computed at multiple distances r with
 * the Schroedinger equation and for given nuclear charge, mass and quantum numbers n and l.
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
  double arg = Z * mu / n;
  double prefac = sqrt(pow(2 * arg, 3) * factorial(n - l - 1) / (2.0 * n * factorial(n + l)));

  if (n < 1 || l < 0 || l >= n) {
    throw "Invalid quantum numbers for hydrogenic wavefunction.";
  }

  // We duplicate the code to avoid recalculating the prefactor at every point
  for (int i = 0; i < r.size(); ++i) {
    R[i] = prefac * r[i] * exp(-arg * r[i]) * pow(2 * arg * r[i], l) *
           genLaguerrePoly(2 * arg * r[i], n - l - 1, 2 * l + 1);
  }

  return R;
}

/**
 * @brief Energy of a given hydrogenic Dirac state
 * @note Returns the known eigenenergy for a hydrogenic atom computed with
 * the Dirac equation and for given nuclear charge, mass and quantum numbers n and k.
 *
 * E_n = mu*c^2*{1+[Z*alpha/(n-|k|+(k^2-Z^2*alpha^2))]^2}^{-1/2}
 *
 * @param Z:    Nuclear charge (default 1)
 * @param mu:   Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param n:    Principal quantum number (default 1)
 * @param k:    Number related to the total orbital momentum j=|l+s|, k = -(j+1/2)*sign(s) (default -1)
 * @param bind: If true, only return the binding energy (substract mc^2)
 * @return:     Energy
 */
double hydrogenicDiracEnergy(double Z, double mu, int n, int k, bool bind) {

  if (n < 0 || abs(k) > n || k == 0 || k == n) {
    throw "Invalid quantum numbers for hydrogenic Dirac energy.";
  }

  double gamma = k * k - pow(Z * Physical::alpha, 2);
  if (gamma < 0) {
    throw "Invalid orbital energy for given Z";
  }
  gamma = sqrt(gamma);

  return mu * pow(Physical::c, 2) / sqrt(1 + pow(Z * Physical::alpha / (n - abs(k) + gamma), 2)) - (bind ? mu * pow(Physical::c, 2) : 0);
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Dirac state
 * @note  Returns the known eigenenfunction for a hydrogenic atom computed at distance r with
 * the Dirac equation and for given nuclear charge, mass and quantum numbers n and k.
 *
 * @param  r:   Distance at which to compute the wavefunction
 * @param  Z:   Nuclear charge (default 1)
 * @param  mu:  Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param  n:   Principal quantum number (default 1)
 * @param k:    Number related to the total orbital momentum j=|l+s|, k = -(j+1/2)*sign(s) (default -1)
 * @retval:     Vector of major and minor component of the Dirac bispinor
 */
vector<double> hydrogenicDiracWavefunction(double r, double Z, double mu, int n, int k) {

  vector<double> pq(2);
  double E, Ek, mc2, C, A, rho, rhodep, lagP, lagM;
  double gamma;

  if (n < 0 || abs(k) > n || k == 0 || k == n) {
    throw "Invalid quantum numbers for hydrogenic Dirac wavefunction.";
  }

  gamma = k * k - pow(Z * Physical::alpha, 2);
  if (gamma < 0) {
    throw "Invalid orbital energy for given Z";
  }
  gamma = sqrt(gamma);

  E = hydrogenicDiracEnergy(Z, mu, n, k);
  mc2 = mu * pow(Physical::c, 2);
  C = sqrt(pow(mc2, 2) - E * E) / Physical::c;
  rho = 2 * C * r;
  rhodep = pow(rho, gamma) * exp(-0.5 * rho);

  // LOG(TRACE) << gamma << " " << C << " "<< rho << " " << rhodep << "\n";

  /* Formulas are from Wikipedia. No reliable reference could be found - Wiki
  references "The Quantum Theory of the Hydrogen Atom" by Felix Nendzig,
  which is not an official publication.

  The formulas in these documents seem correct, but here we add a minus sign
  in front of the small component (the Q, or [1] here), because that's in accord
  with the conventions used elsewhere, including in boundary conditions for
  Coulomb potentials. This should not have a big effect on the final result as
  long as, in case of reconstruction of the full wavefunction, care is taken
  in picking the spherical harmonic/Weyl spinor part.
  */
  if (k == -n) {
    A = sqrt(C / (2 * n * (n + gamma) * gamma * tgamma(2 * gamma)));
    pq[0] = A * (n + gamma) * rhodep;
    pq[1] = -A * Z * Physical::alpha * rhodep;
  } else {
    Ek = E * k / (gamma * mc2);
    A = sqrt(C * factorial(n - abs(k) - 1) / (4 * k * (k - gamma) * (n - abs(k) + gamma) * tgamma(n - abs(k) + 2 * gamma + 1)) *
             (Ek + pow(Ek, 2)));
    lagP = rho * genLaguerrePoly(rho, n - abs(k) - 1, 2 * gamma + 1);
    lagM = (gamma * mc2 - k * E) / (Physical::c * C) * genLaguerrePoly(rho, n - abs(k), 2 * gamma - 1);
    pq[0] = A * rhodep * (Z * Physical::alpha * lagP + (gamma - k) * lagM);
    pq[1] = -A * rhodep * (Z * Physical::alpha * lagM + (gamma - k) * lagP);
  }

  return pq;
}

/**
 * @brief  Radial eigenfunction of a given hydrogenic Dirac state
 * @note  Returns the known eigenenfunction for a hydrogenic atom computed at multiple distances r with
 * the Dirac equation and for given nuclear charge, mass and quantum numbers n and k.
 *
 * @param  r:   Distances at which to compute the wavefunction
 * @param  Z:   Nuclear charge (default 1)
 * @param  mu:  Reduced mass of the system, m*m_N/(m+m_N) (default 1)
 * @param  n:   Principal quantum number (default 1)
 * @param k:    Number related to the total orbital momentum j=|l+s|, k = -(j+1/2)*sign(j-l) (default -1)
 * @retval:     Vector of major and minor component of the Dirac bispinor
 */
vector<vector<double>> hydrogenicDiracWavefunction(vector<double> r, double Z, double mu, int n, int k) {
  vector<vector<double>> PQ(2);
  PQ[0] = vector<double>(r.size());
  PQ[1] = vector<double>(r.size());

  double E, Ek, mc2, C, A, rho, rhodep, lagP, lagM;
  double gamma;

  if (n < 0 || abs(k) > n || k == 0 || k == n) {
    throw "Invalid quantum numbers for hydrogenic Dirac wavefunction.";
  }

  gamma = k * k - pow(Z * Physical::alpha, 2);
  if (gamma < 0) {
    throw "Invalid orbital energy for given Z";
  }
  gamma = sqrt(gamma);

  E = hydrogenicDiracEnergy(Z, mu, n, k);
  mc2 = mu * pow(Physical::c, 2);
  C = sqrt(pow(mc2, 2) - E * E) / Physical::c;

  /* Formulas are from Wikipedia. No reliable reference could be found - Wiki
  references "The Quantum Theory of the Hydrogen Atom" by Felix Nendzig,
  which is not an official publication.

  The formulas in these documents seem correct, but here we add a minus sign
  in front of the small component (the Q, or [1] here), because that's in accord
  with the conventions used elsewhere, including in boundary conditions for
  Coulomb potentials. This should not have a big effect on the final result as
  long as, in case of reconstruction of the full wavefunction, care is taken
  in picking the spherical harmonic/Weyl spinor part.
  */
  if (k == -n) {
    A = sqrt(C / (2 * n * (n + gamma) * gamma * tgamma(2 * gamma)));
    for (int i = 0; i < r.size(); ++i) {
      rho = 2 * C * r[i];
      rhodep = pow(rho, gamma) * exp(-0.5 * rho);
      PQ[0][i] = A * (n + gamma) * rhodep;
      PQ[1][i] = -A * Z * Physical::alpha * rhodep;
    }
  } else {
    Ek = E * k / (gamma * mc2);
    A = sqrt(C * factorial(n - abs(k) - 1) / (4 * k * (k - gamma) * (n - abs(k) + gamma) * tgamma(n - abs(k) + 2 * gamma + 1)) *
             (Ek + pow(Ek, 2)));
    for (int i = 0; i < r.size(); ++i) {
      rho = 2 * C * r[i];
      rhodep = pow(rho, gamma) * exp(-0.5 * rho);
      lagP = rho * genLaguerrePoly(rho, n - abs(k) - 1, 2 * gamma + 1);
      lagM = (gamma * mc2 - k * E) / (Physical::c * C) * genLaguerrePoly(rho, n - abs(k), 2 * gamma - 1);
      PQ[0][i] = A * rhodep * (Z * Physical::alpha * lagP + (gamma - k) * lagM);
      PQ[1][i] = -A * rhodep * (Z * Physical::alpha * lagM + (gamma - k) * lagP);
    }
  }

  return PQ;
}