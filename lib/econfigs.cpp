/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * econfigs.hpp
 *
 * Electronic configuration parsing and building
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "econfigs.hpp"

/**
 * @brief  Create an ElectronicConfiguration instance from a string
 * @note   Create an ElectronicConfiguration instance by parsing a string
 * in which the various shells are separated by spaces;
 * e.g. "1s2 2s2 2p4" would be the standard configuration of
 * oxygen.
 * Can include a single element symbol in square brackets
 * to mean the full configuration of that element ("[Ne] 3s2" for Mg), or even
 * be just a single element symbol ("Mg" for non-ionized Mg).
 * The other parameters define the behaviour used when computing the precise
 * charge density at a given radius.
 *
 * @param  config:      The electronic configuration string
 * @param  Z:           Nuclear charge. If negative, calculate it to neutralise the electrons
 * @param  mu:          Effective mass of the electron
 * @param  shield:      If true, use the method outlined by Tauscher, Z. Phys. A, 1978 to
 *                      account for the effect of shielding of internal electrons. All
 *                      electrons of shell N will thus see a nuclear charge shielded by
 *                      the sum total of electrons in the N-1 shells below
 * @param  dirac:       If true, use Dirac instead of Schroedinger orbitals for the density
 * @retval
 */
ElectronicConfiguration::ElectronicConfiguration(string config, int Z, double mu, bool shield, bool dirac) {
  epop = this->parseConfig(config);

  if (Z < 0) {
    Z = -totQ();
  }

  this->Z = Z;
  this->mu = mu;
  this->shield = shield;
  this->dirac = dirac;

  // Assign the shell Z
  if (shield) {
    Zshell = vector<int>(0);
    int etot = 0;
    for (int n = 1; n <= epop.size(); ++n) {
      if (Z-etot <= 0) {
        LOG(INFO) << "Truncating electronic configuration at n = " << n << " due to the entire nuclear charge being shielded\n";
        epop = vector<vector<int>>(epop.begin(), epop.begin()+n-1);
        break;
      }
      Zshell.push_back(Z-etot);
      LOG(TRACE) << "Shell n = " << n << ", effective Z = " << Zshell[n-1] << "\n";
      for (int l = 0; l < n; ++l) {
        etot += epop[n-1][l];
        LOG(TRACE) << "Channel l = " << l << ", population = " << epop[n-1][l] << "\n";
      }
    }
  } else {
    Zshell = vector<int>(epop.size(), Z);
  }
}

/**
 * @brief  Return the population of a given shell
 * @note   Return the population of a shell given the quantum
 * numbers n and l
 *
 * @param  n:   Principal quantum number
 * @param  l:   Orbital quantum number
 * @retval      Population
 */
int ElectronicConfiguration::getPopulation(int n, int l) {
  if (l >= n || n < 0 || l < 0) {
    throw invalid_argument("Invalid quantum numbers for getPopulation");
  }
  if (n <= epop.size()) {
    return epop[n - 1][l];
  }

  return 0;
}

/**
 * @brief  Maximum principal quantum number of this configuration
 * @note   Maximum principal quantum number of this configuration.
 *
 * @retval Maximum n for this configuration
 */
int ElectronicConfiguration::maxn() {
  return epop.size();
}

/**
 * @brief  Radius of the innermost electronic shell
 * @note   Radius of the innermost electronic shell estimated
 * as
 *
 * a = n/Z
 *
 * for n the smallest shell with population greater than zero.
 *
 * @retval The innermost shell's radius
 */
double ElectronicConfiguration::innerShellRadius() {
  int nmin = 0;

  for (int n = 1; n <= epop.size(); ++n) {
    int pop = 0;
    for (int l = 0; l < epop[n - 1].size(); ++l) {
      pop += epop[n - 1][l];
    }
    if (pop > 0) {
      nmin = n;
      break;
    }
  }

  return nmin * 1.0 / Z;
}

/**
 * @brief  Radius of the outermost electronic shell
 * @note   Radius of the outermost electronic shell estimated empirically as
 *
 * a = n**3/(2Z)
 *
 * for n the largest shell with population greater than zero.
 *
 * @retval The outermost shell's radius
 */
double ElectronicConfiguration::outerShellRadius() {
  int nmax = epop.size();

  if (nmax < 1) {
    return 0.0;
  }

  return pow(nmax, 3.0) * 0.5 / Zshell[nmax-1];
}

/**
 * @brief  Total charge of this configuration
 * @note   Total charge of this configuration
 *
 * @retval Total electronic charge for this configuration
 */
int ElectronicConfiguration::totQ() {
  int Q = 0;
  for (int n = 0; n < epop.size(); ++n) {
    for (int l = 0; l < epop[n].size(); ++l) {
      Q -= epop[n][l];
    }
  }

  return Q;
}

/**
 * @brief  Electronic charge density with hydrogen-like orbitals
 * @note   Electronic charge density at a given distance from the nucleus,
 * assuming that orbitals have an hydrogen-like structure.
 *
 * @param  r:           Radial distance at which to compute the density
 * @retval
 */
double ElectronicConfiguration::hydrogenicChargeDensity(double r) {
  double rho = 0;
  int Zn;

  for (int n = 1; n <= epop.size(); ++n) {
    Zn = Zshell[n-1];

    for (int l = 0; l < n; ++l) {
      if (dirac) {
        /* We need to take care here: for the same l-shell, the j = l-1/2 orbital ( k > 0 )
        is always lower energy than the other, so we fill that first */
        int ku, kd, pu, pd;
        vector<double> PQ;

        if (l > 0) {
          // j = l - 1/2
          kd = l;
          pd = min(epop[n - 1][l], 2 * l);
          // j = l + 1/2
          ku = -l - 1;
          pu = max(epop[n - 1][l] - pd, 0);
        } else {
          kd = 1;
          pd = 0;
          ku = -1;
          pu = epop[n - 1][l];
        }

        if (pd > 0) {
          PQ = hydrogenicDiracWavefunction(r, Zn, mu, n, kd);
          rho += pd * (pow(PQ[0], 2) + pow(PQ[1], 2));
        }
        if (pu > 0) {
          PQ = hydrogenicDiracWavefunction(r, Zn, mu, n, ku);
          rho += pu * (pow(PQ[0], 2) + pow(PQ[1], 2));
        }
      } else {
        rho += epop[n - 1][l] * pow(hydrogenicSchroWavefunction(r, Zn, mu, n, l), 2);
      }
    }
  }

  return -rho;
}

vector<vector<int>> ElectronicConfiguration::parseConfig(string config) {
  vector<vector<int>> pop;

  try {
    string c = econfig_data.at(config);
    LOG(TRACE) << "Configuration for " << config << " identified as " << c << "\n";
    return this->parseConfig(c);
  } catch (out_of_range e) {
    // Nothing; it just means we need to do this the long way
  }

  // Split by space
  vector<string> ctok = splitString(config, " ", true);

  for (int i = 0; i < ctok.size(); ++i) {
    if (ctok[i][0] == '[') {
      // Assume it's a symbol
      vector<vector<int>> subpop = this->parseConfig(stripString(ctok[i], "[]"));

      // Now add them together
      if (subpop.size() > pop.size()) {
        for (int j = pop.size(); j < subpop.size(); ++j) {
          pop.push_back(vector<int>(j + 1, 0));
        }
      }

      for (int n = 1; n <= subpop.size(); ++n) {
        for (int l = 0; l < n; ++l) {
          pop[n - 1][l] += subpop[n - 1][l];
          if (pop[n - 1][l] > 2 * (2 * l + 1)) {
            throw(invalid_argument("Invalid electronic configuration string"));
          }
        }
      }
    } else {
      // Identify the various parts
      vector<string> np = splitString(ctok[i], "spdf");
      string orb = stripString(ctok[i], "0123456789");

      if (np.size() != 2) {
        throw(invalid_argument("Invalid electronic configuration string"));
      }

      int n, p, l;
      try {
        n = stoi(np[0]);
        p = stoi(np[1]);
      } catch (invalid_argument e) {
        throw(invalid_argument("Invalid electronic configuration string"));
      }

      l = string("spdf").find(orb);

      if (n < 0 || l < 0 || l >= n || p < 0 || p > 2 * (2 * l + 1)) {
        throw(invalid_argument("Invalid electronic configuration string"));
      }

      LOG(TRACE) << "Configuration " << ctok[i] << " interpreted as n = " << n << ", l = " << l << ", pop = " << p << "\n";

      // Edit pop as required
      if (pop.size() < n) {
        for (int j = pop.size(); j < n; ++j) {
          pop.push_back(vector<int>(j + 1, 0));
        }
      }

      pop[n - 1][l] += p;

      if (pop[n - 1][l] > 2 * (2 * l + 1)) {
        throw(invalid_argument("Invalid electronic configuration string"));
      }
    }
  }

  return pop;
}