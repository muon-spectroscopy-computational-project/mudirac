/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * econfigs.hpp
 * 
 * Electronic configuration parsing and building
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
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
 * @param  Z:           Nuclear charge
 * @param  mu:          Effective mass of the electron
 * @param  shield:      If true, use the method outlined by Tauscher, Z. Phys. A, 1978 to
 *                      account for the effect of shielding of internal electrons. All 
 *                      electrons of shell N will thus see a nuclear charge shielded by
 *                      the sum total of electrons in the N-1 shells below
 * @param  dirac:       If true, use Dirac instead of Schroedinger orbitals for the density
 * @retval              
 */
ElectronicConfiguration::ElectronicConfiguration(string config, int Z, double mu, bool shield, bool dirac)
{
    epop = this->parseConfig(config);
    this->Z = Z;
    this->mu = mu;
    this->shield = shield;
    this->dirac = dirac;
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
int ElectronicConfiguration::getPopulation(int n, int l)
{
    if (l >= n || n < 0 || l < 0)
    {
        throw invalid_argument("Invalid quantum numbers for getPopulation");
    }
    if (n <= epop.size())
    {
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
int ElectronicConfiguration::maxn()
{
    return epop.size();
}

/**
 * @brief  Electronic charge density with hydrogen-like orbitals
 * @note   Electronic charge density at a given distance from the nucleus,
 * assuming that orbitals have an hydrogen-like structure.
 * 
 * @param  r:           Radial distance at which to compute the density
 * @retval 
 */
double ElectronicConfiguration::hydrogenicChargeDensity(double r)
{
    double Zn = Z;
    double rho = 0;

    LOG(INFO) << "Computing hydrogen-like charge density at r = " << r << " for atom with Z = " << Z << ", mu = " << mu << "\n";
    LOG(INFO) << "Shielding: " << shield << " - Dirac: " << dirac << "\n";

    for (int n = 1; n <= epop.size(); ++n)
    {
        int npop = 0;

        for (int l = 0; l < n; ++l)
        {
            if (dirac)
            {
                /* We need to take care here: for the same l-shell, the j = l-1/2 orbital ( k > 0 )
                is always lower energy than the other, so we fill that first */
                int ku, kd, pu, pd;
                // j = l - 1/2
                kd = l;
                pd = min(epop[n - 1][l], 2 * l);
                // j = l + 1/2
                ku = -l - 1;
                pu = max(epop[n - 1][l], 0);

                vector<double> PQ;

                PQ = hydrogenicDiracWavefunction(r, Zn, mu, n, kd);
                rho += pd * (pow(PQ[0], 2) + pow(PQ[1], 2));
                if (pu > 0)
                {
                    PQ = hydrogenicDiracWavefunction(r, Zn, mu, n, ku);
                    rho += pu * (pow(PQ[0], 2) + pow(PQ[1], 2));
                }
            }
            else
            {
                LOG(TRACE) << "Psi: " << hydrogenicSchroWavefunction(r, Zn, mu, n, l) << "\n";
                rho += epop[n - 1][l] * pow(hydrogenicSchroWavefunction(r, Zn, mu, n, l), 2);
            }

            LOG(TRACE) << pow(hydrogenicSchroWavefunction(r, Zn, mu, n, l), 2) << "\n";
            LOG(TRACE) << "n = " << n << ", l = " << l << ", rhoTot =" << rho << "\n";

            npop += epop[n - 1][l];
        }
        if (shield)
        {
            Zn -= npop;
        }
    }

    return -rho;
}

vector<vector<int>> ElectronicConfiguration::parseConfig(string config)
{
    vector<vector<int>> pop;

    try
    {
        string c = econfig_data.at(config);
        LOG(TRACE) << "Configuration for " << config << " identified as " << c << "\n";
        return this->parseConfig(c);
    }
    catch (out_of_range e)
    {
        // Nothing; it just means we need to do this the long way
    }

    // Split by space
    vector<string> ctok = splitString(config, " ", true);

    for (int i = 0; i < ctok.size(); ++i)
    {
        if (ctok[i][0] == '[')
        {
            // Assume it's a symbol
            vector<vector<int>> subpop = this->parseConfig(stripString(ctok[i], "[]"));

            // Now add them together
            if (subpop.size() > pop.size())
            {
                for (int j = pop.size(); j < subpop.size(); ++j)
                {
                    pop.push_back(vector<int>(j + 1, 0));
                }
            }

            for (int n = 1; n <= subpop.size(); ++n)
            {
                for (int l = 0; l < n; ++l)
                {
                    pop[n - 1][l] += subpop[n - 1][l];
                    if (pop[n - 1][l] > 2 * (2 * l + 1))
                    {
                        throw(invalid_argument("Invalid electronic configuration string"));
                    }
                }
            }
        }
        else
        {
            // Identify the various parts
            vector<string> np = splitString(ctok[i], "spdf");
            string orb = stripString(ctok[i], "0123456789");

            if (np.size() != 2)
            {
                throw(invalid_argument("Invalid electronic configuration string"));
            }

            int n, p, l;
            try
            {
                n = stoi(np[0]);
                p = stoi(np[1]);
            }
            catch (invalid_argument e)
            {
                throw(invalid_argument("Invalid electronic configuration string"));
            }

            l = string("spdf").find(orb);

            if (n < 0 || l < 0 || l >= n || p < 0 || p > 2 * (2 * l + 1))
            {
                throw(invalid_argument("Invalid electronic configuration string"));
            }

            LOG(TRACE) << "Configuration " << ctok[i] << " interpreted as n = " << n << ", l = " << l << ", pop = " << p << "\n";

            // Edit pop as required
            if (pop.size() < n)
            {
                for (int j = pop.size(); j < n; ++j)
                {
                    pop.push_back(vector<int>(j + 1, 0));
                }
            }

            pop[n - 1][l] += p;

            if (pop[n - 1][l] > 2 * (2 * l + 1))
            {
                throw(invalid_argument("Invalid electronic configuration string"));
            }
        }
    }

    return pop;
}