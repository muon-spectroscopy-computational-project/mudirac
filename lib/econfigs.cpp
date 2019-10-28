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
 * 
 * @param  config:      The electronic configuration string
 * @retval              
 */
ElectronicConfiguration::ElectronicConfiguration(string config)
{
    epop = this->parseConfig(config);
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