/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * elements.hpp
 * 
 * Atomic mass and spin data extracted by AME2016 and NUBASE databases.
 * References:
 * 
 * AME2016 (masses):
 *      "The Ame2016 atomic mass evaluation (I)"   by W.J.Huang, G.Audi, M.Wang, F.G.Kondev, S.Naimi and X.Xu
 *         Chinese Physics C41 030002, March 2017.
 *      "The Ame2016 atomic mass evaluation (II)"  by M.Wang, G.Audi, F.G.Kondev, W.J.Huang, S.Naimi and X.Xu
 *         Chinese Physics C41 030003, March 2017.
 * 
 * NUBASE (spins):
 *      "The NUBASE2016 evaluation of nuclear properties"   by G.Audi, F.G.Kondev, M.Wang, W.J.Huang and S.Naimi
 *         Chinese Physics C41 030001, March 2017.
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <map>
#include <cmath>
#include <limits>

using namespace std;

struct isotope
{
    double mass;
    double spin;
};

struct element
{
    int Z;
    map<int, isotope> isotopes;
};

//{HERE GOES THE ACTUAL DATA}//