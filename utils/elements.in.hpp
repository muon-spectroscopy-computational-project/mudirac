/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
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
 * @version 1.0 20/03/2020
 */

#include <map>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>

using namespace std;

#ifndef MUDIRAC_ELEMENTS
#define MUDIRAC_ELEMENTS
struct isotope
{
    const double mass;
    const double spin;
    const double radius;
};

struct element
{
    const int Z;
    const int maxA;
    const map<int, isotope> isotopes;
};

double getIsotopeMass(string symbol, int isotope=-1);
double getIsotopeMass(int Z, int isotope=-1);
double getIsotopeSpin(string symbol, int isotope=-1);
double getIsotopeSpin(int Z, int isotope=-1);
double getIsotopeRadius(string symbol, int isotope=-1);
double getIsotopeRadius(int Z, int isotope=-1);
vector<int> getAllIsotopes(string symbol);
vector<int> getAllIsotopes(int Z);

int getElementZ(string symbol);
string getElementSymbol(int Z);
int getElementMainIsotope(string symbol);
int getElementMainIsotope(int Z);
#endif

#ifndef ATOMIC_DATA
#define ATOMIC_DATA
//{HERE GOES THE ACTUAL DATA}//
#endif