#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <tuple>
#include <limits>
#include "../lib/atom.hpp"

using namespace std;

int main()
{
    double E;
    DiracAtom da = DiracAtom(80, Physical::m_mu, -1);

    for (int n = 1; n < 10; ++n) {
        E = hydrogenicDiracEnergy(da.getZ(), da.getmu(), n, -1);
        da.optimalGrid(E, -1);
    }
}
