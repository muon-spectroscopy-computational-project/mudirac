#include "../lib/elements.hpp"
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{

    for (int Z = 1; Z < 92; ++Z) {
        vector<int> As = getAllIsotopes(Z);
        for (int i = 0; i < As.size(); ++i) {
            int A = As[i];
            int N = A - Z;
            double A3 = pow(A, 1.0/3.0);
            double Z3 = pow(Z, 1.0/3.0);
            cout << Z << '\t' << A << '\t' << 1.223*A3 << '\t' << 1.631*Z3 << '\t';
            cout << sqrt(5.0/3.0*(pow(1.242*Z3, 2.0)+0.64)) << '\t' << 1.269*(1-(0.252*(N-Z))/A)*A3 << '\t';
            cout << 1.235*(1-(0.177*(N-Z))/A+1.960/A)*A3 << '\n';
        }
    }
}