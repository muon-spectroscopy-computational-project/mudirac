#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main(int argc, char **argv)
{
    AixLog::Log::init({make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "calc_Etrans.log"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::special, "calc_Etrans.err")});

    int Z;
    int n1, l1, n2, l2, k1, k2;
    bool s1, s2;
    double A;
    DiracState ds1, ds2;

    cout << "Input atomic number\n";
    cin >> Z;
    cout << "Input isotope mass (if <= 0, default is used)\n";
    cin >> A;

    if (A <= 0)
    {
        A = getIsotopeMass(Z);
    }

    cout << "Input first state\n";
    cin >> n1 >> l1 >> s1;
    cout << "Input second state\n";
    cin >> n2 >> l2 >> s2;

    qnumSchro2Dirac(l1, s1, k1);
    qnumSchro2Dirac(l2, s2, k2);

    // Now do the computation
    DiracAtom da_p = DiracAtom(Z, Physical::m_mu, A, POINT, 1, 0.001);
    DiracAtom da_s = DiracAtom(Z, Physical::m_mu, A, SPHERE, 1, 0.001);
    DiracAtom da_su = DiracAtom(Z, Physical::m_mu, A, SPHERE, 1, 0.001);
    da_su.setUehling(true, 100);

    ds1 = da_p.getState(n1, l1, s1);
    ds2 = da_p.getState(n2, l2, s2);
    double E_p = (ds2.E - ds1.E) / Physical::eV;

    try
    {
        ds1 = da_s.getState(n1, l1, s1);
        ds2 = da_s.getState(n2, l2, s2);
    }
    catch (char const *e)
    {
        cout << e << '\n';
        return -1;
    }

    double E_s = (ds2.E - ds1.E) / Physical::eV;

    ds1 = da_su.getState(n1, l1, s1);
    ds2 = da_su.getState(n2, l2, s2);

    double E_su = (ds2.E - ds1.E) / Physical::eV;

    cout << "E_PT\tdE_FS\taZa\tE_th\n";
    cout << E_p << '\t' << E_s - E_p << '\t' << E_su - E_s << '\t' << E_su << '\n';
}