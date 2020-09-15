#include "../lib/atom.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"
#include "datapath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "../vendor/catch/catch.hpp"

using namespace std;

/*  These tests use reference values found from

E. Borie, G.A. Rinker
"The energy levels of muonic atoms"
Rev. Mod. Phys. 54, 1 (1982)

Table VI
*/

vector<double> energyCorrections(int Z, int n1, int l1, bool s1, int n2, int l2,
                                 bool s2, NuclearRadiusModel rmodel = SPHERE)
{
  vector<double> E_corr;
  DiracAtom da_point = DiracAtom(Z, Physical::m_mu, getElementMainIsotope(Z));
  DiracAtom da_sphere =
      DiracAtom(Z, Physical::m_mu, getElementMainIsotope(Z), rmodel);
  DiracAtom da_uehling =
      DiracAtom(Z, Physical::m_mu, getElementMainIsotope(Z), rmodel);
  da_uehling.setUehling(true, 100);

  double E_point = da_point.getState(n2, l2, s2).E;
  E_point -= da_point.getState(n1, l1, s1).E;
  E_point /= Physical::eV;

  E_corr.push_back(E_point);

  double E_sphere = da_sphere.getState(n2, l2, s2).E;
  E_sphere -= da_sphere.getState(n1, l1, s1).E;
  E_sphere /= Physical::eV;

  E_corr.push_back(E_sphere - E_point);

  double E_uehling = da_uehling.getState(n2, l2, s2).E;
  E_uehling -= da_uehling.getState(n1, l1, s1).E;
  E_uehling /= Physical::eV;

  E_corr.push_back(E_uehling - E_sphere);

  return E_corr;
}

TEST_CASE("Calculated energies and corrections", "[energyCorrections]")
{
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::warning,
                                      AixLog::Type::all);
  SECTION("12Mg")
  {
    // 3d3/2-2p1/2
    vector<double> E_corr = energyCorrections(12, 2, 1, false, 3, 2, false);

    REQUIRE(E_corr[0] == Approx(56213.9).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-2.8).epsilon(1)); // This is very loose
    REQUIRE(E_corr[2] == Approx(179.36).epsilon(1e-3));

    // 3d5/2-2p3/2
    E_corr = energyCorrections(12, 2, 1, true, 3, 2, true);

    REQUIRE(E_corr[0] == Approx(56038.9).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-0.9).epsilon(1));
    REQUIRE(E_corr[2] == Approx(177.55).epsilon(1e-3));
  }

  SECTION("14Si")
  {
    // 3d3/2-2p1/2
    vector<double> E_corr = energyCorrections(14, 2, 1, false, 3, 2, false);

    REQUIRE(E_corr[0] == Approx(76670.7).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-7.4).epsilon(1));
    REQUIRE(E_corr[2] == Approx(276.8).epsilon(1e-3));

    // 3d5/2-2p3/2
    E_corr = energyCorrections(14, 2, 1, true, 3, 2, true);

    REQUIRE(E_corr[0] == Approx(76345.5).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-2.6).epsilon(1));
    REQUIRE(E_corr[2] == Approx(273.2).epsilon(1e-3));
  }

  SECTION("82Pb")
  {
    // 5g7/2-4f5/2
    vector<double> E_corr = energyCorrections(82, 4, 3, false, 5, 4, false);

    REQUIRE(E_corr[0] == Approx(435664.0).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-10.0).epsilon(1));
    REQUIRE(E_corr[2] == Approx(2189.0).epsilon(1e-3));

    // 5g9/2-4f7/2
    E_corr = energyCorrections(82, 4, 3, true, 5, 4, true);

    REQUIRE(E_corr[0] == Approx(429343.0).epsilon(1e-3));
    REQUIRE(E_corr[1] == Approx(-4.0).epsilon(1));
    REQUIRE(E_corr[2] == Approx(2105.0).epsilon(1e-3));

    // 2p3/2-1s1/2
    E_corr = energyCorrections(82, 1, 0, true, 2, 1, true,
                               NuclearRadiusModel::FERMI2);
    double E_tot = E_corr[0] + E_corr[1] + E_corr[2];

    REQUIRE(E_tot == Approx(5.9627e6).epsilon(1e-3));
  }
}