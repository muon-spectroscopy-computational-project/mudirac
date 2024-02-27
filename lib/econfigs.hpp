/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * econfigs.hpp
 *
 * Electronic configuration parsing and building - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <map>
#include <string>
#include <vector>

#include "utils.hpp"
#include "elements.hpp"
#include "hydrogenic.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_ECONFIGS
#define MUDIRAC_ECONFIGS

class ElectronicConfiguration {
 public:
  ElectronicConfiguration(string config="", int Z=-1, double mu=1.0, bool shield=false, bool dirac=false);

  int Z;
  double mu;
  bool shield;
  bool dirac;

  int getPopulation(int n, int l);
  int maxn();
  int totQ();
  double innerShellRadius();
  double outerShellRadius();
  double hydrogenicChargeDensity(double r);

 private:
  vector<vector<int>> epop;
  vector<int> Zshell;
  vector<vector<int>> parseConfig(string config);
};

const map<string, string> econfig_data = {
  {"H", "1s1"},
  {"He", "1s2"},
  {"Li", "[He] 2s1"},
  {"Be", "[He] 2s2"},
  {"B", "[He] 2s2 2p1"},
  {"C", "[He] 2s2 2p2"},
  {"N", "[He] 2s2 2p3"},
  {"O", "[He] 2s2 2p4"},
  {"F", "[He] 2s2 2p5"},
  {"Ne", "[He] 2s2 2p6"},
  {"Na", "[Ne] 3s1"},
  {"Mg", "[Ne] 3s2"},
  {"Al", "[Ne] 3s2 3p1"},
  {"Si", "[Ne] 3s2 3p2"},
  {"P", "[Ne] 3s2 3p3"},
  {"S", "[Ne] 3s2 3p4"},
  {"Cl", "[Ne] 3s2 3p5"},
  {"Ar", "[Ne] 3s2 3p6"},
  {"K", "[Ar] 4s1"},
  {"Ca", "[Ar] 4s2"},
  {"Sc", "[Ar] 3d1 4s2"},
  {"Ti", "[Ar] 3d2 4s2"},
  {"V", "[Ar] 3d3 4s2"},
  {"Cr", "[Ar] 3d5 4s1"},
  {"Mn", "[Ar] 3d5 4s2"},
  {"Fe", "[Ar] 3d6 4s2"},
  {"Co", "[Ar] 3d7 4s2"},
  {"Ni", "[Ar] 3d8 4s2"},
  {"Cu", "[Ar] 3d10 4s1"},
  {"Zn", "[Ar] 3d10 4s2"},
  {"Ga", "[Ar] 3d10 4s2 4p1"},
  {"Ge", "[Ar] 3d10 4s2 4p2"},
  {"As", "[Ar] 3d10 4s2 4p3"},
  {"Se", "[Ar] 3d10 4s2 4p4"},
  {"Br", "[Ar] 3d10 4s2 4p5"},
  {"Kr", "[Ar] 3d10 4s2 4p6"},
  {"Rb", "[Kr] 5s1"},
  {"Sr", "[Kr] 5s2"},
  {"Y", "[Kr] 4d1 5s2"},
  {"Zr", "[Kr] 4d2 5s2"},
  {"Nb", "[Kr] 4d4 5s1"},
  {"Mo", "[Kr] 4d5 5s1"},
  {"Tc", "[Kr] 4d5 5s2"},
  {"Ru", "[Kr] 4d7 5s1"},
  {"Rh", "[Kr] 4d8 5s1"},
  {"Pd", "[Kr] 4d10"},
  {"Ag", "[Kr] 4d10 5s1"},
  {"Cd", "[Kr] 4d10 5s2"},
  {"In", "[Kr] 4d10 5s2 5p1"},
  {"Sn", "[Kr] 4d10 5s2 5p2"},
  {"Sb", "[Kr] 4d10 5s2 5p3"},
  {"Te", "[Kr] 4d10 5s2 5p4"},
  {"I", "[Kr] 4d10 5s2 5p5"},
  {"Xe", "[Kr] 4d10 5s2 5p6"},
  {"Cs", "[Xe] 6s1"},
  {"Ba", "[Xe] 6s2"},
  {"La", "[Xe] 5d1 6s2"},
  {"Ce", "[Xe] 4f1 5d1 6s2"},
  {"Pr", "[Xe] 4f3 6s2"},
  {"Nd", "[Xe] 4f4 6s2"},
  {"Pm", "[Xe] 4f5 6s2"},
  {"Sm", "[Xe] 4f6 6s2"},
  {"Eu", "[Xe] 4f7 6s2"},
  {"Gd", "[Xe] 4f7 5d1 6s2"},
  {"Tb", "[Xe] 4f9 6s2"},
  {"Dy", "[Xe] 4f10 6s2"},
  {"Ho", "[Xe] 4f11 6s2"},
  {"Er", "[Xe] 4f12 6s2"},
  {"Tm", "[Xe] 4f13 6s2"},
  {"Yb", "[Xe] 4f14 6s2"},
  {"Lu", "[Xe] 4f14 5d1 6s2"},
  {"Hf", "[Xe] 4f14 5d2 6s2"},
  {"Ta", "[Xe] 4f14 5d3 6s2"},
  {"W", "[Xe] 4f14 5d4 6s2"},
  {"Re", "[Xe] 4f14 5d5 6s2"},
  {"Os", "[Xe] 4f14 5d6 6s2"},
  {"Ir", "[Xe] 4f14 5d7 6s2"},
  {"Pt", "[Xe] 4f14 5d9 6s1"},
  {"Au", "[Xe] 4f14 5d10 6s1"},
  {"Hg", "[Xe] 4f14 5d10 6s2"},
  {"Tl", "[Xe] 4f14 5d10 6s2 6p1"},
  {"Pb", "[Xe] 4f14 5d10 6s2 6p2"},
  {"Bi", "[Xe] 4f14 5d10 6s2 6p3"},
  {"Po", "[Xe] 4f14 5d10 6s2 6p4"},
  {"At", "[Xe] 4f14 5d10 6s2 6p5"},
  {"Rn", "[Xe] 4f14 5d10 6s2 6p6"},
  {"Fr", "[Rn] 7s1"},
  {"Ra", "[Rn] 7s2"},
  {"Ac", "[Rn] 6d1 7s2"},
  {"Th", "[Rn] 6d2 7s2"},
  {"Pa", "[Rn] 5f2 6d1 7s2"},
  {"U", "[Rn] 5f3 6d1 7s2"},
  {"Np", "[Rn] 5f4 6d1 7s2"},
  {"Pu", "[Rn] 5f6 7s2"},
  {"Am", "[Rn] 5f7 7s2"},
  {"Cm", "[Rn] 5f7 6d1 7s2"},
  {"Bk", "[Rn] 5f9 7s2"},
  {"Cf", "[Rn] 5f10 7s2"},
  {"Es", "[Rn] 5f11 7s2"},
  {"Fm", "[Rn] 5f12 7s2"},
  {"Md", "[Rn] 5f13 7s2"},
  {"No", "[Rn] 5f14 7s2"},
  {"Lr", "[Rn] 5f14 7s2 7p1"},
  {"Rf", "[Rn] 5f14 6d2 7s2"},
  {"Db", "[Rn] 5f14 6d3 7s2"},
  {"Sg", "[Rn] 5f14 6d4 7s2"},
  {"Bh", "[Rn] 5f14 6d5 7s2"},
  {"Hs", "[Rn] 5f14 6d6 7s2"},
  {"Mt", "[Rn] 5f14 6d7 7s2"},
  {"Ds", "[Rn] 5f14 6d9 7s1"},
  {"Rg", "[Rn] 5f14 6d10 7s1"},
  {"Cn", "[Rn] 5f14 6d10 7s2"}
};

#endif