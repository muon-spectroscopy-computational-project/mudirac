/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * config.hpp
 *
 * Functions and classes useful for interpreting configuration files - header
 * file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "atom.hpp"
#include "econfigs.hpp"
#include "elements.hpp"
#include "input.hpp"
#include <map>
#include <vector>

#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_CONFIG
#define MUDIRAC_CONFIG



/**
 * @class MuDiracInputFile
 *
 * @brief  A specialised InputFile class for the config input file in mudirac
 * @note   A specialised InputFile class for the config input file in mudirac.
 * Initialises all the relevant keywords.
 *
 * @retval None
 */
class MuDiracInputFile : public InputFile {
 public:
  MuDiracInputFile(void);
  DiracAtom makeAtom();
  void validateOptimisation(int args, string &coords, string &min_2pF_algo);
  
  /**
   * @brief A function which parses the input xray transition lines and returns the
   * related quantum numbers
   * @note converts the IUPAC notation for the transition lines to find the quantum numbers for
   * the levels in each transition. The transitions are validated so invalid transitions are skipped.
   * 
   * @returns a vector of TransLineSpec which stores n l and s quantum numbers
   * for the initial and final state for each transition.
   */
  vector<TransLineSpec> parseXRLines();

 private:
  map<string, NuclearRadiusModel> nucmodelmap = {
    {"POINT", POINT}, {"SPHERE", SPHERE}, {"FERMI2", FERMI2}
  };
};

#endif