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
 * @class BaseInputFile
 *
 * @brief  A specialised InputFile class for the input files in MuDirac
 * @note   This class initialises only the "xr_lines" keyword and defines a method
 * to parse these lines into transition quantum numbers. These are both required for the
 * dependent classes MuDiracInputFile and ExperimentalResultFile.
 *
 * @retval None
 */
class BaseInputFile : public InputFile {
  public:
    BaseInputFile(void);
  
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
    
};


/**
 * @class MuDiracInputFile
 *
 * @brief  A specialised BaseInputFile class for the config input file in mudirac
 * @note   A specialised BaseInputFile class for the config input file in mudirac.
 * Initialises all the relevant keywords.
 *
 * @retval None
 */
class MuDiracInputFile : public BaseInputFile {
 public:
  MuDiracInputFile(void);
  DiracAtom makeAtom();
  /**
  * @brief  validates the input config file to mudirac
   * @note  checks the number of arguments passed to mudirac and parses the
   * first argument as a config file. the filename also provides the seed
   * which is used for outputfiles in MuDirac.
   *
   * @param  argc: number of arguments to run mudirac
   * @param argv: list of the arguments.
   * @param seed: seed updated by reference
   * @retval None
   */
  void validate(int argc, char *argv[], string & seed);

  /**
   * 
   * @brief: Validates the input file for optimisations settings
   * @note: validates the number of arguments passed to mudirac, the nuclear model,
   * the choice of coordinate system and the optimisation algorithm.
   * 
   * @param args: argc the number of arguments passed to MuDirac
   * @param coords: the coordinate system (corrected to default if invalid)
   * @param min_2pF_algo: the minimisation algorithm (corrected to default if invalid)
   * @retval: None
   */
  void validateOptimisation(int args, string &coords, string &min_2pF_algo);

 private:
  map<string, NuclearRadiusModel> nucmodelmap = {
    {"POINT", POINT}, {"SPHERE", SPHERE}, {"FERMI2", FERMI2}
  };
};

#endif