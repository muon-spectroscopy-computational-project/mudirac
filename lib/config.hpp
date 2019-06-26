/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * config.hpp
 * 
 * Functions and classes useful for interpreting configuration files - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <vector>
#include <map>
#include "input.hpp"
#include "atom.hpp"
#include "elements.hpp"

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
class MuDiracInputFile : public InputFile
{
public:
    MuDiracInputFile(void);
    DiracAtom makeAtom();

private:
    map<string, NuclearRadiusModel> nucmodelmap = {{"POINT", POINT},
                                                   {"SPHERE", SPHERE}};
};

#endif