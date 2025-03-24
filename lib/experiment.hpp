/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * experiment.hpp
 *
 * Functions and classes useful for interpreting experimental results files - header
 * file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "input.hpp"
#include <vector>
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

/**
 * @class ExperimentalResultFile
 *
 * @brief  A specialised InputFile class for the experimental results input file in mudirac
 * @note   A specialised InputFile class for the experimental results input file in mudirac.
 * Initialises all the relevant keywords.
 *
 * @retval None
 */
class ExperimentalResultFile : public InputFile {
 public:
  ExperimentalResultFile(void);
};
