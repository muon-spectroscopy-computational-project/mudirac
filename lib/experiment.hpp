/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * experiment.hpp
 *
 * Functions and classes useful for interpreting experimental results files - header
 * file
 *
 * @author Milan Kumar
 * @version 1.0 09/05/2025
 */

#include "input.hpp"
#include <vector>
#include "../vendor/aixlog/aixlog.hpp"
#include <math.h>
#include <utility>
#include "config.hpp"
using namespace std;


/**
 * @class ExperimentalResultFile
 *
 * @brief  A specialised BaseInputFile class for the experimental results input file in mudirac
 * @note   A specialised BaseInputFile class for the experimental results input file in mudirac.
 * Initialises all the relevant keywords.
 *
 * @retval None
 */
class ExperimentalResultFile : public BaseInputFile {
 public:
  ExperimentalResultFile(void);

  /**
   * @brief: validates experminental xray input file
   * @note: validates the XR input file for correct keywords, types and
   * completeness. Exits the program if invalid.
   * @param xr_infile: name of the file containing the experimental results
   * @retval: none
   */
  void validate(const string xr_infile);
};



