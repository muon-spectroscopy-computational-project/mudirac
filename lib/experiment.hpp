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
#include <math.h>
#include <utility>

using namespace std;

#ifndef OPT_DATA
#define OPT_DATA

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


// a constant in the ellipse equation containing fermi_c, fermi_t and uniform radius R
const double ellipse_const = (7*M_PI*M_PI)/(48*log(3)*log(3));

/**
 * 
 *
 * @brief  A function to transform a parameterisation from polar coordinates to fermi parameters
 * @note   A function which takes the polar coordinates rms_radius, theta and calculates the corresponding
 * parameters for the 2 parameter fermi distribution function.
 * 
 *
 * @retval a pair of doubles fermi_c and fermi_t
 */
pair<double, double> fermiParameters(double rms_radius, double theta);



struct OptimisationData {
    double rms_radius;
    double theta;
    double mse;
    double fermi_c;
    double fermi_t;
  };

#endif