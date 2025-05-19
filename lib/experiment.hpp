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


/**
 *
 *
 * @brief  A function to read the experimental measurements input file
 * and store the data in referenced variables
 * @note   A function which reads the experimental measurements input file
 * into the referenced xr_measurements object. The data of xr_lines, energies errors
 * are copied into referenced vector variables. if successful the references read_success
 * variable is true so the main function can continue with 2 parameter fermi optimisation.
 *
 * @param xr_measurements:    object to contain data from xr input file
 * @param read_success:     switch to control whether 2pF optimisation should occur
 * @param xr_lines:     list of measured xray transitions
 * @param xr_energies:      list of measured xray energies
 * @param xr_errors:      list of measured xray errors
 * @param argv:      MuDirac argument containing the experimental input file
 *
 */
void readXrayMeasurements(ExperimentalResultFile &xr_measurements, bool &read_success, vector<string> &xr_lines, vector<double> &xr_energies, vector<double> &xr_errors, char *argv[]);


#endif