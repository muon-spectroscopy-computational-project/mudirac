/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * parameterise.hpp
 *
 * Functions and classes useful parameterising the 
 * 2 parameter fermi domain into polar coordinates - header
 * file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <math.h>
#include <utility>
using namespace std;
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
    double rms_radius_opt;
    double theta_opt;
    double mse;
    double fermi_c;
    double fermi_t;
  };