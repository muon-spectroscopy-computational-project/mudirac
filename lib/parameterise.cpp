/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * parameterise.cpp
 *
 * Functions and classes useful parameterising the 
 * 2 parameter fermi domain into polar coordinates
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */


#include "parameterise.hpp"

// a constant in the ellipse equation containing fermi_c, fermi_t and uniform radius R
const double ellipse_const = (7*M_PI*M_PI)/(48*log(3)*log(3));


pair<double, double> fermiParameters(double rms_radius, double theta){
    // reusing these calculated trig values
    double var_cos = cos(theta);
    double var_sin = sin(theta);
    double uniform_sphere_radius = sqrt(5.0/3.0) * rms_radius;
    double common_constant = uniform_sphere_radius/sqrt((var_cos*var_cos) + (ellipse_const*var_sin*var_sin));
    pair<double, double> fermi_parameters = {common_constant*var_cos, common_constant*var_sin};
    return fermi_parameters;
}
