/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * experiment.cpp
 *
 * Functions and classes useful for interpreting experimental results files
 * and fitting the 2 parameter fermi distribution to the data.
 *
 * @author Milan Kumar
 * @version 1.0 09/05/2025
 */

#include "experiment.hpp"

ExperimentalResultFile::ExperimentalResultFile() : InputFile() {
  // Vector string keywords
  this->defineStringNode("xr_lines", InputNode<string>(vector<string> {""}, false)); // List of spectral lines experimentally measured

  // Vector double keywords
  this->defineDoubleNode("xr_energy", InputNode<double>(vector<double> {0}, false)); // experimentally measured xray energies
  this->defineDoubleNode("xr_error", InputNode<double>(vector<double> {0}, false)); // energy uncertainty for the measured xrays

}


pair<double, double> fermiParameters(double rms_radius, double theta) {
  // reusing these calculated trig values
  double var_cos = cos(theta);
  double var_sin = sin(theta);
  double uniform_sphere_radius = sqrt(5.0/3.0) * rms_radius;
  double common_constant = uniform_sphere_radius/sqrt((var_cos*var_cos) + (ellipse_const*var_sin*var_sin));
  pair<double, double> fermi_parameters = {common_constant*var_cos, common_constant*var_sin};
  return fermi_parameters;
}

double rmsRadius(double fermi_c, double fermi_t){
    double t_const = (7.0*M_PI*M_PI)/(16.0*log(3.0)*log(3.0));
    double r_sq = (3*fermi_c*fermi_c + t_const*fermi_t*fermi_t)*0.2;
    double rms_radius = sqrt(r_sq);
    return rms_radius;
}
