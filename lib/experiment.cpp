/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * experiment.cpp
 *
 * Functions and classes useful for interpreting experimental results files
 * and fitting the 2 parameter fermi distribution to the data.
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "experiment.hpp"

ExperimentalResultFile::ExperimentalResultFile() : InputFile() {
  // Vector string keywords
  this->defineStringNode("xr_lines", InputNode<string>(vector<string> {"K1-L2"}, false)); // List of spectral lines experimentally measured
  
  // Vector double keywords
  this->defineDoubleNode("xr_energy", InputNode<double>(vector<double>{0}, false)); // experimentally measured xray energies
  this->defineDoubleNode("xr_error", InputNode<double>(vector<double>{0}, false)); // energy uncertainty for the measured xrays

  /* These keywords are reserved for developers and debugging */

}


pair<double, double> fermiParameters(double rms_radius, double theta){
    // reusing these calculated trig values
    double var_cos = cos(theta);
    double var_sin = sin(theta);
    double uniform_sphere_radius = sqrt(5.0/3.0) * rms_radius;
    double common_constant = uniform_sphere_radius/sqrt((var_cos*var_cos) + (ellipse_const*var_sin*var_sin));
    pair<double, double> fermi_parameters = {common_constant*var_cos, common_constant*var_sin};
    return fermi_parameters;
}

void readXrayMeasurements(ExperimentalResultFile &xr_measurements, bool &read_success, vector<string> &xr_lines, vector<double> &xr_energies, vector<double> &xr_errors, char *argv[]){
  xr_measurements.parseFile(argv[2]);
  xr_lines = xr_measurements.getStringValues("xr_lines");
  LOG(DEBUG) << "Reading experimental Xray measurments for transitions: ";
  for (auto transition: xr_lines){
    LOG(DEBUG) << transition << ", ";
  }
  LOG(DEBUG) << "\n";
  xr_energies = xr_measurements.getDoubleValues("xr_energy");

  LOG(DEBUG) << "Reading experimental Xray energies: ";
  for (auto transition_energy: xr_energies){
    LOG(DEBUG) << transition_energy << ", ";
  }
  LOG(DEBUG) << "\n";

  xr_errors = xr_measurements.getDoubleValues("xr_error");

  LOG(DEBUG) << "Reading experimental Xray energy errors: ";
  for (auto transition_energy_error: xr_errors){
    LOG(DEBUG) << transition_energy_error << ", ";
  }
  LOG(DEBUG) << "\n";
  read_success = true;
}

