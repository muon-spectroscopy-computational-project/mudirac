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

ExperimentalResultFile::ExperimentalResultFile() : BaseInputFile() {
  // Vector string keywords
  this->defineStringNode("xr_lines", InputNode<string>(vector<string> {""}, false)); // List of spectral lines experimentally measured

  // Vector double keywords
  this->defineDoubleNode("xr_energy", InputNode<double>(vector<double> {0}, false)); // experimentally measured xray energies
  this->defineDoubleNode("xr_error", InputNode<double>(vector<double> {0}, false)); // energy uncertainty for the measured xrays

}

void ExperimentalResultFile::validate() {
  // read the measured transition lines
  vector<string> xr_lines_measured = getStringValues("xr_lines");
  LOG(DEBUG) << "Reading experimental Xray measurments for transitions: ";
  for (auto transition: xr_lines_measured) {
    LOG(DEBUG) << transition << ", ";
  }
  LOG(DEBUG) << "\n";

  // read the measured transition energies
  vector<double> xr_energies = getDoubleValues("xr_energy");
  LOG(DEBUG) << "Reading experimental Xray energies: ";
  for (auto transition_energy: xr_energies) {
    LOG(DEBUG) << transition_energy << ", ";
  }
  LOG(DEBUG) << "\n";

  // read the measured transition errors
  vector<double> xr_errors = getDoubleValues("xr_error");
  LOG(DEBUG) << "Reading experimental Xray energy errors: ";
  for (auto transition_energy_error: xr_errors) {
    LOG(DEBUG) << transition_energy_error << ", ";
  }
  LOG(DEBUG) << "\n";

  // checking that the file has contents and not the default values
  LOG(DEBUG) << "Validating experimental results input \n";
  if (xr_lines_measured[0] == "") {
    cout << "Experimental results input file is empty\n";
    cout << "Please check the filename of the experimental results input file \n";
    cout << "Quitting...\n";
    exit(0);
  }

  // check that the data provided is complete: all transitions measured have energies and errors
  if (xr_lines_measured.size() == xr_energies.size() && xr_energies.size() == xr_errors.size()) {
    LOG(DEBUG) << "Experimental results valid \n";
  } else {
    cout << "Invalid experimental measurements file: Missing input values\n";
    cout << "please check energies and errors are listed for each xray transition line \n";
    cout << "Quitting...\n";
    exit(0);
  }
}
