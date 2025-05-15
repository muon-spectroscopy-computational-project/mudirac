/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * experiment.cpp
 *
 * Functions and classes useful for interpreting experimental results files
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
