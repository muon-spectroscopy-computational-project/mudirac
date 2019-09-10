/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * config.cpp
 * 
 * Functions and classes useful for interpreting configuration files
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "config.hpp"

MuDiracInputFile::MuDiracInputFile() : InputFile()
{
    // Definition of all input keywords that can be used

    // String keywords
    this->defineStringNode("element", InputNode<string>("H"));                     // Element to compute the spectrum for
    this->defineStringNode("nuclear_model", InputNode<string>("POINT", false));    // Model used for nucleus
    this->defineStringNode("uehling_correction", InputNode<string>("OFF", false)); // Whether to use the Uehling potential correction
    // Double keywords
    this->defineDoubleNode("mass", InputNode<double>(1));                   // Mass of orbiting particle (in muon masses)
    this->defineDoubleNode("energy_tol", InputNode<double>(1e-7));          // Tolerance for electronic convergence
    this->defineDoubleNode("energy_damp", InputNode<double>(0.5));          // "Damping" used in steepest descent energy search
    this->defineDoubleNode("max_dE_ratio", InputNode<double>(0.1));         // Maximum |dE|/E ratio in energy search
    this->defineDoubleNode("node_tol", InputNode<double>(1e-6));            // Tolerance parameter used for counting nodes in wavefunctions
    this->defineDoubleNode("loggrid_step", InputNode<double>(0.005));       // Logarithmic grid step
    this->defineDoubleNode("loggrid_center", InputNode<double>(1.0));       // Logarithmic grid center (in units of 1/(Z*m))
    this->defineDoubleNode("uehling_lowcut", InputNode<double>(0.0));       // Low cutoff parameter for Uehling potential (approximation of r ~ 0)
    this->defineDoubleNode("uehling_highcut", InputNode<double>(INFINITY)); // High cutoff parameter for Uehling potential (approximation of r >> 1/2c)
    // Integer keywords
    this->defineIntNode("isotope", InputNode<int>(-1));         // Isotope to use for element
    this->defineIntNode("max_E_iter", InputNode<int>(100));     // Max iterations in energy search
    this->defineIntNode("max_nodes_iter", InputNode<int>(100)); // Max iterations in nodes search
    this->defineIntNode("max_state_iter", InputNode<int>(100)); // Max iterations in state search
    this->defineIntNode("uehling_steps", InputNode<int>(100));  // Uehling correction integration steps
    this->defineIntNode("verbosity", InputNode<int>(1));        // Verbosity level (1 to 3)
    this->defineIntNode("output", InputNode<int>(1));           // Output level (1 to 3)
    // Vector string keywords
    this->defineStringNode("xr_lines", InputNode<string>(vector<string>{"K1-L2"}, false)); // List of spectral lines to compute
}

DiracAtom MuDiracInputFile::makeAtom()
{
    // Now extract the relevant parameters
    int Z = getElementZ(this->getStringValue("element"));
    double m = this->getDoubleValue("mass") * Physical::m_mu;
    int A = this->getIntValue("isotope");
    if (nucmodelmap.find(this->getStringValue("nuclear_model")) == nucmodelmap.end())
    {
        throw invalid_argument("Invalid nuclear_model parameter in input file");
    }
    NuclearRadiusModel nucmodel = nucmodelmap[this->getStringValue("nuclear_model")];
    double fc = this->getDoubleValue("loggrid_center");
    double dx = this->getDoubleValue("loggrid_step");

    // Prepare the DiracAtom
    DiracAtom da = DiracAtom(Z, m, A, nucmodel, fc, dx);
    da.Etol = this->getDoubleValue("energy_tol");
    da.Edamp = this->getDoubleValue("energy_damp");
    da.max_dE_ratio = this->getDoubleValue("max_dE_ratio");
    da.nodetol = this->getDoubleValue("node_tol");
    da.maxit_E = this->getIntValue("max_E_iter");
    da.maxit_nodes = this->getIntValue("max_nodes_iter");
    da.maxit_state = this->getIntValue("max_state_iter");

    if (this->getStringValue("uehling_correction") == "ON")
    {
        da.setUehling(true, this->getIntValue("uehling_steps"),
                      this->getDoubleValue("uehling_lowcut"),
                      this->getDoubleValue("uehling_highcut"));
    }

    return da;
}