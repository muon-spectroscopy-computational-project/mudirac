/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * config.cpp
 * 
 * Functions and classes useful for interpreting configuration files
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "config.hpp"

MuDiracInputFile::MuDiracInputFile() : InputFile()
{
    // Definition of all input keywords that can be used

    // String keywords
    this->defineStringNode("element", InputNode<string>("H"));                  // Element to compute the spectrum for
    this->defineStringNode("nuclear_model", InputNode<string>("POINT", false)); // Model used for nucleus
    this->defineStringNode("electronic_config", InputNode<string>(""));         // Electronic configuration for background charge
    this->defineStringNode("ideal_atom_minshell", InputNode<string>(""));       // Shell above which to treat the atom as ideal,
                                                                                // and simply use standard hydrogen-like orbitals
    // Boolean keywords
    this->defineBoolNode("uehling_correction", InputNode<bool>(false, false)); // Whether to use the Uehling potential correction
    this->defineBoolNode("write_spec", InputNode<bool>(false, false));         // If true, write a simulated spectrum with the lines found
    this->defineBoolNode("sort_byE", InputNode<bool>(false, false));           // If true, sort output transitions by energy in report
    // Double keywords
    this->defineDoubleNode("mass", InputNode<double>(Physical::m_mu));      // Mass of orbiting particle (default: muon mass)
    this->defineDoubleNode("energy_tol", InputNode<double>(1e-7));          // Tolerance for electronic convergence
    this->defineDoubleNode("energy_damp", InputNode<double>(0.5));          // "Damping" used in steepest descent energy search
    this->defineDoubleNode("max_dE_ratio", InputNode<double>(0.1));         // Maximum |dE|/E ratio in energy search
    this->defineDoubleNode("node_tol", InputNode<double>(1e-6));            // Tolerance parameter used for counting nodes in wavefunctions
    this->defineDoubleNode("loggrid_step", InputNode<double>(0.005));       // Logarithmic grid step
    this->defineDoubleNode("loggrid_center", InputNode<double>(1.0));       // Logarithmic grid center (in units of 1/(Z*m))
    this->defineDoubleNode("uehling_lowcut", InputNode<double>(0.0));       // Low cutoff parameter for Uehling potential (approximation of r ~ 0)
    this->defineDoubleNode("uehling_highcut", InputNode<double>(INFINITY)); // High cutoff parameter for Uehling potential (approximation of r >> 1/2c)
    this->defineDoubleNode("econf_rhoeps", InputNode<double>(1e-4));        // Density threshold at which to truncate the electronic charge background
    this->defineDoubleNode("econf_rin_max", InputNode<double>(-1));         // Upper limit to innermost radius for electronic charge background grid
    this->defineDoubleNode("econf_rout_min", InputNode<double>(-1));        // Lower limit to outermost radius for electronic charge background grid
    this->defineDoubleNode("spec_step", InputNode<double>(1e2));            // Simulated spectrum: energy step (eV)
    this->defineDoubleNode("spec_linewidth", InputNode<double>(1e3));       // Simulated spectrum: width of Gaussian-broadened lines (eV)
    this->defineDoubleNode("spec_expdec", InputNode<double>(-1.0));         // Simulated spectrum: exponential decay factor
                                                                            // (reproduces instrumental sensitivity)
    // Integer keywords
    this->defineIntNode("isotope", InputNode<int>(-1));            // Isotope to use for element
    this->defineIntNode("max_E_iter", InputNode<int>(100));        // Max iterations in energy search
    this->defineIntNode("max_nodes_iter", InputNode<int>(100));    // Max iterations in nodes search
    this->defineIntNode("max_state_iter", InputNode<int>(100));    // Max iterations in state search
    this->defineIntNode("uehling_steps", InputNode<int>(100));     // Uehling correction integration steps
    this->defineIntNode("xr_print_precision", InputNode<int>(-1)); // Number of digits to print out in values in .xr.out file
    this->defineIntNode("verbosity", InputNode<int>(1));           // Verbosity level (1 to 3)
    this->defineIntNode("output", InputNode<int>(1));              // Output level (1 to 3)
    // Vector string keywords
    this->defineStringNode("xr_lines", InputNode<string>(vector<string>{"K1-L2"}, false)); // List of spectral lines to compute

    /* These keywords are reserved for developers and debugging */

    // String keywords
    this->defineStringNode("devel_debug_task", InputNode<string>("")); // Which debugging task to perform

    // Integer keywords
    this->defineIntNode("devel_EdEscan_k", InputNode<int>(-1));      // Value of quantum number k for E->dE scan
    this->defineIntNode("devel_EdEscan_steps", InputNode<int>(100)); // Energy steps for E->dE scan

    // Double keywords
    this->defineDoubleNode("devel_EdEscan_minE", InputNode<double>(-INFINITY)); // Minimum binding energy for E->dE scan
    this->defineDoubleNode("devel_EdEscan_maxE", InputNode<double>(0));         // Maximum binding energy for E->dE scan

    // Boolean keywords
    this->defineBoolNode("devel_EdEscan_log", InputNode<bool>(false, false)); // Make the energy scan logarithmic
}

DiracAtom MuDiracInputFile::makeAtom()
{
    // Now extract the relevant parameters
    int Z = getElementZ(this->getStringValue("element"));
    double m = this->getDoubleValue("mass");
    int A = this->getIntValue("isotope");
    if (A == -1)
    {
        A = getElementMainIsotope(Z);
    }
    if (nucmodelmap.find(this->getStringValue("nuclear_model")) == nucmodelmap.end())
    {
        throw invalid_argument("Invalid nuclear_model parameter in input file");
    }
    NuclearRadiusModel nucmodel = nucmodelmap[this->getStringValue("nuclear_model")];
    double fc = this->getDoubleValue("loggrid_center");
    double dx = this->getDoubleValue("loggrid_step");

    int idshell = -1;
    string idshell_str = this->getStringValue("ideal_atom_minshell");
    if (idshell_str.length() > 1)
    {
        throw invalid_argument("Invalid string for ideal_atom_minshell");
    }
    else if (idshell_str.length() == 1)
    {
        idshell = this->getStringValue("ideal_atom_minshell")[0] - 'J';
    }

    // Prepare the DiracAtom
    DiracAtom da;
    da = DiracAtom(Z, m, A, nucmodel, fc, dx, idshell);
    da.Etol = this->getDoubleValue("energy_tol");
    da.Edamp = this->getDoubleValue("energy_damp");
    da.max_dE_ratio = this->getDoubleValue("max_dE_ratio");
    da.nodetol = this->getDoubleValue("node_tol");
    da.maxit_E = this->getIntValue("max_E_iter");
    da.maxit_nodes = this->getIntValue("max_nodes_iter");
    da.maxit_state = this->getIntValue("max_state_iter");

    if (this->getBoolValue("uehling_correction"))
    {
        da.setUehling(true, this->getIntValue("uehling_steps"),
                      this->getDoubleValue("uehling_lowcut"),
                      this->getDoubleValue("uehling_highcut"));
    }

    if (this->getStringValue("electronic_config") != "")
    {
        double e_mu = effectiveMass(1.0, da.getM() * Physical::amu);
        LOG(TRACE) << "Electronic effective mass: " << e_mu << "\n";
        ElectronicConfiguration econf(this->getStringValue("electronic_config"), da.getZ() - 1, e_mu, true, true);
        da.setElectBkgConfig(true, econf, this->getDoubleValue("econf_rhoeps"),
                             this->getDoubleValue("econf_rin_max"),
                             this->getDoubleValue("econf_rout_min"));
    }

    return da;
}