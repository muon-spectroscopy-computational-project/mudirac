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


BaseInputFile::BaseInputFile() : InputFile () {
  // String keywords
  this->defineStringNode("xr_lines", InputNode<string>(vector<string> {"K1-L2"}, false)); // List of spectral lines to compute

}

MuDiracInputFile::MuDiracInputFile() : BaseInputFile() {
  // Definition of all input keywords that can be used

  // String keywords
  this->defineStringNode("element", InputNode<string>("H"));                  // Element to compute the spectrum for
  this->defineStringNode("nuclear_model", InputNode<string>("POINT", false)); // Model used for nucleus
  this->defineStringNode("electronic_config", InputNode<string>(""));         // Electronic configuration for background charge
  this->defineStringNode("ideal_atom_minshell", InputNode<string>(""));       // Shell above which to treat the atom as ideal, and simply use standard hydrogen-like orbitals

  // Boolean keywords
  this->defineBoolNode("uehling_correction", InputNode<bool>(false, false));        // Whether to use the Uehling potential correction
  this->defineBoolNode("write_spec", InputNode<bool>(false, false));                // If true, write a simulated spectrum with the lines found
  this->defineBoolNode("sort_by_energy", InputNode<bool>(false, false));                  // If true, sort output transitions by energy in report
  this->defineBoolNode("optimise_fermi_parameters", InputNode<bool>(false, false)); // If true, perform least squares optimisation using user-provided experimental energies for the fermi model

  // Double keywords
  this->defineDoubleNode("mass", InputNode<double>(Physical::m_mu));      // Mass of orbiting particle (default: muon mass)
  this->defineDoubleNode("radius", InputNode<double>(-1));                // Solid sphere equivalent radius
  this->defineDoubleNode("fermi_t", InputNode<double>(-1));               // Skin thickness for Fermi model
  this->defineDoubleNode("fermi_c", InputNode<double>(-1));               // Half radius parameter c for Fermi-2 distribution
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
  this->defineDoubleNode("spec_expdec", InputNode<double>(-1.0));         // Simulated spectrum: exponential decay factor (reproduces instrumental sensitivity)

  // Integer keywords
  this->defineIntNode("isotope", InputNode<int>(-1));            // Isotope to use for element
  this->defineIntNode("max_E_iter", InputNode<int>(100));        // Max iterations in energy search
  this->defineIntNode("max_nodes_iter", InputNode<int>(100));    // Max iterations in nodes search
  this->defineIntNode("max_state_iter", InputNode<int>(100));    // Max iterations in state search
  this->defineIntNode("uehling_steps", InputNode<int>(100));     // Uehling correction integration steps
  this->defineIntNode("xr_print_precision", InputNode<int>(-1)); // Number of digits to print out in values in .xr.out file
  this->defineIntNode("state_print_precision", InputNode<int>(-1)); // Number of digits to print out in values in Dirac state output .{state_name}.out files
  this->defineIntNode("verbosity", InputNode<int>(1));           // Verbosity level (1 to 3)
  this->defineIntNode("output", InputNode<int>(1));              // Output level (1 to 3)
  // Vector string keywords
  this->defineStringNode("xr_lines", InputNode<string>(vector<string> {"K1-L2"}, false)); // List of spectral lines to compute

  /* These keywords are reserved for developers and debugging */

  // String keywords
  this->defineStringNode("devel_debug_task", InputNode<string>("")); // Which debugging task to perform
  this->defineStringNode("2pF_coords", InputNode<string>("polar", true)); // 2pF optimisation coordinate system
  this->defineStringNode("min_2pF_algorithm", InputNode<string>("bfgs", true)); // 2pF optimisation coordinate system


  // Integer keywords
  this->defineIntNode("devel_EdEscan_k", InputNode<int>(-1));      // Value of quantum number k for E->dE scan
  this->defineIntNode("devel_EdEscan_steps", InputNode<int>(100)); // Energy steps for E->dE scan
  this->defineIntNode("rms_radius_decimals", InputNode<int>(2));      //N decimal places for rms radius iterations in fermi parameter optimisation

  // Double keywords
  this->defineDoubleNode("devel_EdEscan_minE", InputNode<double>(-INFINITY)); // Minimum binding energy for E->dE scan
  this->defineDoubleNode("devel_EdEscan_maxE", InputNode<double>(0));         // Maximum binding energy for E->dE scan

  // Boolean keywords
  this->defineBoolNode("devel_EdEscan_log", InputNode<bool>(false, false)); // Make the energy scan logarithmic
  this->defineBoolNode("reduced_mass", InputNode<bool>(true, false)); // Use the reduced mass
}


vector<TransLineSpec> BaseInputFile::parseXRLines() {
  LOG(DEBUG) << "parsing transition lines into quantum numbers for start and end states in transitions\n";
  // First we unravel the user specified string
  vector<string> xr_lines = getStringValues("xr_lines");
  LOG(DEBUG) << "xr Lines: " << xr_lines.size() << "\n";

  // Convert the user specified strings into quantum numbers for the start
  // and end states in the transition
  vector<TransLineSpec> transqnums;

  for (int i = 0; i < xr_lines.size(); ++i) {
    vector<string> ranges = splitString(xr_lines[i], "-");

    // initialise princial (n) orbital (l) and spin (s) quantum number ranges
    vector<int> n1_range, n2_range, l1_range, l2_range;
    vector<bool> s1_range, s2_range;

    LOG(DEBUG) << "Parsing XR line specification " << xr_lines[i] << "\n";

    // transitions can only be from 1 initial to 1 final state
    if (ranges.size() != 2) {
      LOG(ERROR) << SPECIAL << "Line " << xr_lines[i] << " can not be interpreted properly\n";
      throw invalid_argument("Invalid spectral line in input file");
    }

    vector<int> nr, lr;
    vector<bool> sr;

    parseIupacRange(ranges[0], nr, lr, sr);
    n1_range.insert(n1_range.end(), nr.begin(), nr.end());
    l1_range.insert(l1_range.end(), lr.begin(), lr.end());
    s1_range.insert(s1_range.end(), sr.begin(), sr.end());

    parseIupacRange(ranges[1], nr, lr, sr);
    n2_range.insert(n2_range.end(), nr.begin(), nr.end());
    l2_range.insert(l2_range.end(), lr.begin(), lr.end());
    s2_range.insert(s2_range.end(), sr.begin(), sr.end());

    for (int j = 0; j < n1_range.size(); ++j) {
      for (int k = 0; k < n2_range.size(); ++k) {
        TransLineSpec tnums;
        tnums.n1 = n1_range[j];
        tnums.l1 = l1_range[j];
        tnums.s1 = s1_range[j];

        tnums.n2 = n2_range[k];
        tnums.l2 = l2_range[k];
        tnums.s2 = s2_range[k];

        // validate transitions
        if (tnums.n2 < tnums.n1 || abs(tnums.l2 - tnums.l1) != 1) {
          continue;
        }

        transqnums.push_back(tnums);

        LOG(TRACE) << "Identified transition: " << tnums.n1 << ", " << tnums.l1 << ", " << tnums.s1 << "\t";
        LOG(TRACE) << tnums.n2 << ", " << tnums.l2 << ", " << tnums.s2 << "\n";
      }
    }
  }
  return transqnums;
}

DiracAtom MuDiracInputFile::makeAtom() {
  // Now extract the relevant parameters
  int Z = getElementZ(this->getStringValue("element"));
  double radius = this->getDoubleValue("radius");
  double t = this->getDoubleValue("fermi_t");
  double c_param = this->getDoubleValue("fermi_c");
  double m = this->getDoubleValue("mass");
  int A = this->getIntValue("isotope");
  bool reduced_mass = this->getBoolValue("reduced_mass");

  if (A == -1) {
    A = getElementMainIsotope(Z);
  }
  if (nucmodelmap.find(this->getStringValue("nuclear_model")) == nucmodelmap.end()) {
    throw invalid_argument("Invalid nuclear_model parameter in input file");
  }
  NuclearRadiusModel nucmodel = nucmodelmap[this->getStringValue("nuclear_model")];
  double fc = this->getDoubleValue("loggrid_center");
  double dx = this->getDoubleValue("loggrid_step");

  int idshell = -1;
  string idshell_str = this->getStringValue("ideal_atom_minshell");
  if (idshell_str.length() > 1) {
    throw invalid_argument("Invalid string for ideal_atom_minshell");
  } else if (idshell_str.length() == 1) {
    idshell = this->getStringValue("ideal_atom_minshell")[0] - 'J';
  }

  // Prepare the DiracAtom
  DiracAtom da;
  da = DiracAtom(Z, m, A, nucmodel, radius, fc, dx, idshell, reduced_mass);
  da.Etol = this->getDoubleValue("energy_tol");
  da.Edamp = this->getDoubleValue("energy_damp");
  da.max_dE_ratio = this->getDoubleValue("max_dE_ratio");
  da.nodetol = this->getDoubleValue("node_tol");
  da.maxit_E = this->getIntValue("max_E_iter");
  da.maxit_nodes = this->getIntValue("max_nodes_iter");
  da.maxit_state = this->getIntValue("max_state_iter");

  if (this->getBoolValue("uehling_correction")) {
    da.setUehling(true, this->getIntValue("uehling_steps"),
                  this->getDoubleValue("uehling_lowcut"),
                  this->getDoubleValue("uehling_highcut"));
  }

  if (this->getStringValue("electronic_config") != "") {
    double e_mu = effectiveMass(1.0, da.getM() * Physical::amu);
    LOG(TRACE) << "Electronic effective mass: " << e_mu << "\n";
    ElectronicConfiguration econf(this->getStringValue("electronic_config"), da.getZ() - 1, e_mu, true, true);
    da.setElectBkgConfig(true, econf, this->getDoubleValue("econf_rhoeps"),
                         this->getDoubleValue("econf_rin_max"),
                         this->getDoubleValue("econf_rout_min"));
  }

  if (t ==-1) {
    t = Physical::fermi2_thickness;
  }

  if (c_param == -1) {
    da.setFermi2(t*Physical::fm);
  } else {
    da.setFermi2(t*Physical::fm, c_param*Physical::fm);
  }
  LOG(INFO) << "t = " << da.fermi2.t << "and c = " << da.fermi2.c <<  "\n";

  return da;
}

void MuDiracInputFile::validate(int argc, char *argv[], string & seed) {
  if (argc < 2) {
    cout << "Input file missing\n";
    cout << "Please use the program as `mudirac <input_file>`\n";
    cout << "Quitting...\n";
    exit(-1);
  }

  seed = splitString(argv[1], ".")[0];
  try {
    this->parseFile(argv[1]);
  } catch (runtime_error e) {
    cout << "Invalid configuration file:\n";
    cout << e.what() << "\n";
    exit(-1);
  }
}

void MuDiracInputFile::validateOptimisation(int args, string &coords, string &min_2pF_algo) {
  // check the experimental results input file is provided
  if (args < 3) {
    LOG(ERROR) << "Experimental results input file missing\n";
    LOG(ERROR) << "When optimise_fermi_parameters is True, an additional experimental results file is expected\n";
    LOG(ERROR) << "Please use the program as `mudirac <input_file> <experimental_results_input_file>`\n";
    LOG(ERROR) << "If experimental results cannot be provided, optimise_fermi_parameters should be set to False \n";
    LOG(ERROR) << "Quitting...\n";
    exit(0);
  }
  LOG(DEBUG) << "Minimum file arguments for 2pF optimisation provided \n";

  // check the nuclear model is suitable for optimisation
  if(getStringValue("nuclear_model") != "FERMI2") {
    LOG(ERROR) << "nuclear model parameters can only be optimised for the 2 parameter Fermi model\n";
    LOG(ERROR) << "Please add the line `nuclear_model: FERMI2` to your first input file\n";
    LOG(ERROR) << "Quitting...\n";
    exit(0);
  }
  LOG(DEBUG) << "nuclear model correctly set to \"FERMI2\" \n";

  // check the coordinate system is valid
  coords = this->getStringValue("2pF_coords");

  // get the 2pF optimsation coordinate system (dev)


  if (!((coords == "ct")||(coords == "polar"))) {
    LOG(WARNING)<< "Invalid 2pF coordinate system choice for minimsation\n";
    LOG(WARNING)<< "please use \"ct\" or \"polar\" (default is \"polar\") \n";
    LOG(WARNING)<< "You used: \""<<coords <<"\" \n";
    LOG(INFO) << "Using default polar coordinate system\n";
    coords = "polar";
  }

  if (this->getIntValue("isotope") < 5) {
    LOG(INFO) << "using ct coordinate system as polar parameterisation no longer holds for A < 5 \n";
    coords = "ct";
  }

  LOG(DEBUG) << "optimisation coordinate system valid\n";

  // check the algorithm is valid
  min_2pF_algo = this->getStringValue("min_2pF_algorithm");

  if (!((min_2pF_algo == "bfgs")||(min_2pF_algo == "global")||(min_2pF_algo == "trust"))) {
    LOG(WARNING)<< "Invalid 2pF algorithm for minimsation\n";
    LOG(WARNING)<< "please use \"bfgs\" or \"global\" (default is \"trust\") \n";
    LOG(WARNING)<< "You used: \""<<min_2pF_algo<<"\" \n";
    LOG(INFO) << "Using default optimisation algorithm bfgs\n";
    min_2pF_algo = "bfgs";
  }
  LOG(DEBUG) << "Settings for optimisation are valid\n";
}