.. mudirac- List of input keywords documentation master file, created by
   sphinx-quickstart on Thu Feb 15 09:03:01 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _section_mudirac_input_keywords:

Input keywords
===========================================================

:literal:`mudirac` takes a single input file, containing multiple lines with the format :literal:`<keyword>: <value>`. Here, we list all the currently available keywords, divided by type, their purpose, and default values.

String keywords
~~~~~~~~~~~~~~~~~
These keywords take a string as value; invalid strings (e.g. a chemical symbol that doesn't correspond to a known element) will give rise to errors.

* :literal:`element`: symbol of the element for the calculation. Determines the nuclear charge. Can be any symbol in the periodic table up to Z=111, Roentgenium (Rg). Default is H.
* :literal:`nuclear_model`: model used to describe the nucleus. Can be POINT (point charge), SPHERE (finite size, uniformly charged spherical nucleus) or FERMI2 (Fermi 2-term charge distribution). Default is POINT.
* :literal:`electronic_config`: electronic configuration to use in order to describe the negative charge background. Can be a full string describing the configuration (e.g. ``1s2 2s2 2p2``), an element symbol to represent the default configuration of that atom when neutral (e.g. ``C``) or a mix of the two (e.g. ``[He] 2s2 2p2``). Default is the empty string (no electrons).
* :literal:`ideal_atom_minshell`: for this shell, and all above it, treat the atom as a simple hydrogen-like point charge Dirac atom, using the known analytical solution and discarding all corrections. Mostly useful for debugging, or when very high shell states have difficulty to converge. The shell must use IUPAC notation (:math:`K \Rightarrow n=1`, :math:`L \Rightarrow n=2`, etc.). Default is the empty string (no ideal solutions used).
* :literal:`xr_lines`: the transition or transitions for which energy and rates are desired. Each line must be expressed using the conventional IUPAC notation [Jenkins et al., 1991]. Multiple lines can be separated by commas. For example:
	
  ::
      
      xr_lines: K1-L2,K1-L3
	
  In addition, colons can be used to indicate ranges of lines. The notation :literal:`K1:L3-M1` would compute the lines K1-M1, L1-M1, L2-M1 and L3-M1. Note that if some of these lines are forbidden by selection rules, they will simply be skipped. A double colon, like :literal:`K1:L3-K1:L3` would loop on both sides, and not count all repeated lines. 

Boolean keywords
~~~~~~~~~~~~~~~~~
These keywords can only have a value of TRUE or FALSE. In order to set them true, either the word 'TRUE' or the letter 'T' (regardless of case) work.

* :literal:`uehling_correction`: whether to turn on the Uehling correction or not. Default is FALSE.
* :literal:`write_spec`:  if true, write a spectrum file using the transition lines found broadened with Gaussian functions. Other :ref:`floating_point_keywords` starting with :literal:`spec_` can then be specified. Default is FALSE.
* :literal:`sort_byE`: if true, print out the transitions sorted by energy instead than by shell. Default is FALSE.

.. _floating_point_keywords:

Floating point keywords
~~~~~~~~~~~~~~~~~~~~~~~~
These keywords accept a non-integer number. It can be written normally (e.g. 105.3) or in scientific notation (e.g. 1.053E2).

* :literal:`mass`:  mass of the particle in atomic units (1 = mass of the electron). By default it's the mass of the muon, 206.7683.
* :literal:`energy_tol`: absolute tolerance for energy convergence when searching for eigenvalues. Iterations will stop once the energy change is smaller than this number, in atomic units. Default is 1E-7.
* :literal:`energy_damp`: a damping parameter used in steepest descent energy search to ease convergence. Used to multiply the suggested step :math:`\delta E` and make it smaller. Helps avoiding overshooting; fine-tuning it might help to converge difficult calculations, while making it bigger might make convergence faster in simple ones. Default is 0.5.
* :literal:`max_dE_ratio`: maximum ratio between energy step, :math:`\delta E`, and current energy :math:`E` in convergence search. If the suggested step exceeds this ratio times the guessed energy, it will be rescaled. This also serves as a measure to avoid overshooting and can be tweaked to get around cases of bad convergence. Default is 0.1.
* :literal:`node_tol`: tolerance parameter used to identify and count nodes in wavefunctions. Very unlikely to need any tweaking. Default is 1E-6.
* :literal:`loggrid_step`: step of the logarithmic grid. Default is 0.005.
* :literal:`loggrid_center`: center of the logarithmic grid in units of :math:`a_0 = 1/(Zm)`. Default is 1.
* :literal:`uehling_lowcut`: low cutoff for Uehling potential, under which the radius will be considered 0. Default is 0.
* :literal:`uehling_highcut`: high cutoff for Uehling potential, over which the radius will be considered :math:`r >> c`. Default is INFINITY.
* :literal:`econf_rhoeps`: charge density threshold under which the electronic charge background will be truncated and treated as zero. Default is 1E-4.
* :literal:`econf_rin_max`: upper limit for the innermost radius of the electronic charge background grid. Default is -1 (no limit).
* :literal:`econf_rout_min`: lower limit for the outermost radius of the electronic charge background grid. Default is -1 (no limit)
* :literal:`spec_step`: energy step for the simulated spectrum, in eV. Only has effect if :literal:`write\_spec = TRUE`. Default is 1E2 eV.
* :literal:`spec_linewidth`: Gaussian broadening width for the simulated spectrum, in eV. Only has effect if :literal:`write\_spec = TRUE`. Default is 1E3 eV.
* :literal:`spec_expdec`: exponential decay parameter :math:`E_{\text{dec}}` for a sensitivity function for the simulated spectrum, in eV. Multiplies the entire spectrum by a function :math:`\exp(-E/E_{\text{dec}})`. Only has effect if :literal:`write\_spec = TRUE`. Default is -1 (no decay).

Integer keywords
~~~~~~~~~~~~~~~~~
Keywords that take an integer number as value.

* :literal:`isotope`: which isotope of the element to consider. Important to determine the mass of the nucleus and its size. Default is -1, which means the most common isotope for the element will be used.
* :literal:`max_E_iter`: maximum number of iterations to perform when searching for the energy of a state. If exceeded, convergence will fail. Increase this value for slow convergences that are however progressing. Default is 100.
* :literal:`max_nodes_iter`: maximum number of iterations to perform when searching for a starting energy value that gives a state the expected number of nodes. If exceeded, convergence will fail. Should generally not need to be adjusted. Default is 100.
* :literal:`max_state_iter`: maximum number of iterations to perform when searching for a state. This loop encloses both node-based and energy-based search. Once a state is converged, the program checks again that it has the correct number of nodes. If it does not, the state is stored for future use and to provide an upper or lower limit to the energy of the searches and then the process is repeated. This number represents how much can the process be repeated before failing. Should not generally need to be adjusted. Default is 100.
* :literal:`uehling_steps`: integration steps for the Uehling potential. Higher numbers will make the Uehling energy more precise but increase computation times. Default is 100.
* :literal:`xr_print_precision`: number of digits after the point to use when printing out energies and transition rates in the :literal:`.xr.out` file. Default is -1 (print as many as possible).
* :literal:`verbosity`: verbosity level. Going from 1 to 3 will increase the amount of information printed to the log file. Default is 1.
* :literal:`output`: output level. Going from 1 to 3 will increase the amount of files produced. Specifically:
   1. will print out only the transition energies and rates in the :literal:`.xr.out` file;
   2. will print out also each of the states in a separate ASCII file as well as the transition matrices for each line;
   3. is reserved for future uses and currently has the same effect as 2.


