Example of MuDirac Usage
========================
Modelling of frequencies and probabilities of transition between energy levels of muonic atoms
-------------------------------------------------------------------------------------------------

To learn how to use MuDirac for this task, let's try a simple example. Open a text editor and write the following:

.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3
    write_spec: T

Save this as :literal:`Au_basic.in` and then pass it to MuDirac. The simulation should be really fast and it should produce the
files :literal:`Au_basic.xr.out`, :literal:`Au_basic.log`, :literal:`Au_basic.err` and :literal:`Au_basic.spec.dat`.
The :literal:`.log` and :literal:`.err` files are just a log of the program's calculations and a file where any errors are 
stored; they are not important unless you're trying to figure out what went wrong in a failed calculation. The :literal:`.xr.out` 
file contains a text summary of the result, and the :literal:`.spec.dat` file contains tabulated data for a simulated spectrum. 
Let's look at the input file and at what each line does.
::

    element: Au 
    isotope: 197

This specifies that we're interested in studying gold, specifically the 197-Au isotope.
::

    xr_lines: K1-L2,K1-L3

This indicates which X-ray transitions we want to know about. The notation is the IUPAC 
standard notation for X-ray spectrometry. These would be the transitions connecting the 1s 
shell (K1) to the 2p1/2 and 2p3/2 shells (L2, L3). Remember that because these orbitals are 
relativistic, spin-orbit coupling is built into them, and orbitals with different total angular 
momentum (orbital + spin) have different energies.
::

   write_spec: T

Finally, this tells to the program to write also a :literal:`.spec.dat` file. Without this line, 
it wouldn't be created. T here stands for True.

Now open the :literal:`Au_basic.xr.out` file. The contents should look something like this:

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    Line    DeltaE (eV)     W_12 (s^-1)
    K1-L2   1.43693e+07             4.94871e+18
    K1-L3   1.48315e+07             4.63201e+18

The first line is a header that records the context of the calculation - the element's atomic 
number and atomic mass, and the mass of the particle in atomic units (for a muon this will always 
be 206.768 au). The next lines show for each line the transition energy in electron volts and the 
transition rate in 1/seconds, which connects to the relative intensity of the line in the spectrum. 
Generally speaking, lines with higher transition rates will be stronger, though the connection isn't 
perfect as there are other factors at play.

Now, this result is achieved with the default settings, that are in fact insufficient to simulate 
accurately an atom with a large Z like gold; as a general rule, the higher the charge of a nucleus, 
the more important all the additional terms. Copy :literal:`Au_basic.in` as :literal:`Au.in` and edit 
it to add lines so that it looks like this:

.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3
    write_spec: T
    nuclear_model: FERMI2
    uehling_correction: T
    electronic_config: Au

This adds three more lines:

:literal:`nuclear_model: FERMI2`: This sets the nucleus to be modelled not as a point charge, but as a Fermi 2-parameter 
(2pF) charge distribution, which is far more accurate to reality. A 2pF charge distribution is defined by a homogeneously 
charged sphere of radius *c* that has a sort of skin around it of thickness *t*, which smooths the charge transition into vacuum. 
MuDirac can automatically access databases containing parameters for this distribution for the isotope being simulated. 
Alternativelly, the *c* and *t* can be manually entered. This nuclear model can account for the finite size of the nucleus, and the 
overlap of the muon orbitals with it.

:literal:`uehling_correction: T`: This accounts for the Uehling correction, a quantum field effect relevant to electrostatics 
at these high energies. It can be understood as accounting for the vacuum itself acting as a polarizable medium; because virtual 
electron-positron pairs can be generated in quantum field theory, these partially shield the charges and lower the traditional 
Coulomb force. This is an important term especially for very massive nuclei like Au or Pb and orbitals close to the nucleus.

:literal:`electronic_config: Au`: This term includes approximatively the effect of the other electrons orbiting the nucleus. 
It does not solve the equations for them, rather it just places them in fixed idealised orbitals and builds a negative charge 
background from them. The result is an additional correction to the energy, that is however tiny compared to the previous two 
terms, and often easily ignored.

Try running again MuDirac with this input. The calculation should take longer, and this time the output in :literal:`Au.xr.out` should be:

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    Line    DeltaE (eV)     W_12 (s^-1)
    K1-L2   5.5936e+06              1.62308e+18
    K1-L3   5.76294e+06             1.76987e+18

Note the significant changes - the energies are almost three times smaller than previously! You can try removing each of the new terms 
individually, or commenting them out by adding :literal:`#` at the beginning of a line, and re-running to see their effects. Now to 
familiarize yourself you can try a few more things:

1. try adding more :literal:`xr_lines`, for example L1-M2 and L1-M3;
2. try adding a range of lines; this can be written as K1:M5-K1:M5. It will compute all transitions within the given ranges that 
   obey the selection rules to be allowed;
3. try plotting the spectra in the :literal:`.spec.dat files`, using Gnuplot or importing them in software like Excel or Origin.

You can also vary the size and shape of the nucleus using the *radius* and *fermi_t* keywords. For example, you can evaluate what 
happens if you increase the radius by 1%, or change the skin thickness parameter in the Fermi model to 2.1 fm 
(rather than the conventional 2.3 fm). Increasing the radius and decreasing the skin thickness should both reduce the transition energies.


.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3
    write_spec: T
    nuclear_model: FERMI2
    uehling_correction: T
    electronic_config: Au
    radius: 7.09
    fermi_t: 2.5


Which changes the output to the following

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    Line    DeltaE (eV)     W_12 (s^-1)
    K1-L2   5.54813e+06             1.60217e+18
    K1-L3   5.71535e+06             1.74712e+18
    # Z = 79, A = 197 amu, m = 206.768 au

The transition energies have indeed shifted down. Even such a minor change has induced a 45 keV shift.

Modelling ground state nuclear properties using muonic x-ray measurements
--------------------------------------------------------------------------

Due to the growing need in the negative muon community, working on nuclear physics, 
for software tools that can be used to estimate the 2pF model's *c* and *t* parameters from experimental 
muonic x-ray measurements and use them to study nuclear propeties, MuDirac has been adapted 
to reverse-engineer the Dirac equation.  It can efficiently obtain the 2pF model parameters *c* and *t* 
using, at least, two experimentally measured muonic x-ray transition energies.

For running the optimisation of the 2pF model's parameters, new keywords have been added to MuDirac's main
input file.  In this example, we will continue working with the 197-Au isotope. Let's define the 
standard MuDirac input file :literal:`Au_nuclear.in` plus an additional input file, :literal:`Au_nuclear.xr.in`,
which will provide the muonic x-rays experimental measurements. Examples of the new 
keywords needed for the main input file are shown in the 197-Au example below:

.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3,L2-M4,L3-M5
    write_spec: T
    nuclear_model: FERMI2
    uehling_correction: T
    electronic_config: Au
    nuclear_model: FERMI2
    fermi_c: 6.6
    fermi_t: 2.0
    optimise_fermi_parameters: True
    rms_radius_decimals: 4

In this input file the optimisation is turned on the by :literal:`optimise_fermi_parameters: T`
keyword. The Fermi model for the nuclei with variable *c* and *t* is chosen by the
:literal:`nuclear_model: FERMI2` parameter. The Fermi parameters provided are used as the
starting point for the optimisation. If they are not provided the optimisation will
start from the default values from nuclear data tables. The :literal:`rms_radius_decimals: 4` 
keyword controls the print precision of the optimised values.  If :literal:`optimise_fermi_parameters: F`, 
the Fermi parameters can be changed manually, like in the example presented in the previous section. 

For :literal:`Au_nuclear.xr.in`, a placeholder text file is being used temporarily to
contain the x-ray measurement data, as there is currently no file standard for this
experimental data. Muonic x-ray energies and their corresponding experimental
uncertainties should be provided in eV. An example :literal:`Au_nuclear.xr.in` file is shown below:

.. code-block:: bash

     # Au 197 input file Au_nuclear.xr.in
     xr_lines:  K1-L2     K1-L3   L2-M4   L3-M5
     xr_energy: 5592800.0 5762500 2474400 2343100.0 
     xr_error:  5000.0    5000.0  2000.0  2500.0

Note that the keyword :literal:`xr_lines` appears in both input files, and its value must
be the same in both files. This is because MuDirac must simulate the same transitions 
contained in the file with the experimental data.

From the command line, MuDirac should be run with both input files as arguments as
follows:

.. code-block:: bash

   $\path\to\mudirac Au_nuclear.in Au_nuclear.xr.in

An additional output file , :literal:`Au_nuclear.fermi_parameters.out` is generated when
:literal:`optimise_fermi_parameters: T`. This file contains the optimised *c* and *t* parameters,
the calculated rms radius, the calculation time and the *c* and *t* parameters in polar coordinates. 
For a more detailed discussion of the meaning of these parameters see:

- L. Liborio *et. al.*, MuDirac 1.3.0: A Sustainable Software Tool for Calculating Ground State Nuclear
  Properties Using Muonic X-Ray Measurements. https://arxiv.org/abs/2605.00554

The result from the 197Au example is shown below:

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    fermi_c fermi_t rms_radius[fm] theta  chi_sq time[s]
    6.6744  1.9265  5.4207         0.2810 0.1560 59.1980