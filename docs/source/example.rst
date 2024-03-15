Example of MuDirac Usage
========================
To learn how to use mudirac, let's try a simple example. Open a text editor and write the following:

.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3
    write_spec: T

Save this as :literal:`Au_basic.in` and then pass it to MuDirac. The simulation should be really fast and it should produce the files :literal:`Au_basic.xr.out`, :literal:`Au_basic.log`, :literal:`Au_basic.err` and :literal:`Au_basic.spec.dat`. The :literal:`.log` and :literal:`.err` files are just a log of the program's calculations and a file where any errors are stored; they are not important unless you're trying to figure out what went wrong in a failed calculation. The :literal:`.xr.out` file contains a text summary of the result, and the :literal:`.spec.dat` file contains tabulated data for a simulated spectrum. Let's look at the input file and at what each line does.
::

    element: Au; isotope: 197

This specifies that we're interested in studying gold, specifically the 197-Au isotope.
::

    xr_lines: K1-L2,K1-L3

This indicates which X-ray transitions we want to know about. The notation is the IUPAC standard notation for X-ray spectrometry. These would be the transitions connecting the 1s shell (K1) to the 2p1/2 and 2p3/2 shells (L2, L3). Remember that because these orbitals are relativistic, spin-orbit coupling is built into them, and orbitals with different total angular momentum (orbital + spin) have different energies.
::

   write_spec: T

Finally, this tells to the program to write also a :literal:`.spec.dat` file. Without this line, it wouldn't be created. T here stands for True.

Now open the :literal:`Au_basic.xr.out` file. The contents should look something like this:

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    Line    DeltaE (eV)     W_12 (s^-1)
    K1-L2   1.43693e+07             4.94871e+18
    K1-L3   1.48315e+07             4.63201e+18

The first line is a header that records the context of the calculation - the element's atomic number and atomic mass, and the mass of the particle in atomic units (for a muon this will always be 206.768 au). The next lines show for each line the transition energy in electron volts and the transition rate in 1/seconds, which connects to the relative intensity of the line in the spectrum. Generally speaking, lines with higher transition rates will be stronger, though the connection isn't perfect as there are other factors at play.

Now, this result is achieved with the default settings, that are in fact insufficient to simulate accurately an atom with a large Z like gold; as a general rule, the higher the charge of a nucleus, the more important all the additional terms. Copy :literal:`Au_basic.in` as :literal:`Au.in` and edit it to add lines so that it looks like this:

.. code-block:: bash

    element: Au
    isotope: 197
    xr_lines: K1-L2,K1-L3
    write_spec: T
    nuclear_model: FERMI2
    uehling_correction: T
    electronic_config: Au

This adds three more lines:

:literal:`nuclear_model: FERMI2`: This sets the nucleus to be modelled not as a point charge, but as a Fermi 2-term charge distribution, which is far more accurate to reality. The program contains parameters for this distribution for all isotopes of interest in the periodic table. This will account for the finite size of the nucleus, and the overlap of the muon orbitals with it.

:literal:`uehling_correction: T`: This accounts for the Uehling correction, a quantum field effect relevant to electrostatics at these high energies. It can be understood as accounting for the vacuum itself acting as a polarizable medium; because virtual electron-positron pairs can be generated in quantum field theory, these partially shield the charges and lower the traditional Coulomb force. This is an important term especially for very massive nuclei like Au or Pb and orbitals close to the nucleus.

:literal:`electronic_config: Au`: This term includes approximatively the effect of the other electrons orbiting the nucleus. It does not solve the equations for them, rather it just places them in fixed idealised orbitals and builds a negative charge background from them. The result is an additional correction to the energy, that is however tiny compared to the previous two terms, and often easily ignored.

Try running again MuDirac with this input. The calculation should take longer, and this time the output in :literal:`Au.xr.out` should be:

.. code-block:: bash

    # Z = 79, A = 197 amu, m = 206.768 au
    Line    DeltaE (eV)     W_12 (s^-1)
    K1-L2   5.5936e+06              1.62308e+18
    K1-L3   5.76294e+06             1.76987e+18

Note the significant changes - the energies are almost three times smaller than previously! You can try removing each of the new terms individually, or commenting them out by adding :literal:`#` at the beginning of a line, and re-running to see their effects. Now to familiarize yourself you can try a few more things:

1. try adding more :literal:`xr_lines`, for example L1-M2 and L1-M3;
2. try adding a range of lines; this can be written as K1:M5-K1:M5. It will compute all transitions within the given ranges that obey the selection rules to be allowed;
3. try plotting the spectra in the :literal:`.spec.dat files`, using Gnuplot or importing them in software like Excel or Origin.