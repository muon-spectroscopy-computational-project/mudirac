# mudirac
MuDirac is a simulation software that integrates the Dirac equation for muonic atoms to compute their X-Ray transition energies.

## Contact
leandro.liborio@stfc.ac.uk

## Compiling

To compile, create a directory called `build`; then, within that directory, execute the commands

    cmake ..
    make mudirac

In order to run the test suite, within the same directory, run

    make tests
    make test

and wait for a few seconds for the tests to complete. If you want `mudirac` to be accessible from any folder in your computer, add the resulting `bin` directory to your system `PATH` environment variable.

## Running MuDirac on the command line

Simulations can be run simply with the command

   mudirac input.in

where the `.in` file can have any name one prefers. Detailed documentation for MuDirac can be found [here](https://muon-spectroscopy-computational-project.github.io/mudirac/index.html).  The documentation contains an explanation of the theory behind MuDirac, a full list of keywords and their meaning, and some examples of MuDirac usage.

## How to Contribute to MuDirac

For contributing to the Muon Spectroscopy Computational Project's software tools, please refer to the [contributing guidelines](https://muon-spectroscopy-computational-project.github.io/contributions.html).

## Citing MuDirac

For the theoretical background on the software and examples of its applications, see the published paper:

*Sturniolo, S, Hillier, A.* Mudirac: A Dirac equation solver for elemental analysis with muonic X‐rays. *X‐Ray Spectrom. ***2020***; 1– 17.* [https://doi.org/10.1002/xrs.3212](https://doi.org/10.1002/xrs.3212)

Cite the above paper if you make use of the software in your work too.

## Acknowledgments 

Written with funding from the Ada Lovelace Centre, in collaboration with the ISIS Muon Group.
