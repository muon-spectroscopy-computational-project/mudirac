# MuDirac
A muonic atom Dirac equation solver.

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

## Running

Simulations can be run simply with the command

   mudirac input.in

where the `.in` file can have any name one prefers. A full list of keywords employable in the `.in` file and their meaning can be found in `docs/Keywords.pdf`.

## Mudirac on Galaxy
Mudirac can be run without any installaiton in Galaxy.
Mudirac is available as a Galaxy tool on [MaterialsGalaxy](https://materialsgalaxy.stfc.ac.uk/) in the `Other Muon Tools` section of the tool panel.

## Citing MuDirac

For the theoretical background on the software and examples of its applications, see the published paper:

*Sturniolo, S, Hillier, A.* Mudirac: A Dirac equation solver for elemental analysis with muonic X‐rays. *X‐Ray Spectrom. ***2020***; 1– 17.* [https://doi.org/10.1002/xrs.3212](https://doi.org/10.1002/xrs.3212)

Cite the above paper if you make use of the software in your work too.

## Contributing

Please follow the muon spectroscopy computational project [contributing guidelines](https://muon-spectroscopy-computational-project.github.io/contributions.html).

## Acknowledgments 

Written with funding from the Ada Lovelace Centre, in collaboration with the ISIS Muon Group.
