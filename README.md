# mudirac
A muonic atom Dirac equation solver.

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

## Citing MuDirac

For the theoretical background on the software and examples of its applications, see the published paper on the arXiv:

[S. Sturniolo, "mudirac: a Dirac equation solver for elemental analysis with muonic X-rays"](https://arxiv.org/abs/2004.11876)

Cite the above paper if you make use of the software in your work too.

## Acknowledgments 

Written with funding from the Ada Lovelace Centre, in collaboration with the ISIS Muon Group.
