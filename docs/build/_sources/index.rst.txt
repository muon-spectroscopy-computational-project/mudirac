.. mudirac- homepage master file

Welcome to mudirac's documentation
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   mudirac- List of input keywords's documentation
   ...

mudirac is a muonic atom Dirac equation solver.

Compiling
==========
To compile, create a directory called :literal:`mudirac`; then, within that directory, execute the commands

.. code-block:: bash

   cmake ..
   make mudirac

In order to run the test suite. within the same directory, run

.. code-block:: bash

   make tests
   make test

and wait for a few seconds for the tests to complete. If you want :literal:`mudirac` to be accessible from any folder in your computer, add the resulting :literal:`bin` directory to your system :literal:`PATH` environment variable.

Running
========
Simulations can be run simply with the command

.. code-block:: bash

   mudirac input.in

where :literal:`.in` file can have any name one prefers. A full list of keywords employable in the :literal:`.in` file and their meaning can be found in :ref:`section_mudirac_input_keywords`.

Citing mudirac
===============
For the theoretical background on the software and examples of its applications, see the published paper:

Sturniolo, S, Hillier, A. Mudirac: A Dirac equation solver for elemental analysis with muonic X‐rays. X‐Ray Spectrom. 2020; 1– 17. https://doi.org/10.1002/xrs.3212

Cite the above paper if you make use of the software in your work too.

Contact
========
leandro.liborio@stfc.ac.uk

Acknowledgments
=================
Written with funding from the Ada Lovelace Centre, in collaboration with the ISIS Muon Group.

