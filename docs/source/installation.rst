Installation of MuDirac
========================
MuDirac uses CMake as a build system, and requires a C++ compiler. In order to compile it and prepare it to be executed on a Linux, Unix, or MacOS system with a working C++ compiler installed, follow these steps:

#. Download and unpack (or :literal:`git clone`) the repository on your local system.
#. Within the main folder of the repository (the one containing the :literal:`READ.md` file), create a subfolder called :literal:`build`.
#. Within the :literal:`build` folder, run the following commands:

.. code-block:: bash

   cmake ..
   make mudirac

In order to run the test suite, within the same directory run:

.. code-block:: bash

   make tests
   make test

and wait for a few seconds for the tests to complete. If you want :literal:`mudirac` to be accessible from any folder in your computer, add the resulting :literal:`bin` directory to your system :literal:`PATH` environment variable.

Usage
--------
MuDirac works simply by running it with an input file:

.. code-block:: bash

   mudirac input.in

where :literal:`.in` file can have any name one prefers. The input file is a text file containing rows of the form :literal:`keyword: value`. A full list of keywords employable in the :literal:`.in` file and their meaning can be found in :ref:`section_mudirac_input_keywords`.
