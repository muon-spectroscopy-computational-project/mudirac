Installation of MuDirac
========================

Dependencies
-------------

MuDirac requires:

* A C++ compiler
* CMake 3.25 or later
* An internet connection on the first build (to fetch Ceres and its dependencies automatically)

`Ceres Solver <http://ceres-solver.org/>`_ and its dependencies (Eigen, glog) are fetched and built automatically by CMake the first time you build MuDirac — no separate installation is needed. If Ceres is already installed on your system, CMake will use it directly and skip the download.

Building MuDirac
-----------------

MuDirac uses CMake as a build system. In order to compile it and prepare it to be executed on a Linux, Unix, or MacOS system with a working C++ compiler installed, follow these steps:

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

Install MuDirac using conda
---------------------------
MuDirac is available as `conda <https://docs.conda.io/projects/conda/en/stable/>`_ package.
To install mudirac using conda run the following

.. code-block:: bash

   conda install conda-forge::mudirac

Usage
--------
To simulate the x-ray transition energies of a muonic atom, MuDirac works simply by running it with an input file:

.. code-block:: bash

   mudirac input.in

where :literal:`.in` file can have any name one prefers. The input file is a text file with rows of the form
:literal:`keyword: value`, which contains all the necessary information about the parameters of the Dirac equation, 
and the muonic atom and transitions being simulated. A full list of keywords employable in the :literal:`.in` file and
their meaning can be found in :ref:`section_mudirac_input_keywords`.
