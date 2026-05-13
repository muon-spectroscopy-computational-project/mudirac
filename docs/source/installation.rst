Installation of MuDirac
========================

Dependencies
-------------

MuDirac requires:

* A C++ compiler
* CMake 3.25 or later
* `Ceres Solver <http://ceres-solver.org/>`_ (required for 2pF nuclear model optimisation)

Installing Ceres Solver
~~~~~~~~~~~~~~~~~~~~~~~~

Ceres can be installed via a package manager on most platforms:

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install libceres-dev

**macOS (Homebrew):**

.. code-block:: bash

   brew install ceres-solver

**Conda:**

.. code-block:: bash

   conda install -c conda-forge ceres-solver

If you prefer to build Ceres from source, first install its dependencies (including glog), then clone and build:

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install libgoogle-glog-dev libgflags-dev libeigen3-dev

**macOS (Homebrew):**

.. code-block:: bash

   brew install glog gflags eigen

Then build and install Ceres. MuDirac only uses Ceres's dense solver (``DENSE_QR``), so SuiteSparse and LAPACK are not needed and can be disabled for a faster, leaner build:

.. code-block:: bash

   git clone https://github.com/ceres-solver/ceres-solver.git
   cd ceres-solver && git checkout tags/2.2.0
   mkdir build && cd build
   cmake -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DSUITESPARSE=OFF -DLAPACK=OFF ..
   make -j4
   sudo make install

.. note::

   Do **not** use :literal:`-DMINIGLOG=ON` when building Ceres for use with MuDirac. MuDirac uses the `aixlog <https://github.com/badaix/aixlog>`_ logging library, which defines the same severity symbols (``INFO``, ``WARNING``, ``ERROR``, ``FATAL``) as Ceres's miniglog substitute. Using miniglog will cause a compilation error due to this conflict. Always build Ceres with the full glog library.

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

Usage
--------
MuDirac works simply by running it with an input file:

.. code-block:: bash

   mudirac input.in

where :literal:`.in` file can have any name one prefers. The input file is a text file containing rows of the form :literal:`keyword: value`. A full list of keywords employable in the :literal:`.in` file and their meaning can be found in :ref:`section_mudirac_input_keywords`.
