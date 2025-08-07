/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * mudirac.hpp
 *
 * Main binary of the program (header)
 *
 * @author Simone Sturniolo
 * @version 1.0.1 20/03/2020
 */

#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <dlib/optimization.h>
#include "../lib/optimisation.hpp"
#include "../lib/config.hpp"
#include "../lib/atom.hpp"
#include "../lib/output.hpp"
#include "../lib/utils.hpp"
#include "../lib/elements.hpp"
#include "../lib/constants.hpp"
#include "../lib/debugtasks.hpp"

#include "../vendor/aixlog/aixlog.hpp"

#define PROJECT_VERSION "1.0.1"

using namespace std;


/**
 * @brief  Sets up aixlog logging
 * @note    sets up aixlogging with the provided verbosity level.
 * Logs are saved to a file name based on the seed. The initial log message is
 *  also output to the file.
 *
 * @param  verbosity:  verbosity level from the config file
 * @param  seed:       seed of the output log file name.
 * @retval None
 */
void setupLogging(const int verbosity, const string & seed);

/**
 * @brief  Runs debug tasks if requested
 * @note   runs debug tasks when present in input config file.
 *
 * @param  config:  config object used to make a dirac atom and run an EdEscan
 * @retval None
 */
void debugTasks(MuDiracInputFile & config);

/**
 * @brief  Writes the potential from electronic cofiguration to file.
 * @note   finds the electric potential in a Dirac atom caused by the
 * electronic configuration. The result is output to a file name based
 * on the seed with verbosity, out_verbosity.
 *
 * @param  da:  dirac atom containing the electronic potential.
 * @param seed: seed name for the output file.
 * @param out_verbosity: verbosity level for the output file.
 * @retval None
 */
void writeEConfV(DiracAtom & da, const string & seed, const int & out_verbosity);

