/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * constants.hpp
 * 
 * Useful physical constants, atomic units
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

// All values from CODATA 2014 unless otherwise specified

#ifndef PHYSICAL_DEF
namespace Physical {
    const double alpha = 7.2973525664e-3;   // Fine structure constant
    const double c = 137.035999139;         // Speed of light

    // Particle constants
    const double m_e = 1.0;                 // Electron mass
    const double m_mu = 1.0/4.83633170e-3;  // Muon mass
    const double m_p = 1.0/5.44617021352e-4;// Proton mass

    // Units
    const double eV = 1.0/27.211385;        // Electron volt
}
#endif
#define PHYSICAL_DEF