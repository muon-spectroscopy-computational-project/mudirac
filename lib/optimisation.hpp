
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include "atom.hpp"


typedef dlib::matrix<double,0,1> column_vector;

/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib bgfs optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. The starting rms radius for each atom is taken from the
 * default values for the Dirac atom. The starting theta value is set at 0.3 radians.
 *
 *
 * @param config:     config object for MuDirac updated in every iteration
 * @param coord_system: "ct" or "polar" coordinate system the optimisation routine optimises over to configure the nuclear model of each iteration.
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void optimizeFermiParameters(DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time);



/**
 *
 *
 * @brief  A function which finds the derivative with respect to the fermi parameters of the mean square error of MuDirac energies compared to experimental energies for selecte transitions. 
 * @note   This function uses the dlib derivative method to numerically calculate the derivate of the calculateMSE function. The derivative value is returned as a dlib column vector and can be used in minimisation routines.
 *
 * @param m: fermi parameters in ct or polar coordinates
* @param da: Dirac atom used to calculate all the energies.
*
 *
 * @retval derivative: dlib::column_vector (length 2)
 *
 */
column_vector MSE_2pF_derivative(const column_vector &m, DiracAtom & da);


/**
 *
 *
 * @brief  A function which finds the hessian with respect to the fermi parameters of the mean square error of MuDirac energies compared to experimental energies for selecte transitions. 
 * @note   This function uses a central differences finite difference method on the MSE_2pF_derivative to find the hessian with respect to the fermi parameters.
 * The Hessian is left un symmetrised as this yielded better optimisation results. 
 *
 * @param m: fermi parameters in ct or polar coordinates
 * @param coord_system: "ct" or "polar" determines the coordinate system the hessian function will use.
 * @param config:     config object for the initial dirac atom
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
*
 *
 * @retval hessian matrix: dlib::matrix (2x2)
 *
 */
dlib::matrix<double> MSE_2pF_hessian(const column_vector & m, DiracAtom & da);

/**
 * 
 * @brief functor class which allows the means square error of the energies simulated by mudirac compared with experimental energies to be optimised with respect to the fermi parameters.
 * @note This class uses a numerically defined derivative and hessian. This functor exists to be used in dlib::find_min_trust_function.
 */
class opt_2pF_model
{
  /*!
    This object is a "function model" which can be used with the
    find_min_trust_region() routine.  
  !*/

  public:
    // common parameters for the objective, derivative and hessian functions
    typedef ::column_vector column_vector;
    typedef dlib::matrix<double> general_matrix;
    DiracAtom da;

    // constructor
    opt_2pF_model(DiracAtom d_a){
      da = d_a;
    }

    double operator() (
      const column_vector& x
    ) const {
      DiracAtom da_1 = da;
      return da_1.calculateMSE(x(0), x(1));
    }

    // function for the dlib minisation routine to get the derivative and hessian
    void get_derivative_and_hessian (
      const column_vector& x,
      column_vector& der,
      general_matrix & hess
    ) const
    {
      DiracAtom da_1 = da;
      der = MSE_2pF_derivative(x, da_1);
      hess = MSE_2pF_hessian(x, da_1);
    }
};

/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib trust region optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. In polar coordinates, the starting rms radius for each atom is taken from the
 * default values for the Dirac atom. The starting theta value is set at 0.3 radians.
 *@param opt_obj: object defined for dlib optmisation routines which has the experimental values as attributes already defined.
 The operator function of this functor is the objective function for minimisation and the functor also has functions for the derivative and hessian.
 * @param coord_system:   either "ct" or "polar" sets the coordinate system ised by the minisation algorithm.
 * @param config:     config object for MuDirac updated in every iteration
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void optimizeFermiParameters(const opt_2pF_model &opt_obj, DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time);


/**
 *
 *
 * @brief  A function which finds the optimal 2pF parameters for given experimentally measured muonic xray transitions.
 * @note   This function uses a dlib global optimization method to minimise the Mean Square Error between simulated energies
 * calculated using 2pF paarameters and experimentally measured energies. The domain is bounded close the rms radius values near those
 * given by default. The optimization is limited to 10 minutes.
 *
 *
 * @param config:     config object for MuDirac updated in every iteration
 * @param coord_system: "ct" or "polar" coordinate system the optimisation routine optimises over to configure the nuclear model of each iteration.
 * @param da:     Dirac atom with 2pF model updated in every iteration
 * @param transqnums:     transition quantum numbers required to index the transition energies
 * @param xr_lines_measured:      vector of measured muonic transitions
 * @param xr_energies:      vector of measured muonic transition energies
 * @param xr_errors:      vector of measured muonic transition energy errors
 * @param fermi_parameters: structure to contain the optimal conventional and polar fermi parameters as well as the minimised MSE.
 *
 * @retval: None
 *
 */
void globalOptimizeFermiParameters(DiracAtom & da, OptimisationData &fermi_parameters, double & opt_time);

void optFermi2(DiracAtom & da, const string algo, OptimisationData & best_fermi_parameters, double & opt_time);
