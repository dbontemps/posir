#ifndef C_N_TRAJ_SIMU_H
#define C_N_TRAJ_SIMU_H

#include <string>
#include <sstream>
#include <RcppArmadillo.h>
using namespace Rcpp;

/// Log and raise an exception (stop)
/**
 * The message is logged at level ERROR, then stop(message) is called.
 * An optional additional message (details) can be given,
 * to be logged at level DEBUG.
 *
 * Calls to R functions from package logger.
 *
 * See R version [posir::log_n_stop()] for details.
 *
 * @param message text (string) to be logged and passed to stop().
 * @param m_details optional additional details to be log at DEBUG level.
 *
 */
void log_n_stop_C(std::string message, std::string m_details="");

/// Random generation following a centered or symmetrized Pareto distribution
/**
 * Simulate following a Pareto distribution, then symmetrize it or center it.
 *
 * If the symmetrized parameter is true,
 * * takes n Pareto random variables;
 * * translates them so that their support is the nonnegative real numbers;
 * * then symmetrizes them (using a Bernoulli variable).
 *
 * Else, takes n Pareto random variables and centers them.
 *
 * @param n number of generated i.i.d. random variables.
 * @param location location parameter of the Pareto distribution
 * (before centering or symetrizing), with default value 1.
 * @param shape Shape parameter of the Pareto distribution
 * (before centering or symetrizing), with default value 3.
 * @param symmetrized boolean(s) indicating whether the Pareto variables
 * are to be symmetrized or centered, with default value false (centered).
 * @seealso R documentation [posir::random_C_or_R()]
 *
 * @return A vector of simulated centered or symmetrized Pareto random variables.
 */
arma::vec rCenteredPareto_C(int n,
                            double location=1,
                            double shape=3,
                            bool symmetrized=false);

/// Random generation using Rcpp if possible, else a call to R
/**
 * Random generation according to a distribution.
 * Calls the Rcpp function if rDisName is among the list
 * {"rnorm", "runif", "rlnorm", "rCenterdPareto"},
 * else calls the R Function with name rDisName.
 *
 * @param n number of generated i.i.d. random variables.
 * @param rDisName name of the distribution function.
 * @return A pseudo-random vector of length n.
 */
arma::vec random_C_or_R(int n, std::string rDisName);

// /// Random generation using Rcpp if possible, else a call to R, inline version
// /**
//  * See random_C_or_R() for details.
//  */
// inline arma::vec random_C_or_R_inline(int const & n, std::string const & rDisName) {
//   // Voir https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
//   if(rDisName=="rnorm") {
//     return as<arma::vec>(Rcpp::rnorm(n));
//   }
//   if(rDisName=="runif") {
//     return as<arma::vec>(Rcpp::runif(n));
//   }
//   if(rDisName=="rlnorm") {
//     return as<arma::vec>(Rcpp::rlnorm(n));
//   }
//   if(rDisName=="rCenteredPareto") {
//     return rCenteredPareto_C(n);
//   }
//   // default case
//   Function aux(rDisName);
//   return as<arma::vec>(aux(n)) ;
// }

/// Pointwise multiplication of a matrix by a colvec
/**
 * Multiply the rows of a matrix by the coordinates of a column vector.
 *
 * The number of rows in the matrix has to be the same as the length of the vector.
 *
 * @param M a matrix.
 * @param v a vector.
 * @return A matrix where each row of M has been multiplied
 * by the corresponding coordinate of v.
 */
arma::mat RowMult(arma::mat const & M, arma::vec const & v);

/// Pointwise division of a matrix by a colvec
/**
 * Divide the rows of a matrix by the coordinates of a column vector.
 *
 * The number of rows in the matrix has to be the same as the length of the vector.
 *
 * @param M a matrix.
 * @param v a vector.
 * @return A matrix where each row of M has been divided
 * by the corresponding coordinate of v.
 */
arma::mat RowDiv(arma::mat const & M, arma::vec const & v);

/**
 * @param rDisName name of the distribution function used to simulate
 * the i.i.d. random variables that are summed in the POSIR process.
 */

/// Random simulation of trajectories of the 1D POSIR process
/**
 * See n_traj_simu_C() below for details.
 */
arma::mat n_traj_simu_1D_C(int n, int Ndis, arma::ivec Intdgrid,
                           std::string rDisName, bool is_standard);

/// Random simulation of trajectories of the 2D POSIR process
/**
 * See n_traj_simu_C() below for details.
 */
arma::mat n_traj_simu_2D_C(int n, int Ndis, arma::ivec Intdgrid,
                           std::string rDisName, bool is_standard);

/// Random simulation of trajectories of the 1D or 2D POSIR process
/**
 * @param n number of trajectories to be simulated.
 * @param Ndis discretisation parameter of the POSIR process.
 * @param rDisName name of the distribution function of the
 * i.i.d. random variables that are summed in the POSIR process.
 * @param is_standard whether the POSIR process is to be standardized
 * (in particular if the standard deviation of rdistrib is not 1).
 * @param d dimension parameter of the POSIR process (1 or 2).
 * @seealso R documentation [posir::n_traj_simu()] and [posir::n_traj_simu_C()].
 */
arma::mat n_traj_simu_C(int n, int Ndis, arma::ivec Intdgrid,
                        std::string rDisName, bool is_standard, int d);



#endif // C_N_TRAJ_SIMU_H
