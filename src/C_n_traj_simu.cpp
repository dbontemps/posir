#include "C_n_traj_simu.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

/**
 * C++ doc is in .h files.
 *
 * Here is only the R doc processed by Roxygen.
 */

/**
 * R doc for Roxygen
 */
//' Log and raise an exception (stop) from C++
//'
//' A message is logged at level ERROR, then stop(message) is called.
//' An optional additional message (details) can be given,
//' to be logged at level DEBUG.
//'
//' @inheritParams log_n_stop
//' @family log and stop from R or C++
//'
//' @examples
//' \dontrun{
//'   tryCatch(log_n_stop_C("essai de log_n_top_C", "via tryCatch()"),
//'            error= function(e){log_info(paste("erreur :",e))},
//'            finally=function(e){log_info(paste("finally :",e))})
//' }
//' @keywords internal
// [[Rcpp::export]]
void log_n_stop_C(std::string message, std::string m_details) {
  Function log_er("log_error"), log_de("log_debug");
  log_er(message);
  if(!m_details.empty()) log_de(m_details);
  stop(message) ;
  return;
}

arma::vec rCenteredPareto_C(int n, double location, double shape,
                            bool symmetrized) {
  arma::vec Z = as<arma::vec>(Rcpp::rexp(n, shape));
  if(symmetrized) {
    Z = (location-1)*exp(Z);
    arma::ivec B = as<arma::ivec>(Rcpp::rbinom(n, 1, .5));
    int i;
    for(i=0; i<n; i++) if(B[i]) Z[i] = -Z[i];
  }
  else {
    Z = location*exp(Z) - (shape*location/(shape-1)) ;
  }
  return Z;
}

/**
 * R doc for Roxygen
 */
//' Random generation using Rcpp if possible, else a call to R
//'
//' Random generation according to a distribution.
//' Calls the [Rcpp][Rcpp::Rcpp-package] function if rDisName is among the list
//' {"rnorm", "runif", "rlnorm", "rCenterdPareto"},
//' else calls the R Function with name rDisName.
//'
//' @param n number of generated i.i.d. random variables.
//' @param rDisName name (string) of the distribution function.
//'
//' @return A pseudo-random vector of length n.
//' @export
//'
//' @examples
//' random_C_or_R(10,"rnorm")
//' myrdistrib=stats::rnorm
//' random_C_or_R(10,"myrdistrib")
// [[Rcpp::export]]
arma::vec random_C_or_R(int n, std::string rDisName) {
  // Voir https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html
  if(rDisName=="rnorm") {
    return as<arma::vec>(Rcpp::rnorm(n));
  }
  if(rDisName=="runif") {
    return as<arma::vec>(Rcpp::runif(n));
  }
  if(rDisName=="rlnorm") {
    return as<arma::vec>(Rcpp::rlnorm(n));
  }
  if(rDisName=="rCenteredPareto") {
    return rCenteredPareto_C(n);
  }
  // default case
  Function aux(rDisName);
  return as<arma::vec>(aux(n)) ;
}

arma::mat RowMult(arma::mat const & M, arma::vec const & v) {
  arma::uword k = M.n_rows, l = M.n_cols, i, j ;
  arma::mat res(k,l);
  if(k!=v.n_elem)
    log_n_stop_C("incompatible dimensions", "in posir/Rcpp RowMult().");
  for(i=0; i<k; i++) for(j=0; j<l; j++) res(i,j)=M(i,j)*v[i];
  return res;
}

arma::mat RowDiv(arma::mat const & M, arma::vec const & v) {
  arma::uword i, l=v.n_elem;
  arma::vec vinv(l) ;
  for(i=0; i<l; i++) vinv(i)=1/v(i);
  return RowMult(M, vinv);
}

arma::mat n_traj_simu_1D_C(int n, int Ndis, arma::ivec Intdgrid,
                           std::string rDisName, bool is_standard) {
  int l=Intdgrid.n_elem, k, i, pos=0 ;
  int q=Intdgrid[l-1] ;
  arma::mat Y(n, l), X(n, Ndis+1) ;
  arma::colvec Ycur(n), Yprec(n) ;
  arma::mat Z = random_C_or_R(n*Ndis, rDisName);
  Z.reshape(n, Ndis);
  if(!is_standard){
    Ycur = stddev(Z,0,1);
    X.submat(0, 1, n-1, Ndis) = cumsum(RowDiv(Z, Ycur), 1);
  }
  else{
    X.submat(0, 1, n-1, Ndis) = cumsum(Z, 1);
  }
  X.submat(0, 0, n-1, 0).zeros();
  Yprec.zeros();
  for(k=Ndis ; k>=q ; k--){
    Ycur=max(abs(X.submat(0, k, n-1, Ndis)-X.submat(0, 0, n-1, Ndis-k))/sqrt(k),1);
    for(i=0;i<n;i++) {
      Ycur[i]=std::max(Ycur[i], Yprec[i]);
    }
    if(k==Intdgrid[pos]) {
      Y.submat(0, pos, n-1, pos) = Ycur;
      pos++;
    }
    Yprec=Ycur;
  }
  return Y ;
}

double max2Dwindow(arma::mat const & cummat, int N, int i, int j) {
  arma::mat dY ;
  dY = abs(cummat.submat(i,j,N,N)
             + cummat.submat(0,0,N-i,N-j)
             - cummat.submat(i,0,N,N-j)
             - cummat.submat(0,j,N-i,N)
             )/sqrt(i*j) ;
  dY.reshape((N+1-i)*(N+1-j),1);
  return as_scalar(max(dY)) ;
}

arma::mat n_traj_simu_2D_C(int n, int Ndis, arma::ivec Intdgrid,
                           std::string rDisName, bool is_standard) {
  int l=Intdgrid.n_elem, k, a, i, j, prev_pos=Ndis, next_pos ;
  double s;
  arma::mat Z, res(n, l), Y(Ndis+1, Ndis+1) ;
  //Function gc("gc");
  for(a=0; a<n; a++) {
    Z.set_size(0,0); // liberate previously used memory.
    Z = random_C_or_R(Ndis*Ndis, rDisName);
    // gc(); // clean R memory after call to rdistrib // trop lent !
    if(!is_standard) {
      s = as_scalar(stddev(Z,0,0)) ;
      Z = Z/s;
    }
    Z.reshape(Ndis, Ndis);
    Y.submat(0, 0, 0, Ndis).zeros();
    Y.submat(1, 0, Ndis, 0).zeros();
    Y.submat(1,1, Ndis, Ndis) = cumsum(cumsum(Z, 0), 1);
    s = std::abs(Y(Ndis,Ndis))/Ndis ;
    for(i=0; i<l; i++) {
      next_pos=Intdgrid[i];
      for(j=prev_pos-1; j>=next_pos; j--){
        for(k=Ndis; k>=next_pos; k--) {
          s = std::max(s, max2Dwindow(Y,Ndis,j,k)) ;
        }
      }
      for(j=Ndis; j>=prev_pos; j--){
        for(k=prev_pos-1; k>=next_pos; k--) {
          s = std::max(s, max2Dwindow(Y,Ndis,j,k)) ;
        }
      }
      res(a,i) = s ;
      prev_pos=next_pos ;
    }
  }
  return res ;
}

/**
 * R doc for Roxygen
 */
//' Random simulation of trajectories of the 1D or 2D POSIR process
//'
//' Simulate n trajectories (as \eqn{\delta} decreases)
//' of the (discretized) 1D or 2D POSIR process.
//'
//' @param Intdgrid decreasing vector of the integers \eqn{1\le k\le Ndis}
//' such that \eqn{k/Ndis} is in the grid for \eqn{\delta}.
//' Generally obtained as Ndis*deltagrid,
//' where the vector dgrid is the desired grid for \eqn{\delta}.
//' (Beware: no checking is done in this internal function.)
//' @inheritParams random_C_or_R
//' @inheritParams n_traj_simu
//' @inheritParams simulationDir
//' @inheritParams check_grid
//' @seealso [run_simu()], [batch_filename()]
//' @family POSIR process generators.
//'
//' @return A matrix with n lines and length(grille) columns.
//' @keywords internal
// [[Rcpp::export]]
arma::mat n_traj_simu_C(int n, int Ndis, arma::ivec Intdgrid,
                        std::string rDisName, bool is_standard, int d){
  int l=Intdgrid.n_elem ;
  if(d==1) return n_traj_simu_1D_C(n, Ndis, Intdgrid, rDisName, is_standard) ;
  if(d==2) return n_traj_simu_2D_C(n, Ndis, Intdgrid, rDisName, is_standard) ;
  else
    log_n_stop_C("unsupported dimensions", "in posir/Rcpp n_traj_simu_C().");
  return arma::mat(n, l) ;
}

