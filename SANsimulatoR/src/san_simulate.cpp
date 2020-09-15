#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>

using namespace Rcpp;

/** SAN model simulation function (internal)
 *
 * Simulates the SAN model starting from the initial counts in \code{C_init}
 * for the number(s) of steps specified in \code{steps}. 
 */ 
// [[Rcpp::export]]
List san_timediscrete_c(List C_init,
                        NumericVector p_S, NumericVector p_0, NumericVector p_R,
                        NumericVector p_A, NumericVector p_N, NumericVector p_D, 
                        IntegerVector steps)
{
  /* Extract & validate parameters */
  IntegerVector C_curr_S(clone(as<IntegerVector>(C_init["S"])));
  IntegerVector C_curr_A(clone(as<IntegerVector>(C_init["A"])));
  IntegerVector C_curr_N(clone(as<IntegerVector>(C_init["N"])));
  const int L_ = IntegerVector(C_curr_S).size();
  if ((C_curr_A.size() != L_) || (C_curr_N.size() != L_))
    Rcpp::stop("list S_init must contain initial count vectors S, A, N of the same length");
  if ((p_S.size() != 1) || (p_0.size() != 1) || (p_R.size() != 1) ||
      (p_A.size() != 1) || (p_N.size() != 1) || (p_D.size() != 1))
    Rcpp::stop("p_S, p_0, p_A, p_N, p_D must have length 1");
  const double p_S_ = DoubleVector(p_S)[0];
  const double p_0_ = DoubleVector(p_0)[0];
  const double p_R_ = DoubleVector(p_R)[0];
  const double p_A_ = DoubleVector(p_A)[0];
  const double p_N_ = DoubleVector(p_N)[0];
  const double p_D_ = DoubleVector(p_D)[0];

  /* Prepare result */
  List result;
  
  /* Run simulation over all requested time intervals */
  int* const C_curr_S_ = INTEGER(C_curr_S);
  int* const C_curr_A_ = INTEGER(C_curr_A);
  int* const C_curr_N_ = INTEGER(C_curr_N);
  for(int i=0; i < steps.size(); ++i) {
    /* Within each time intervals, lineages are process in batches of L_B lineages.
     * This facilitates parallelization, and even on a single core one can hope
     * that the overhead of copying to and from the batch buffers is made up for
     * by the more efficient cache use.
     */
    const int L_B = 64;

    #pragma omp parallel
    {
      /* Per-thread output buffers */
      int C_batch_S[L_B], C_batch_A[L_B], C_batch_N[L_B];
      
      /* Process batches in parallel if possible */
      #pragma omp for schedule(static)
      for(int j_min=0; j_min < L_; j_min += L_B) {
        /* Current batch contains the k lineages j_min, ..., j_max-1 */
        const int j_max = std::min(j_min + L_B, L_);
        const int k = j_max - j_min;

        /* Copy current state into thread-local buffer */
        for(int b=0; b < k; ++b) {
          C_batch_S[b] = C_curr_S_[j_min + b];
          C_batch_A[b] = C_curr_A_[j_min + b];
          C_batch_N[b] = C_curr_N_[j_min + b];
        }
        
        /* Simulate until the end of the current interval */
        for(int l=0, l_end=steps[i]; l < l_end; ++l) {
          for(int b=0; b < k; ++b) {
            /* Process b-th lineage in batch with has the index j_min + b */
            
            /* Prevent overflows */
            if ((C_batch_S[b] >= 1000000000) || (C_batch_A[b] >= 1000000000) || (C_batch_N[b] >= 1000000000))
              Rcpp::stop("S, A or N cell count overflowed (>= 1,000,000,000)");
            
            /* Determine the number of cells affected by each type of event */
            const int dS = (p_S_ > 0) ? R::rbinom(C_batch_S[b], p_S_) : 0;
            const int d0 = (p_0_ > 0) ? R::rbinom(C_batch_S[b], p_0_) : 0;
            const int dR = (p_R_ > 0) ? R::rbinom(C_batch_S[b], p_R_) : 0;
            const int dA = (p_A_ > 0) ? R::rbinom(C_batch_S[b], p_A_) : 0;
            const int dN = (p_N_ > 0) ? R::rbinom(C_batch_A[b], p_N_) : 0;
            const int dD = (p_D_ > 0) ? R::rbinom(C_batch_A[b], p_D_) : 0;
            
            /* Update cell counts */
            C_batch_S[b] = std::max(C_batch_S[b] + dS - d0 - dR - dA          , 0);
            C_batch_A[b] = std::max(C_batch_A[b]                + dA     - dD , 0);
            C_batch_N[b] = std::max(C_batch_N[b]           + dR      +dN + dD , 0);
          }
        }
        
        /* Copy final state back */
        for(int b=0; b < k; ++b) {
          C_curr_S_[j_min + b] = C_batch_S[b];
          C_curr_A_[j_min + b] = C_batch_A[b];
          C_curr_N_[j_min + b] = C_batch_N[b];
        }
      }
    } /* parallel end */

    /* Append current per-lineage cell counts to the result list */
    List C_curr = List::create(
      Named("i") = IntegerVector(1, i+1),
      Named("S") = clone(C_curr_S),
      Named("A") = clone(C_curr_A),
      Named("N") = clone(C_curr_N)
    );
    result.push_back(C_curr);
    
    /* Allow interruptions */
    Rcpp::checkUserInterrupt();
  }
  
  /* Return result */
  return result;
}
