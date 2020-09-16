#include <algorithm>
#include <exception>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>

// [[Rcpp::depends(dqrng, BH) ]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <R_randgen.h>
#include <xoshiro.h>
#include <boost/random/binomial_distribution.hpp>

#ifdef _OPENMP
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#endif

using namespace Rcpp;

/** The RNG we use.
 *
 * The R RNG cannot be used from multiple threads, so we use one of the RNGs
 * provided by the dqrng package instead. xoroshiro128plus should be the fastest
 * one. Other RNGs have slightly better statistical properties, but given that
 * our requirements are low, this one should suffice.
 */
typedef dqrng::xoroshiro128plus RNGType;

/** Thread-save version of `R::rbinom`
 *
 */
int rbinom_threadsafe(RNGType& rng, int n, double p) {
  return boost::random::binomial_distribution<int>(n, p)(rng);
}

/** SAN model simulation function (internal)
 *
 * Simulates the SAN model starting from the initial counts in \code{C_init}
 * for the number(s) of steps specified in \code{steps}.
 */
// [[Rcpp::export]]
List san_timediscrete_c(const List C_init,
                        const NumericVector p_S, const NumericVector p_0, const NumericVector p_R,
                        const NumericVector p_A, const NumericVector p_N, const NumericVector p_D,
                        const IntegerVector steps)
{
  /* Declared *first* to ensure it's destroyed last. Once the result object is gone,
   * result is unprotected from the GC, so no allocations must occur that could trigger
   * a GC run until the caller has protected it again.
   */
  List result;

  /* Extract & validate parameters */
  const IntegerVector C_init_S(as<IntegerVector>(C_init["S"]));
  const IntegerVector C_init_A(as<IntegerVector>(C_init["A"]));
  const IntegerVector C_init_N(as<IntegerVector>(C_init["N"]));
  const int L = IntegerVector(C_init_S).size();
  if ((C_init_A.size() != L) || (C_init_N.size() != L))
    Rcpp::stop("list C_init must contain initial count vectors S, A, N of the same length");
  const int* const C_init_S_ = INTEGER(C_init_S);
  const int* const C_init_A_ = INTEGER(C_init_A);
  const int* const C_init_N_ = INTEGER(C_init_N);

  if ((p_S.size() != 1) || (p_0.size() != 1) || (p_R.size() != 1) ||
      (p_A.size() != 1) || (p_N.size() != 1) || (p_D.size() != 1))
    Rcpp::stop("p_S, p_0, p_A, p_N, p_D must have length 1");
  const double p_S_ = DoubleVector(p_S)[0];
  const double p_0_ = DoubleVector(p_0)[0];
  const double p_R_ = DoubleVector(p_R)[0];
  const double p_A_ = DoubleVector(p_A)[0];
  const double p_N_ = DoubleVector(p_N)[0];
  const double p_D_ = DoubleVector(p_D)[0];

  if (any(steps < 0).is_true())
    stop("steps vector must be non-negative");
  const int T = steps.size();
  const int* const steps_ = INTEGER(steps);

  /* Thead-local RNG template, seeded from R's RNG */
  const RNGType RNG_tpl(((uint64_t)dqrng::R_random_u32() << 32) | (uint64_t)dqrng::R_random_u32());

  /* Prepare result and pointers to the per-interval S, A and N vectors */
  std::vector<int*> result_S_, result_A_, result_N_;
  for(int i=0; i < T; ++i) {
    List ri = List::create(
      Named("i") = IntegerVector(1, i+1),
      Named("S") = IntegerVector(L, NA_INTEGER),
      Named("A") = IntegerVector(L, NA_INTEGER),
      Named("N") = IntegerVector(L, NA_INTEGER)
    );
    result_S_.emplace_back(INTEGER(ri["S"]));
    result_A_.emplace_back(INTEGER(ri["A"]));
    result_N_.emplace_back(INTEGER(ri["N"]));
    result.push_back(ri);
  }

  /* Parallelisation happens on the level of lineages in batches of size L_B */
  static const int L_B = 256;

#if _OPENMP
  /* Stores the first error that occurs in any of the parallel tasks */
  std::exception_ptr caught_exception;
#endif

  #pragma omp parallel default(none)        \
    shared(result_S_, result_A_, result_N_, caught_exception) \
    firstprivate(C_init_S_, C_init_A_, C_init_N_, L, steps_, T, \
                 p_S_, p_0_, p_R_, p_A_, p_N_, p_D_, RNG_tpl)
  {
      /* Per-thread RNG */
      RNGType RNG(RNG_tpl);

#if _OPENMP
      /* Use thread-specific sub-sequence of the RNG */
      RNG.long_jump(omp_get_thread_num());
      const bool is_master_thread = (omp_get_thread_num() == 0);
#else
      const bool is_master_thread = true;
#endif

    /* Per-thread buffers */
    int C_batch_S[L_B], C_batch_A[L_B], C_batch_N[L_B];

    /* Process batches of lineages in parallel if possible */
    #pragma omp for schedule(static)
    for(int j_min=0; j_min < L; j_min += L_B) {
      try {
        /* Current batch contains the k lineages j_min, ..., j_max-1 */
        const int j_max = std::min(j_min + L_B, L);
        const int k = j_max - j_min;

        /* Copy initial state into thread-local buffer */
        for(int b=0; b < k; ++b) {
          C_batch_S[b] = C_init_S_[j_min + b];
          C_batch_A[b] = C_init_A_[j_min + b];
          C_batch_N[b] = C_init_N_[j_min + b];
        }

        /* Run simulation over all requested time intervals */
        for(int i=0; i < T; ++i) {
          /* Simulate until the end of the current interval */
          for(int l=0, l_end = steps_[i]; l < l_end; ++l) {
            for(int b=0; b < k; ++b) {
              /* Process b-th lineage in batch with has the index j_min + b */

              /* Prevent overflows */
              if ((C_batch_S[b] >= 1000000000) || (C_batch_A[b] >= 1000000000) || (C_batch_N[b] >= 1000000000))
                Rcpp::stop("S, A or N cell count overflowed (>= 1,000,000,000)");

              /* Determine the number of cells affected by each type of event */
              const int dS = (p_S_ > 0) ? rbinom_threadsafe(RNG, C_batch_S[b], p_S_) : 0;
              const int d0 = (p_0_ > 0) ? rbinom_threadsafe(RNG, C_batch_S[b], p_0_) : 0;
              const int dR = (p_R_ > 0) ? rbinom_threadsafe(RNG, C_batch_S[b], p_R_) : 0;
              const int dA = (p_A_ > 0) ? rbinom_threadsafe(RNG, C_batch_S[b], p_A_) : 0;
              const int dN = (p_N_ > 0) ? rbinom_threadsafe(RNG, C_batch_A[b], p_N_) : 0;
              const int dD = (p_D_ > 0) ? rbinom_threadsafe(RNG, C_batch_A[b], p_D_) : 0;

              /* Update cell counts */
              C_batch_S[b] = std::max(C_batch_S[b] + dS - d0 - dR - dA          , 0);
              C_batch_A[b] = std::max(C_batch_A[b]                + dA     - dD , 0);
              C_batch_N[b] = std::max(C_batch_N[b]           + dR      +dN + dD , 0);
            }

            /* Check if there was a cancel request */
            #pragma omp cancellation point for
          }

          /* Copy result from thread-local buffer to result vectors */
          for(int b=0; b < k; ++b) {
            result_S_[i][j_min + b] = C_batch_S[b];
            result_A_[i][j_min + b] = C_batch_A[b];
            result_N_[i][j_min + b] = C_batch_N[b];
          }

          /* Allow the simulation to be interrupted. This will throw in case of an interrupt */
          if (is_master_thread)
            checkUserInterrupt();
        } /* for(int i... */
      } /* try ... */
      catch (...) {
#if _OPENMP
        /* Save first error, cancel parallel tasks and re-throw on main thread */
        #pragma omp critical (caught_exception)
        {
          if (!caught_exception)
            caught_exception = std::current_exception();
        }
        #pragma omp cancel for
#else
        throw;
#endif
      } /* try ... catch ... */
    } /* parallel for(int j_min... */
  }  /* omp parallel */

#if _OPENMP
  /* Throw possible pending error now that the parallel section has ended */
  if (caught_exception)
    std::rethrow_exception(caught_exception);
#endif

  return result;
}
