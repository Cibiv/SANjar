library(deSolve)

test_deterministic_san <- function(x0, tend, rates, expect.maxerr, expect.regime) {
  times <- seq(from=0, to=tend, by=0.01)
  
  # Solve numerically
  M <- with(rates, matrix(c(r_S-r_0-r_R-r_A,        0, 0,
                            r_A,                 -r_D, 0,
                            r_R,              r_N+r_D, 0), byrow=TRUE, ncol=3))
  r.num <- ode(y=x0, times=times, func=function(t, y, params, ...) {
    list(M %*% y)
  }, method="rk4", parm=list())
  colnames(r.num) <- c("t", "S", "A", "N")
  
  # Evaluate analytical solution
  r.ana <- san_deterministic_eval_fixedrates(x0, times, rates)
  
  # Compute errors
  err <- cbind(t=times,
               S=as.vector(abs(r.num[,'S'] - r.ana[,'S'])),
               A=as.vector(abs(r.num[,'A'] - r.ana[,'A'])),
               N=as.vector(abs(r.num[,'N'] - r.ana[,'N'])))
  max.err <- max(err[,c('S', 'A', 'N')])
  
  # Report
  message("Tested regime ", attr(r.ana, "regime"), ", found max.err ", max.err)
  
  # Check
  expect_equal(attr(r.ana, "regime"), expect.regime)
  expect_lt(max.err, expect.maxerr)
}

test_that('solution accuracy', {
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=0.5, r_R=0.1, r_A=0.3, r_N=0.8, r_D=0.4),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_fin.general")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=1.5, r_R=0.1, r_A=0.3, r_N=0.8, r_D=0.4),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_fin.general")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=0.5, r_R=0.6, r_A=0.3, r_N=0.8, r_D=0.4),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_fin.g_eq_mrD")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=1.5, r_R=0.6, r_A=0.3, r_N=0.8, r_D=1.4),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_fin.g_eq_mrD")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=0.5, r_R=0.1, r_A=0.3, r_N=0.8, r_D=0),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_inf")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=1.5, r_R=0.1, r_A=0.3, r_N=0.8, r_D=0),
                         expect.maxerr=1e-6, expect.regime="s_noneq.a_inf")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=0.5, r_R=0.1, r_A=0.4, r_N=0.8, r_D=0.3),
                         expect.maxerr=1e-6, expect.regime="s_equil.a_fin")
  test_deterministic_san(c(S=1,A=2,N=3), tend=100, rates=list(r_S=1, r_0=0.5, r_R=0.1, r_A=0.4, r_N=0.8, r_D=0),
                         expect.maxerr=1e-6, expect.regime="s_equil.a_inf")
})
