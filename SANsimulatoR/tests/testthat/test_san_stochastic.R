test_stochastic_san <- function(M, s0, rates, p_errprob) {
  # Simulate model and compute expectations using the deterministic model
  Tmax <- max(rates$Tmax)
  r.exp <- san_deterministic(s0, rates)
  r.sim <- san_stochastic(M, s0, rates, p_cutoff=1e-6)

  # Simple one-sample t-test
  t.test.onesample <- function(x0, x.avg, x.sd, nsample) {
    if (x.sd > 0) {
      t <- (x.avg - x0) / (x.sd / sqrt(nsample))
      2*pt(abs(t), lower.tail=FALSE, df=nsample-1)
    }
    else if (abs(x.avg - x0) < 1/M)
      1
    else
      0
  } 

  # Join simulation results to expected value at every time point,
  # and compute sample mean and standard deviation
  r <- r.exp[r.sim, on="t"][, list(
    t=t,
    S.sim=i.S, A.sim=i.A, N.sim=i.N,
    S.exp=S, A.exp=A, N.exp=N
  )][, list(
    S.exp=S.exp[1], S.sim.avg=mean(S.sim), S.sim.sd=sd(S.sim),
    A.exp=A.exp[1], A.sim.avg=mean(A.sim), A.sim.sd=sd(A.sim),
    N.exp=N.exp[1], N.sim.avg=mean(N.sim), N.sim.sd=sd(N.sim),
    nsample=.N
  ), by="t"]

  # Pperform a one-sample t-test to compare the sample mean to
  # the expected value
  p <- r[, list(
    S.pval=t.test.onesample(S.exp, S.sim.avg, S.sim.sd, nsample),
    A.pval=t.test.onesample(A.exp, A.sim.avg, A.sim.sd, nsample),
    N.pval=t.test.onesample(N.exp, N.sim.avg, N.sim.sd, nsample)
  ), by="t"]

  # Shouldn't get too small p-values
  expect_gt(min(p$S.pval, p$A.pval, p$N.pval), 1e-3)
}

test_that('solution accuracy_equilibrium', {
  test_stochastic_san(M=1e3, s=1e3, data.table(Tmax=10, r_S=1, r_0=0.2, r_R=0.3, r_A=0.5, r_N=0.15, r_D=0.25), p_errprob=1e-3)
})

test_that('solution accuracy_growing', {
  test_stochastic_san(M=1e3, s=1e3, data.table(Tmax=10, r_S=1.5, r_0=0.2, r_R=0.3, r_A=0.5, r_N=0.15, r_D=0.25), p_errprob=1e-3)
})

test_that('solution accuracy_shrinking', {
  test_stochastic_san(M=1e3, s=1e3, data.table(Tmax=10, r_S=0.5, r_0=0.2, r_R=0.3, r_A=0.5, r_N=0.15, r_D=0.25), p_errprob=1e-3)
})
