test_stochastic_san <- function(s0, rates, p_errprob) {
  Tmax <- max(rates$Tmax)
  r.exp <- san_deterministic(s0, rates)
  r.sim <- san_stochastic(1e4, s0, rates, p_cutoff=1e-6)
  
  print(r.exp[r.sim, on="t"][, list(
    S.exp=S[1],
    S.sim.avg=mean(i.S),
    S.sim.sd=sd(i.S),
    S.sim.se=sd(i.S) / sqrt(.N),
    S.relerr=mean(i.S)/S[1]
  ), by="t"])
  
  r.exp[r.sim, on="t"][, {
    stopifnot(all(S==S[1]))
    if (sd(i.S) > 0) {
      t <- (mean(i.S) - S[1]) / (sd(i.S) / sqrt(.N))
      2*pt(abs(t), lower.tail=FALSE, df=.N-1)
    }
    else if (abs(mean(i.S) - S[1]) < 1e-6)
      1
    else
      0
  }, by="t"]
}

#test_that('solution accuracy', {
#  test_stochastic_san(1e2, data.table(Tmax=10, r_S=1.3, r_0=0.1, r_R=0.2, r_A=0.3, r_N=1.4, r_D=1.5))
#})


