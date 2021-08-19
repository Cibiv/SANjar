message("*** Loading LT47 MCMC results from lt47.mcmcabc.1000chains1000accepts.noreeval.initnorep.minrprop02.maxcorr095.nohistory.rd")
load("data-raw/lt47.mcmcabc.1000chains1000accepts.noreeval.initnorep.minrprop02.maxcorr095.nohistory.rd")

message("*** Synthesizing a SANjar-compatible SANMCMC object")
# The parameter values were taken from lt47.mcmcabc.1000chains1000accepts.noreeval.initnorep.minrprop02.maxcorr095.nohistory.R,
# augmented by the default values from mcmcmabc.R from the SVN repository. Parameter names were then changed to match the
# current mcmc function (e.g. initial.samples -> candidate.samples)
# The likelihood function is set to NULL, including this would require including the whole san_posterior() implemention
# from the SVN repository.
lt47.mcmc <- structure(
  list(
    variables=list(`11S`=c(0, 4), `11R`=c(0, 4), `40S`=c(0, 4), `40A`=c(0, 4), `40N`=c(0, 4), `40D`=c(0, 4)),
    llfun=NULL,
    arguments=list(
      steps=NA, accepts=1000, chains=1000, candidate.samples=10e6,
      initial.samplewithreplacement=FALSE,
      reevaluate.likelihood=FALSE, acceptance.target=0.234,
      initial.Rproposal=1.0, minimal.Rproposal=0.2,
      initial.Vproposal=NULL, minimal.Vproposal.CV=0.1, maximal.Vproposal.corr=0.95,
      proposal.update.rate=0.2, delta.ll.cutoff=15, initial.ll.cutoff=-100),
    initial=LT47.MCMCABC$initial,
    final=LT47.MCMCABC$final),
  class="SANMCMC")

message("*** Saving MCMC results lt47.mcmc to data/lt47.mcmc.RData")
save(lt47.mcmc, file="data/lt47.mcmc.RData")
