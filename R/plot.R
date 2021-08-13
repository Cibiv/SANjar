#' Plot log-normal confidence intervals defining the posterior overlayed with model predictions
#'
#' @export
plot.SANPosterior <- function(posterior, params=NULL, llout=NULL) {
  if (is.null(llout)) {
    if (is.null(params))
      stop("either a parameter vector or the result of the posterior distribution loglikelihood function must be specified")
    message("llout not specified, evaluating model for parameters ",
            paste0(names(params), "=", signif(params, 3), collapse=", "))
    llmeta <- p$loglikelihood(params)
  }
  
  # Convert output of loglikelihood function into a data.table  
  llout.tab <- as.data.table(llout)
  llout.tab[, index := 1:.N ]
  
  # Nicer log-scales
  make_scale_log10 <- function(scale, ...) scale(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    oob=scales::oob_keep,
    ...
  )
  
  # Transform used for likelihood color scale
  log_ll <- scales::trans_new("log_ll", function(x) log10(-x), function(y) -10^y)
  midpoint <- function(x) sum(range(log_ll$transform(x)))/2
  
  # Extract per-day cell counts as a table suitable for plotting
  llout.cc <- melt(llout.tab, id.vars=c("index", "ll_cc"),
                   measure.vars=paste0("CC", posterior$parametrization$cc_days))[, day := {
                     as.integer(sub("^CC([0-9]+)$", "\\1", variable))
                   }]
  
  # Plot total organoid sizes
  p.cc <- ggplot2::ggplot() +
    ggplot2::geom_errorbar(data=posterior$data$cc, aes(x=day, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
    ggplot2::geom_line(data=llout.cc, aes(x=day, y=pmax(value, posterior$arguments$min.size),
                                          color=ll_cc, group=index)) +
    make_scale_log10(scale_y_log10) +
    ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                   low="red", mid="green", midpoint=midpoint(llout.cc$ll_cc), high="blue") +
    ggplot2::labs(y="cells")
  
  # Plot rank-sizes
  p.rs <- lapply(posterior$parametrization$rs_days, function(rs_day) {
    # Extract per-day cell counts as a table suitable for plotting
    llout.rs <- melt(llout.tab, id.vars=c("index", paste0("ll_rs_", rs_day)),
                     measure.vars=paste0("RS", rs_day, "_", posterior$parametrization$rs_ranks))
    llout.rs <- llout.rs[, c("ll_rs", "day", "rank") := {
      list(eval(as.name(paste0("ll_rs_", rs_day))),
           rs_day,
           as.integer(sub(paste0("^RS", rs_day, "_([0-9]+)$"), "\\1", variable)))
    }]
    
    # Plot
    ggplot() +
      ggplot2::geom_errorbar(data=posterior$data$rs[day==rs_day],
                             aes(x=rank, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
      ggplot2::geom_line(data=llout.rs, aes(x=rank, y=pmax(value, posterior$arguments$min.size),
                                            color=ll_rs, group=index)) +
      annotate("text", x=Inf, y=Inf, hjust=1, vjust=1, col="black", label=paste0("day ", rs_day)) +
      make_scale_log10(scale_x_log10) +
      make_scale_log10(scale_y_log10) +
      ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                     low="red", mid="green", midpoint=midpoint(llout.rs$ll_rs), high="blue") +
      ggplot2::labs(y=posterior$data$unit)
  })
  
  # Return list of plots
  return(c(list(p.cc), p.rs))
}