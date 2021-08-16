#' Plot log-normal confidence intervals defining the posterior overlayed with model predictions
#'
#' @export
plot.SANPosterior <- function(posterior, params=NULL, ll=NULL, plotlist=FALSE, ...) {
  if (is.null(ll)) {
    if (is.null(params))
      stop("either a parameter vector or the result of the posterior distribution loglikelihood function must be specified")
    message("ll not specified, evaluating model for parameters ",
            paste0(colnames(params), "=", signif(params, 3), collapse=", "))
    ll <- posterior$loglikelihood(params)
  }
  
  # Convert output of loglikelihood function into a data.table  
  ll.tab <- as.data.table(ll)
  ll.tab[, index := 1:.N ]
  
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
  ll.cc <- melt(ll.tab, id.vars=c("index", "ll_cc"),
                measure.vars=paste0("CC", posterior$parametrization$cc_days))[, day := {
                  as.integer(sub("^CC([0-9]+)$", "\\1", variable))
                }]
  
  # Plot total organoid sizes
  p.cc <- ggplot2::ggplot() +
    ggplot2::geom_errorbar(data=posterior$data$cc, ggplot2::aes(x=day, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
    ggplot2::geom_line(data=ll.cc, ggplot2::aes(x=day, y=pmax(value, posterior$arguments$min.size),
                                                color=ll_cc, group=index)) +
    make_scale_log10(scale_y_log10) +
    ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                   low="red", mid="green", midpoint=midpoint(ll.cc$ll_cc), high="blue") +
    ggplot2::labs(y="cells")
  
  # Plot rank-sizes
  rs_days <- posterior$parametrization$rs_days
  if (length(rs_days) > 0)
    names(rs_days) <- paste0("day", rs_days)
  p.rs <- lapply(rs_days, function(rs_day) {
    # Extract per-day cell counts as a table suitable for plotting
    ll.rs <- melt(ll.tab, id.vars=c("index", paste0("ll_rs_", rs_day)),
                  measure.vars=paste0("RS", rs_day, "_", posterior$parametrization$rs_ranks))
    ll.rs <- ll.rs[, c("ll_rs", "day", "rank") := {
      list(eval(as.name(paste0("ll_rs_", rs_day))),
           rs_day,
           as.integer(sub(paste0("^RS", rs_day, "_([0-9]+)$"), "\\1", variable)))
    }]
    
    # Plot
    ggplot() +
      ggplot2::geom_errorbar(data=posterior$data$rs[day==rs_day],
                             ggplot2::aes(x=rank, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
      ggplot2::geom_line(data=ll.rs, ggplot2::aes(x=rank, y=pmax(value, posterior$arguments$min.size),
                                                  color=ll_rs, group=index)) +
      annotate("text", x=Inf, y=Inf, hjust=1, vjust=1, col="black", label=paste0("day ", rs_day)) +
      make_scale_log10(scale_x_log10) +
      make_scale_log10(scale_y_log10) +
      ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                     low="red", mid="green", midpoint=midpoint(ll.rs$ll_rs), high="blue") +
      ggplot2::labs(y=posterior$data$unit)
  })
  
  # Either return a list of plots, or arrange them in a grid with ggarrange
  p <- c(list(cc=p.cc), p.rs)
  if (!plotlist)
    return(ggpubr::ggarrange(plotlist=p, align="hv", ...))
  else
    return(p)
}

#' Combine the plots from all components of a combined posterior distribution
#' 
#' @export
plot.SANCombinedPosterior <-  function(posterior, params=NULL, ll=NULL, plotlist=FALSE, ...) {
  ls <- names(posterior$components)
  names(ls) <- ls
  p <- do.call(c, lapply(ls, function(l) {
    # If only parameters, no loglikelihood output was specified, simply pass the params
    # on to the individual plotting method
    if (is.null(ll))
      return(plot(posterior$components[[l]], params=params, plotlist=TRUE))
    
    # Otherwise, extract the columns beloning to the current component from the loglikelihood
    # output (format is <componentlabel>.<colname>), and remove the prefix
    l.prefix <- paste0(l, ".")
    ll.cols <- colnames(ll)
    ll.cols.comp <- ll.cols[startsWith(ll.cols, l.prefix)]
    ll.comp <- ll[, ll.cols.comp, with=FALSE]
    colnames(ll.comp) <- substr(ll.cols.comp, start=nchar(l.prefix)+1, stop=nchar(ll.cols.comp))
    # Call the component's plotting function
    plot(posterior$components[[l]], ll=ll.comp, plotlist=TRUE)
  }))
  
  # Either return a list of plots, or arrange them in a grid with ggarrange
  if (!plotlist)
    return(ggpubr::ggarrange(plotlist=p, align="hv", ...))
  else
    return(p)
}
