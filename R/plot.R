plot.data <- function(posterior, data) {
  pars <- names(posterior$parameters)
  aux <-  c(posterior$auxiliary$ll, posterior$auxiliary$cc, posterior$auxiliary$sc, posterior$auxiliary$rs)
  if ("SANMCMC" %in% class(data)) {
    # Collect auxiliary data columns
    data$final[, c("ll", aux), with=FALSE]
  } else if ("data.table" %in% class(data)) {
    if (all(c("ll_tot", aux) %in% colnames(data)))
      # Collect auxiliary data columns
      data[, c("ll_tot", aux), with=FALSE]
    else if (all(pars %in% colnames(data)))
      # Collect parameter columns and re-simulate model
      posterior$loglikelihood(data[, pars, with=FALSE])
    else
      stop("data must contains either all parameters or all auxiliary data columns from the posterior")
  } else {
    # Make sure that x has two dimensions,
    data <- rbind(data)
    if (all(c("ll_tot", aux) %in% colnames(data)))
      # Collect auxiliary data columns
      data[, c("ll_tot", aux)]
    else if (all(pars %in% colnames(data)))
      # Collect parameter columns and re-simulate model
      posterior$loglikelihood(data[, pars])
    else
      stop("data must contains either all parameters or all auxiliary data columns from the posterior")
  }
}

#' Plot log-normal confidence intervals defining the posterior overlayed with model predictions
#'
#' @export
plot.SANPosterior <- function(posterior, data, params=NULL, ll=NULL, plotlist=FALSE, ...) {
  # Convert input data (MCMC results, parameter vector or result of posterior$loglikelihood) into usable form
  data <- as.data.table(plot.data(posterior, data))
  data[, index := 1:.N ]

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
  ll.cc <- melt(data, id.vars=c("index", "ll_cc"),
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
    ll.rs <- melt(data, id.vars=c("index", paste0("ll_rs_", rs_day)),
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
plot.SANCombinedPosterior <-  function(posterior, data, plotlist=FALSE, ...) {
  # Convert input data (MCMC results, parameter vector or result of posterior$loglikelihood) into usable form
  data <- as.data.table(plot.data(posterior, data))

  # Plot individual components
  ls <- names(posterior$components)
  names(ls) <- ls
  p <- do.call(c, lapply(ls, function(l) {
    # Extract the columns belonging to the current component from the data.
    prefix <- paste0(l, ".")
    cols <- colnames(data)
    cols.comp <- cols[startsWith(cols, prefix)]
    data.comp <- data[, cols.comp, with=FALSE]
    colnames(data.comp) <- substr(cols.comp, start=nchar(prefix)+1, stop=nchar(cols.comp))
    # Call the component's plotting function
    plot(posterior$components[[l]], data.comp, plotlist=TRUE)
  }))
  
  # Either return a list of plots, or arrange them in a grid with ggarrange
  if (!plotlist)
    return(ggpubr::ggarrange(plotlist=p, align="hv", ...))
  else
    return(p)
}
