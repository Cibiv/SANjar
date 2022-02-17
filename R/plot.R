plot.data <- function(posterior, data, include=character()) {
  pars <- names(posterior$parameters)
  aux <-  c(posterior$auxiliary$ll, posterior$auxiliary$cc, posterior$auxiliary$sc, posterior$auxiliary$rs)
  include <- setdiff(include, aux)
  if ("SANMCMC" %in% class(data)) {
    # Collect auxiliary data columns, rename ll to ll_tot
    d <- data$final[, c("ll", aux), with=FALSE]
    colnames(d)[1] <- "ll_tot"
    d
  } else if ("data.table" %in% class(data)) {
    if (all(c("ll_tot", aux) %in% colnames(data)))
      # Collect auxiliary data columns
      data[, c("ll_tot", aux, include), with=FALSE]
    else if (all(c("ll", aux) %in% colnames(data))) {
      # Collect auxiliary data columns, rename ll to ll_tot
      d <- data[, c("ll", aux), with=FALSE]
      colnames(d)[1] <- "ll_tot"
      d
    }  
    else if (all(pars %in% colnames(data)))
      # Collect parameter columns and re-simulate model
      posterior$loglikelihood(data[, pars, with=FALSE])
    else
      stop("data must contains either all parameters or all auxiliary data columns from the posterior ",
           "for of the former is missing ", paste0(setdiff(pars, colnames(data)), sep=", "), " and of ",
           "the latter is missing ", paste0(setdiff(c("ll_tot", aux), colnames(data))))
  } else {
    # Make sure that x has two dimensions,
    data <- rbind(data)
    if (all(c("ll_tot", aux) %in% colnames(data)))
      # Collect auxiliary data columns
      data[, c("ll_tot", aux, include)]
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
plot.SANPosterior <- function(posterior, data, ll="ll_tot", plotlist=FALSE, rs_days, ...) {
  # For combined posterior, the leading dot indicates the field within each component is
  # to be used to color trajectories. Since there are no further nesting levels, simply
  # strip any leadind dot left in place
  if (startsWith(ll, "."))
    ll <- substring(ll, 2, nchar(ll))
  # Convert input data (MCMC results, parameter vector or result of posterior$loglikelihood) into usable form
  data <- as.data.table(plot.data(posterior, data, include=as.character(ll)))
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
  llfield.cc <- if (!is.null(ll)) as.name(ll) else as.name("ll_cc")
  ll.cc <- melt(data, id.vars=c("index", ll, "ll_cc"),
                measure.vars=paste0("CC", posterior$parametrization$cc_days))[, day := {
                  as.integer(sub("^CC([0-9]+)$", "\\1", variable))
                }]
  # Index groups in order of increasing likelihood to plot likelist trajectories on top
  i <- 0
  ll.cc[ll.cc[order(eval(llfield.cc)), list(index=unique(index))]
        , index_ll := {i <<- i + 1}
        , on=.(index), by=.EACHI]

  # Plot total organoid sizes
  p.cc <- ggplot2::ggplot() +
    ggplot2::geom_line(data=ll.cc, ggplot2::aes(x=day, y=pmax(value, posterior$arguments$min.size),
                                                color=!!llfield.cc, group=index_ll)) +
    ggplot2::geom_errorbar(data=posterior$data$cc, ggplot2::aes(x=day, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
    make_scale_log10(ggplot2::scale_y_log10) +
    ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                   low="red", mid="green", midpoint=midpoint(ll.cc[, eval(llfield.cc)]), high="blue") +
    ggplot2::labs(y="cells")
  
  # Plot rank-sizes
  if (is.null(rs_days))
    rs_days <- posterior$parametrization$rs_days
  if (length(rs_days) > 0)
    names(rs_days) <- paste0("day", rs_days)
  p.rs <- lapply(rs_days, function(rs_day) {
    llfield.rs <- if (!is.null(ll)) as.name(ll) else  as.name(paste0("ll_rs_", rs_day))

    # Extract per-day cell counts as a table suitable for plotting
    ll.rs <- melt(data, id.vars=c("index", ll, paste0("ll_rs_", rs_day)),
                  measure.vars=paste0("RS", rs_day, "_", posterior$parametrization$rs_ranks))
    ll.rs <- ll.rs[, c("ll_rs", "day", "rank") := {
      list(eval(llfield.rs),
           rs_day,
           as.integer(sub(paste0("^RS", rs_day, "_([0-9]+)$"), "\\1", variable)))
    }]
    # Index groups in order of increasing likelihood to plot likelist trajectories on top
    i <- 0
    ll.rs[ll.rs[order(ll_rs), list(index=unique(index))]
          , index_ll := {i <<- i + 1}
          , on=.(index), by=.EACHI]

    # Plot
    ggplot2::ggplot() +
      ggplot2::geom_line(data=ll.rs, ggplot2::aes(x=rank, y=pmax(value, posterior$arguments$min.size),
                                                  color=!!llfield.rs, group=index_ll)) +
      ggplot2::geom_errorbar(data=posterior$data$rs[day==rs_day],
                             ggplot2::aes(x=rank, ymin=10^(logmean-logsd), ymax=10^(logmean+logsd))) +
      ggplot2::annotate("text", x=Inf, y=Inf, hjust=1, vjust=1, col="black", label=paste0("day ", rs_day)) +
      make_scale_log10(ggplot2::scale_x_log10) +
      make_scale_log10(ggplot2::scale_y_log10) +
      ggplot2::scale_color_gradient2(name=expression(log~L), trans=log_ll,
                                     low="red", mid="green", midpoint=midpoint(ll.rs[, eval(llfield.rs)]), high="blue") +
      ggplot2::labs(y=posterior$data$unit)
  })
  
  # Either return a list of plots, or arrange them in a grid with ggarrange
  p <- c(list(cc=p.cc), p.rs)
  if (!plotlist) {
    args <- list(align="hv", common.legend=(!is.null(ll) && (ll=="ll")))
    args <- modifyList(args, list(...))
    return(do.call(ggpubr::ggarrange, c(list(plotlist=p), args)))
  } else
    return(p)
}

#' Combine the plots from all components of a combined posterior distribution
#' 
#' @export
plot.SANCombinedPosterior <-  function(posterior, data, ll=".ll_tot", plotlist=FALSE, rs_days=NULL, ...) {
  # Convert input data (MCMC results, parameter vector or result of posterior$loglikelihood) into usable form
  ll.include <- if (!startsWith(ll, ".")) ll else character()
  data <- as.data.table(plot.data(posterior, data, include=ll.include))

  # Plot individual components
  plots <- list()
  labels <- list()
  i <- 1
  for(cn in names(posterior$components)) {
    # Extract the columns belonging to the current component from the data.
    prefix <- paste0(cn, ".")
    cols <- colnames(data)
    cols.comp <- cols[startsWith(cols, prefix)]
    data.comp <- data[, c(ll.include, cols.comp), with=FALSE]
    colnames(data.comp) <- c(ll.include, substr(cols.comp, start=nchar(prefix)+1, stop=nchar(cols.comp)))
    # Call the component's plotting function
    p <- lapply(plot(posterior$components[[cn]], data.comp, ll=ll, rs_days=rs_days, plotlist=TRUE), function(p) {
      p + ggplot2::theme(plot.margin = margin(t=20, r=5)) 
    })
    j <- i + length(p)
    plots[i:(j-1)] <- p
    labels[i:(j-1)] <- rep(list(cn), j-i)
    i <- j
  }

  # Either return a list of plots, or arrange them in a grid with ggarrange
  if (!plotlist) {
    args <- list(align="hv",
                 labels=labels, label.x=0.5, label.y=1, hjust=0.5, vjust=1.5,
                 font.label=list(size=11, face="bold"),
                 common.legend=(!is.null(ll) && (ll=="ll_grandtot")))
    args <- modifyList(args, list(...))
    return(do.call(ggpubr::ggarrange, c(list(plotlist=plots), args)))
  } else
    return(plots)
}

#' Plot posterior distribution
#'
#' @export
plot.SANMCMC <- function(sanmcmc, type="violins", which="final", point=NULL, expressions=names(sanmcmc$variables), extra=character(), ...) {
  type <- match.arg(type, c("violins", "parcord"))
  switch(type,
         `violins`=plot.SANMCMC.violins(sanmcmc, which=which, point=point, expressions=expressions, extra=extra, ...),
         `parcord`=plot.SANMCMC.parcord(sanmcmc, which=which, point=point, expressions=expressions, extra=extra, ...))
}

plot.SANMCMC.violins <- function(sanmcmc, which="final", point, expressions, extra) {
  data <- sanmcmc.evaluate(sanmcmc, which=which, expressions=expressions, extra=extra)
  p <- ggplot2::ggplot(data=melt(data, id.vars="chain", measure.vars=attr(data, "expressions"))) +
    ggplot2::geom_violin(aes(x=variable, y=value), scale="width")
  if (!is.null(point))
    p <- p + ggplot2::geom_segment(data=melt(as.data.table(as.list(point)), measure.vars=attr(data, "expressions")),
                                   aes(x=as.integer(variable)-0.45, xend=as.integer(variable)+0.45,
                                       y=value, yend=value),
                                   linetype="dashed")
  return(p)
}

plot.SANMCMC.parcord <- function(sanmcmc, which, point, expressions, extra, linecolor="#00000010", pointcolor="darkblue") {
  data <- sanmcmc.evaluate(sanmcmc, which=which, expressions=expressions, extra=extra)
  p <- ggplot2::ggplot(data=melt(data, id.vars="chain", measure.vars=attr(data, "expressions"))) +
    ggplot2::geom_line(aes(x=variable, y=value, group=chain), color=linecolor)
  if (!is.null(point))
    p <- p + ggplot2::geom_point(data=melt(as.data.table(as.list(point)), measure.vars=attr(data, "expressions")),
                                 aes(x=as.integer(variable), y=value), color=pointcolor)
  return(p)
}
