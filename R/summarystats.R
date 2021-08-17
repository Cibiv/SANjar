#' Locates the distribution's mode(s) using the mean-shift algorithm
#'
#' @export
locate.modes <- function(x, ...) UseMethod("locate.modes")

#' Computes various summary statistics
#' 
#' @export
summarystats <- function(x, ...) UseMethod("summarystats")

#' Locates the posterior distribution's mode(s) using the mean-shift algorithm
#'
#' @export
locate.modes.SANMCMC <- function(sanmcmc, tolerance=0.1, adjust=1.0) {
  x <- as.matrix(sanmcmc$final[, names(sanmcmc$variables), with=FALSE])
  H <- ks::Hpi(x, deriv.order=1, nstage=2-(ncol(x)>2))
  ks::kms(x, H=H*(adjust^2), min.clust.size=0.1*nrow(x), tol.clust=tolerance)
}

#' Computes various summary statistics of the posterior distribution
#' 
#' @export
summarystats.SANMCMC <- function(sanmcmc, modes, expressions=names(sanmcmc$variables)) {
  stats <- c("mean", "std. dev.", "median", "mad*1.48", "mode", "ML")
  
  # Run mean-shift algorithm if necessary and find mode
  if (missing(modes))
    modes <- locate.modes(sanmcmc)
  i.mode <- which.max(modes$nclust.table)
  mode <- modes$mode[i.mode,]
  
  sanmcmc$final[, {
    stats <- lapply(expressions, function(expr) {
      # Translate strings into expressions. For added convenience, column
      # names are quoted automatically, this makes e.g. writing "40S-40A" work.
      if (is.character(expr)) {
        expr <- gsub(paste0("(", paste0("(?:", colnames(.SD), ")", collapse="|"), ")"), "`\\1`", expr)
        expr <- parse(text=expr)
      }
      if (!is.expression(expr) && !is.name(expr))
        stop("invalid expression of type ", class(expr))
      
      # Evaluate expression
      v <- eval(expr)
      
      # Determine mode, either for mean-shift results or through KDE
      expr.str <- as.character(expr)
      m <- if (expr.str %in% names(mode))
        mode[expr.str]
      else {
        d <- density(v)
        d$x[which.max(d$y)]
      }
      
      # Compute summary statistics 
      c(`mean`=mean(v),
        `std. dev.`=sd(v),
        `median`=median(v),
        `mad*1.48`=mad(v),
        `mode`=m,
        `ML`=v[which.max(ll)])
    })
    
    # Prepend "stat" column which indicates which row represents
    # which statistic and name the columns
    names(stats) <- as.character(expressions)
    c(list(statistic=names(stats[[1]])), stats)
  }]
}
