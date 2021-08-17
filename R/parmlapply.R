#' Parallel version of mapply that operates on lists
#'
#' @export
parMLapply <- function(cluster, FUN, ..., .batchcount=10*length(cluster)) {
  args <- list(...)
  if (length(args) == 0)
    stop("parMLapply needs at least one argument")
  
  # Arguments must be lists
  if (!Reduce(`&&`, Map(is.list, args), TRUE))
    stop("arguments to parMLapply must be lists")
  
  # Arguments must have the same length or length one. Arguments of length 1
  # are recycled as many times as required.
  ls <- as.vector(sapply(args, length))
  l <- ls[1]
  if (any((ls != l[1]) & (ls != 1)))
    stop("arguments to parMLapply must all have the same length or length 1")
  if (l > 1)
    args <- Map(function(a) { if (length(a) == 1) rep(a, l) else a }, args)
  
  # Split work into at most the specified number of batches (.batchcount=B)
  # If the number of work units (l) is not an integral multiple of the
  # number of batches, the first 1...f batches will have the
  # full size (b), and batches (f+1)...B contain one unit less.
  stopifnot(.batchcount > 0)
  b <- ceiling(l / .batchcount)
  f <- l - (b-1)*.batchcount
  k <- 1
  batches <- lapply(1:min(.batchcount, l), function(i) {
    s <- if (i <= f) b else b - 1
    k <<- k + s
    if (s > 0) (k-s):(k-1) else NULL
  })
  stopifnot(unlist(batches) == 1:l)
  args.batches <- Map(function(a) { Map(function(b) a[b], batches) }, args)
  
  # Distribute the batches over the nodes in the cluster,
  # and concatenate results back together.
  fun.batch <- function(...) {
    mapply(FUN, ..., SIMPLIFY=FALSE)
  }
  environment(fun.batch) <- list2env(list(FUN=FUN), parent=baseenv())
  r <- do.call(c, do.call(parallel::clusterMap, c(list(cl=cluster, fun=fun.batch), args.batches,
                                                  list(SIMPLIFY=FALSE, .scheduling="dynamic"))))
  return(r)
}
