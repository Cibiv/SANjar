#' Runs the SAN model explorer app
#'
#' @export
san_model_explorer <- function(..., location.models=".") {
  app_dir <- system.file("SANModelExplorer", package = "SANjar", mustWork=TRUE)

  # Parameters:
  #  LOCATION.MODELS
  #  DEFAULT.MODEL
  #  MODELS
  #  DATASET

  app_env <- new.env()
  
  app_env$LOCATION.MODELS <- location.models

  ps_expr <- substitute(list(...))
  ps <- list(...)
  ps_labs <- names(ps)
  ps_lab <- if (is.null(ps_labs)) rep(FALSE, length(ps)) else (ps_labs != "")
  names(ps)[!ps_lab] <- as.character(ps_expr[2:length(ps_expr)])[!ps_lab]
  app_env$DATASETS <- Filter(function(p) class(p)=="LTData", ps)
  app_env$MODELS <- Filter(function(p) class(p)=="SANModel", ps)
  if (length(app_env$MODELS) == 0)
    app_env$MODELS <- list(`lt47`=lt47.model)
  app_env$DEFAULT.MODEL <- paste0("preset:", names(app_env$MODELS)[1])
  
  
  app_globalR <- file.path(app_dir, "global.R")
  if (file.exists(app_globalR))
    local(source(app_globalR, local=TRUE, echo=FALSE, keep.source=TRUE), envir=app_env)
  app_ui <- source(file.path(app_dir, "ui.R"), local=new.env(parent=app_env),
                   echo=FALSE, keep.source=TRUE)$value
  app_server <- source(file.path(app_dir, "server.R"), local=new.env(parent=app_env),
                       echo=FALSE, keep.source=TRUE)$value
  
  shinyApp(app_ui, app_server)
}
