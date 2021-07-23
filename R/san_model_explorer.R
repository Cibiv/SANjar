#' Runs the SAN model explorer app
#'
#' @export
san_model_explorer <- function(dataset, location.parametersets=".") {
  app_dir <- system.file("SANModelExplorer", package = "SANjar", mustWork=TRUE)

  # Parameters:
  #  LOCATION.PARAMETERSETS
  #  DEFAULT.PARAMETERSET
  #  DATASET

  app_env <- new.env()
  app_env$LOCATION.PARAMETERSETS <- location.parametersets
  app_env$DEFAULT.PARAMETERSET <- lt47.ps
  app_env$DATASET <- dataset

  app_globalR <- file.path(app_dir, "global.R")
  if (file.exists(app_globalR))
    local(source(app_globalR, local=TRUE, echo=FALSE, keep.source=TRUE), envir=app_env)
  app_ui <- source(file.path(app_dir, "ui.R"), local=new.env(parent=app_env),
                   echo=FALSE, keep.source=TRUE)$value
  app_server <- source(file.path(app_dir, "server.R"), local=new.env(parent=app_env),
                       echo=FALSE, keep.source=TRUE)$value
  
  shinyApp(app_ui, app_server)
}
