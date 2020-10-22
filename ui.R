library(readr)
library(shiny)
library(rhandsontable)
library(shinyjs)
library(shinycssloaders)
library(bsplus)

CONTAINER_REVISION <- (if (file.exists("CONTAINER_REVISION"))
    trimws(read_file("CONTAINER_REVISION"))
else
    NA_character_)

fluidPage(
    # Application title
    titlePanel(paste0("SAN Model Explorer",
                      if (!is.na(CONTAINER_REVISION)) paste0(" (", CONTAINER_REVISION ,")") else "")),
    
    # Enable ShinyJS
    useShinyjs(),

    # Tabs
    tabsetPanel(
        type="tabs",
        tabPanel("Rates & Deterministic Model",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Rates"),
                                 rHandsontableOutput("rates", width="100%"),
                                 helpText(HTML("At each time <i>t</i> and for each rate (i.e. column), the first ",
                                          "non-empty value whose <i>T</i><sub>max</sub> is larger than <i>t</i> applies, ",
                                          "or zero if no such value exists. Empty cells thus effectively contain the ",
                                          "first non-empty value <i>below</i> the empty cell, or zero if all following ",
                                          "values are empty as well."))),
                          column(width=6,
                                 numericInput("s0", HTML("Initial number of cells (s<sub>0</sub>)"), min=1, max=3e4, value=0))),
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Cell counts, model vs. data"),
                                 plotOutput("deterministic_cellcounts", height="500px")),
                          column(width=6,
                                 tags$label(class="control-label", "Organoid composition, model vs. data"),
                                 plotOutput("deterministic_celltypes", height="500px"))),
                 fluidRow(column(width=3,
                                 checkboxInput("deterministic_cellcounts_incsim", width="100%",
                                               label = "also show simulated cell counts", value=FALSE)))
        ),
        tabPanel("Lineage Size Distribution",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Lineage size distribution, model vs. data"),
                                 withSpinner(plotOutput("stochastic_lsd", width="100%"), hide.ui=FALSE)),
                          column(width=6,
                                 tags$label(class="control-label", "Number of lineages, model vs. data"),
                                 withSpinner(plotOutput("stochastic_nlineages", width="100%"), hide.ui=FALSE))),
                 fluidRow(column(width=6,
                                 checkboxInput("stochastic_lsd_incpuremodel", label="show pure SAN model without PCR+Sequencing", value=FALSE, width="100%"),
                                 checkboxInput("stochastic_lsd_incpowerlaw", label="show powerlaw model", value=FALSE, width="100%")),
                          column(width=6,
                                 checkboxInput("stochastic_nlineages_logy", label="logarithmic y-axis", value=TRUE, width="100%"),
                                 checkboxInput("stochastic_nlineages_incpowerlaw", label="show powerlaw model", value=FALSE, width="100%"))),
                 fluidRow(
                          column(width=3,
                                 selectInput("day_lsd", label = "Day to show",
                                             selected=character(0),
                                             choices=character(0))),
                          column(width=3,
                                 selectInput("stochastic_lsd_plottype", label="Plot lineage size distribution as",
                                             selected="log-rank vs. log-lineagesize",
                                             choices=c("log-rank vs. log-lineagesize",
                                                       "density of log-lineagesize"))),
                          column(width=6)),
                 fluidRow(column(width=3,
                                 bs_button("Show/hide simulation settings") %>%
                                     bs_attach_collapse("simulation_settings")),
                          column(width=3,
                                 downloadButton("stochastic_lsd_download", label="Download Simulation Results"))),
                 bs_collapse("simulation_settings",
                             fluidRow(column(width=3,
                                             selectInput("p_cutoff", "Discretization err.", width="100%",
                                                         choices=c("1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"), selected="1e-3"),
                                             htmlOutput("p_cutoff_message")),
                                      column(width=3,
                                             numericInput("pcr_efficiency", "PCR efficiency", width="100%",
                                                          min=0, max=1, value=NA_real_, step=0.05),
                                             htmlOutput("pcr_efficiency_auto_message")),
                                      column(width=3,
                                             numericInput("library_size", "Library size", width="100%",
                                                          min=1, max=Inf, value=NA_real_, step=1e6),
                                             htmlOutput("library_size_auto_message")),
                                      column(width=3,
                                             numericInput("phantom_threshold", "Phantom threshold", width="100%",
                                                          min=1, max=Inf, value=NA_integer_, step=1),
                                             htmlOutput("phantom_threshold_auto_message")))
                 )
        ),
        tabPanel("Pareto Equality Indices",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Estimated Pareto equalities indices over time"),
                                 withSpinner(plotOutput("pareto_alpha", width="100%"), hide.ui=FALSE))),
                 fluidRow(column(width=6,
                                 checkboxInput("pareto_alpha_incpuremodel", label="show pure SAN model without PCR+Sequencing", value=FALSE, width="100%")))),
        tabPanel("S-cell extinction dynamics",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Lineage size on day 40 vs. S-cell extinction time"),
                                 withSpinner(plotOutput("stochastic_scellext_vs_lineagesize", width="100%", height="250px"), hide.ui=FALSE),
                                 tags$label(class="control-label", "Fraction of lineages with extant S-cell population"),
                                 withSpinner(plotOutput("stochastic_scellext_vs_time", width="100%", height="250px"), hide.ui=FALSE)),
                          column(width=6,
                                 tags$label(class="control-label",  htmlOutput("stochastic_scellext_trajectory_label")),
                                 withSpinner(plotOutput("stochastic_scellext_trajectory", width="100%"), hide.ui=FALSE),
                                 sliderInput("day_scellext", "S-cell extinction time", width="100%",
                                             min=0, max=0, step=1, value=0),
                                 tags$span(
                                     tags$div(selectInput("stochastic_scellext_trajectory_video_fps", choices=c(1, 6, 12, 24), label="FPS", width="100px"), style="display:inline-block; vertical-align: middle"),
                                     tags$div(actionButton("stochastic_scellext_trajectory_video_generate", "Generate video", width="150px"), style="display:inline-block; vertical-align: middle"),
                                     tags$div(downloadButton("stochastic_scellext_trajectory_video_download", "Download video", width="150px"), style="display:inline-block; vertical-align: middle")
                                 ))),
                 fluidRow(column(width=6,
                                 bs_button("Show/hide simulation settings") %>%
                                     bs_attach_collapse("scellext_simulation_settings"))),
                 bs_collapse("scellext_simulation_settings",
                             fluidRow(column(width=3,
                                             selectInput("scellext_nlineages", "Number of lineages to simulate", width="100%",
                                                         choices=c("1e4", "1e5", "1e6", "1e7"), selected="1e5")))
                 )
        ),
        tabPanel("Parameter Sets",
                 fluidRow(column(width=6,
                                 selectInput("loadfrom", label="Use parameter set", width="100%",
                                             selected=character(0), choices=c("select set to load"="", character(0))))),
                 fluidRow(column(width=6,
                          textInput("saveas", label="Save parameter set as", width="100%",
                                    placeholder="filename"),
                          actionButton("save", label="Save")))
        )
    )
)
