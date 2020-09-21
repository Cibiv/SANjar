library(readr)
library(shiny)
library(rhandsontable)
library(shinycssloaders)
library(bsplus)

CONTAINER_REVISION <- (if (file.exists("CONTAINER_REVISION"))
    read_file("CONTAINER_REVISION")
else
    NA_character_)

fluidPage(
    # Application title
    titlePanel(paste0("SAN Model Explorer",
                      if (!is.na(CONTAINER_REVISION)) paste0(" (", CONTAINER_REVISION ,")") else "")),
    
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
                                 tags$label(class="control-label", "Cell type composition, model vs. data"),
                                 plotOutput("deterministic_celltypes", height="500px"))),
                 fluidRow(column(width=6,
                                 checkboxInput("deterministic_cellcounts_incsim", width="100%",
                                               label = "show simulated cell counts in additional to theoretical counts", value=FALSE)))
        ),
        tabPanel("Lineage Size Distribution",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Lineage size distribution, model vs. data"),
                                 withSpinner(plotOutput("stochastic_lsd", width="100%")),
                                 selectInput("day_lsd", label = "Day to show",
                                             selected=character(0),
                                             choices=character(0))),
                          column(width=6,
                                 tags$label(class="control-label", "Number of lineages, model vs. data"),
                                 withSpinner(plotOutput("stochastic_nlineages", width="100%")),
                                 checkboxInput("stochastic_nlineages_logy", label="logarithmic y-axis", value=TRUE))),
                 fluidRow(column(width=12,
                                 bs_button("Show/hide simulation settings") %>%
                                     bs_attach_collapse("simulation_settings"))),
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
