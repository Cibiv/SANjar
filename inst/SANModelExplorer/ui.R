library(shiny)
library(rhandsontable)
library(shinyjs)
library(shinycssloaders)
library(bsplus)

fluidPage(
    # Application title
    titlePanel(paste0("SAN Model Explorer",
                      if (exists("CONTAINER_REVISION")) paste0(" (", CONTAINER_REVISION ,")") else "")),
    
    # Enable ShinyJS
    useShinyjs(),

    # Tabs
    tabsetPanel(
        type="tabs",
        tabPanel("Rates & Deterministic Model",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Rates"),
                                 rHandsontableOutput("rates", width="99%"),
                                 helpText(HTML("At each time <i>t</i> and for each rate (i.e. column), the first ",
                                          "non-empty value whose <i>T</i><sub>max</sub> is larger than <i>t</i> applies, ",
                                          "or zero if no such value exists. Empty cells thus effectively contain the ",
                                          "first non-empty value <i>below</i> the empty cell, or zero if all following ",
                                          "values are empty as well."))),
                          column(width=6,
                                 numericInput("s0", HTML("Initial number of cells (s<sub>0</sub>)"), min=1, max=Inf, value=0),
                                 checkboxInput("s0_scale", label=HTML("scale s<sub>0</sub> for partial organoids"), value=TRUE, width="100%"),
                                 htmlOutput("s0_message"))),
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Cell counts, model vs. data"),
                                 plotOutput("deterministic_cellcounts", height="500px"))),
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
                                 checkboxInput("stochastic_lsd_normlibsize", label="normalize library sizes", value=TRUE, width="100%")),
                          column(width=6,
                                 checkboxInput("stochastic_nlineages_logy", label="logarithmic y-axis", value=TRUE, width="100%")),
                          column(width=6,
                                 checkboxInput("stochastic_lsd_incpuremodel", label="show pure SAN model without PCR+Sequencing", value=FALSE, width="100%"))),
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
        tabPanel("Data & Parameter Sets",
                 fluidRow(column(width=6,
                                 selectInput("datagroup", label="Group", width="100%",
                                             selected=character(0), choices=character(0)))),
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
