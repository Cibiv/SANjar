library(shiny)
library(rhandsontable)
library(shinycssloaders)

fluidPage(
    withMathJax(),
    HTML(paste0("<script type='text/x-mathjax-config'>",
                "MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});",
                "</script>")),
    
    # Application title
    titlePanel("SAN Model Explorer"),
    
    # Tabs
    tabsetPanel(
        type="tabs",
        tabPanel("Rates & Deterministic Model",
                 fluidRow(column(width=6,
                                 tags$label(class="control-label", "Rates"),
                                 rHandsontableOutput("rates", width="450px"),
                                 helpText(HTML("At each point in time, the non-empty rate specified for the ",
                                          "smallest T<sub>max</sub> that still lies in the future applies, ",
                                          "or zero if no such non-empty rate exists. Empty cells thus effectively ",
                                          "contain the value from the <i>next</i> non-empty row, or zero if ",
                                          "they are the last empty cell in the respective column."))),
                          column(width=6,
                                 numericInput("s0", HTML("Initial number of cells (s<sub>0</sub>)"), min=1, max=3e4, value=0))),
                 fluidRow(column(width=12,
                                 plotOutput("deterministic")))
        ),
        tabPanel("Lineage Size Distribution",
                 fluidRow(column(width=2,
                                 numericInput("p_cutoff", "Discretization err.", width="150px",
                                              min=1e-6, max=1e-1, value=1e-2)),
                          column(width=2,
                                 numericInput("pcr_efficiency", "PCR efficiency", width="150px",
                                              min=0, max=1, value=NA_real_, step=0.05),
                                 htmlOutput("pcr_efficiency_auto_message")),
                          column(width=2,
                                 numericInput("library_size", "Library size", width="150px",
                                              min=1, max=Inf, value=NA_real_, step=1e6),
                                 htmlOutput("library_size_auto_message")),
                          column(width=2,
                                 numericInput("phantom_threshold", "Phantom threshold", width="150px",
                                              min=1, max=Inf, value=NA_integer_, step=1),
                                 htmlOutput("phantom_threshold_auto_message"))),
                 fluidRow(column(width=5,
                                 withSpinner(plotOutput("stochastic_lsd", width="500px")),
                                 selectInput("day_lsd", label = "Day to show", width="150px",
                                             selected=character(0),
                                             choices=character(0))),
                          column(width=5,
                                 withSpinner(plotOutput("stochastic_nlineages", width="500px")),
                                 checkboxInput("stochastic_nlineages_logy", label="logarithmic y-axis", value=TRUE)))
        ),
        tabPanel("Parameter Sets",
                 fluidRow(column(width=6,
                                 selectInput("loadfrom", label="Use parameter set", width="500px",
                                             selected=character(0), choices=c("select set to load"="", character(0))))),
                 fluidRow(column(width=6,
                          textInput("saveas", label="Save parameter set as", width="500px",
                                    placeholder="filename"),
                          actionButton("save", label="Save")))
        )
    )
)
