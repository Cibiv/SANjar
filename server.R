library(data.table)
library(Hmisc)
library(gwpcR)
library(SANsimulatoR)
library(ggplot2)
library(scales)
library(gganimate)
library(av)           # mp4 output for gganimate
library(transformr)   # required by gganimate
library(shiny)
library(shinyjs)
library(rhandsontable)
source("rank_size_utilities.R")

# Pattern defining valid "saveas" filenames for parameter sets
FILENAME.PATTERN <- "^[A-Za-z0-9. -]*$"

# Location to store the parameter sets in
LOCATION.PARAMETERSETS <- "parametersets"

# Line width
LWD <- 0.8
LWD2 <- 1.5

# Empty rates table
RATES.EMPTY <- as.data.table({
    x <- rep(list(numeric(0)), 1+length(SAN.RATENAMES))
    names(x) <- c("Tmax", SAN.RATENAMES)
    x
})

# ggplot2 setup
theme_set(theme_grey(base_size = 20))
theme_update(legend.position="bottom",
             legend.key.width = unit(15,"mm"))

# Load data
load("data/organoidsizes.rd")
load("data/celltypes.rd")
load("data/lt47.processed.rd")

# Return a list of saved parameter sets
parametersets.list <- function() {
    files <- list.files(file.path(LOCATION.PARAMETERSETS))
    p <- "\\.rd$"
    sub(x=files[grepl(p, files)], pattern=p, replacement="")
}

# Convert a parameterset name to a full path
parameterset.fullpath <- function(name) {
    file.path(LOCATION.PARAMETERSETS, paste0(name, ".rd"))
}

# Fill in empty cells in the rate table
parameterset.fillempty <- function(value) {
    for(c in SAN.RATENAMES) {
        # Fetch row corresponding to rate x
        x <- value[[c]]
        x.set <- !is.na(x)
        # Set empty entries to the next non-empty entry or zero if no such entry exists
        value[[c]] <- c(x[x.set], 0)[c(1, head(cumsum(x.set) + 1, n=-1))]
    }
    return(value)
}

# Default parameter set is specified by parametersets/load_by_default.txt
DEFAULT.PARAMETERSET <- local({
    link_fn <- file.path(LOCATION.PARAMETERSETS, "load_by_default.txt")
    default_name <- trimws(readChar(link_fn, file.info(link_fn)$size))
    default <- parameterset.fullpath(default_name)
    if (file.exists(default))
        default_name
    else
        character(0)
})

# Better ggplot2 logscale
my_scale_log10 <- function(scale, ...) scale(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    ...
)

my_scale_log10_pretransformed <- function(scale, ...) scale(
    labels = scales::math_format(10^.x),
    ...
)

#' Hack for deterministic_celltypes to show the correct legend
make_legend_key_twocolor <- function(legend, row, column, color1, color2, pattern="dashed") {
    idx <- grep(paste0("key-", row, "-", column*5 - 3), legend$grobs[[1]]$layout$name)
    for(i in 2:(length(idx)-2)) {
        # Dunno why there are so many overlapping grobs that all draw the same line
        # here, but to avoid visual artifacts we set all their widths to zero.
        legend$grobs[[1]]$grobs[[idx[i]]]$gp$lwd <- 0
    }
    # Make second to last grob a solid line of color2
    legend$grobs[[1]]$grobs[[idx[length(idx)-1]]]$gp$col <- color2
    legend$grobs[[1]]$grobs[[idx[length(idx)-1]]]$gp$lty <- 'solid'
    # Make last grob a dashed line of color1
    legend$grobs[[1]]$grobs[[idx[length(idx)]]]$gp$col <- color1
    legend$grobs[[1]]$grobs[[idx[length(idx)]]]$gp$lty <- pattern
    legend
}

#' Return true if x and y have the same value, including both being NA
is.samevalue <- function(x, y) {
    if (is.atomic(x) && is.atomic(y)) {
        xv <- !is.na(x)
        yv <- !is.na(y)
        return(!((length(x) != length(y)) || any(xv != yv) || any(x[xv] != y[xv])))
    }
    else
        stop("non-atomic values are currently unsupported")
}

#' Return true if the two rates tables are identical
are.rates.identical <- function(r1, r2) {
    if (nrow(r1) != nrow(r2))
        return(FALSE)
    if (any(colnames(r1) != colnames(r2)))
        return(FALSE)
    if (any(unlist(lapply(colnames(r1), function(c) { !is.samevalue(r1[[c]], r2[[c]]) }))))
        return(FALSE)
    return(TRUE)
}

# Server logic
function(input, output, session) {
    # The rates table set programmatically
    # It seems that even though output$rates consumes output_rates(), updating
    # output_rates doesn't always trigger an update of the on-screen table. We
    # thus have a separate trigger which we simply increment whenever output_rates
    # is updated, and which is consumed by output$rates as well.
    output_rates <- reactiveVal(RATES.EMPTY)
    output_rates_trigger <- reactiveVal(0)
    
    # The current rates table, includes updates performed via the UI
    input_rates <- reactiveVal(RATES.EMPTY)
    updateRateInput <- function(session, value, alwaysUpdateOutputAsWell=FALSE) {
        input.update <- function(v) {
            if (!are.rates.identical(isolate(input_rates()), v))
                input_rates(v)
        }
        output.updated <- FALSE
        output.update.maybe.once <- function(v, override=alwaysUpdateOutputAsWell) {
            if (!override || output.updated)
                return()
            message("Updating rates table programatically")
            output_rates(v)
            output_rates_trigger(output_rates_trigger() + 1)
            output.updated <<- TRUE
        }
        if (nrow(value) == 0) {
            input.update(RATES.EMPTY)
            output.update.maybe.once(RATES.EMPTY)
            return()
        }
        if (any(is.na(value$Tmax))) {
            # Ignore any row where Tmax is NA, and if such rows are present, don't reorder
            # Must update the output before removing rows
            output.update.maybe.once(value)
            value <- value[!is.na(Tmax)]
        } else if (is.unsorted(value$Tmax)) {
            # If all rows specify Tmax but they are out of order, reorder them
            # We always update the output in this case, otherwise users would have to reorder manually
            value <- value[order(Tmax)]
            output.update.maybe.once(value, override=TRUE)
        }
        # Update the output before filling in missing values
        output.update.maybe.once(value)
        # Fill empty cells with surrounding values
        value <- parameterset.fillempty(value)
        input.update(value)
    }
    observeEvent(input$rates, {
        if (is.null(input$rates))
            return(RATES.EMPTY)
        message("Rates table was updated via the UI")
        updateRateInput(session, as.data.table(hot_to_r(input$rates)), alwaysUpdateOutputAsWell=FALSE)
    })

    # The final time, which is the Tmax in the last row of the rates table
    Tfinal <- reactive({
        if (nrow(input_rates()) > 0)
            max(input_rates()$Tmax, na.rm=TRUE)
        else
            0
    })
    
    # Rate matrix as a list that can be fed to evaluate_deterministic_san and simulate_san
    input_rates_list <- reactive({
        if (nrow(input_rates()) > 0)
            lapply(split(input_rates(), f=1:nrow(input_rates())), as.list)
        else
            list()
    })

    # Update the possible times for which the LT47 is plotted if Tfinal changes
    observeEvent(Tfinal(), {
        # Update day_lsd
        sel <- input$day_lsd
        choices <- sort(unique(c(LT47$day, Tfinal())))
        if ((length(sel) < 1) || is.na(sel) || (sel == ""))
            sel <- max(LT47$day, na.rm=TRUE)
        if (!(sel %in% choices))
            sel <- choices[length(choices)]
        updateSelectInput(session, "day_lsd", choices=choices, selected=sel)

        # Update day_scellext
        updateSliderInput(session, "day_scellext", min=0, max=Tfinal()+1, value=Tfinal()+1)
    })
    
    san_deterministic_results <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            r <- san_deterministic(s0=input$s0, rates=input_rates(), samples_per_day=10)[, list(
                sid=-1, day=t,
                S, A, N,
                C=S+A+N
            )]
            setkey(r, sid, day)
            r
        } else NULL
    })

    # Simulating the stochastic model is slow, so we cache the result
    # The result is a table with columns day, lid, S, A, N and C=S+A+N
    san_stochastic_results <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            message("Simulating the SAN model")
            r <- san_stochastic(L=input$s0, rates=input_rates(), samples_per_day=1,
                                p_cutoff=as.numeric(input$p_cutoff))[, list(
                sid=-1,
                day=t,
                dt, lid, S, A, N,
                C=S+A+N
            )]
            setkey(r, sid, day, lid)
            r
        } else NULL
    })
    san_stochastic_ranksize <- reactive({
        if (!is.null(san_stochastic_results())) {
            rank_size(san_stochastic_results()[, list(sid, day, lsize=C)])
        } else NULL
    })
    san_stochastic_powerlawfit <- reactive({
        if (!is.null(san_stochastic_ranksize())) {
            ranksize <- san_stochastic_ranksize()[size > 0]
            pareto <- fit_pareto(ranksize)
            limit.alpha <- pareto[day >= LT47.LIMIT.ALPHA.STARTDAY, signif(mean(alpha), digits=2)]
            r <- fit_powerlaw_model(ranksize, alpha=limit.alpha)
            r[, pareto.alpha.sid.day := pareto[r, alpha] ]
            r
        } else NULL
    })

    san_stochastic_extinction_trajectories <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            L=as.integer(input$scellext_nlineages)
            message("Simulating ", L, " lineages under the SAN model to obtain average extinction trajectories")
            r <- san_stochastic(L=L, rates=input_rates(), samples_per_day=1,
                                p_cutoff=as.numeric(input$p_cutoff))[, list(
                day=t,
                dt, lid, S, A, N,
                C=S+A+N
            )]
            r[, Text.S := {
                t.is.ext.S <- day[S == 0]
                if (length(t.is.ext.S) > 0) min(t.is.ext.S) else Inf
            }, by=lid]
            r[, list(
                C=mean(C), C.sd=sd(C),
                S=mean(S), S.sd=sd(S),
                A=mean(A), A.sd=sd(A),
                N=mean(N), N.sd=sd(N),
                nlineages=.N
            ), keyby=.(day, Text.S)]
        } else NULL
    })

    p_cutoff_message <- reactive({
        # Compute expected values for S, A, N using the deterministic model
        mod.exp <- san_deterministic(s0=input$s0, rates=input_rates())
        setkey(mod.exp, "t")
        # Compute the total S, A, N cell counts from the simulation results
        mod.sim <- san_stochastic_results()
        mod.sim.agg <- mod.sim[, list(
            S=sum(S), S.sd=sd(S)*sqrt(.N),
            A=sum(A), A.sd=sd(A)*sqrt(.N),
            N=sum(N), N.sd=sd(N)*sqrt(.N),
            dt=if(any(!is.na(dt))) min(dt, na.rm=TRUE) else NA_real_
        ), keyby=day]
        # Compute the z-scores for the deviation of the total counts from the expectation
        mod.sim.zscore <- mod.sim.agg[mod.exp, list(
            S=(S - i.S) / S.sd,
            A=(A - i.A) / A.sd,
            N=(N - i.N) / N.sd
        )]
        mod.sim.zscore.mat <- as.matrix(mod.sim.zscore)
        max.zscore <- mod.sim.zscore.mat[which.max(abs(mod.sim.zscore.mat))]
        # Compute relative deviations of the total counts compared to the expectation
        mod.sim.relerr <- mod.sim.agg[mod.exp, list(
            S=(S - i.S) / i.S,
            A=(A - i.A) / i.A,
            N=(N - i.N) / i.N
        )]
        mod.sim.relerr.mat <- as.matrix(mod.sim.relerr)
        max.relerr <- mod.sim.relerr.mat[which.max(abs(mod.sim.relerr.mat))]
        # Range of time steps
        min.dt <- mod.sim.agg[, min(dt, na.rm=TRUE)]
        max.dt <- mod.sim.agg[, max(dt, na.rm=TRUE)]
        # Report
        paste0(round(1/max.dt), "-", round(1/min.dt), " steps/day, ",
               "max dev=", signif(100*max.relerr, digits=2), "%, ",
               "max z=", signif(max.zscore, digits=2), "")
    })
    
    gwpcr_parameters_interpolate <- function(p) {
        if (nrow(p[day==Tfinal()]) > 0)
            return(p)
        pext <- rbind(p, p[day <= Tfinal()][day==max(day)][, day := Tfinal()])
        attr(pext, "auto") <- attr(p, "auto")
        pext
    }
    
    # PCR efficiency parameter
    pcr_efficiency_manual_or_auto <- reactive({
        if (!is.na(input$pcr_efficiency) && (input$pcr_efficiency >= 0) && (input$pcr_efficiency <= 1)) {
            eff <- LT47.GWPCR.PARAMETERS[, list(pcr_efficiency=input$pcr_efficiency), keyby=day]
            attr(eff, "auto") <- FALSE
        } else {
            eff <- LT47.GWPCR.PARAMETERS[, list(pcr_efficiency), keyby=day]
            attr(eff, "auto") <- TRUE
        }
        gwpcr_parameters_interpolate(eff)
    })
    pcr_efficiency_auto_message <- reactive({
        eff <- pcr_efficiency_manual_or_auto()
        if (!attr(eff, "auto")) return("")
        paste0("est. from day 0 data, =",eff[day==input$day_lsd, pcr_efficiency], " uniformly")
    })

    # Sequencing library size parameter
    library_size_manual_or_auto <- reactive({
        if (!is.na(input$library_size) && (input$library_size > 0)) {
            ls <- LT47.GWPCR.PARAMETERS[, list(library_size=(input$library_size)), keyby=day]
            attr(ls, "auto") <- FALSE
        } else {
            ls <- LT47.GWPCR.PARAMETERS[, list(library_size), keyby=day]
            attr(ls, "auto") <- TRUE
        }
        gwpcr_parameters_interpolate(ls)
    })
    library_size_auto_message <- reactive({
        ls <- library_size_manual_or_auto()
        if (!attr(ls, "auto")) return("")
        ls_daylsd <- ls[day==input$day_lsd, library_size]
        ls_daylsd_exp <- if (ls_daylsd >= 1) floor(log10(ls_daylsd)) else 0
        ls_daylsd_mant <- signif(ls_daylsd / 10^ls_daylsd_exp, 2)
        paste0("est. from data, =", ls_daylsd_mant,
               if (ls_daylsd_exp > 0) paste0("&middot;10<sup>", ls_daylsd_exp, "</sup>") else "",
               " for day ", input$day_lsd)
    })

    # Sequencing read-count threshold parameter
    phantom_threshold_manual_or_auto <- reactive({
        if (!is.na(input$phantom_threshold) && (input$phantom_threshold > 0)) {
            th <- LT47.GWPCR.PARAMETERS[, list(phantom_threshold=(input$phantom_threshold)), keyby=day]
            attr(th, "auto") <- FALSE
        } else {
            th <-  LT47.GWPCR.PARAMETERS[, list(phantom_threshold), keyby=day]
            attr(th, "auto") <- TRUE
        }
        gwpcr_parameters_interpolate(th)
    })
    phantom_threshold_auto_message <- reactive({
        th <- phantom_threshold_manual_or_auto()
        if (!attr(th, "auto")) return("")
        paste0("est. from data, =", th[day==input$day_lsd, phantom_threshold],
               " for day ", input$day_lsd)
    })
    
    # Simulate PCR+Sequencing
    san_stochastic_results_with_pcr <- reactive({
        if (!is.null(san_stochastic_results())) {
            message("Simulating PCR and sequencing")
            # Days for which experimental data is available
            days <- LT47[, list(), keyby=day]
            # Compute total number of cells in each simulated sample
            total_sizes <- san_stochastic_results()[, list(C.total=sum(as.numeric(C))), keyby=day]
            # Simulate read counts under the gwpcR PCR model
            # The sequencing depth is computed from the total number of cells,
            # and the library size. The resulting table has the additional column
            # "R" containing the simulated read count.
            r <- (san_stochastic_results()
                [total_sizes, on="day"]
                [pcr_efficiency_manual_or_auto(), on="day"]
                [library_size_manual_or_auto(), on="day"]
                [, list(
                    dt, lid, S, A, N,
                    R=if (C > 0)
                        rgwpcrpois(.N, efficiency=pcr_efficiency[1], threshold=0,
                                   molecules=C, lambda0=as.numeric(C)*library_size[1]/C.total[1])
                    else 0
                ), by=.(sid, day, C)])
            setkey(r, sid, day, lid)
            r
        } else NULL
    })
    san_stochastic_ranksize_with_pcr <- reactive({
        if (!is.null(san_stochastic_results_with_pcr())) {
            rank_size(san_stochastic_results_with_pcr()[, list(sid, day, lsize=R)])
        } else NULL
    })
    
    # Simulate read-count thresholding
    san_stochastic_results_with_pcr_filtered <- reactive({
        if (!is.null(san_stochastic_results_with_pcr())) {
            sim <- san_stochastic_results_with_pcr()
            sim[phantom_threshold_manual_or_auto()[sim, phantom_threshold, on="day"] <= R]
        }
    })
    san_stochastic_ranksize_with_pcr_filtered <- reactive({
        if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            rank_size(san_stochastic_results_with_pcr_filtered()[, list(sid, day, lsize=R)])
        } else NULL
    })
    san_stochastic_powerlawfit_with_pcr <- reactive({
        if (!is.null(san_stochastic_ranksize_with_pcr_filtered())) {
            ranksize_pcr <- san_stochastic_ranksize_with_pcr_filtered()
            pareto <- fit_pareto(ranksize_pcr)
            limit.alpha <- pareto[day >= LT47.LIMIT.ALPHA.STARTDAY, signif(mean(alpha), digits=2)]
            r <- fit_powerlaw_model(ranksize_pcr, alpha=limit.alpha)
            r[, pareto.alpha.sid.day := pareto[r, alpha] ]
            r
        } else NULL
    })

    
    # The currently loaded parameterset
    # Note that this reflects the values of the parameters at loading time
    parameterset.loaded <- NULL
    
    # Load parameters and update UI
    observeEvent(input$loadfrom, {
        # Return if no parameter set is selected
        if ((length(input$loadfrom) != 1) || (input$loadfrom == "")) {
            parameterset.loaded <<- NULL
            return()
        }
        # Return if the selected file is already loaded
        if (!is.null(parameterset.loaded) && (parameterset.loaded$name == input$loadfrom))
            return()
        # Load file
        message("Loading parameter set ", input$loadfrom, " from ", parameterset.fullpath(input$loadfrom))
        parameterset.loaded <<- local({
            load(file=parameterset.fullpath(input$loadfrom))
            c(SAN.PARAMS, list(name=input$loadfrom))
        })
        # Update parameter values
        if (!is.null(parameterset.loaded$s0))
            updateNumericInput(session, "s0", value=parameterset.loaded$s0)
        if (!is.null(parameterset.loaded$rates))
            updateRateInput(session, value=parameterset.loaded$rates, alwaysUpdateOutputAsWell=TRUE)
        if (!is.null(parameterset.loaded$pcr_efficiency))
            updateNumericInput(session, "pcr_efficiency", value=parameterset.loaded$pcr_efficiency)
        if (!is.null(parameterset.loaded$library_size))
            updateNumericInput(session, "library_size", value=parameterset.loaded$library_size)
        if (!is.null(parameterset.loaded$phantom_threshold))
            updateNumericInput(session, "phantom_threshold", value=parameterset.loaded$phantom_threshold)
    })
    # Load parameter set "default" on startup
    updateSelectInput(session, "loadfrom", selected=DEFAULT.PARAMETERSET,
                      choices=c("select set to load"="", parametersets.list()))
    
    # If a loaded parameter set is changed, don't show it as selected anymore
    observeEvent(input$s0, {
        if (is.null(parameterset.loaded$s0) || is.samevalue(parameterset.loaded$s0, input$s0))
            return()
        message("Parameter set has been modified (s0 changed), no longer matches ", parameterset.loaded$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
    })
    observeEvent(input$rates, {
        if (is.null(parameterset.loaded$rates) || are.rates.identical(parameterset.loaded$rates, as.data.table(hot_to_r(input$rates))))
            return()
        message("Parameter set has been modified (rates changed), no longer matches ", parameterset.loaded$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
    })
    observeEvent(input$pcr_efficiency, {
        if (is.null(parameterset.loaded$pcr_efficiency) || is.samevalue(parameterset.loaded$pcr_efficiency, input$pcr_efficiency))
            return()
        message("Parameter set has been modified (pcr_efficiency changed), no longer matches ", parameterset.loaded$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
    })
    observeEvent(input$library_size, {
        if (is.null(parameterset.loaded$library_size) || is.samevalue(parameterset.loaded$library_size, input$library_size))
            return()
        message("Parameter set has been modified (library_size changed), no longer matches ", parameterset.loaded$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
    })
    observeEvent(input$phantom_threshold, {
        if (is.null(parameterset.loaded$phantom_threshold) || is.samevalue(parameterset.loaded$phantom_threshold, input$phantom_threshold))
            return()
        message("Parameter set has been modified (phantom_threshold changed), no longer matches ", parameterset.loaded$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
    })

    # Save the current set of parameters
    parameterset.saveas <- function() {
        SAN.PARAMS <- list(rates=as.data.table(hot_to_r(input$rates)),
                           s0=input$s0,
                           pcr_efficiency=input$pcr_efficiency,
                           library_size=input$library_size,
                           phantom_threshold=input$phantom_threshold)
        tryCatch({
            if (!grepl(FILENAME.PATTERN, input$saveas))
                stop("invalid filename")
            save(SAN.PARAMS, file=parameterset.fullpath(input$saveas))
        }, error=function(e) {
            # Show error
            showModal(modalDialog(
                title=paste0("Error saving ", input$saveas),
                conditionMessage(e),
                footer = modalButton("OK")
            ))
        })
        parameterset.loaded <<- c(SAN.PARAMS, list(name=input$saveas))
        # Clear saveas filename after saving
        updateTextInput(session, "saveas", value="")
        # Indicate the just-saved parameter set as being currently loaded
        updateSelectInput(session, "loadfrom", choices=parametersets.list(), selected=input$saveas)
    }
    
    # Save parameters to file
    observeEvent(input$save, {
        if ((length(input$saveas) != 1) || (input$saveas == "") || !grepl(FILENAME.PATTERN, input$saveas)) {
            # Invalid filename
            showModal(modalDialog(
                title=paste0("Invalid filename '", input$saveas, "'"),
                "File names must be non-empty and consist only of letters, numbers, spaces, dots and dashes",
                footer = modalButton("OK")
            ))
        } else if (file.exists(parameterset.fullpath(input$saveas))) {
            # Ask whether to overwrite an existing file
            showModal(modalDialog(
                title=paste0("Overwrite ", input$saveas, "?"),
                paste0(input$saveas, " already exists, shoult it be overwritten?"),
                footer = tagList(actionButton("saveOverwrite", "Overwrite"),
                                 modalButton("Cancel")
                )
            ))
        } else {
            # New file
            parameterset.saveas()
        }
    })
    observeEvent(input$saveOverwrite, {
        removeModal()
        # Overwrite existing file
        message("Saving current parameter set as ", input$saveas, " to ", parameterset.fullpath(input$saveas))
        parameterset.saveas()
    })
    
    # Rate table
    output$rates <- renderRHandsontable({
        message("Rendering rates table")
        # Consume trigger, this allows forcing an update by incrementing output_rates_trigger.
        output_rates_trigger()
        # Construct nice header names
        headers <- colnames(output_rates())
        events <- list(Tmax="", r_S="S &#8594; S S", r_0="S &#8594; &empty;", r_R="S &#8594; N",
                       r_A="S &#8594; A", r_N="A &#8594; A N", r_D="A &#8594; N")[headers]
        headers <- sub("Tmax", "T<sub>max</sub>", headers, fixed=TRUE)
        headers <- sub("^r_([A-Z0-9])$", "g<sub>\\1</sub>", headers)
        headers <- paste0(ifelse(events!="",
                                 paste0("<span style='white-space: nowrap;'>(",
                                        events,
                                        ")</span>"),
                                 ""),
                          "<br/>", headers)
        # Output table
        rhandsontable(output_rates(), useTypes=TRUE, search=FALSE, colHeaders=headers)
    })

    # Deterministic model dynamics
    output$deterministic_cellcounts <- renderPlot({
        message("Rendering deterministic cell counts plot")
        if (is.null(san_deterministic_results()))
            return()
        
        # Evaluate deterministic SAN model
        # We cut cell counts off at 0.1 to avoid plotting problems
        r <- san_deterministic_results()[, list(
            day,
            C=pmax(C, 1),
            S=pmax(S, 1),
            A=pmax(A, 1),
            N=pmax(N, 1)
        )]

        # Show stochastic simulation results for comparison if requested
        r.sim <- if (input$deterministic_cellcounts_incsim && !is.null(san_stochastic_results())) {
            san_stochastic_results()[, list(
                C=pmax(sum(C), 1),
                S=pmax(sum(S), 1),
                A=pmax(sum(A), 1),
                N=pmax(sum(N), 1)
            ), keyby=day]
        } else data.table()

        # Plot results
        p <- ggplot(data=r, aes(x=day)) +
            stat_summary(data=ORGANOIDSIZES, aes(y=count, col='exp.C'), fun=mean, geom="line", linetype="dashed", size=LWD) +
            stat_summary(data=ORGANOIDSIZES, aes(y=count, col='exp.C'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(aes(y=C, col='mod.C', linetype='mod.ana'), size=LWD2) +
            geom_line(aes(y=S, col='mod.S', linetype='mod.ana'), size=LWD) +
            geom_line(aes(y=A, col='mod.A', linetype='mod.ana'), size=LWD) +
            geom_line(aes(y=N, col='mod.N', linetype='mod.ana'), size=LWD) +
            my_scale_log10(scale_y_log10, limits=c(1e2, 1e7), oob=oob_keep) +
            annotation_logticks(sides="l") +
            scale_color_manual(breaks=c('exp.C', 'mod.C', 'mod.S', 'mod.A', 'mod.N'),
                               labels=c('total (data)', 'total (model)',
                                        'S (model)', 'A (model)', 'N (model)'),
                               values=c('maroon', 'black', 'cornflowerblue', 'darkgoldenrod1', 'darkolivegreen4'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of cells") +
            guides(col=guide_legend(override.aes=list(linetype=c('dashed', 'solid', 'solid', 'solid', 'solid'),
                                                      size=c(LWD, LWD2, LWD, LWD, LWD)),
                                    ncol=2, byrow=FALSE),
                   linetype=if (nrow(r.sim) > 0)
                       guide_legend(override.aes=list(size=c(LWD)),
                                    ncol=1, byrow=FALSE)
                   else
                       guide_none())

        if (nrow(r.sim) > 0) {
            p <- p + geom_line(data=r.sim, aes(y=C, col='mod.C', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(y=S, col='mod.S', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(y=A, col='mod.A', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(y=N, col='mod.N', linetype='mod.sim'), size=LWD) +
                scale_linetype_manual(breaks=c('mod.ana', 'mod.sim'),
                                      labels=c('theoretical', 'simulation'),
                                      values=c('solid', 'dotted'),
                                      name=NULL)
        }

        p
    })

    output$deterministic_celltypes <- renderPlot({
        message("Rendering deterministic organoid composition plot ")
        if (is.null(san_deterministic_results()))
            return()
        
        # Evaluate deterministic SAN model
        # We cut cell counts off at 0.1 to avoid plotting problems
        r <- san_deterministic_results()[,  list(
            day,
            S=100 * pmax(S, 1) / C,
            A=100 * pmax(A, 1) / C,
            N=100 * pmax(N, 1) / C
        )]

        CT.CXRC4 <- CELLTYPES[antibody=="CXCR4"]
        CT.NCAM <- CELLTYPES[antibody=="NCAM"]
        CT.TRA160 <- CELLTYPES[antibody=="TRA1-60"]

        # Plot results
        p <- ggplot(data=as.data.table(r), aes(x=day, y=percent)) +
            stat_summary(data=CT.CXRC4, aes(col='exp.CXRC4', linetype='exp.CXRC4'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.CXRC4, aes(col='exp.CXRC4'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            stat_summary(data=CT.NCAM, aes(col='exp.NCAM', linetype='exp.NCAM'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.NCAM, aes(col='exp.NCAM'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            stat_summary(data=CT.TRA160, aes(col='exp.TRA160', linetype='exp.TRA160'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.TRA160, aes(col='exp.TRA160'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(aes(y=S, col='mod.S', linetype='mod.S'), size=LWD) +
            geom_line(aes(y=A, col='mod.A', linetype='mod.A'), size=LWD) +
            geom_line(aes(y=N, col='mod.N', linetype='mod.N'), size=LWD) +
            geom_line(aes(y=S+A, col='mod.S+A'), size=LWD) +
            geom_line(aes(y=S+A, col='mod.A', linetype='mod.S+A'), size=LWD) +
            geom_line(aes(y=A+N, col='mod.A+N'), size=LWD) +
            geom_line(aes(y=A+N, col='mod.N', linetype='mod.A+N'), size=LWD) +
            scale_color_manual(breaks=c('exp.CXRC4', 'exp.NCAM', 'exp.TRA160',
                                        'mod.S', 'mod.A', 'mod.N',
                                        'mod.S+A', 'mod.A+N'),
                               labels=c('CXRC4+ (data)', 'NCAM+ (data)', 'TRA1-60+ (data)',
                                        'S (model)', 'A (model)', 'N (model)',
                                        'S+A (model)', 'A+N (model)'),
                               values=c('brown2', 'darkcyan',  'violet',
                                        'cornflowerblue', 'darkgoldenrod1', 'darkolivegreen4',
                                        'cornflowerblue', 'darkgoldenrod1'),
                               name="") +
            scale_linetype_manual(breaks=c('exp.CXRC4', 'exp.NCAM', 'exp.TRA160',
                                           'mod.S', 'mod.A', 'mod.N',
                                           'mod.S+A', 'mod.A+N'),
                                 values=c('dashed', 'dashed', 'dashed',
                                          'solid', 'solid', 'solid',
                                          'FF', 'FF'),
                                 guide=NULL) +
            xlab("time [days]") +
            ylab("fraction of organoid [%]") +
            guides(col=guide_legend(override.aes=list(linetype=c('dashed', 'dashed', 'dashed',
                                                                 'solid', 'solid', 'solid',
                                                                 'solid', 'solid'),
                                                      size=c(LWD, LWD, LWD,
                                                             LWD, LWD, LWD,
                                                             LWD, LWD)),
                                    ncol=2, byrow=FALSE))

        # Patch legend by rendering the plot as a grob and patching the legens in there
        g <- ggplotGrob(p)
        l <- g$grobs[[which(g$layout$name == "guide-box")]]
        l <- make_legend_key_twocolor(l, 3, 2, 'cornflowerblue', 'darkgoldenrod1', pattern="FF")
        l <- make_legend_key_twocolor(l, 4, 2, 'darkgoldenrod1', 'darkolivegreen4', pattern="FF")
        g$grobs[[which(g$layout$name == "guide-box")]] <- l

        # Plot modified grob
        plot(g)
    })
    
    # Message about simulation accuracy
    output$p_cutoff_message <- renderUI({
        helpText(HTML(p_cutoff_message()))
    })

    # Message about automatically selected parameters 
    output$pcr_efficiency_auto_message <- renderUI({
        helpText(HTML(pcr_efficiency_auto_message()))
    })
    output$library_size_auto_message <- renderUI({
        helpText(HTML(library_size_auto_message()))
    })
    output$phantom_threshold_auto_message <- renderUI({
        helpText(HTML(phantom_threshold_auto_message()))
    })

    # Lineage size distribution
    plot_stochastic_lsd_logrank_loglineagesize <- reactive({
        message("Rendering rank-abundance plot of the lineage size distribution ")

        # Simulate stochastic SAN model
        rank_size.model <- if (input$stochastic_lsd_incpuremodel && !is.null(san_stochastic_ranksize())) {
            san_stochastic_ranksize()[(day==input$day_lsd) & (size > 0)]
        } else data.table()

        # Simulate PCR+Sequencing
        rank_size.model_pcr <- if (!is.null(san_stochastic_ranksize_with_pcr_filtered())) {
            san_stochastic_ranksize_with_pcr_filtered()[day==input$day_lsd]
        } else data.table()
        
        # Do the same (model-free!) powerlaw fit that data/lt47.processed.R applies to exp. data
        powerlaw <- if (input$stochastic_lsd_incpowerlaw && !is.null(san_stochastic_powerlawfit_with_pcr())) {
            san_stochastic_powerlawfit_with_pcr()[day==input$day_lsd]
        } else data.table()

        # Fetch experimental data
        rank_size.experiment <- LT47.RANKSIZE[day==input$day_lsd]
        
        # Abort if no data is available to plot
        if (nrow(rank_size.model) + nrow(rank_size.model_pcr) + nrow(rank_size.experiment) == 0)
            return(NULL)

        # Rank-size normalization and bounds to use when plotting curves
        rank.max <- max(rank_size.experiment$rank, rank_size.model$rank, rank_size.model_pcr$rank)
        size.min <- min(rank_size.experiment$size, rank_size.model$size, rank_size.model_pcr$size)
        size.max <- max(rank_size.experiment$size, rank_size.model$size, rank_size.model_pcr$size)
        xmax <- 10^(ceiling(10*log10(rank.max))/10)
        ymin <- 10^(  floor(log10(size.min)))
        ymax <- 10^(ceiling(log10(size.max)))

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid", `powerlaw`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2, `powerlaw`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10(scale_x_log10, limits=c(1, xmax), oob=oob_keep) +
            my_scale_log10(scale_y_log10, limits=c(ymin, ymax), oob=oob_keep) +
            annotation_logticks(sides="bl") +
            xlab("rank") +
            ylab("lineage size [norm. reads]") +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.', 'powerlaw'),
                               values=c('maroon', 'violet', 'black', 'darkcyan'),
                               name=NULL)
        if (input$stochastic_lsd_incpowerlaw && (nrow(powerlaw) > 0)) {
            model_pcr.rm <- max(rank_size.model_pcr$rank)
            model_pcr.sm <- min(rank_size.model_pcr$size)
            model_pcr.r.scale <- rank.max / model_pcr.rm
            model_pcr.s.scale <- size.min / model_pcr.sm
            p <- p + geom_vline(data=powerlaw, aes(xintercept=(zipf.rank.min-1)*model_pcr.r.scale + 1),
                                col='darkcyan', linetype="dotted", size=LWD)
        }
        
        # Draw rank-size curves for the experimental lineage sizes in all replicates
        plot_rank_size <- function(rank_size, col, ...) {
            align_rank_size(rank_size, rank.max=rank.max, size.min=size.min)[, {
                p <<- p + geom_line(data=copy(.SD), aes(x=rank.aligned, y=size.aligned, col=col), ...)
                NULL
            }, by=sid]
        }
        if (nrow(rank_size.experiment) > 0) {
            groups[length(groups)+1] <- "data"
            plot_rank_size(rank_size.experiment, col="data", size=LWD, linetype="dashed")
        }
        if (input$stochastic_lsd_incpuremodel && (nrow(rank_size.model) > 0)) {
            groups[length(groups)+1] <- "model"
            plot_rank_size(rank_size.model, col="model", size=LWD)
        }
        if (nrow(rank_size.model_pcr) > 0) {
            groups[length(groups)+1] <- "model+seq."
            plot_rank_size(rank_size.model_pcr, col="model+seq.", size=LWD2)
        }
        if (input$stochastic_lsd_incpowerlaw && (nrow(powerlaw) > 0)) {
            groups[length(groups)+1] <- "powerlaw"
            stopifnot(nrow(powerlaw) == 1)
            d <- with(powerlaw, {
                transform <- function(r.scale, s.scale, f) { function(r) { f(r/r.scale)*s.scale } }
                f.nonzipf <- transform(model_pcr.r.scale, model_pcr.s.scale,
                                       function(r) { 10^nonzipf.d * r^nonzipf.k })
                f.zipf <- transform(model_pcr.r.scale, model_pcr.s.scale,
                                    function(r) { 10^zipf.d * r^zipf.k })
                x <- c(1, zipf.rank.min*model_pcr.r.scale, zipf.rank.min*model_pcr.r.scale, model_pcr.rm*model_pcr.r.scale)
                y <- c(f.nonzipf(x[1]), f.nonzipf(x[2]), f.zipf(x[3]), f.zipf(x[4]))
                list(x=x, y=y)
            })
            p <- (p +
                geom_segment(data=powerlaw, aes(x=d$x[1], y=d$y[1], xend=d$x[2], yend=d$y[2], col='powerlaw'),
                              linetype='solid', size=LWD2) +
                geom_segment(data=powerlaw, aes(x=d$x[3], y=d$y[3], xend=d$x[4], yend=d$y[4], col='powerlaw'),
                              linetype='solid', size=LWD2) +
                annotate("text", x=d$x[2], hjust=1.03, y=Inf, vjust=1.5, col='darkcyan',
                         label=" <- roughly uniform size", size=5)) +
                annotate("text", x=d$x[3], hjust=-0.03, y=Inf, vjust=1.5, col='darkcyan',
                         label=" powerlaw tail ->", size=5)
        }
        
        #  If nothing was plottet, return empty
        if (length(groups) == 0)
            return()

        # Finish and return plot
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=3, byrow=TRUE))
    })

    plot_stochastic_lsd_density_loglineagesize <- reactive({
        message("Rendering density plot of the logarithmic lineage size distribution ")
        if (is.null(san_stochastic_results()))
            return()

        # Simulate stochastic SAN model
        log_lsize.model <- if (input$stochastic_lsd_incpuremodel) {
            san_stochastic_results()[(day==input$day_lsd) & (C > 0), list(
                sid, log_lsize=log10(C)
            )]
        } else data.table()

        # Simulate PCR+Sequencing
        log_lsize.model_pcr <-  if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            san_stochastic_results_with_pcr_filtered()[day==input$day_lsd, list(
                sid, log_lsize=log10(R)
            )]
        } else data.table()

        # Do the same (model-free!) powerlaw fit that data/lt47.processed.R applies to exp. data
        log_powerlaw <- if (input$stochastic_lsd_incpowerlaw && (!is.null(san_stochastic_powerlawfit_with_pcr()))) {
            san_stochastic_powerlawfit_with_pcr()[day==input$day_lsd][, list(
                sid, log_zipf.rank.min=log10(zipf.rank.min)
            )]
        } else data.table()

        # Fetch experimental data
        log_lsize.experiment <- LT47[day==input$day_lsd, list(
            sid, log_lsize=log10(lsize)
        )]

        # Size normalization
        log_lsize.min <- min(log_lsize.experiment$log_lsize, log_lsize.model_pcr$log_lsize, log_lsize.model$log_lsize)

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10_pretransformed(scale_x_continuous) +
            annotation_logticks(sides="b") +
            xlab("lineage size [norm. reads]") +
            ylab("density") +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL)
        if (input$stochastic_lsd_incpowerlaw && (nrow(log_powerlaw) > 0)) {
            p <- p + geom_vline(data=log_powerlaw, aes(xintercept=log_zipf.rank.min), col='darkcyan') +
                annotate("text", x=log_powerlaw$log_zipf.rank.min, hjust=-0.03, y=Inf, vjust=1.5, col='darkcyan',
                         label=" powertail tail ->", size=5) +
                annotate("text", x=log_powerlaw$log_zipf.rank.min, hjust=1.03, y=Inf, vjust=1.5, col='darkcyan',
                         label=" <- roughtly uniform size ", size=5)
        }

        # Draw densities
        if (nrow(log_lsize.experiment) > 0) {
            groups[length(groups)+1] <- "data"
            log_lsize.experiment[, {
                p <<- p + stat_density(data=copy(.SD), aes(x=log_lsize - min(log_lsize) + log_lsize.min, col='data'),
                                       geom="line", size=LWD, linetype="dashed")
                NULL
            }, by=sid]
        }
        if (input$stochastic_lsd_incpuremodel && (nrow(log_lsize.model) > 0)) {
            groups[length(groups)+1] <- "model"
            p <- p + stat_density(data=log_lsize.model, aes(x=log_lsize - min(log_lsize) + log_lsize.min, col='model'),
                                  geom="line", size=LWD)
        }
        if (nrow(log_lsize.model_pcr) > 0) {
            groups[length(groups)+1] <- "model+seq."
            p <- p + stat_density(data=log_lsize.model_pcr, aes(x=log_lsize - min(log_lsize) + log_lsize.min, col='model+seq.'),
                                  geom="line", size=LWD2)
        }

        #  If nothing was plotted, return empty
        if (length(groups) == 0)
            return()

        # Finish plot
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=3, byrow=TRUE))
    })
    
    output$stochastic_lsd <- renderPlot({
        switch(input$stochastic_lsd_plottype,
               `log-rank vs. log-lineagesize`=plot_stochastic_lsd_logrank_loglineagesize(),
               `density of log-lineagesize`=plot_stochastic_lsd_density_loglineagesize(),
               stop("unknown plot type ", input$stochastic_lsd_plottype))
    })

    # Download LSD simulation results
    output$stochastic_lsd_download <- downloadHandler(
        filename="san_results.tab.gz",
        contentType="text/tab-separated-values+gzip",
        content=function(file) {
            message("Saving simulation results to ", file)
            file_gz <- gzfile(file, open="wb+")
            write_tsv(san_stochastic_results_with_pcr()[, list(day, lid, S, A, N, C, R)],
                      file=file_gz,
                      col_names=TRUE, na="")
            flush(file_gz)
            close(file_gz)
            return(file)
        }
    )

    # Number of surviving lineages
    output$stochastic_nlineages <- renderPlot({
        message("Rendering number-of-lineages plot ")
        if (is.null(san_stochastic_results()))
            return()
        
        # Compute the numbef of lineages and the number of lineages with extant S-cells
        results <- san_stochastic_results()
        lsd <- results[, list(nlineages=sum(C > 0)), by=day]
        lsd_scells <- results[, list(nlineages=sum(S > 0)), by=day]

        # Compute the number of lineages observed through sequencing
        lsd_pcr <- if (!is.null(san_stochastic_ranksize_with_pcr_filtered())) {
            san_stochastic_ranksize_with_pcr_filtered()[, list(nlineages=.N), by=day]
        } else data.table()

        # Do the same (model-free!) powerlaw fit that data/lt47.processed.R applies to exp. data
        lsd_powerlaw <- if (input$stochastic_nlineages_incpowerlaw && (!is.null(san_stochastic_powerlawfit_with_pcr()))) {
            san_stochastic_powerlawfit_with_pcr()
        } else data.table()
        
        # Legend overrides
        legend.override <- list(linetype=c(`exp.obs`="dashed", `mod.all`="solid",  `mod.obs`="solid",
                                           `mod.S`="solid", `exp.nonzipf`="dashed", `mod.nonzipf`="solid"),
                                size=c(`exp.obs`=LWD, `mod.all`=LWD, `mod.obs`=LWD2,
                                       `mod.S`=LWD, `exp.nonzipf`=LWD, `mod.nonzipf`=LWD2))
        groups <- list()

        # Draw plot
        groups[length(groups)+1] <- "exp.obs"
        groups[length(groups)+1] <- "mod.all"
        p <- ggplot(data=NULL, aes(x=day, y=nlineages)) +
            stat_summary(data=LT47.NLINEAGES, aes(col='exp.obs'), fun=mean, geom="line", size=LWD, linetype="dashed") +
            stat_summary(data=LT47.NLINEAGES, aes(col='exp.obs'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(data=lsd, aes(col='mod.all'), size=LWD)
        if (input$stochastic_nlineages_logy)
            p <- p +
                my_scale_log10(scale_y_log10) +
                annotation_logticks(sides="l")
        if (nrow(lsd_pcr) > 0) {
            groups[length(groups)+1] <- "mod.obs"
            p <- p + geom_line(data=lsd_pcr, aes(col='mod.obs'), size=LWD2)
        }
        if (input$stochastic_nlineages_incpowerlaw && (nrow(lsd_powerlaw) > 0)) {
            groups[length(groups)+1] <- "mod.nonzipf"
            groups[length(groups)+1] <- "exp.nonzipf"
            p <- (p +
                geom_line(data=lsd_powerlaw, aes(y=zipf.rank.min, col='mod.nonzipf'), size=LWD2) +
                stat_summary(data=LT47.POWERLAW, aes(y=zipf.rank.min, col='exp.nonzipf'), fun=mean, geom="line", size=LWD, linetype="dashed") +
                stat_summary(data=LT47.POWERLAW, aes(y=zipf.rank.min, col='exp.nonzipf'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD))
        }

        # Add legend
        groups[length(groups)+1] <- "mod.S"
        groups.all <- c('exp.obs', 'exp.nonzipf', 'mod.obs', 'mod.all', 'mod.S', 'mod.nonzipf')
        groups <- groups.all[groups.all %in% groups]
        p <- p +
            geom_line(data=lsd_scells, aes(col='mod.S'), size=LWD) +
            scale_color_manual(breaks=groups.all,
                               labels=c('obs. (data)', 'outside powerlaw tail (data)', 'obs. (model+seq.)',
                                        'total (model)', 'extant S-cells (model)', 'outside powerlaw tail (model+seq.)'),
                               values=c('maroon', 'darkcyan', 'black', 'violet', 'cornflowerblue', 'darkcyan'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of lineages")

        # Finish plot
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=2, byrow=FALSE))
    })

    output$pareto_alpha <- renderPlot({
        message("Rendering Pareto-alpha plot ")
        if (is.null(san_stochastic_results()))
            return()

        # Compute the alpha estimated from the simulation results
        powerlaw_model <- if (input$pareto_alpha_incpuremodel && !is.null(san_stochastic_powerlawfit())) {
            san_stochastic_powerlawfit()
        } else data.table()

        # Compute the alpha estimated from the simulation results
        powerlaw_model_pcr <- if (!is.null(san_stochastic_powerlawfit_with_pcr())) {
            san_stochastic_powerlawfit_with_pcr()
        } else data.table()

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Draw plot
        groups[length(groups)+1] <- "data"
        p <- ggplot(data=NULL, aes(x=day, y=pareto.alpha.sid.day)) +
            stat_summary(data=LT47.POWERLAW, aes(col='data'), fun=mean, geom="line", size=LWD, linetype="dashed") +
            stat_summary(data=LT47.POWERLAW, aes(col='data'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD)
        if (nrow(powerlaw_model) > 0) {
            groups[length(groups)+1] <- "model"
            p <- p + geom_line(data=powerlaw_model, aes(col='model'), size=LWD2)
        }
        if (nrow(powerlaw_model_pcr) > 0) {
            groups[length(groups)+1] <- "model+seq."
            p <- p + geom_line(data=powerlaw_model_pcr, aes(col='model+seq.'), size=LWD2)
        }

        # Add legend
        groups.all <- c('data', 'model', 'model+seq.')
        groups <- groups.all[groups.all %in% groups]
        p <- p +
            scale_color_manual(breaks=groups.all,
                               labels=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL) +
            xlab("time [days]") +
            ylab(expression("Pareto equality index (" ~ alpha ~ ")"))

        # Finish plot
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=2, byrow=FALSE))
    })

    # S-cell extinction time vs. lineage size
    output$stochastic_scellext_vs_lineagesize <- renderPlot({
        message("Rendering S-cell extinction time vs. lineage size plot")
        if (is.null(san_stochastic_extinction_trajectories()))
            return()

        ext_traj <- san_stochastic_extinction_trajectories()
        ymax <- ext_traj[, 10^ceiling(log10(max(C)))]

        day.max <- ext_traj[, max(day)]
        ggplot(data=ext_traj[day==day.max],
               aes(x=ifelse(is.finite(Text.S), Text.S, day.max+1),
                   y=pmax(C, 0.1),
                   ymin=pmax(C-C.sd, 0.1),
                   ymax=pmax(C+C.sd, 0.1))) +
            geom_ribbon(alpha=0.15) +
            geom_line(size=LWD2) +
            lims(x=c(0, Tfinal()+1)) +
            my_scale_log10(scale_y_log10, limits=c(1e0, ymax), oob=oob_keep) +
            annotation_logticks(sides="l") +
            xlab("S-cell extinction time [days]") +
            ylab(paste0("lineage size [cells]"))
    })

    # S-cell-containing lineages
    output$stochastic_scellext_vs_time <- renderPlot({
        message("Rendering surving fraction of lineages plot")
        if (is.null(san_stochastic_extinction_trajectories()))
            return()

        ext_traj <- san_stochastic_extinction_trajectories()
        data <- rbind(ext_traj[, list(survived=sum(nlineages[Text.S > day]) / sum(nlineages))
                               , by=day],
                      data.table(day=Tfinal()+1,
                                 survived=ext_traj[day==Tfinal(), sum(nlineages[is.infinite(Text.S)]) / sum(nlineages)]))
        ymin <- data[, min(survived)]
        ggplot(data=data, aes(x=day, y=100*survived)) +
            geom_line(col='cornflowerblue') +
            lims(x=c(0, Tfinal()+1)) +
            scale_y_log10(limits=c(100*ymin, 100), oob=oob_keep) +
            xlab("time [days]") +
            ylab("lineages [%]")
    })

    # S-cell extinction trajectories
    stochastic_scellext_trajectory_plot <- function(data, ymax) {
        ggplot(data=data, aes(x=day)) +
            geom_ribbon(aes(ymin=pmax(C-C.sd, 0.1), ymax=pmax(C+C.sd, 0.1), fill='mod.C'), alpha=0.2) +
            geom_ribbon(aes(ymin=pmax(S-S.sd, 0.1), ymax=pmax(S+S.sd, 0.1), fill='mod.S'), alpha=0.15) +
            geom_ribbon(aes(ymin=pmax(A-A.sd, 0.1), ymax=pmax(A+A.sd, 0.1), fill='mod.A'), alpha=0.15) +
            geom_ribbon(aes(ymin=pmax(N-N.sd, 0.1), ymax=pmax(N+N.sd, 0.1), fill='mod.N'), alpha=0.15) +
            geom_line(aes(y=pmax(C, 0.1), col='mod.C'), size=LWD2) +
            geom_line(aes(y=pmax(S, 0.1), col='mod.S'), size=LWD) +
            geom_line(aes(y=pmax(A, 0.1), col='mod.A'), size=LWD) +
            geom_line(aes(y=pmax(N, 0.1), col='mod.N'), size=LWD) +
            my_scale_log10(scale_y_log10, limits=c(1e0, ymax), oob=oob_keep) +
            annotation_logticks(sides="l") +
            scale_color_manual(breaks=c('mod.C', 'mod.S', 'mod.A', 'mod.N'),
                               labels=c('total (model)', 'S (model)', 'A (model)', 'N (model)'),
                               values=c('black', 'cornflowerblue', 'darkgoldenrod1', 'darkolivegreen4'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("cells per lineage") +
            guides(col=guide_legend(override.aes=list(linetype=c('solid', 'solid', 'solid', 'solid'),
                                                      size=c(LWD2, LWD, LWD, LWD)),
                                    ncol=2, byrow=FALSE),
                   fill=guide_none())
    }
    stochastic_scellext_trajectory_label <- reactiveVal()
    output$stochastic_scellext_trajectory_label <- renderUI({
        HTML(stochastic_scellext_trajectory_label())
    })
    output$stochastic_scellext_trajectory <- renderPlot({
        message("Rendering S-cell trajectory plot")
        if (is.null(san_stochastic_extinction_trajectories()))
            return()

        # Fetch average trajectories
        ext_traj <- san_stochastic_extinction_trajectories()
        ymax <- ext_traj[, 10^ceiling(log10(max(C)))]

        # Restrict to those with the selected S-cell extinctin time
        if (input$day_scellext <= Tfinal()) {
            ext_traj_extday <- ext_traj[Text.S == input$day_scellext]
            stochastic_scellext_trajectory_label(paste0("SAN trajectories assuming S-cell extinction on day ", input$day_scellext))
            if (nrow(ext_traj_extday) < 1)
                stop("simulation produced no lineages with S-cell extinction on day ", input$day_scellext)
        } else {
            ext_traj_extday <- san_stochastic_extinction_trajectories()[is.infinite(Text.S)]
            stochastic_scellext_trajectory_label(paste0("SAN trajectories assuming S-cell extinction <b>after</b> day ", Tfinal()))
            if (nrow(ext_traj_extday) < 1)
                stop("simulation produced no lineages with S-cells surving past day ", Tfinal())
        }

        # Plot results
        stochastic_scellext_trajectory_plot(ext_traj_extday, ymax)
    })

    # S-cell extinction trajectory video
    stochastic_scellext_trajectory_video_file <- reactiveVal(NULL)
    observeEvent({
        # All inputs that the video depends on
        input_rates()
        input$s0
        input$scellext_nlineages
        input$stochastic_scellext_trajectory_video_fps
    }, {
        # Remove outdated video if the simulated trajectories change
        if (!is.null(stochastic_scellext_trajectory_video_file())) {
            file.remove(stochastic_scellext_trajectory_video_file())
            stochastic_scellext_trajectory_video_file(NULL)
        }
    })
    observeEvent(input$stochastic_scellext_trajectory_video_generate, {
        # Generate video upon explicit request and only if its outdated (i.e. NULL)
        if (is.null(san_stochastic_extinction_trajectories()) ||
            !is.null(stochastic_scellext_trajectory_video_file()))
            return()
        file <- tempfile(fileext=".mp4")
        onStop(function() { file.remove(file) })
        message("Rendering S-cell trajectory video to temporary file ", file)

        # Fetch trajectories and setup plot
        ext_traj <- san_stochastic_extinction_trajectories()
        ymax <- ext_traj[, 10^ceiling(log10(max(C)))]
        p <- stochastic_scellext_trajectory_plot(ext_traj, ymax) +
            transition_states(Text.S, transition_length=1, state_length=0, wrap=FALSE) +
            ease_aes(y="linear")

        # Create animation
        withProgress(min=0, max=1, value=0, {
            progress.last <- 0
            progress.pattern <- "^.*\\[(=*)(>-*)\\].*"
            incProgress(0, message="preparing animation")

            # Render animation as h264-encoded mp4 file
            withCallingHandlers({
                animate(p, duration=Tfinal(), fps=as.numeric(input$stochastic_scellext_trajectory_video_fps), rewind=FALSE,
                        width=9, height=6, units="in", res=200,
                        renderer=av_renderer(codec="libx264", file=file))
            }, message=function(m) {
                # Update shiny progress indicator if gganimate::animate outputs a progress message
                # Note: These messages are output only if r-lib/progres believes that we're in some
                # kind of interactive session, which is the case if we're running in RStudio *or* if
                # stderr is a tty.
                if (grepl(progress.pattern, m)) {
                    # Message looks like a progress message from gganimate::animate
                    progress.done <- nchar(gsub(progress.pattern, "\\1", m))
                    progress.todo <-  nchar(gsub(progress.pattern, "\\2", m))
                    progress <- progress.done / (progress.done + progress.todo)
                    incProgress(progress - progress.last,
                                message=if(progress < 1) "rendering animation" else "encoding video")
                    progress.last <<- progress
                }
            }, error=function(e) {
                message("Removing incomplete file ", file)
                file.remove(file)
                file <<- NULL
            })
        })

        # Update video file
        stochastic_scellext_trajectory_video_file(file)
    })
    observeEvent(stochastic_scellext_trajectory_video_file(), ignoreNULL=FALSE, {
        # Toggle video generate/download button visibility depending on whether the file is up-to-date
        if (is.null(stochastic_scellext_trajectory_video_file())) {
            message("Enabling S-cell trajectory video generation")
            show(id="stochastic_scellext_trajectory_video_generate")
            hide(id="stochastic_scellext_trajectory_video_download")
        } else {
            message("Enabling S-cell trajectory video download")
            hide(id="stochastic_scellext_trajectory_video_generate")
            show(id="stochastic_scellext_trajectory_video_download")
        }
    })
    output$stochastic_scellext_trajectory_video_download <- downloadHandler(
        filename="scell_extinction.mp4",
        contentType="video/mp4",
        content=function(file) { file.copy(stochastic_scellext_trajectory_video_file(), file) }
    )

}
