library(data.table)
library(Hmisc)
library(gwpcR)
library(SANsimulatoR)
library(ggplot2)
library(scales)
library(shiny)
library(rhandsontable)

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

# Hake for deterministic_celltypes to show the correct legend
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

#' Compute rank-size table from a table with observed lineage sizes
#'
#' @param subset a `data.table` containing columns `sid` (sample id) and `lsize` (lineage size)
#' @return a `data.table`  with columns `sid` (sample id), `rank` (lineage size rank) and `size` (lineage size)
rank_size <- function(subset) {
    subset[, {
        .SD[order(lsize, decreasing=TRUE)][, list(rank=1:.N, size=lsize)]
    }, by="sid"]
}

#' Align rank-size curves
align_rank_size <- function(ranked) {
    rank.max.all <- max(ranked$rank)
    size.min.all <- min(ranked$size)
    ranked[, {
        rank.scale <- rank.max.all / (max(rank)-1)
        size.scale <- size.min.all / min(size)
        sm <- min(size)
        list(rank, size,
             rank.aligned=(rank-1) * rank.scale + 1,
             size.aligned=size * size.scale,
             rank.scale, size.scale)
    }, keyby="sid"]
}

# Return true if x and y have the same value, including both being NA
is.samevalue <- function(x, y) {
    if (is.atomic(x) && is.atomic(y)) {
        xv <- !is.na(x)
        yv <- !is.na(y)
        return(!((length(x) != length(y)) || any(xv != yv) || any(x[xv] != y[xv])))
    }
    else
        stop("non-atomic values are currently unsupported")
}

# Return true if the two rates tables are identical
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
            r <- san_deterministic(s0=input$s0, rates=input_rates(), samples_per_day=10)
            r[, C := S + A + N]
            r
        } else NULL
    })

    # Simulating the stochastic model is slow, so we cache the result
    # The result is a table with columns t, S, A, N and C, where C
    # is the total organoid size, i.e. S+A+N.
    san_stochastic_results <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            message("Simulating the SAN model")
            r <- san_stochastic(L=input$s0, rates=input_rates(), samples_per_day=1,
                                p_cutoff=as.numeric(input$p_cutoff))
            r[, C := S + A + N]
            r
        } else NULL
    })
    
    san_stochastic_extinction_trajectories <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            L=as.integer(input$scellext_nlineages)
            message("Simulating ", L, " lineages under the SAN model to obtain average extinction trajectories")
            r <- san_stochastic(L=L, rates=input_rates(), samples_per_day=1,
                                p_cutoff=as.numeric(input$p_cutoff))
            r[, C := S + A + N]
            r[, Text.S := {
                t.is.ext.S <- t[S == 0]
                if (length(t.is.ext.S) > 0) min(t.is.ext.S) else Inf
            }, by="lid"][, list(
                C=mean(C), C.sd=sd(C),
                S=mean(S), S.sd=sd(S),
                A=mean(A), A.sd=sd(A),
                N=mean(N), N.sd=sd(N),
                nlineages=.N
            ), keyby=c("t", "Text.S")]
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
        ), keyby="t"]
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
            eff <- LT47.GWPCR.PARAMETERS[, list(pcr_efficiency=input$pcr_efficiency), keyby="day"]
            attr(eff, "auto") <- FALSE
        } else {
            eff <- LT47.GWPCR.PARAMETERS[, list(pcr_efficiency), keyby="day"]
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
            ls <- LT47.GWPCR.PARAMETERS[, list(library_size=(input$library_size)), keyby="day"]
            attr(ls, "auto") <- FALSE
        } else {
            ls <- LT47.GWPCR.PARAMETERS[, list(library_size), keyby="day"]
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
            th <- LT47.GWPCR.PARAMETERS[, list(phantom_threshold=(input$phantom_threshold)), keyby="day"]
            attr(th, "auto") <- FALSE
        } else {
            th <-  LT47.GWPCR.PARAMETERS[, list(phantom_threshold), keyby="day"]
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
            days <- LT47[, list(), keyby="day"]
            # Compute total number of cells in each simulated sample
            total_sizes <- san_stochastic_results()[, list(C.total=sum(as.numeric(C))), keyby="t"]
            # Simulate read counts under the gwpcR PCR model
            # The sequencing depth is computed from the total number of cells,
            # and the library size. The resulting table has the additional column
            # "R" containing the simulated read count.
            (san_stochastic_results()
                [total_sizes, on="t"]
                [pcr_efficiency_manual_or_auto(), on=c(t="day")]
                [library_size_manual_or_auto(), on=c(t="day")]
                [, list(
                    t, dt, lid, S, A, N, C,
                    R=if (C > 0)
                        rgwpcrpois(.N, efficiency=pcr_efficiency[1], threshold=0,
                                   molecules=C, lambda0=as.numeric(C)*library_size[1]/C.total[1])
                    else 0
                ), by=c("t", "C")])
        } else NULL
    })
    
    # Simulate read-count thresholding
    san_stochastic_results_with_pcr_filtered <- reactive({
        if (!is.null(san_stochastic_results_with_pcr())) {
            sim <- san_stochastic_results_with_pcr()
            sim[phantom_threshold_manual_or_auto()[sim, phantom_threshold, on=c(day="t")] <= R]
        }
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
            t,
            C=pmax(C, 1),
            S=pmax(S, 1),
            A=pmax(A, 1),
            N=pmax(N, 1)
        )]

        # Plot results
        p <- ggplot(data=r) +
            stat_summary(data=ORGANOIDSIZES, aes(x=day, y=count, col='exp.C'), fun=mean, geom="line", linetype="dashed", size=LWD) +
            stat_summary(data=ORGANOIDSIZES, aes(x=day, y=count, col='exp.C'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(aes(x=t, y=C, col='mod.C', linetype='mod.ana'), size=LWD2) +
            geom_line(aes(x=t, y=S, col='mod.S', linetype='mod.ana'), size=LWD) +
            geom_line(aes(x=t, y=A, col='mod.A', linetype='mod.ana'), size=LWD) +
            geom_line(aes(x=t, y=N, col='mod.N', linetype='mod.ana'), size=LWD) +
            my_scale_log10(scale_y_log10, limits=c(1e2, 1e7), oob=oob_keep) +
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
                   linetype=if (input$deterministic_cellcounts_incsim && !is.null(san_stochastic_results()))
                       guide_legend(override.aes=list(size=c(LWD)),
                                    ncol=1, byrow=FALSE)
                   else
                       guide_none())

        if (input$deterministic_cellcounts_incsim && !is.null(san_stochastic_results())) {
            r.sim <- san_stochastic_results()[, list(
                C=pmax(sum(C), 1),
                S=pmax(sum(S), 1),
                A=pmax(sum(A), 1),
                N=pmax(sum(N), 1)
            ), keyby="t"]
            p <- p + geom_line(data=r.sim, aes(x=t, y=C, col='mod.C', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(x=t, y=S, col='mod.S', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(x=t, y=A, col='mod.A', linetype='mod.sim'), size=LWD) +
                geom_line(data=r.sim, aes(x=t, y=N, col='mod.N', linetype='mod.sim'), size=LWD) +
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
            t,
            S=100 * pmax(S, 1) / C,
            A=100 * pmax(A, 1) / C,
            N=100 * pmax(N, 1) / C
        )]

        CT.CXRC4 <- CELLTYPES[antibody=="CXCR4"]
        CT.NCAM <- CELLTYPES[antibody=="NCAM"]
        CT.TRA160 <- CELLTYPES[antibody=="TRA1-60"]

        # Plot results
        #            my_scale_log10(scale_y_log10, limits=c(0.1, 100), oob=oob_keep) +
        p <- ggplot(data=as.data.table(r)) +
            stat_summary(data=CT.CXRC4, aes(x=day, y=percent, col='exp.CXRC4', linetype='exp.CXRC4'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.CXRC4, aes(x=day, y=percent, col='exp.CXRC4'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            stat_summary(data=CT.NCAM, aes(x=day, y=percent, col='exp.NCAM', linetype='exp.NCAM'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.NCAM, aes(x=day, y=percent, col='exp.NCAM'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            stat_summary(data=CT.TRA160, aes(x=day, y=percent, col='exp.TRA160', linetype='exp.TRA160'), fun=mean, geom="line", size=LWD) +
            stat_summary(data=CT.TRA160, aes(x=day, y=percent, col='exp.TRA160'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(aes(x=t, y=S, col='mod.S', linetype='mod.S'), size=LWD) +
            geom_line(aes(x=t, y=A, col='mod.A', linetype='mod.A'), size=LWD) +
            geom_line(aes(x=t, y=N, col='mod.N', linetype='mod.N'), size=LWD) +
            geom_line(aes(x=t, y=S+A, col='mod.S+A'), size=LWD) +
            geom_line(aes(x=t, y=S+A, col='mod.A', linetype='mod.S+A'), size=LWD) +
            geom_line(aes(x=t, y=A+N, col='mod.A+N'), size=LWD) +
            geom_line(aes(x=t, y=A+N, col='mod.N', linetype='mod.A+N'), size=LWD) +
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
        if (is.null(san_stochastic_results()))
            return()

        # Simulate stochastic SAN model
        rank_size.model <- if (input$stochastic_lsd_incpuremodel) {
            rank_size(san_stochastic_results()[t==input$day_lsd, list(sid=-1, lsize=C)][lsize > 0])
        } else data.table()

        # Simulate PCR+Sequencing
        rank_size.model_pcr <- if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            rank_size(san_stochastic_results_with_pcr_filtered()[t==input$day_lsd, list(sid=-1, lsize=R)])
        } else data.table()

        
        # Fetch experimental data
        rank_size.experiment <- LT47.RANKSIZE[day==input$day_lsd]
        
        # X coordinates to use when plotting curves
        curves.xmax <- max(rank_size.experiment$rank, rank_size.model$rank, rank_size.model_pcr$rank)
        curves.ymin <- min(rank_size.experiment$size, rank_size.model$size, rank_size.model_pcr$size)
        curves.ymax <- max(rank_size.experiment$size, rank_size.model$size, rank_size.model_pcr$size)
        xmax <- 10^(ceiling(10*log10(curves.xmax))/10)
        ymin <- 10^(  floor(log10(curves.ymin)))
        ymax <- 10^(ceiling(log10(curves.ymax)))

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10(scale_x_log10, limits=c(1, xmax)) +
            xlab("rank") +
            my_scale_log10(scale_y_log10, limits=c(ymin, ymax)) +
            ylab("lineage size [norm. reads]") +
            annotation_logticks() +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL)
        
        # Draw rank-size curves for the experimental lineage sizes in all replicates
        plot_rank_size <- function(rank_size, col, ...) {
            rank_size[, {
                # Plot a single curve, align the rightmost datapoints of all curves
                rm <- max(rank)
                sm <- min(size)
                curve.aligned <- .SD[, list(sid, rank.shift=(curves.xmax*rank/rm),
                                            size.norm=curves.ymin*size/sm)]
                p <<- p + geom_line(data=curve.aligned, aes(x=rank.shift, y=size.norm, col=col), ...)
                NULL
            }, by="sid"]
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
        
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=3, byrow=TRUE))
    })

    plot_stochastic_lsd_density_loglineagesize <- reactive({
        message("Rendering density plot of the logarithmic lineage size distribution ")
        if (is.null(san_stochastic_results()))
            return()

        # Simulate stochastic SAN model
        log_lsize.model <- if (input$stochastic_lsd_incpuremodel) {
            san_stochastic_results()[
                (t==input$day_lsd) & (C > 0),
                list(sid=-1, log_lsize=log10(C))]
        } else NULL

        # Fetch experimental data
        log_lsize.experiment <- LT47[
            day==input$day_lsd,
            list(sid, log_lsize=log10(lsize))]

        log_lsize.model_pcr <-  if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            san_stochastic_results_with_pcr_filtered()[
                (t==input$day_lsd) & (R > 0),
                list(sid=-1, log_lsize=log10(R))]
        } else NULL

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10_pretransformed(scale_x_continuous) +
            xlab("lineage size [norm. reads]") +
            ylab("density") +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL)

        # Draw densities
        if (nrow(log_lsize.experiment) > 0) {
            groups[length(groups)+1] <- "data"
            log_lsize.experiment[, {
                p <<- p + stat_density(data=copy(.SD), aes(x=log_lsize, col='data'),
                                       geom="line", size=LWD, linetype="dashed")
                NULL
            }, by="sid"]
        }
        if (input$stochastic_lsd_incpuremodel && !is.null(log_lsize.model)) {
            groups[length(groups)+1] <- "model"
            p <- p + stat_density(data=log_lsize.model, aes(x=log_lsize, col='model'),
                                  geom="line", size=LWD)
        }
        if (!is.null(log_lsize.model_pcr)) {
            groups[length(groups)+1] <- "model+seq."
            p <- p + stat_density(data=log_lsize.model_pcr, aes(x=log_lsize, col='model+seq.'),
                                  geom="line", size=LWD2)
        }

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

    # Number of surviving lineages
    output$stochastic_nlineages <- renderPlot({
        message("Rendering number-of-lineages plot ")
        if (is.null(san_stochastic_results()))
            return()
        
        lsd <- san_stochastic_results()[, list(nlineages=sum(C > 0)), by="t"]
        lsd_scells <- san_stochastic_results()[, list(nlineages=sum(S > 0)), by="t"]
        
        # Setup plot
        p <- ggplot() +
            stat_summary(data=LT47.NLINEAGES, aes(x=day, y=nlineages, col='exp.obs'), fun=mean, geom="line", size=LWD, linetype="dashed") +
            stat_summary(data=LT47.NLINEAGES, aes(x=day, y=nlineages, col='exp.obs'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(data=lsd, aes(x=t, y=nlineages, col='mod.all'), size=LWD)
        if (input$stochastic_nlineages_logy)
            p <- p + my_scale_log10(scale_y_log10)

        # Draw lines
        if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            lsd_pcr <- san_stochastic_results_with_pcr_filtered()[, list(nlineages=.N), by="t"]
            p <- p + geom_line(data=lsd_pcr, aes(x=t, y=nlineages, col='mod.obs'), size=LWD2)
        }
        p <- p +
            geom_line(data=lsd_scells, aes(x=t, y=nlineages, col='mod.S'), size=LWD) +
            scale_color_manual(breaks=c('exp.obs', 'mod.obs', 'mod.all', 'mod.S'),
                               labels=c('obs. (data)', 'obs. (model+seq.)',
                                        'total (model)', 'w/ S-cells (model)'),
                               values=c('maroon', 'black', 'violet', 'cornflowerblue'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of lineages") +
            guides(col=guide_legend(override.aes=list(linetype=c("dashed", "solid", "solid", "solid"),
                                                      size=c(LWD, LWD2, LWD, LWD)),
                                    ncol=2, byrow=FALSE))

        # Finish plot
        p
    })

    # S-cell extinction dynamics
    output$stochastic_scellext_vs_lineagesize <- renderPlot({
        ymax <- san_stochastic_extinction_trajectories()[, 10^ceiling(log10(max(C)))]

                t.max <- san_stochastic_extinction_trajectories()[, max(t)]
        ggplot(data=san_stochastic_extinction_trajectories()[t==t.max],
               aes(x=ifelse(is.finite(Text.S), Text.S, max(t)+1), y=pmax(C, 0.1), ymin=pmax(C-C.sd, 0.1), ymax=pmax(C+C.sd, 0.1))) +
            geom_ribbon(alpha=0.15) +
            geom_line(size=LWD2) +
            lims(x=c(0, Tfinal()+1)) +
            my_scale_log10(scale_y_log10, limits=c(1e0, ymax), oob=oob_keep) +
            xlab("S-cell extinction time [days]") +
            ylab(paste0("lineage size [cells]"))
    })

    output$stochastic_scellext_vs_time <- renderPlot({
        ext_traj <- san_stochastic_extinction_trajectories()
        data <- rbind(ext_traj[, list(survived=sum(nlineages[Text.S > t]) / sum(nlineages)), by="t"],
                      data.table(t=Tfinal()+1, survived=ext_traj[t==Tfinal(), sum(nlineages[is.infinite(Text.S)]) / sum(nlineages)]))
        ymin <- data[, min(survived)]
        ggplot(data=data, aes(x=t, y=100*survived)) +
            geom_line(col='cornflowerblue') +
            lims(x=c(0, Tfinal()+1)) +
            scale_y_log10(limits=c(100*ymin, 100), oob=oob_keep) +
            xlab("time [days]") +
            ylab("lineages [%]")
    })

    stochastic_scellext_trajectory_label <- reactiveVal()
    output$stochastic_scellext_trajectory_label <- renderUI({
        HTML(stochastic_scellext_trajectory_label())
    })
    output$stochastic_scellext_trajectory <- renderPlot({
        ymax <- san_stochastic_extinction_trajectories()[, 10^ceiling(log10(max(C)))]

        # Plot results
        if (input$day_scellext <= Tfinal()) {
            data <- san_stochastic_extinction_trajectories()[Text.S == input$day_scellext]
            stochastic_scellext_trajectory_label(paste0("SAN trajectories assuming S-cell extinction on day ", input$day_scellext))
            if (nrow(data) < 1)
                stop("simulation produced no lineages with S-cell extinction on day ", input$day_scellext)
        } else {
            data <- san_stochastic_extinction_trajectories()[is.infinite(Text.S)]
            stochastic_scellext_trajectory_label(paste0("SAN trajectories assuming S-cell extinction <b>after</b> day ", Tfinal()))
            if (nrow(data) < 1)
                stop("simulation produced no lineages with S-cells surving past day ", Tfinal())
        }

        ggplot(data=data, aes(x=t)) +
            geom_ribbon(aes(ymin=pmax(C-C.sd, 0.1), ymax=pmax(C+C.sd, 0.1), fill='mod.C'), alpha=0.2) +
            geom_ribbon(aes(ymin=pmax(S-S.sd, 0.1), ymax=pmax(S+S.sd, 0.1), fill='mod.S'), alpha=0.15) +
            geom_ribbon(aes(ymin=pmax(A-A.sd, 0.1), ymax=pmax(A+A.sd, 0.1), fill='mod.A'), alpha=0.15) +
            geom_ribbon(aes(ymin=pmax(N-N.sd, 0.1), ymax=pmax(N+N.sd, 0.1), fill='mod.N'), alpha=0.15) +
            geom_line(aes(y=pmax(C, 0.1), col='mod.C'), size=LWD2) +
            geom_line(aes(y=pmax(S, 0.1), col='mod.S'), size=LWD) +
            geom_line(aes(y=pmax(A, 0.1), col='mod.A'), size=LWD) +
            geom_line(aes(y=pmax(N, 0.1), col='mod.N'), size=LWD) +
            my_scale_log10(scale_y_log10, limits=c(1e0, ymax), oob=oob_keep) +
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
    })
}
