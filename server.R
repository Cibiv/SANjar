library(data.table)
library(gwpcR)
library(ggplot2)
library(shiny)
library(rhandsontable)

source("san_deterministic.R")
#source("san_simulate.R")

# Pattern defining valid "saveas" filenames for parameter sets
FILENAME.PATTERN <- "^[A-Za-z0-9. -]*$"

# Location to store the parameter sets in
LOCATION.PARAMETERSETS <- "parametersets"

# Better ggplot2 logscale
my_scale_log10 <- function(scale, ...) scale(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    ...
)

# Transform LT47 to rank-size form
rank_size <- function(subset) {
    subset[, {
        .SD[order(lsize, decreasing=TRUE)][, list(rank=1:.N, size=lsize)]
    }, by="sid"]
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

# Empty rates table
RATES.EMPTY <- as.data.table({
    x <- rep(list(numeric(0)), 1+length(SAN.RATENAMES))
    names(x) <- c("Tmax", SAN.RATENAMES)
    x
})

# Load cell count data
ORGANOIDSIZES <- data.table(read.table("data/organoidsizes.tab", header=TRUE))

# Load lineage size data
load("data/lt47.rd")

# Compute number of observed lineages at each time point
NLINEAGES <- LT47[, list(nlineages=sum(lsize > 0)), by=c("day", "sid")] 

# Estimate PCR and sequencing parameters
GWPCR.PARAMETERS <- LT47[,
    .SD[, list(
        library_size=sum(lsize),
        phantom_threshold=min(lsize)
    ), by="sid"][, list(
        library_size=signif(round(median(library_size)), 2),
        phantom_threshold=round(median(phantom_threshold))
    )]
, by="day"]
GWPCR.PARAMETERS[, pcr_efficiency := {
    LT47[day==0][, {
        gwpcrpois.est(lsize, threshold=min(lsize), molecules=1)
    }, by="sid"][, list(
        pcr_efficiency=signif(median(efficiency), 2)
    )]
}]

# ggplot2 setup
theme_set(
    theme_grey(base_size = 20)
)

# Server logic
function(input, output, session) {
    # The rates table set programmatically
    output_rates <- reactiveVal(RATES.EMPTY)
    
    # The current rates table, includes updates performed via the UI
    input_rates <- reactive({
        if (is.null(input$rates))
            return(RATES.EMPTY)
        message("Rates table was updated via the UI")
        r <- as.data.table(hot_to_r(input$rates))
        if (nrow(r) == 0)
            return(RATES.EMPTY)
        if (any(is.na(r$Tmax))) {
            # Ignore any row where Tmax is NA, and if such rows are present, don't reorder
            r <- r[!is.na(Tmax)]
        } else if (is.unsorted(r$Tmax)) {
            # If all rows specify Tmax but they are out of order, reorder them
            r <- r[order(Tmax)]
            output_rates(r)
        }
        # Fill empty cells with surrounding values
        for(c in SAN.RATENAMES) {
            # Fetch row corresponding to rate r
            v <- r[[c]]
            v.set <- !is.na(v)
            # Set empty entries to the next non-empty entry or zero if no such entry exists
            r[[c]] <- c(v[v.set], 0)[c(1, head(cumsum(v.set) + 1, n=-1))]
        }
        r
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
        sel <- input$day_lsd
        choices <- sort(unique(c(LT47$day, Tfinal())))
        if ((length(sel) < 1) || is.na(sel) || (sel == ""))
            sel <- max(LT47$day, na.rm=TRUE)
        if (!(sel %in% choices))
            sel <- choices[length(choices)]
        updateSelectInput(session, "day_lsd", choices=choices, selected=sel)
    })
    
    # Simulating the stochastic model is slow, so we cache the result
    # The result is a table with columns t, S, A, N and C, where C
    # is the total organoid size, i.e. S+A+N.
    san_simulation <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (input$s0 > 0)) {
            message("Simulating the SAN model")
            r <- simulate_san(L=input$s0, rates=input_rates_list(), samples_per_day=1,
                              p_cutoff=input$p_cutoff)
            r[, C := S + A + N]
            r
        } else NULL
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
            eff <- GWPCR.PARAMETERS[, list(pcr_efficiency=input$pcr_efficiency), keyby="day"]
            attr(eff, "auto") <- FALSE
        } else {
            eff <- GWPCR.PARAMETERS[, list(pcr_efficiency), keyby="day"]
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
            ls <- GWPCR.PARAMETERS[, list(library_size=(input$library_size)), keyby="day"]
            attr(ls, "auto") <- FALSE
        } else {
            ls <- GWPCR.PARAMETERS[, list(library_size), keyby="day"]
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
            th <- GWPCR.PARAMETERS[, list(phantom_threshold=(input$phantom_threshold)), keyby="day"]
            attr(th, "auto") <- FALSE
        } else {
            th <-  GWPCR.PARAMETERS[, list(phantom_threshold), keyby="day"]
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
    san_simulation_with_pcr <- reactive({
        if (!is.null(san_simulation())) {
            message("Simulating PCR and sequencing")
            # Days for which experimental data is available
            days <- LT47[, list(), keyby="day"]
            # Compute total number of cells in each simulated sample
            total_sizes <- san_simulation()[, list(C.total=sum(as.numeric(C))), keyby="t"]
            # Simulate read counts under the gwpcR PCR model
            # The sequencing depth is computed from the total number of cells,
            # and the library size. The resulting table has the additional column
            # "R" containing the simulated read count.
            (san_simulation()
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
    san_simulation_with_pcr_filtered <- reactive({
        if (!is.null(san_simulation_with_pcr())) {
            sim <- san_simulation_with_pcr()
            sim[phantom_threshold_manual_or_auto()[sim, phantom_threshold, on=c(day="t")] <= R]
        }
    })
    
    # Return a list of saved parameter sets
    parametersets.list <- function() {
        sub(x=list.files(file.path(LOCATION.PARAMETERSETS)),
            pattern="\\.rd", replacement="")
    }
    
    # Convert a parameterset name to a full path
    parameterset.fullpath <- function(name) {
        file.path(LOCATION.PARAMETERSETS, paste0(name, ".rd"))
    }
    
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
        if (!is.null(parameterset.loaded$s0) && !is.samevalue(parameterset.loaded$s0, input$s0))
            updateNumericInput(session, "s0", value=parameterset.loaded$s0)
        if (!is.null(parameterset.loaded$rates) && (!are.rates.identical(parameterset.loaded$rates, as.data.table(hot_to_r(input$rates)))))
            output_rates(parameterset.loaded$rates)
        if (!is.null(parameterset.loaded$pcr_efficiency) && !is.samevalue(parameterset.loaded$pcr_efficiency, input$pcr_efficiency))
            updateNumericInput(session, "pcr_efficiency", value=parameterset.loaded$pcr_efficiency)
        if (!is.null(parameterset.loaded$library_size) && !is.samevalue(parameterset.loaded$library_size, input$library_size))
            updateNumericInput(session, "library_size", value=parameterset.loaded$library_size)
        if (!is.null(parameterset.loaded$phantom_threshold) && !is.samevalue(parameterset.loaded$phantom_threshold, input$phantom_threshold))
            updateNumericInput(session, "phantom_threshold", value=parameterset.loaded$phantom_threshold)
    })
    # Load parameter set "default" on startup
    updateSelectInput(session, "loadfrom", selected="default",
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
        headers <- colnames(output_rates())
        headers <- sub("Tmax", "T<sub>max</sub>", headers, fixed=TRUE)
        headers <- sub("^r_([A-Z0-9])$", "g<sub>\\1</sub>", headers)
        rhandsontable(output_rates(), useTypes=TRUE, search=FALSE, colHeaders=headers)
        # %>% hot_validate_numeric(cols=colnames(output_rates()), min=0, max=Inf)
    })

    # Deterministic model dynamics
    output$deterministic <- renderPlot({
        # Do nothing if the rates table is empty
        if (length(input_rates_list()) == 0)
            return()
        
        # Evaluate deterministic SAN model
        r <- evaluate_deterministic_san(s0=input$s0, rates=input_rates_list(), samples_per_day=10)
        
        # Plot results
        ggplot(data=as.data.table(r)) +
            geom_line(aes(x=t, y=S+A+N, col='total size in model')) +
            geom_line(aes(x=t, y=S, col='S-cells in model')) +
            geom_line(aes(x=t, y=A, col='A-cells in model')) +
            geom_line(aes(x=t, y=N, col='N-cells in model')) +
            stat_summary(data=ORGANOIDSIZES, aes(x=day, y=count, col=I('total size in experiment')), fun=mean, geom="point") +
            stat_summary(data=ORGANOIDSIZES, aes(x=day, y=count), fun.data=mean_sdl, geom="errorbar", width=0.5, col="purple") +
            my_scale_log10(scale_y_log10, limits=c(1e2, 1e7)) +
            scale_color_manual(breaks=c('total size in experiment', 'total size in model',
                                        'S-cells in model', 'A-cells in model', 'N-cells in model'),
                               values=c('purple', 'black', 'blue', 'orange', 'green'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of cells")
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
    output$stochastic_lsd <- renderPlot({
        if (is.null(san_simulation()))
            return()

        # Simulate stochastic SAN model
        rank_size.model <- rank_size(san_simulation()[t==input$day_lsd, list(sid=-1, lsize=C)][lsize > 0])
        
        # Fetch experimental data
        rank_size.experiment <- rank_size(LT47[day==input$day_lsd])
        
        # X coordinates to use when plotting curves
        curves.xmax <- max(rank_size.experiment$rank, rank_size.model$rank)
        curves.ymin <- min(rank_size.experiment$size, rank_size.model$size)
        curves.ymax <- max(rank_size.experiment$size, rank_size.model$size)
        xmax <- 10^(ceiling(10*log10(curves.xmax))/10)
        ymin <- 10^(  floor(log10(curves.ymin)))
        ymax <- 10^(ceiling(log10(curves.ymax)))

        # Setup plot
        p <- ggplot() +
            my_scale_log10(scale_x_log10, limits=c(1, xmax)) +
            xlab("rank") +
            my_scale_log10(scale_y_log10, limits=c(ymin, ymax)) +
            ylab("lineage size [norm. reads]") +
            annotation_logticks() +
            scale_color_manual(breaks=c('experiment', 'model', 'model+PCR+seq'),
                               values=c('purple', 'black', 'darkred'),
                               name=NULL) +
            theme(legend.position="bottom") +
            guides(col=guide_legend(ncol=3,byrow=TRUE))
        
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
        if (nrow(rank_size.experiment) > 0)
            plot_rank_size(rank_size.experiment, col="experiment")
        if (nrow(rank_size.model) > 0)
            plot_rank_size(rank_size.model, col="model")
        if (!is.null(san_simulation_with_pcr_filtered())) {
            rank_size.model_pcr <- rank_size(san_simulation_with_pcr_filtered()[t==input$day_lsd, list(sid=-1, lsize=R)])
            plot_rank_size(rank_size.model_pcr, col="model+PCR+seq")
        }
        
        p
    })
    
    output$stochastic_nlineages <- renderPlot({
        if (is.null(san_simulation()))
            return()
        
        lsd <- san_simulation()[, list(nlineages=sum(C > 0)), by="t"]
        lsd_scells <- san_simulation()[, list(nlineages=sum(S > 0)), by="t"]
        
        # Plot results
        p <- ggplot() +
            geom_line(data=lsd, aes(x=t, y=nlineages, col='all lineages in model')) +
            geom_line(data=lsd_scells, aes(x=t, y=nlineages, col='lineages w/ S-cells in model')) +
            stat_summary(data=NLINEAGES, aes(x=day, y=nlineages, col=I('all lineages obs. in experiment')), fun=mean, geom="point") +
            stat_summary(data=NLINEAGES, aes(x=day, y=nlineages), fun.data=mean_sdl, geom="errorbar", width=0.5, col="purple") +
            scale_color_manual(breaks=c('all lineages obs. in experiment',
                                        'all lineages in model', 'lineages w/ S-cells in model',
                                        'lineages obs. in model+PCR+seq'),
                               values=c('purple', 'black', 'blue', 'darkred'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of lineages") +
            theme(legend.position="bottom") +
            guides(col=guide_legend(ncol=1,byrow=FALSE))
        if (!is.null(san_simulation_with_pcr_filtered())) {
            lsd_pcr <- san_simulation_with_pcr_filtered()[, list(nlineages=.N), by="t"]
            p <- p + geom_line(data=lsd_pcr, aes(x=t, y=nlineages, col='lineages obs. in model+PCR+seq'))
        }
        if (input$stochastic_nlineages_logy)
            p <- p + my_scale_log10(scale_y_log10)

        
        p
    })
}
