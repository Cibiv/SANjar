library(data.table)
library(shiny)
library(rhandsontable)
library(ggplot2)
library(scales)
library(extraDistr)

# Parameters:
#  LOCATION.PARAMETERSETS
#  DEFAULT.PARAMETERSET
#  DATASETS

# Pattern defining valid "saveas" filenames for parameter sets
FILENAME.PATTERN <- "^[A-Za-z0-9. -]*$"

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

# Equation to solve to determine the lambda parameter
# for a zero-truncated Poisson distribution with mean m
tpois.lambda.eq <- deriv(expression((m - l/(1-exp(-l)))^2), namevec=c("l"),
                         function.arg=c("l", "m"))

# Compute lambda parameter for a zero-truncated Poisson distribution with  mean m
tpois.lambda <- function(mean) {
    if (mean <= 1)
        return(NA_real_)
    r <- nlm(tpois.lambda.eq, p=mean, m=mean)
    return(r$estimate)
}

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
    # Datasets
    if (length(DATASETS) > 1) {
        enable("dataset")
        updateSelectInput(session, "dataset", choices=names(DATASETS))
    } else {
        disable("dataset")
        updateSelectInput(session, "dataset", choices=character())
    }
    
    # Dataset (all groups)
    dataset <- reactive({
        if (length(DATASETS) > 1) {
            if ((length(input$dataset) == 1) && !is.null(input$dataset) && (input$dataset != "")) {
                message("Switching to dataset ", input$dataset)
                DATASETS[[input$dataset]]
            } else {
                message("Switching to first dataset by default, selected dataset ", input$dataset, " does not exist")
                DATASETS[[1]]
            }
        } else {
            message("Switching to first dataset by default, no dataset selected")
            DATASETS[[1]]
        }
    })
    
    # Datagroups
    # Datasets can contain multiple groups of samples, distinguished by
    # the group key column(s) dataset()$groups. dataset_groups() lists the
    # available groups and their label (last column and key).
    dataset_groups <- reactiveVal()
    observeEvent(dataset_groups(), {
        dg <- dataset_groups()
        if (is.null(dg) || (nrow(dg) < 1)) {
            updateSelectInput(session, "datagroup", choices=character())
            disable("datagroup")
        } else {
            enable("datagroup")
            labels <- dg[, key(dg), with=FALSE]$V1
            updateSelectInput(session, "datagroup", choices=labels)
        }
    })
    observeEvent(dataset(), {
        message("Dataset changed, updating data group list")
        ds <- dataset()
        dg <- if (length(ds$groups) > 0) {
            # Create table listing groups and their label (first column and key)
            dg <- ds$lineagesizes[, list(paste0(names(.BY), "=", lapply(.BY, as.character),
                                                collapse=" "))
                                  , by=eval(ds$groups) ]
            k <- length(ds$groups)
            dg[, colnames(dg)[c(k+1, 1:k)], with=FALSE]
        } else {
            # Dataset has no groups
            dg <- data.table(label=character())
        }
        # Update group list
        setkeyv(dg, colnames(dg)[1])
        dataset_groups(dg)
    })

    # Lineagesize unit
    dataset_lineagesize_expr <- reactiveVal()
    observeEvent(dataset(), {
        ds <- dataset()
        if (ds$unit == "cells") {
            disable("stochastic_lsd_normlibsize")
            dataset_lineagesize_expr(expression(cells))
        } else if (ds$unit == "reads") {
            enable("stochastic_lsd_normlibsize")
            dataset_lineagesize_expr(expression(reads))
        } else
            stop("lineagesizes must contain column 'cells' or 'reads'")
    })

    # Dataset (after filtering down to a single group)
    dataset_group <- reactive({
        # It's a waste to respond to dataset changes here, since changing
        # the dataset will also update dataset_groups() 
        ds <- isolate(dataset())
        # Select group
        group <- if (nrow(dataset_groups()) > 1) {
            # Translate label to group key
            # If the translation fails, select the first group in the dataset
            groupkey <- NULL
            if ((length(input$datagroup) == 1) && !is.null(input$datagroup) && (input$datagroup != "")) {
                groupkey <- dataset_groups()[list(input$datagroup), !1]
                if (nrow(groupkey) != 1)
                    groupkey <- NULL
            }
            # Validate the group key, if invalid select first group
            if (is.null(groupkey) || !all(names(groupkey) %in% ds$groups))
                groupkey <- dataset_groups()[1, !1]
            # Filter dataset down to only the selected group
            message("Selecting data group ",
                    paste0(names(groupkey), "=", lapply(groupkey, as.character), collapse=", "))
            do.call(subset, c(list(ds), groupkey))
        } else {
            message("Selecting full dataset")
            ds
        }
        # Normalize library sizes if lineage sizes are relative and normalization was requested
        if ((ds$unit == "reads") && input$stochastic_lsd_normlibsize)
            group <- SANjar::normalize_library_sizes(group, method="scale")
        group
    })
    dataset_tfinal <- reactive(({
        dataset_group()$lineagesizes[, max(day)]
    }))
    dataset_s0scale <- reactive({
        f <- dataset_group()$organoidsizes[day==0, median(fraction)]
        if (is.finite(f) && (f > 0)) f else 1.0
    })
    dataset_group_nlineages <- reactive({
        dataset_group()$lineagesizes[, list(
            nlineages=sum(eval(dataset_lineagesize_expr()) > 0)
        ), keyby=.(day, sid)]
    })
    dataset_group_ranksize <- reactive({
        rank_size(dataset_group()$lineagesizes[, list(
            sid, day, size=eval(dataset_lineagesize_expr())
        )])
    })
    dataset_group_logsize <- reactive({
        dataset_group()$lineagesizes[, list(
            sid, log_size=log10(eval(dataset_lineagesize_expr()))
        ), keyby=.(day, sid)]
    })  

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
            max(input_rates()$Tmax, dataset_tfinal(), na.rm=TRUE)
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
        choices <- sort(unique(c(dataset_group()$lineagesizes$day, Tfinal())))
        if ((length(sel) < 1) || is.na(sel) || (sel == ""))
            sel <- max(dataset_group()$lineagesizes$day, na.rm=TRUE)
        if (!(sel %in% choices))
            sel <- choices[length(choices)]
        updateSelectInput(session, "day_lsd", choices=choices, selected=sel)

        # Update day_scellext
        updateSliderInput(session, "day_scellext", min=0, max=Tfinal()+1, value=Tfinal()+1)
    })
    
    s0_message <- reactiveVal()
    s0_actual <- reactive({
        if (!input$s0_scale) {
            s0_message("")
            return(input$s0)
        }
        s0 <- ceiling(input$s0 * dataset_s0scale())
        s0_message(paste0("Using s0=", s0, " (", signif(100*dataset_s0scale(), digits=2), "% of ", input$s0, ")"))
        return(s0)
    })
    
    san_deterministic_results <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (s0_actual() > 0)) {
            r <- san_deterministic(s0=s0_actual(), rates=input_rates(), samples_per_day=10)[, list(
                sid=-1, day=t,
                S, A, N,
                C=S+A+N
            )]
            setkey(r, day, sid)
            r
        } else NULL
    })

    # Simulating the stochastic model is slow, so we cache the result
    # The result is a table with columns day, lid, S, A, N and C=S+A+N
    san_stochastic_results <- reactive({
        # Do nothing if the rates table is empty or s0 is zero
        if ((length(input_rates_list()) > 0) && (s0_actual() > 0)) {
            message("Simulating the SAN model")
            rt <- copy(input_rates())
            # If the last Tmax in the rates table precedes Tfinal(), adjust
            # the rates table accordingly
            if (rt[nrow(rt), Tmax] < Tfinal())
                rt[nrow(rt), Tmax := Tfinal()]
            # Simulate model
            r <- san_stochastic(L=s0_actual(), rates=rt, samples_per_day=1,
                                p_cutoff=as.numeric(input$p_cutoff))[, list(
                sid=-1,
                day=t,
                dt, lid, S, A, N,
                C=S+A+N
            )]
            setkey(r, day, sid, lid)
            r
        } else NULL
    })

    p_cutoff_message <- reactive({
        # Compute expected values for S, A, N using the deterministic model
        mod.exp <- san_deterministic(s0=s0_actual(), rates=input_rates())
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
            eff <- dataset_group()$sequencing[, list(pcr_efficiency=input$pcr_efficiency), keyby=day]
            attr(eff, "auto") <- FALSE
        } else {
            eff <- dataset_group()$sequencing[, list(pcr_efficiency=median(pcr_efficiency)), keyby=day]
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
            ls <- dataset_group()$sequencing[, list(library_size=as.numeric(input$library_size)), keyby=day]
            attr(ls, "auto") <- FALSE
        } else {
            ls <- dataset_group()$sequencing[, list(library_size=round(median(library_size))), keyby=day]
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
            th <- dataset_group()$sequencing[, list(phantom_threshold=(input$phantom_threshold)), keyby=day]
            attr(th, "auto") <- FALSE
        } else {
            th <-  dataset_group()$sequencing[, list(phantom_threshold=round(median(phantom_threshold))), keyby=day]
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
    
    # Lineage aliases parameter
    lineage_aliases_manual_or_auto <- reactive({
        if (!is.na(input$lineage_aliases) && (input$lineage_aliases > 0)) {
            a <- dataset_group()$sequencing[, list(lineage_aliases=(input$lineage_aliases)), keyby=day]
            auto <- FALSE
        } else {
            a <- dataset_group()$sequencing[, list(lineage_aliases=median(lineage_aliases, na.rm=TRUE)), keyby=day]
            a[is.na(lineage_aliases), lineage_aliases := 1]
            auto <- TRUE
        }
        # Compute lambda parameter of a zero-truncated Poisson distribution that
        # yields the requested average number of lineage_aliases per cell
        a <- a[, list(day, lineage_aliases_lambda=tpois.lambda(lineage_aliases)), by=.(lineage_aliases)][, list(
            lineage_aliases, lineage_aliases_lambda
        ), keyby=day]
        attr(a, "auto") <- auto
        gwpcr_parameters_interpolate(a)
    })
    lineage_aliases_auto_message <- reactive({
        lineage_aliases <- lineage_aliases_manual_or_auto()
        if (!attr(lineage_aliases, "auto")) return("")
        paste0("est. from data, =", signif(lineage_aliases[day==input$day_lsd, lineage_aliases], 3),
               " uniformly")
    })
    
    san_stochastic_results_aliases <- reactive({
        if (!is.null(san_stochastic_results())) {
            message("Simulating lineage aliasing")
            # Simulate lineage aliasing, i.e. lineages which have multiple labels
            (san_stochastic_results()
             [lineage_aliases_manual_or_auto(), on="day", nomatch=NULL][, {
                 if (is.na(lineage_aliases_lambda[1]))
                     list(day=day, sid=sid, lid=lid, dt=dt, S=S, A=A, N=N, C=C)
                 else {
                     r <- rtpois(length(day), lambda=lineage_aliases_lambda[1], a=0)
                     list(day=rep(day, r), sid=rep(sid, r), lid=rep(lid, r), dt=rep(dt, r),
                          S=rep(S, r), A=rep(A, r), N=rep(N, r), C=rep(C, r))
                 }
             }])
        } else NULL
    })
    
    # Scale cell counts 
    san_stochastic_results_reads <- reactive({
        if (!is.null(san_stochastic_results_aliases())) {
            message("Scaling cell counts to make them comparable with reads")
            # Translate cell counts to reads according to library size 
            # In contrast to san_stochastic_results_with_pcr(), we don't
            # actually introduce any stochasticity here, but simply adjust
            # the scale. C.obs is therefore identical to C here.
            r <- (san_stochastic_results_aliases()
                  [pcr_efficiency_manual_or_auto(), on="day", nomatch=NULL]
                  [library_size_manual_or_auto(), on="day", nomatch=NULL]
                  [, list(
                      dt, lid, S, A, N, C,
                      R=as.numeric(C)*(library_size[1]/sum(C)),
                      reads_per_cell=library_size[1]/sum(C)
                  ), by=.(day, sid)])
            r[, C.obs := C]
            setkey(r, day, sid, lid)
            r
        } else NULL
    })
    
    san_stochastic_ranksize <- reactive({
        if (!is.null(san_stochastic_results_reads())) {
            rank_size(san_stochastic_results_reads()[, list(
                sid, day, size=eval(san_stochastic_results_lineagesize_expr()))])
        } else NULL
    })
    
    # Simulate PCR+Sequencing
    san_stochastic_results_with_pcr <- reactive({
        if (!is.null(san_stochastic_results_aliases())) {
            message("Simulating PCR and sequencing")
            # Simulate read counts under the gwpcR PCR model
            # The sequencing depth is computed from the total number of cells,
            # and the library size. The resulting table has the additional column
            # "R" containing the simulated read count.
            r <- (san_stochastic_results_aliases()
                [pcr_efficiency_manual_or_auto(), on="day", nomatch=NULL]
                [library_size_manual_or_auto(), on="day", nomatch=NULL]
                [, list(
                    dt, lid, S, A, N, C,
                    R=seqsim(C, reads.target=library_size[1], efficiency=pcr_efficiency[1]),
                    reads_per_cell=library_size[1]/sum(C)
                ), by=.(day, sid)])
            r[, C.obs := R / reads_per_cell]
            setkey(r, day, sid, lid)
            r
        } else NULL
    })
    san_stochastic_results_lineagesize_expr <- reactive({
        ds <- dataset()
        if (ds$unit == "cells")
            expression(C.obs)
        else if (ds$unit == "reads")
            expression(R)
        else
            stop("lineagesizes must contain column 'cells' or 'reads'")
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
            rank_size(san_stochastic_results_with_pcr_filtered()[, list(
                sid, day, size=eval(san_stochastic_results_lineagesize_expr()))])
        } else NULL
    })

    # The currently loaded parameterset
    # Note that this reflects the values of the parameters at loading time
    parameterset.loaded <- reactiveVal()
    
    # Load parameters and update UI
    observeEvent(input$loadfrom, {
        # Return if no parameter set is selected
        if ((length(input$loadfrom) != 1) || (input$loadfrom == "")) {
            return()
        }
        # Return if the selected file is already loaded
        ps <- isolate(parameterset.loaded())
        if (!is.null(ps) && is.character(ps$name) && (ps$name == input$loadfrom)) {
            message("Parameter set ", input$loadfrom, " is already loaded")
            return()
        }
        # Load file
        message("Loading parameter set ", input$loadfrom, " from ", parameterset.fullpath(input$loadfrom))
        parameterset.loaded(local({
            load(file=parameterset.fullpath(input$loadfrom))
            c(SAN.PARAMS, list(name=input$loadfrom))
        }))
    })
    observeEvent(parameterset.loaded(), {
        ps <- parameterset.loaded()
        if (is.null(ps))
            return()
        # Update parameter values
        message("Parameter set loaded, updating UI")
        if (!is.null(ps$s0))
            updateNumericInput(session, "s0", value=ps$s0)
        if (!is.null(ps$rates))
            updateRateInput(session, value=ps$rates, alwaysUpdateOutputAsWell=TRUE)
    })
    # Load parameter set "default" on startup
    if (is.character(DEFAULT.PARAMETERSET))
        updateSelectInput(session, "loadfrom", selected=DEFAULT.PARAMETERSET,
                          choices=c("select set to load"="", parametersets.list()))
    else if (is.list(DEFAULT.PARAMETERSET)) {
        updateSelectInput(session, "loadfrom", selected=character(0),
                          choices=c("select set to load"="", parametersets.list()))
        parameterset.loaded(DEFAULT.PARAMETERSET)
    } else {
        stop("DEFAULT.PARAMETERSET be a parameterset object or the name of a parameter set to load")
    }
    
    # If a loaded parameter set is changed, don't show it as selected anymore
    observeEvent(input$s0, {
        ps <- isolate(parameterset.loaded())
        if (is.null(ps$name))
            return()
        if (is.null(ps$s0) || is.samevalue(ps$s0, input$s0))
            return()
        message("Parameter set has been modified (s0 changed), no longer matches ", parameterset.loaded()$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
        parameterset.loaded(NULL)
    })
    observeEvent(input$rates, {
        ps <- isolate(parameterset.loaded())
        if (is.null(ps$name))
            return()
        if (is.null(ps$rates) || are.rates.identical(ps$rates, as.data.table(hot_to_r(input$rates))))
            return()
        message("Parameter set has been modified (rates changed), no longer matches ", parameterset.loaded()$name)
        updateSelectInput(session, "loadfrom", selected=character(0))
        parameterset.loaded(NULL)
    })

    # Save the current set of parameters
    parameterset.saveas <- function() {
        SAN.PARAMS <- list(rates=as.data.table(hot_to_r(input$rates)),
                           s0=input$s0)
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
        parameterset.loaded(c(SAN.PARAMS, list(name=input$saveas)))
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

    # Message about s0 scaling 
    output$s0_message <- renderUI({
        helpText(HTML(s0_message()))
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
            stat_summary(data=dataset_group()$organoidsizes, aes(y=cells, col='exp.C'), fun=mean, geom="line", linetype="dashed", size=LWD) +
            stat_summary(data=dataset_group()$organoidsizes, aes(y=cells, col='exp.C'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
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
    output$lineage_aliases_auto_message <- renderUI({
        helpText(HTML(lineage_aliases_auto_message()))
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
        
        # Fetch experimental data
        rank_size.experiment <- dataset_group_ranksize()[day==input$day_lsd]
            
        # Abort if no data is available to plot
        if (nrow(rank_size.model) + nrow(rank_size.model_pcr) + nrow(rank_size.experiment) == 0)
            return(NULL)

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10(scale_x_log10, oob=oob_keep) +
            my_scale_log10(scale_y_log10, oob=oob_keep) +
            annotation_logticks(sides="bl") +
            xlab("rank") +
            ylab(paste0("lineage size [", dataset()$unit, "]")) +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL)

        # Draw rank-size curves for the experimental lineage sizes in all replicates
        plot_rank_size <- function(rank_size, col, ...) {
            rank_size[, {
                p <<- p + geom_line(data=copy(.SD), aes(x=rank, y=size, col=col), ...)
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
        log_size.model <- if (input$stochastic_lsd_incpuremodel) {
            san_stochastic_results_reads()[(day==input$day_lsd) & (C > 0), list(
                sid, log_size=log10(eval(san_stochastic_results_lineagesize_expr()))
            )]
        } else data.table()

        # Simulate PCR+Sequencing
        log_size.model_pcr <-  if (!is.null(san_stochastic_results_with_pcr_filtered())) {
            san_stochastic_results_with_pcr_filtered()[day==input$day_lsd, list(
                sid, log_size=log10(eval(san_stochastic_results_lineagesize_expr()))
            )]
        } else data.table()

        # Fetch experimental data
        log_size.experiment <- dataset_group_logsize()[day==input$day_lsd]

        # Legend overrides
        legend.override <- list(linetype=c(data="dashed", model="solid", `model+seq.`="solid"),
                                size=c(data=LWD, model=LWD, `model+seq.`=LWD2))
        groups <- list()

        # Setup plot
        p <- ggplot() +
            my_scale_log10_pretransformed(scale_x_continuous) +
            annotation_logticks(sides="b") +
            xlab(paste0("lineage size [", dataset()$unit, "]")) +
            ylab("density") +
            scale_color_manual(breaks=c('data', 'model', 'model+seq.'),
                               values=c('maroon', 'violet', 'black'),
                               name=NULL)

        # Draw densities
        if (nrow(log_size.experiment) > 0) {
            groups[length(groups)+1] <- "data"
            log_size.experiment[, {
                p <<- p + stat_density(data=copy(.SD), aes(x=log_size, col='data'),
                                       geom="line", size=LWD, linetype="dashed")
                NULL
            }, by=sid]
        }
        if (input$stochastic_lsd_incpuremodel && (nrow(log_size.model) > 0)) {
            groups[length(groups)+1] <- "model"
            p <- p + stat_density(data=log_size.model, aes(x=log_size, col='model'),
                                  geom="line", size=LWD)
        }
        if (nrow(log_size.model_pcr) > 0) {
            groups[length(groups)+1] <- "model+seq."
            p <- p + stat_density(data=log_size.model_pcr, aes(x=log_size, col='model+seq.'),
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
            stat_summary(data=dataset_group_nlineages(), aes(col='exp.obs'), fun=mean, geom="line", size=LWD, linetype="dashed") +
            stat_summary(data=dataset_group_nlineages(), aes(col='exp.obs'), fun.data=mean_sdl, geom="errorbar", width=0.8, size=LWD) +
            geom_line(data=lsd, aes(col='mod.all'), size=LWD)
        if (input$stochastic_nlineages_logy)
            p <- p +
                my_scale_log10(scale_y_log10) +
                annotation_logticks(sides="l")
        if (nrow(lsd_pcr) > 0) {
            groups[length(groups)+1] <- "mod.obs"
            p <- p + geom_line(data=lsd_pcr, aes(col='mod.obs'), size=LWD2)
        }

        # Add legend
        groups[length(groups)+1] <- "mod.S"
        groups.all <- c('exp.obs', 'mod.obs', 'mod.all', 'mod.S')
        groups <- groups.all[groups.all %in% groups]
        p <- p +
            geom_line(data=lsd_scells, aes(col='mod.S'), size=LWD) +
            scale_color_manual(breaks=groups.all,
                               labels=c('obs. (data)', 'obs. (model+seq.)',
                                        'total (model)', 'extant S-cells (model)'),
                               values=c('maroon', 'black', 'violet', 'cornflowerblue'),
                               name=NULL) +
            xlab("time [days]") +
            ylab("number of lineages")

        # Finish plot
        p + guides(col=guide_legend(override.aes=lapply(legend.override, function(o) { o[unlist(groups)] }),
                                    ncol=2, byrow=FALSE))
    })
}
