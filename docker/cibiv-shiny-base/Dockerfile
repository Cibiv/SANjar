FROM openanalytics/r-base:4.0.2

MAINTAINER Florian G. Pflug <florian.pflug@univie.ac.at>

VOLUME /state

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    cargo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0 \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libavfilter-dev \
    libavcodec-dev \
    libavformat-dev \
    libavutil-dev \
    libavresample-dev \
    libswresample-dev \
    libswscale-dev \
    libpostproc-dev

# Install devtools plus basic libraries usefull for packages containing C/C++ code
RUN R -e "install.packages(c( \
        'devtools', 'inline', \
        'Rcpp', 'RcppParallel', 'RcppEigen', 'RcppArmadillo', 'BH', \
        'bsplus', 'dqrng' \
      ), repos='https://cloud.r-project.org/')"

# Install basic R packages
RUN R -e "install.packages(c( \
        'data.table', 'MASS', 'VGAM', 'Deriv', 'lamW' \
      ), repos='https://cloud.r-project.org/')"

# Install shiny, ggplot, tidyverse and related packages
RUN R -e "install.packages(c( \
        'shiny', 'shinycssloaders', 'DT', 'rhandsontable', \
        'tidyverse', \
        'ggplot2', 'scales', 'cowplot', 'gridExtra', 'Hmisc', \
        'gganimate', 'gifski', 'av', 'transformr' \
      ), repos='https://cloud.r-project.org/')"
