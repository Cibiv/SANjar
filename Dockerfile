FROM openanalytics/r-base:4.0.2

MAINTAINER Florian G. Pflug <florian.pflug@univie.ac.at>

VOLUME /state

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0 \
    libxml2-dev

# Install inline package and anything that might be usefull or is a common dependency
RUN R -e "install.packages(c( \
        'inline', 'Rcpp', 'RcppParallel', 'RcppEigen', 'RcppArmadillo', 'BH' \
      ), repos='https://cloud.r-project.org/')"

# Install basic R packages
RUN R -e "install.packages(c( \
        'data.table', 'MASS', 'VGAM', 'Deriv', 'lamW' \
      ), repos='https://cloud.r-project.org/')"

# Install R packages data.table, inline, shiny, ggplot2 and extensions
# Note that Hmisc is listed here because it depends on ggplot2
RUN R -e "install.packages(c( \
        'shinycssloaders', 'DT', 'rhandsontable', \
        'ggplot2', 'scales', 'cowplot', 'gridExtra', 'Hmisc' \
      ), repos='https://cloud.r-project.org/')"

# Install R packages not on CRAN via devtools
RUN R -e "install.packages(c( \
        'devtools' \
      ), repos='https://cloud.r-project.org/')" && \
    R -e "library(devtools); install_github('schaffman5/shinyTree', ref='e41e200437d828faf203c22b98f9f6d6cac58206')" && \
    R -e "library(devtools); install_github('Cibiv/gwpcR', ref='v0.9.10')"

# Install additional packages
RUN R -e "install.packages(c( \
        'bsplus' \
      ), repos='https://cloud.r-project.org/')"

# R seettings
COPY Rprofile.site /usr/lib/R/etc/

# Setup directories
RUN mkdir -p /sanmodelexplorer && \
    ln -sf /state/parametersets /sanmodelexplorer/

# Copy data
COPY data/organoidsizes.tab     /sanmodelexplorer/data/
COPY data/lt47.rd               /sanmodelexplorer/data/

# Run App
EXPOSE 3838
ENTRYPOINT ["/entrypoint.sh"]
CMD ["R", "-e", "shiny::runApp('/sanmodelexplorer')"]

# Copy App
COPY entrypoint.sh              /
COPY *.R                        /sanmodelexplorer/

# Copy upstream parametersets
COPY parametersets/*            /sanmodelexplorer/parametersets.dist/
