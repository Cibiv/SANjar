FROM cibiv-shiny-base-r40:2020-10-02

MAINTAINER Florian G. Pflug <florian.pflug@univie.ac.at>

VOLUME /state

# Install R packages not on CRAN via devtools
RUN R -e "library(devtools); install_github('Cibiv/gwpcR', ref='v0.9.10')"

# R seettings
COPY Rprofile.site /usr/lib/R/etc/

ARG CONTAINER_REVISION

# Setup directories
RUN mkdir -p /src && \
    mkdir -p /sanmodelexplorer && \
    echo "$CONTAINER_REVISION" > /sanmodelexplorer/CONTAINER_REVISION && \
    ln -sf /state/parametersets /sanmodelexplorer/

# Copy SANsimulatoR source and install package
COPY SANsimulatoR /src/SANsimulatoR
RUN R CMD INSTALL --preclean --clean /src/SANsimulatoR

# Copy data
COPY data/lt47.processed.rd     /sanmodelexplorer/data/
COPY data/organoidsizes.rd      /sanmodelexplorer/data/
COPY data/celltypes.rd          /sanmodelexplorer/data/

# Run App
EXPOSE 3838
ENTRYPOINT ["/entrypoint.sh"]
CMD ["R", "-e", "shiny::runApp('/sanmodelexplorer')"]

# Copy App
COPY entrypoint.sh              /
COPY *.R                        /sanmodelexplorer/

# Copy upstream parametersets
COPY parametersets/*            /sanmodelexplorer/parametersets.dist/
