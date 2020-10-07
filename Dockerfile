FROM cibiv-shiny-base-r40:2020-10-07

MAINTAINER Florian G. Pflug <florian.pflug@univie.ac.at>

# Install R packages not on CRAN via devtools
RUN R -e "library(devtools); install_github('Cibiv/gwpcR', ref='v0.9.10')"

# Persistent data
VOLUME /state

# Setup directories, make sure stored parametersets persist
RUN mkdir -p /src && \
    mkdir -p /sanmodelexplorer && \
    ln -sf /state/parametersets /sanmodelexplorer/

# Run App
EXPOSE 3838
COPY entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]
CMD ["R", "-e", "shiny::runApp('/sanmodelexplorer')"]

# R seettings
COPY Rprofile.site /usr/lib/R/etc/

# Copy data
COPY data/lt47.processed.rd     /sanmodelexplorer/data/
COPY data/organoidsizes.rd      /sanmodelexplorer/data/
COPY data/celltypes.rd          /sanmodelexplorer/data/

# Copy SANsimulatoR source and install package
COPY SANsimulatoR /src/SANsimulatoR
RUN R CMD INSTALL --preclean --clean /src/SANsimulatoR

# Copy App
COPY *.R                        /sanmodelexplorer/

# Copy upstream parametersets
COPY parametersets/*            /sanmodelexplorer/parametersets.dist/

# Set container revision
ARG CONTAINER_REVISION
RUN echo "$CONTAINER_REVISION" > /sanmodelexplorer/CONTAINER_REVISION
