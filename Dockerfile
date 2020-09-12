FROM ubuntu:20.04

MAINTAINER Florian G. Pflug <florian.pflug@univie.ac.at>

# Volumes
VOLUME /state

# Versions
ENV VER_SHINY_SERVER=1.5.14.948

# Install Debian/Ubuntu packages
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -qq && \
    apt-get install -qq -y \
      -o Dpkg::Options::="--force-confdef" \
      -o Dpkg::Options::="--force-confold" \
      --no-install-recommends  \
      curl ca-certificates lsb-release gnupg2 libnss-wrapper && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    (echo "\n# CRAN R packages\ndeb http://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" \
      >> /etc/apt/sources.list) && \
    apt-get update -qq && \
    apt-get install -qq -y \
      -o Dpkg::Options::="--force-confdef" \
      -o Dpkg::Options::="--force-confold" \
      --no-install-recommends  \
      r-base-core && \
    curl -# -O https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-${VER_SHINY_SERVER}-amd64.deb && \
    dpkg --install shiny-server-${VER_SHINY_SERVER}-amd64.deb && \
    rm -r /var/lib/shiny-server /var/log/shiny-server /srv/shiny-server \
      shiny-server-${VER_SHINY_SERVER}-amd64.deb && \
    apt-get clean && \
    apt-get purge -qq -y curl && \
    apt autoremove -qq -y

# Install R packages
# Note that we leave the -dev packages installed, because the inline package
# might need them, and it will definitely need the compiler and binutils.
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -qq && \
    apt-get install -qq -y \
      -o Dpkg::Options::="--force-confdef" \
      -o Dpkg::Options::="--force-confold" \
      --no-install-recommends  \
      r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev && \
    R -e "install.packages(c( \
        'devtools', 'inline', \
        'data.table', 'ggplot2', 'Deriv', 'lamW', \
        'shiny', 'shinycssloaders', 'DT', 'rhandsontable' \
      ), repos='https://cran.rstudio.com/')" && \
    R -e "library(devtools); install_github('schaffman5/shinyTree', ref='e41e200437d828faf203c22b98f9f6d6cac58206')" && \
    R -e "library(devtools); install_github('Cibiv/gwpcR', ref='v0.9.10')" && \
    apt-get clean && \
    apt autoremove -qq -y

# Remover all users but root and nobody
RUN find ./ -mount -a \( \( ! -uid 0 \) -o \( ! -gid 0 \) \) -print0 | \
       xargs -0 chown root.root && \
    ((echo 'root:x:0:0:root:/root:/bin/bash'; \
      echo 'nobody:x:65534:65534:nobody:/nonexistent:/usr/sbin/nologin') > /etc/passwd) && \
    ((echo 'root:*:17211:0:99999:7:::'; \
      echo 'nobody:*:17211:0:99999:7:::') > /etc/shadow) && \
    ((echo 'root:x:0:'; \
      echo 'nogroup:x:65534:') > /etc/group)

# Copy SANModelExplorer
RUN mkdir -p /srv/sanmodelexplorer && ln -sf /state/parametersets /srv/sanmodelexplorer/
COPY entrypoint.sh /
COPY *.R /srv/sanmodelexplorer/
COPY data/*.tab /srv/sanmodelexplorer/data
COPY data/*.rd /srv/sanmodelexplorer/data/
COPY shiny-server.conf /etc/shiny-server

# Run App
EXPOSE 3838
ENTRYPOINT ["/entrypoint.sh"]
CMD ["/opt/shiny-server/bin/shiny-server"]
