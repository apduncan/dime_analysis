FROM bioconductor/bioconductor_docker:latest
# Required libraries
RUN apt update && apt install -y \
    libglpk-dev 

# Using bioconductor dockers, as standard install too time consuming
# bioconductor based on rocker/rstudio, so missing some of the tidyverse
# packages I like. Took the below from the rocker instal_tidyverse.sh
ADD install_tidyverse.sh /tmp/
RUN bash /tmp/install_tidyverse.sh

# Bioconductor installs
RUN R -e 'BiocManager::install(ask = FALSE)' && \
    R -e 'BiocManager::install(c(\
    "microbiome" \
    ), ask = FALSE)'
RUN R -e 'install.packages("littler", repos="https://cran.r-project.org")'

# Challenging installs of NetCoMi and SpiecEasi
RUN R -e 'BiocManager::install("SpiecEasi")'
RUN install2.r --error \
    devtools \
    remotes
RUN R -e 'pak::pak("GraceYoon/SPRING")'
RUN R -e 'devtools::install_github("stefpeschel/NetCoMi", \
repos = c("https://cloud.r-project.org/", \
BiocManager::repositories()))'

# Install anything possible from CRAN
RUN install2.r --error \
    targets \
    vegan \
    quarto \
    tarchetypes \
    devtools \
    patchwork \
    readxl \
    ggrepel \
    ggpubr \
    pheatmap \
    ggplotify \
    usedist \
    reticulate \
    lsa \
    rtk \
    logger
    # openxslx \
    # colorpsace

# Install anything needed which has been removed from CRAN or needs a pinned
# version in the Dockerfile
RUN R -e 'library(devtools);\
    install_version("httpgd", version="2.0.2", ask = FALSE)'

# Install some packages from custom repos, often need for VSCode extentions
RUN install2.r --error --repos https://p3m.dev/cran/__linux__/noble/latest \
    languageserver && \
    install2.r --error --repos https://manuelhentschel.r-universe.dev \
    vscDebugger

RUN install2.r --error tidyverse