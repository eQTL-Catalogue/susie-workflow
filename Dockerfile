FROM bioconductor/release_core2
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for finemapping analysis with susieR"

RUN R -e "BiocManager::install(c('GDSArray','devtools', 'SNPRelate', 'dplyr','purrr', 'readr', 'optparse','SummarizedExperiment','expm','matrixStats','lumi','limma','tidyr','cqn','plotly'))"
RUN R -e "devtools::install_github('stephenslab/susieR@0.8.0')"
RUN R -e "devtools::install_github('kauralasoo/eQTLUtils')"
