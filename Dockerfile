FROM nfcore/base
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for the SuSiE fine-mapping pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/finemapping-1.0dev/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN R -e "devtools::install_github('stephenslab/susieR@0.8.0')"
RUN R -e "devtools::install_github('kauralasoo/eQTLUtils')"