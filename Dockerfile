################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="rnaseq-nf"
LABEL software.version="2.4"
LABEL about.summary="Container image containing all requirements for rnaseq-nf"
LABEL about.home="http://github.com/IARCbioinfo/RNAseq-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/RNAseq-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/RNAseq-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda config --set channel_priority strict
RUN conda env create -n rnaseq-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rnaseq-nf/bin:$PATH
