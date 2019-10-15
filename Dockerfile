################## BASE IMAGE #####################
FROM nfcore/base


################## METADATA #######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="rnaseq-nf"
LABEL software.version="2.1"
LABEL about.summary="Container image containing all requirements for rnaseq-nf"
LABEL about.home="http://github.com/IARCbioinfo/RNAseq-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/RNAseq-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/RNAseq-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>


#RUN mkdir -p /var/cache/apt/archives/partial && \
#	touch /var/cache/apt/archives/lock && \
#	chmod 640 /var/cache/apt/archives/lock && \
#	apt-get update -y &&\
#	apt-get install -y gnupg2
#	RUN	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys F76221572C52609D && \
#	apt-get clean && \
#	apt-get update -y && \

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env create -n rnaseq-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rnaseq-nf/bin:$PATH
#RUN echo ". /opt/conda/etc/profile.d/conda.sh"  >> ~/.bashrc 
