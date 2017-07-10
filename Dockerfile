# Set the base image to Debian
FROM debian:9.0

# File Author / Maintainer
MAINTAINER **nalcala** <**alcalan@fellows.iarc.fr**>

RUN mkdir -p /var/cache/apt/archives/partial && \
	touch /var/cache/apt/archives/lock && \
	chmod 640 /var/cache/apt/archives/lock && \
	apt-get update -y &&\
	apt-get install -y gnupg2
	
RUN	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys F76221572C52609D && \
	apt-get clean && \
	apt-get update -y && \


  # Install dependences
  DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  make \
  g++ \
  perl \
  default-jre \
  zlib1g-dev \
  libncurses5-dev \
  libncurses5 \
  git \
  wget \
  ca-certificates \
  python-dev \
  python-pip \
  bzip2 \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libfreetype6-dev \
  libpng-dev \
  unzip \
  r-base \
  r-cran-ggplot2 \
  r-cran-gplots \
  r-cran-reshape && \
  cp /usr/include/freetype2/*.h /usr/include/. && \

  Rscript -e 'install.packages("gsalib",repos="http://cran.us.r-project.org")' && \
  
  # Install samtools specific version manually
  wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar -jxf samtools-1.3.1.tar.bz2 && \
  cd samtools-1.3.1 && \
  make && \
  make install && \
  cd .. && \
  rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2 && \

  # Install FastQC
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
  unzip fastqc_v0.11.5.zip && \
  chmod 755 FastQC/fastqc && \
  cp -r FastQC /usr/local/bin/. && \
  ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/ && \
  rm -rf fastqc_v0.11.5.zip FastQC && \

  # Install cutadapt
  pip install cutadapt && \

  # Install trim_galore
  wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz && \
  tar xvzf 0.4.3.tar.gz && \
  mv TrimGalore-0.4.3/trim_galore /usr/bin && \
  rm -rf TrimGalore-0.4.3 0.4.3.tar.gz && \

  # Install hisat2

  # Install htseq
  pip install numpy && \
  pip install setuptools && \
  pip install HTSeq && \

  # Install multiqc
  pip install --upgrade --force-reinstall git+https://github.com/nalcala/MultiQC.git && \

  # Install STAR specific version manually
  wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
  tar -xzf 2.5.3a.tar.gz && \
  cp STAR-2.5.3a/bin/Linux_x86_64_static/STAR /usr/local/bin/. && \
  rm -rf 2.5.3a.tar.gz STAR-2.5.3a && \

  # Install hisat2
  wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
  unzip hisat2-2.1.0-Linux_x86_64.zip && \
  cp -r hisat2-2.1.0/. /usr/local/bin/. && \
  rm -rf hisat2-2.1.0-Linux_x86_64.zip hisat2-2.1.0 && \

  # Install RSeQC
  pip install RSeQC && \

  # Install samblaster specific version manually
  wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz && \
  tar -xzf samblaster-v.0.1.24.tar.gz && \
  cd samblaster-v.0.1.24 && \
  make && \
  cp samblaster /usr/local/bin/. && \
  cd .. && \
  rm -rf samblaster-v.0.1.24.tar.gz samblaster-v.0.1.24 && \

  # Install sambamba specific version manually
  wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && \
  tar -jxf sambamba_v0.6.6_linux.tar.bz2 && \
  cp sambamba_v0.6.6 /usr/local/bin/sambamba && \
  rm -rf sambamba_v0.6.6_linux.tar.bz2 && \

  # Remove unnecessary dependences
  DEBIAN_FRONTEND=noninteractive apt-get remove -y \
  make \
  g++ \
  wget \
  bzip2 \
  git \
  zlib1g-dev \
  libncurses5-dev && \

  # Clean
  DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \
  apt-get clean
