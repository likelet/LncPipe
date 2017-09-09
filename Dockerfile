FROM ubuntu

LABEL authors="zhaoqi@sysucc.org.cn,sun_yu@mail.nankai.edu.cn" \
	description="Docker image containing all requirements for LncPipe"

# pigz will decompress all archieves to this directory
WORKDIR /LncPipeDB/

# Once the container started, decompress all databases
ENTRYPOINT pigz -d /LncPipeDB/*.gz && \
	mkdir /LncPipeDB/hg38 /LncPipeDB/hg19 && \
	mv /LncPipeDB/*hg38*.gtf /LncPipeDB/hg38 && \
	mv /LncPipeDB/gencode.v27.annotation.gtf /LncPipeDB/hg38 && \
	mv /LncPipeDB/*hg19*.gtf /LncPipeDB/hg19 && \
	mv /LncPipeDB/gencode.v27lift37.annotation.gtf /LncPipeDB/hg19 && \
	/bin/bash

# Update OS
# Relieve the dependence of readline perl library by prohibiting interactive frontend first
RUN export DEBIAN_FRONTEND=noninteractive && \
	apt-get -qq update &&\
	apt-get -qq install -y --no-install-recommends \
	# For Nextflow (run groovy)
	default-jre \
	# For decompress GitHub archieve
	unzip \
	# pbzip2 \
	pigz \
	aria2 \
	# Below two is needed for CPAT and PLEK compiling
	gcc \
	g++ \
	# Needed when compiling R package hexbin (required by plotly)
	gfortran \
	# For exec makefile of libsvm-3.0 used by CNCI
	make \
	# Provide head file like Python.h for CPAT compiling
	python-dev \
	# Must install cython HERE, DO NOT use pip, which will cause missing .h files
	cython \
	# For CPAT compiling dependency
	zlib1g-dev \
	# For samtools compiling dependency
	# libncurses5-dev \
	# Required by R package openssl
	libssl-dev \
	# Required by R package curl
	libcurl4-openssl-dev \
	# Required by cpanm, or will get "Can't locate PerlIO.pm in @INC" error
	# FindBin module in
	perl

# Download databases
# ADD CANNOT download FTP links and CANNOT resume from break point
RUN aria2c ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz -q && \
	aria2c ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz -q

# COPY local databases (we don't need auto-decompression, so DO NOT use ADD)
COPY *.gz /LncPipeDB/

# Set working directory back to /
WORKDIR /

# Install latest pip WITHOUT setuptools and wheel
# DO NOT use apt-get python-pip in ubuntu to prevent from complicated related tools and libraries
# Keep the image size down
RUN aria2c https://bootstrap.pypa.io/get-pip.py -q -o /opt/get-pip.py && \
	python /opt/get-pip.py --no-setuptools --no-wheel && \
	rm /opt/get-pip.py

# Install required python packages
RUN pip -qqq install numpy

# Install nextflow
RUN aria2c https://github.com/nextflow-io/nextflow/releases/download/v0.25.6/nextflow -q -o /opt/nextflow && \
	chmod 755 /opt/nextflow && \
	ln -s /opt/nextflow /usr/local/bin

# Install STAR
RUN aria2c https://raw.githubusercontent.com/alexdobin/STAR/master/bin/Linux_x86_64/STAR -q -o /opt/STAR && \
	chmod 755 /opt/STAR && \
	ln -s /opt/STAR /usr/local/bin

# Install cufflinks
RUN aria2c https://github.com/bioinformatist/cufflinks/releases/download/v2.2.1/cufflinks-2.2.1.Linux_x86_64.tar.gz -q -o /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
	tar xf /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz --use-compress-prog=pigz -C /opt/ && \
	rm /opt/cufflinks-2.2.1.Linux_x86_64/README && \
	ln -s /opt/cufflinks-2.2.1.Linux_x86_64/* /usr/local/bin/ && \
	rm /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz

# Install CPAT
RUN aria2c https://jaist.dl.sourceforge.net/project/rna-cpat/v1.2.3/CPAT-1.2.3.tar.gz -q -o /opt/CPAT-1.2.3.tar.gz && \
	tar xf /opt/CPAT-1.2.3.tar.gz --use-compress-prog=pigz -C /opt/ && \
	# DO NOT use absolute path here, changing directory is necessary, python interpreter will check current directory for dependencies
	cd /opt/CPAT-1.2.3/ && \
	mv dat/* /LncPipeDB/ && \
	python setup.py install > /dev/null 2>&1 && \
	rm -rf /opt/CPAT*

# Install PLEK
RUN aria2c https://nchc.dl.sourceforge.net/project/plek/PLEK.1.2.tar.gz -q -o /opt/PLEK.1.2.tar.gz && \
	tar xf /opt/PLEK.1.2.tar.gz --use-compress-prog=pigz -C /opt/ && \
	cd /opt/PLEK.1.2/ && \
	python PLEK_setup.py || : && \
	# Remove documents, demo files, source files, object files and R scripts
	rm *.pdf *.txt *.h *.c *.model *.range *.fa *.cpp *.o *.R *.doc PLEK_setup.py && \
	chmod 755 * && \
	# dos2unix in perl one-liner: remove BOM head and deal with \r problem
	perl -CD -pi -e'tr/\x{feff}//d && s/[\r\n]+/\n/' *.py && \
	ln -s /opt/PLEK.1.2/* /usr/local/bin/ && \
	rm /opt/PLEK.1.2.tar.gz

# Install CNCI
# Use bash instead of sh for shopt only works with bash
SHELL ["/bin/bash", "-c"]
RUN aria2c https://codeload.github.com/www-bioinfo-org/CNCI/zip/master -q -o /opt/CNCI-master.zip && \
	unzip -qq /opt/CNCI-master.zip -d /opt/ && \
	rm /opt/CNCI-master.zip && \
	unzip -qq /opt/CNCI-master/libsvm-3.0.zip -d /opt/CNCI-master/ && \
	rm /opt/CNCI-master/libsvm-3.0.zip && \
	cd /opt/CNCI-master/libsvm-3.0 && \
	make > /dev/null 2>&1 && \
	# enable the extglob shell option
	shopt -s extglob && \
	# Parentheses and the pipe symbol should be escaped
	rm -rfv !\("svm-predict"\|"svm-scale"\) && \
	cd .. && \
	rm draw_class_pie.R LICENSE README.md && \
	chmod -R 755 * && \
	ln -s /opt/CNCI-master/*.py /usr/local/bin/

# Set back to default shell
# SHELL ["/bin/sh", "-c"]

# Install StringTie
RUN aria2c http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz -q -o /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz && \
	tar xf /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz --use-compress-prog=pigz -C /opt/ && \
	rm /opt/stringtie-1.3.3b.Linux_x86_64/README && \
	ln -s /opt/stringtie-1.3.3b.Linux_x86_64/stringtie /usr/local/bin/stringtie && \
	rm /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz

# Install Hisat2
RUN aria2c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -q -o /opt/hisat2-2.1.0-Linux_x86_64.zip && \
	unzip -qq /opt/hisat2-2.1.0-Linux_x86_64.zip -d /opt/ && \
	rm /opt/hisat2-2.1.0-Linux_x86_64.zip && \
	cd /opt/hisat2-2.1.0 && \
	rm -rf doc example *debug MANUAL* NEWS TUTORIAL && \
	ln -s /opt/hisat2-2.1.0/hisat2* /usr/local/bin/ && \
	ln -sf /opt/hisat2-2.1.0/*.py /usr/local/bin/

# Install Kallisto
RUN aria2c https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz -q -o  /opt/kallisto_linux-v0.43.1.tar.gz && \
	tar xf /opt/kallisto_linux-v0.43.1.tar.gz --use-compress-prog=pigz -C /opt/ && \
	# There's some trashy pointers in Kallisto tarball
	cd /opt && \
	rm ._* kallisto_linux-v0.43.1.tar.gz && \
	cd kallisto_linux-v0.43.1 && \
	rm -rf ._* 	README.md test && \
	ln -s /opt/kallisto_linux-v0.43.1/kallisto /usr/local/bin/

# Install Microsoft-R-Open with MKL
RUN aria2c https://mran.microsoft.com/install/mro/3.4.0/microsoft-r-open-3.4.0.tar.gz -q -o /opt/microsoft-r-open-3.4.0.tar.gz && \
	tar xf /opt/microsoft-r-open-3.4.0.tar.gz --use-compress-prog=pigz -C /opt/ && \
	cd /opt/microsoft-r-open && \
	./install.sh -as && \
	rm -rf /opt/microsoft-r*

# Cleaning up the apt cache helps keep the image size down (must be placed here, since MRO installation need the cache)
RUN rm -rf /var/lib/apt/lists/*

# Install R packages (only for LncPipe-Reporter at current stage)
RUN echo 'install.packages("devtools")' > /opt/packages.R && \
	echo 'install.packages(c("curl", "httr"))' >> /opt/packages.R && \
	echo 'install.packages("data.table")' >> /opt/packages.R && \
	echo 'install.packages("cowplot")' >> /opt/packages.R && \
	echo 'install.packages("DT")' >> /opt/packages.R && \
	echo 'devtools::install_github("ramnathv/htmlwidgets")' >> /opt/packages.R && \
	# Plotly always has too many bugs fixed, so keep using develop branch version :)
	echo 'devtools::install_github("ropensci/plotly")' >> /opt/packages.R && \
	echo 'devtools::install_github("vqv/ggbiplot")' >> /opt/packages.R && \
	echo 'source("https://bioconductor.org/biocLite.R")' >> /opt/packages.R && \
	echo 'biocLite()' >> /opt/packages.R && \
	echo 'biocLite("edgeR")' >> /opt/packages.R && \
	# Heatmaply and plotly may need different version of ggplot2 in dev-branch, keep this statement last to avoid compatibility problems
	printf 'install.packages("heatmaply")' >> /opt/packages.R && \
	Rscript /opt/packages.R && \
	rm /opt/packages.R

# Install cpanminus
#RUN aria2c https://cpanmin.us/ -q -o /opt/cpanm && \
#	chmod +x /opt/cpanm && \
#	ln -s /opt/cpanm /usr/local/bin/

# Install Perl module FindBin, which is required by FastQC
# RUN cpanm FindBin

# Install FastQC
RUN aria2c https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -q -o /opt/fastqc_v0.11.5.zip && \
	unzip -qq /opt/fastqc_v0.11.5.zip -d /opt/ && \
	rm /opt/fastqc_v0.11.5.zip && \
	cd /opt/FastQC && \
	shopt -s extglob && \
	rm -rfv !\("fastqc"\|*.jar\) && \
	chmod 755 * && \
	ln -s /opt/FastQC/fastqc /usr/local/bin/

# Lines below maybe used later
# Install BWA
#RUN bash -c 'aria2c https://codeload.github.com/lh3/bwa/zip/master -q -o /opt/bwa-master.zip && \
#	unzip -qq /opt/bwa-master.zip -d /opt/ && \
#	rm /opt/bwa-master.zip && \
#	cd /opt/bwa-master && \
#	make > /dev/null 2>&1 && \
#	shopt -s extglob && \
#	rm -rfv !\("bwa"\|"qualfa2fq.pl"\|"xa2multi.pl"\|"COPYING"\) && \
#	ln -s /opt/bwa-master/bwa /usr/local/bin/ && \
#	ln -s /opt/bwa-master/*.pl /usr/local/bin/'

# Install SAMtools (incomplete)
#RUN aria2c https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -q -o /opt/samtools-1.5.tar.bz2 && \
#	tar xf /opt/samtools-1.5.tar.bz2 --use-compress-prog=pbzip2 -C /opt/ && \
#	cd /opt/samtools-1.5 && \
#	make && \
#	make install && \
#	rm /opt/samtools-1.5.tar.bz2

# Install BCFtools
# https://github.com/samtools/bcftools/releases/download/1.5/bcftools-1.5.tar.bz2
