
# Dependencies for run lncPipe Locally

Prerequisites install command (required when docker image is not favored, you should execute them via root)

* [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)


		aria2c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -q -o /opt/hisat2-2.1.0-Linux_x86_64.zip && \
		unzip -qq /opt/hisat2-2.1.0-Linux_x86_64.zip -d /opt/ && \
		rm /opt/hisat2-2.1.0-Linux_x86_64.zip && \
		cd /opt/hisat2-2.1.0 && \
		rm -rf doc example *debug MANUAL* NEWS TUTORIAL && \
		ln -s /opt/hisat2-2.1.0/hisat2* /usr/local/bin/ && \
		ln -sf /opt/hisat2-2.1.0/*.py /usr/local/bin/


* [StringTie](http://www.ccb.jhu.edu/software/stringtie/)


		aria2c http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz -q -o /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz && \
		tar xf /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz --use-compress-prog=pigz -C /opt/ && \
		rm /opt/stringtie-1.3.3b.Linux_x86_64/README && \
		ln -s /opt/stringtie-1.3.3b.Linux_x86_64/stringtie /usr/local/bin/stringtie && \
		rm /opt/stringtie-1.3.3b.Linux_x86_64.tar.gz


* [gffcompare](http://www.ccb.jhu.edu/software/stringtie/gff.shtml#gffcompare)


		aria2c https://github.com/gpertea/gffcompare/archive/master.zip -q -o /opt/gffcompare-master.zip && \
		aria2c https://github.com/gpertea/gclib/archive/master.zip -q -o /opt/gclib-master.zip && \
		unzip -qq /opt/gffcompare-master.zip -d /opt/ && \
		unzip -qq /opt/gclib-master.zip -d /opt/ && \
		rm /opt/gffcompare-master.zip /opt/gclib-master.zip && \
		cd /opt/gffcompare-master && \
		make release


* [Bedops](http://bedops.readthedocs.io/en/latest/):


		aria2c https://github.com/bedops/bedops/releases/download/v2.4.29/bedops_linux_x86_64-v2.4.29.tar.bz2 -q -o /opt/bedops_linux_x86_64-v2.4.29.tar.bz2 && \
		tar xf /opt/bedops_linux_x86_64-v2.4.29.tar.bz2 --use-compress-prog=pbzip2 -C /opt/ && \
		ln -s /opt/bin/* /usr/local/bin/ && \
		rm /opt/bedops_linux_x86_64-v2.4.29.tar.bz2


* [PLEK](www.ibiomedical.net):


		aria2c https://nchc.dl.sourceforge.net/project/plek/PLEK.1.2.tar.gz -q -o /opt/PLEK.1.2.tar.gz && \
		tar xf /opt/PLEK.1.2.tar.gz --use-compress-prog=pigz -C /opt/ && \
		cd /opt/PLEK.1.2/ && \
		python PLEK_setup.py || : && \
		rm *.pdf *.txt *.h *.c *.fa *.cpp *.o *.R *.doc PLEK_setup.py && \
		chmod 755 * && \
		perl -CD -pi -e'tr/\x{feff}//d && s/[\r\n]+/\n/' *.py && \
		ln -s /opt/PLEK.1.2/* /usr/local/bin/ && \
		rm /opt/PLEK.1.2.tar.gz


* [CNCI](https://github.com/www-bioinfo-org/CNCI):


		aria2c https://codeload.github.com/www-bioinfo-org/CNCI/zip/master -q -o /opt/CNCI-master.zip && \
		unzip -qq /opt/CNCI-master.zip -d /opt/ && \
		rm /opt/CNCI-master.zip && \
		unzip -qq /opt/CNCI-master/libsvm-3.0.zip -d /opt/CNCI-master/ && \
		rm /opt/CNCI-master/libsvm-3.0.zip && \
		cd /opt/CNCI-master/libsvm-3.0 && \
		make > /dev/null 2>&1 && \
		shopt -s extglob && \
		rm -rfv !\("svm-predict"\|"svm-scale"\) && \
		cd .. && \
		rm draw_class_pie.R LICENSE README.md && \
		chmod -R 755 * && \
		ln -s /opt/CNCI-master/*.py /usr/local/bin/


* [CPAT](http://rna-cpat.sourceforge.net):[Citation](https://academic.oup.com/nar/article/41/6/e74/2902455/CPAT-Coding-Potential-Assessment-Tool-using-an)


		aria2c https://jaist.dl.sourceforge.net/project/rna-cpat/v1.2.3/CPAT-1.2.3.tar.gz -q -o /opt/CPAT-1.2.3.tar.gz && \
		tar xf /opt/CPAT-1.2.3.tar.gz --use-compress-prog=pigz -C /opt/ && \
		cd /opt/CPAT-1.2.3/ && \
		mv dat/* /LncPipeDB/ && \
		python setup.py install > /dev/null 2>&1 && \
		rm -rf /opt/CPAT*


* [fastp](https://github.com/OpenGene/fastp)

        RUN aria2c http://opengene.org/fastp/fastp -q -o /usr/local/bin/fastp && \
            chmod a+x /usr/local/bin/fastp


* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)


		aria2c https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -q -o /opt/fastqc_v0.11.5.zip && \
		unzip -qq /opt/fastqc_v0.11.5.zip -d /opt/ && \
		rm /opt/fastqc_v0.11.5.zip && \
		cd /opt/FastQC && \
		shopt -s extglob && \
		rm -rfv !\("fastqc"\|*.jar\) && \
		chmod 755 * && \
		ln -s /opt/FastQC/fastqc /usr/local/bin/


* or [AfterQC](https://github.com/OpenGene/AfterQC)


		aria2c https://github.com/OpenGene/AfterQC/archive/v0.9.7.tar.gz -q -o /opt/AfterQC-0.9.7.tar.gz && \
		tar xf /opt/AfterQC-0.9.7.tar.gz --use-compress-prog=pigz -C /opt/ && \
		cd /opt/AfterQC-0.9.7 && \
		make && \
		perl -i -lape's/python/pypy/ if $. == 1' after.py && \
		rm -rf Dockerfile Makefile README.md testdata report_sample && \
		rm editdistance/*.cpp editdistance/*.h && \
		ln -s /opt/AfterQC-0.9.7/*.py /usr/local/bin/ && \
		rm /opt/AfterQC-0.9.7.tar.gz


When using afterQC, we recommend that users install `pypy` in their operation system, which can accelerate about 3X speed for raw reads processing, as [suggested]((https://github.com/OpenGene/AfterQC#pypy-suggestion)) by the author of AfterQC.

* [LncPipeReporter](https://github.com/bioinformatist/LncPipe-Reporter)

		Install [pandoc](https://pandoc.org/installing.html) first. Then run commands:

		Rscript -e "install.packages('devtools'); devtools::install_github('bioinformatist/LncPipeReporter')"

		For detailed usage of LncPipeReporter in case you are going to run it separately, plz refers to [README](https://github.com/bioinformatist/LncPipeReporter#lncpipereporter) of LncPipeReporter.

* [kallisto](https://github.com/pachterlab/kallisto)


		aria2c https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz -q -o  /opt/kallisto_linux-v0.43.1.tar.gz && \
		tar xf /opt/kallisto_linux-v0.43.1.tar.gz --use-compress-prog=pigz -C /opt/ && \
		cd /opt && \
		rm ._* kallisto_linux-v0.43.1.tar.gz && \
		cd kallisto_linux-v0.43.1 && \
		rm -rf ._* 	README.md test && \
		ln -s /opt/kallisto_linux-v0.43.1/kallisto /usr/local/bin/


* [sambamba](http://lomereiter.github.io/sambamba/)


        aria2c https://github.com/biod/sambamba/releases/download/v0.6.7/sambamba_v0.6.7_linux.tar.bz2 -q -o /opt/sambamba_v0.6.7_linux.tar.bz2 && \
        tar xf /opt/sambamba_v0.6.7_linux.tar.bz2 --use-compress-prog=pbzip2 -C /opt/ && \
        ln -s /opt/sambamba /usr/local/bin/ && \
        rm /opt/sambamba_v0.6.7_linux.tar.bz2


**Alternatively, when you are going to using STAR-Cufflinks in your analysis, the corresponding install cmd are as follows:**

* [STAR](https://github.com/alexdobin/STAR)


		aria2c https://raw.githubusercontent.com/alexdobin/STAR/master/bin/Linux_x86_64/STAR -q -o /opt/STAR && \
		chmod 755 /opt/STAR && \
		ln -s /opt/STAR /usr/local/bin


* [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks)


		aria2c https://github.com/bioinformatist/cufflinks/releases/download/v2.2.1/cufflinks-2.2.1.Linux_x86_64.tar.gz -q -o /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
		tar xf /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz --use-compress-prog=pigz -C /opt/ && \
		rm /opt/cufflinks-2.2.1.Linux_x86_64/README && \
		ln -s /opt/cufflinks-2.2.1.Linux_x86_64/* /usr/local/bin/ && \
		rm /opt/cufflinks-2.2.1.Linux_x86_64.tar.gz


> The `gffcompare` utility share the same function as `cuffcompare`, therefore, in STAR-cufflinks analysis pipe, `gffcompare` is not required.
