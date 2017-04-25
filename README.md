# Workflow
## 1. Description
A collection of self-developed bioinformatics analysis pipeline coding by a DSL language, nextflow.
To run the pipeline files with nf suffixed, plz just type the following command.

```
nextflow <your nf file> -c nextflow.config -with-trace
```

If the pipeline fails at any point and you fix the issue, you can revoke the processes with job avoidance using this command:
```
nextflow <your nf file> -c nextflow.config -with-trace -resume
```

All those pipelines were written in Nextflow commands. For more details, please see the following link:
https://www.nextflow.io/



## 2. Installation of dependencies.
### Install NextFlow
To run our pipelines. NextFlow should be preinstalled at  POSIX compatible system (Linux, Solaris, OS X, etc), It requires BASH and Java 7 or higher to be installed. We do not recommend running the pipes in the Windows since most of bioinformatic tools do not supported.
Here, we show the step by step installation of NextFlow in linux system as an example, which adapted from [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html).

1. Download the executable package by copying and pasting the following command in your terminal window: ```wget -qO- get.nextflow.io | bash```. It will create the nextflow main executable file in the current directory.
2. Optionally, move the nextflow file in a directory accessible by your `$PATH` variable (this is only required to avoid to remember and type the Nextflow full path each time you need to run it).
3. Download the lastest binary verion of NextFlow from the https://github.com/nextflow-io/nextflow/releases and add the path into your system environment.
### Install third-party software and database required by each pipe.
#### Pipe 1
##### introduction
##### References, index and annotation files（required）
1. STAR index (hg38 genome index etc.)
2. Genome reference (genome fasta file with suffix .fa , UCSC etc).
3. GENCODE gene annotation file in GTF format:
      [gencode.v26.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz)
4. LNCipedia gene annotation file in GTF format:
      [lncipedia_4_0_hc_hg38.gtf](http://www.lncipedia.org/downloads/lncipedia_4_0_hc_hg38.gtf)
##### software and tools (required)
1. [STAR](https://github.com/alexdobin/STAR), Reference https://www.ncbi.nlm.nih.gov/pubmed/23104886
2. [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks), Reference https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146043/
3. [Bedops](http://bedops.readthedocs.io/en/latest/), Reference https://www.ncbi.nlm.nih.gov/pubmed/22576172/
4. [PLEK]
      ```
      wget https://sourceforge.net/projects/plek/files/PLEK.1.2.tar.gz/download
      tar -zvxf PLEK.1.2.tar.gz 
      cd PLEK.1.2
      python PLEK_setup.py 
      ```
5. [CNCI]
      ```
      git clone git@github.com:www-bioinfo-org/CNCI.git
      cd CNCI
      unzip libsvm-3.0.zip
      cd libsvm-3.0
      make
      cd ..
      ```
6. [CPAT]
      ```
      wget https://sourceforge.net/projects/rna-cpat/files/?source=navbar
      tar zxf CPAT-VERSION.tar.gz
      cd CPAT-VERSION
      python setup.py install
      python setup.py install --root=/home/user/CPAT
      export PYTHONPATH=/home/user/CPAT/usr/local/lib/python2.7/site-packages:$PYTHONPATH.
      export PATH=/home/user/CPAT/usr/local/bin:$PATH #setup PATH
      ```


