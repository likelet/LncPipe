# [LncPipe](https://github.com/likelet/LncPipe) 
[![AUR](https://img.shields.io/aur/license/yaourt.svg)](https://github.com/likelet/LncPipe/blob/master/LICENSE)
 [![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
## Overall
Recently, long noncoding RNA molecules (lncRNA) captured widespread attentions for its critical 
roles in diverse biological process and important implications in variety of human diseases and 
cancers. Identification and profiling of lncRNAs is a fundamental step to advance our knowledge 
on their function and regulatory mechanisms. However, RNA sequencing based lncRNA discovery is 
limited due to complicated operations and implementation. Therefore, we presented a one-stop 
pipeline called [LncPipe](https://github.com/likelet/LncPipe) focused on characterizing lncRNAs from raw transcriptome sequencing 
data. The pipeline was developed based on a popular workflow framework [Nextflow](https://github.com/nextflow-io/nextflow), composed of 
four core procedures including reads alignment, assembly, identification and quantification. 
It contains various unique features such as well-designed lncRNAs annotation strategy, optimized 
calculating efficiency, diversified classification and interactive analysis report. Additionally, 
[LncPipe](https://github.com/likelet/LncPipe)  allows users cancel pipeline, reset parameters from command or modifying main script 
directly and resume analysis from continues checkpoint. 

## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Schematic diagram](#schematic-diagram)
- [Installation](#installation)
- [Run Docker](#run-docker)
- [Interactive reports](#interactive-reports)
- [Parameters](#parameters)
- [FAQ](#faq)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)
- [License](#license)

## Schematic diagram


## Installation 
[Nextflow](https://github.com/nextflow-io/nextflow)  
LncPipe is implemented with Nextflow pipeline manage system. To run our pipelines. [Nextflow](https://github.com/nextflow-io/nextflow) should be preinstalled at  POSIX compatible system (Linux, Solaris, OS X, etc), It requires BASH and Java 7 or higher to be installed. We do not recommend running the pipes in the Windows since most of bioinformatic tools do not supported.
Here, we show the step by step installation of [Nextflow](https://github.com/nextflow-io/nextflow) in linux system as an example, which adapted from [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html).

1. Download the executable package by copying and pasting the following command in your terminal window:


        wget -qO- get.nextflow.io | bash  
        
    
> It will create the [Nextflow](https://github.com/nextflow-io/nextflow) main executable file in the current directory.

2. Optionally, move the nextflow file in a directory accessible by your `$PATH` variable (this is only required to avoid to remember and type the Nextflow full path each time you need to run it).
Of course you can download the lastest binary version of NextFlow by yourself from the https://github.com/nextflow-io/nextflow/releases and add the path into your system environment.
All those pipelines were written in [Nextflow](https://github.com/nextflow-io/nextflow) commands. For more details, please see [here](https://www.nextflow.io).

3. A type command for run nextflow  

nextflow run LncRNAanalysisPipe.nf <parameters>
### Prepare input files 
#### References, index and annotation files（required）

1. [Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml) index (e.g. human index can be downloaded from ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz) or [STAR](https://github.com/alexdobin/STAR) index (hg38 genome index etc.) according aligner your are going to use. 
    > Building index of hisat relatively require large amount of memory, thus we sugguested that users downloaded it directly from the hisat website.  
    2. Genome reference (genome fasta file with suffix `.fa` , `UCSC` etc.).   
    3. GENCODE gene annotation file in GTF format  
    4. LNCipedia gene annotation file in GTF format.(set null if not available for your species)  
    5. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed  
    
#### Supported species

    >Currently, LncPipe was designed for human only, but it also support other species which required users providing "known_protein_coding.gtf" and  "known_lncRNA.gtf "
    the detail usage for non-human species could be found here.  

* human 

    1. hisat index: ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz
    2. Genome reference: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz
    3. GENCODE gene annotation: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
    4. LNCipedia gene annotation: https://lncipedia.org/downloads/lncipedia_5_0_hc_hg38.gtf
    5. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed

* mouse 

    1. hisat index: ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz
    2. Genome reference: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz
    3. GENCODE gene annotation: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
    4. LNCipedia gene annotation: null
    5. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed

## Run Docker 

1. Prepare input files.
2. Modify the `docker.config` in `mandatory` section.
3. Install docker
4. Command


        nextflow -c docker.config run LncRNAanalysisPipe.nf

   docker image can be downloaded from https://hub.docker.com/r/bioinformatist/lncpipe/tags/ with the latest tag. 
> Alternatively, nextflow can automatically pull image from docker.io. `Dockerfile` recorded  that what we have done with the image. For docker pull image from local china, [we suggest users using mirror site instead](https://github.com/likelet/Blogs_tips/blob/master/README.md#setting-docker-download-mirror-site).

## Dependencies 

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
		  
		
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)

		
		aria2c https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -q -o /opt/fastqc_v0.11.5.zip && \
		unzip -qq /opt/fastqc_v0.11.5.zip -d /opt/ && \
		rm /opt/fastqc_v0.11.5.zip && \
		cd /opt/FastQC && \
		shopt -s extglob && \
		rm -rfv !\("fastqc"\|*.jar\) && \
		chmod 755 * && \
		ln -s /opt/FastQC/fastqc /usr/local/bin/
		  
			
or [AfterQC](https://github.com/OpenGene/AfterQC)
			
		
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

## Interactive reports
LncPipe output was well-summarized in an interactive manner, which was carried out by our novel-developing R package
[LncPipeReporter](https://github.com/bioinformatist/LncPipeReporter).

## Configuration
As a nextflow-based analysis pipeline, LncPipe allow users edit configure file `nextflow.config` to set the index files and default file path parameters instead of typing in command.
We strongly recommend users using config file rather than command input to adjust their parameters.
To configure, plz go to `params` line, and set the following information of your operation system and environment
groovy

        params {
        /*
            User setting options (mandatory)
             */
        // input file and genome reference()
            fastq_ext = '*_{1,2}.clean.fq.gz'
            fasta_ref = '/data/database/hg38/genome.fa'
            design = 'design.file'
            hisat2_index = '/data/database/hg38/hisatIndex/grch38_snp_tran/genome_snp_tran'
            gencode_annotation_gtf = "/data/database/hg38/Annotation/gencode.v24.annotation.gtf"
            lncipedia_gtf = "/data/database/hg38/Annotation/lncipedia_4_0_hg38.gtf"
            cpatpath = '/home/zhaoqi/software/CPAT/CPAT-1.2.2'
        
        /*
            User setting options (optional)
             */
            star_idex = ''//set if star used
            bowtie2_index = ''//set if tophat used
            aligner = "hisat" // or "star","tophat"
            sam_processor="sambamba"//or "samtools"
            qctools = "fastqc" // or "afterqc","fastp"
            singleEnd = false
            unstrand = false
            skip_combine = false
            lncRep_Output = 'reporter.html'
            lncRep_theme = 'npg'
            lncRep_cdf_percent = 10
            lncRep_max_lnc_len = 10000
            lncRep_min_expressed_sample = 50
            mem=60//Gb
            cpu=30
        }



## Parameters 
> Those parameters would cover the setting from `nextflow.config` file
* Mandatory(plz configure those in *nextflow.config* file)

| Name | Example/Default value | Description |
|-----------|--------------:|-------------|
|--input_folder | `.` | input folder |
|--species | `human` | Your species, mouse, fly and zebra fish are also supported |
|--fastq_ext | `*_{1,2}.fastq.gz` | input raw paired reads |
|--out_folder |  `.` | output folder |
|--design     | `FALSE` | a txt file that stored experimental design information, plz see details from `--design` section below |

* References

| Name | Required | Description |
|-----------|--------------|-------------|
|--star_index/--bowtie2_index/--hisat2_index  | -| Path to STAR／bowtie2/hisat2(mutually exclusive) index(required if not set in config file)  |
|--fasta  | `-` | Path to Fasta reference(required if not set in config file)|
|--gencode_annotation_gtf  | `-` | An annotation file from GENCODE database for annotating lncRNAs(required if not set in config file). e.g. [gencode.v26.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz) |
|--lncipedia_gtf  | `-` | An annotation file from LNCipedia database for annotating lncRNAs(required if not set in config file) e.g. [lncipedia_4_0_hc_hg38.gtf](http://www.lncipedia.org/downloads/lncipedia_4_0_hc_hg38.gtf) |

* software path (should not setting when using docker )


| Name | Required | Description |  
|-----------|--------------|-------------|  
|--cpatpath  | `-` | Home folder of cpat installed location |

> since cpat may call model data from its home path, users should specified where the model file is located in. Especially users install cpat by themselves without our install code.  

* Optional

| Name | Default value | Description |
|-----------|--------------|-------------|
|--singleEnd  | `FALSE` | specify that the reads are single ended  |
|--merged_gtf | `FALSE` | Skip mapping and assembly step by directly providing assembled merged gtf files|
|--unstrand     | `FALSE` | specify that library is unstrand specific  |
|--aligner |  `star` | Aligner for reads mapping (optional), STAR is default and supported only at present,*star*/*tophat*/*hisat2*|
|--qctools |  `afterqc` | Tools for assess raw reads quality or filtered by AfterQC, *fastqc* or *afterqc*|

* LncPipeReporter options

| Name | Default value | Description |
|-----------|--------------|-------------|
|--lncRep_Output  | `reporter.html` | Specify report file name. |
|--lncRep_theme  | `npg` | Plot theme setting in interactive plot. Values from [ggsci](https://github.com/road2stat/ggsci) |
|--lncRep_min_expressed_sample  | `50` |Minimum expressed gene allowed in each sample, 50 default. Samples not passed were filtered from analysis|

`--fastq_ext`
> Raw fastq files were required for denovo analysis.This parameters should be set according to your paired or singled reads file names.
Suppose your paired end sequence files are compressed with `.gz` suffixed.
For example:

        Sample1_1.fq.gz
        Sample1_2.fq.gz
        Sample2_1.fq.gz
        Sample2_2.fq.gz

Then you can input pattern `*_{1,2}.fq.gz` to make the all paired end file recognized by [LncPipe](https://github.com/likelet/LncPipe) .

For singled reads file, file pattern should be feed with `--singleEnd` specified.


`--star_idex／--bowtie2_index/--hisat2_index`

> This parameter is *required* when not configured in nextflow.config file. It specify the star/tophat/hisat2(mutually exclusive) index folder built before running [LncPipe](https://github.com/likelet/LncPipe) .
If you don't know what it is，You can use `--fasta` to specify the reference sequence data. The index file would be built by [LncPipe](https://github.com/likelet/LncPipe)  automatically.


`--design`
> Experimental design file matrix for differential expression analysis. Default: `null`
Format:

    WT:Sample1,Sample2,Sample3
    KO:Sample1,Sample2,Sample3

While `KO/WT` represents the two experimental condition, and sample1, sample2, sample3 are replicates which should be comma-delimited in the same line .

For sample names, it should be the sample as the prefix of fastq files which was trimmed by `--fastq_ext`.

For example:

 if fastq file names are `Sample1_1.fq.gz, Sample1_2.fq.gz` that comes from one sample and your `--fastq_ext` is set as `*_{1,2}.fq.gz`, the sample name
should be Sample1.

## Output
While the whole pipeline is finished properly, there is `Result` folder under current path(default) or output_folder set by user. The basic structure of Result is follows:

        Result/
        ├── QC
        │   ├── N1141_1.clean_fastqc.html
        │   ├── N1141_2.clean_fastqc.html
        │   ├── N1177_1.clean_fastqc.html
        │   └── N1177_2.clean_fastqc.html
        ├── Identified_lncRNA
        │   ├── all_lncRNA_for_classifier.gtf
        │   ├── final_all.fa
        │   ├── final_all.gtf
        │   ├── lncRNA.fa
        │   ├── protein_coding.fa
        │   └── protein_coding.final.gtf
        ├── LncReporter
        │   ├── Differential_Expression_analysis.csv
        │   └── Report.html
        ├── Quantification
        │   ├── kallisto.count.txt
        │   └── kallisto.tpm.txt
        └── Star_alignment
            ├── STAR_N1141
            │   ├── N1141Aligned.sortedByCoord.out.bam
            │   ├── N1141Log.final.out
            │   ├── N1141Log.out
            │   ├── N1141Log.progress.out
            │   └── N1141SJ.out.tab
            └── STAR_N1177
                ├── N1177Aligned.sortedByCoord.out.bam
                ├── N1177Log.final.out
                ├── N1177Log.out
                ├── N1177Log.progress.out
                └── N1177SJ.out.tab


* `QC` stored the Quality control output generated by FastQC or AfterQC software.<br>
* `Identified_lncRNA` contains all assembled lncRNA and their sequences. *all_lncRNA_for_classifier.gtf* includes both novel and known lncRNA features in [GTF format](http://www.ensembl.org/info/website/upload/gff.html);
*lncRNA.fa* is all lncRNA sequences in fasta format. *protein_coding.final.gtf* and *protein_coding.fa* are protein coding information extracted from gencode annotation. *final_all.gtf* and *final_all.fa* are combined files for further analysis.<br>
* `Alignment` are hisat/tophat/STAR aligner standard output<br>
* `Quantification` are estimated abundance using kallisto. *kallisto.count.txt* stored reads count matrix and *kallisto.tpm.txt* are tpm(Transcripts Per Kilobase Million) matrix.
* `LncReporter` stored the interactive report file and differential expression matrix generated by LncPipeReporter which wrapped EdgeR.

## Tips	 

* :blush:Plz keep the consistency of your genome sequence, index library and annotation files: genome version, chromosome format, gtf coordinated e.g. The third-party software may stop for any of the above reasons. 
* :confused:Setting your analysis parameters always in config file, differ project should  corresponding to differ configurations for reproductive analysis. To rerun a project, you can just specify -c `your.config` in your command, which can also help you to record analysis parameters.
* :open_mouth:Run analysis on docker container, no much to say.
* :grimacing:Always use the latest version to be away from the known bugs. 

## Acknowledgement  

 Thanks to the author of [AfterQC](https://github.com/OpenGene/AfterQC), Shifu Chen, for his help on providing a gzip output support to meet the require of LncPipe.  Thanks to the internal test by Hongwan Zhang and Yan Wang from SYSUCC Cancer bioinformatics platform.

## FAQ  

* *1. PLEK throws an error "/data/software/PLEK.1.2/PLEK.py:line12: $'\r': can not find command", how to fix?*
>A: using the follow command as suggested in the installation section. 
    
        perl -CD -pi -e'tr/\x{feff}//d && s/[\r\n]+/\n/' *.py 
    
* *2. IOError: [Errno 2] No such file or directory: '/opt/CPAT-1.2.3/dat/Human_Hexamer.tsv'?*
>A: The cpat command required  the `Human_Hexamer.tsv` to predict lncRNA coding potential, plz check your `cpatpath` parameters. 
* *3. When using htseq to quantify transicript, it throws "Error occured when reading beginning of SAM/BAM file. 'csamtools.AlignedRead' object has no attribute 'reference_start' "*
>A: It's a version conflict caused by htseq and hisat generated bamfile, a possible solution for this is to install the old version of htseq


    
## Contact

For implementation:  
* [Qi Zhao](https://github.com/likelet) zhaoqi@sysucc.org.cn, Sun Yat-sen University Cancer Center;  
* [Yu Sun](http://icannotendure.space)  sun_yu@mail.nankai.edu.cn, Nan kai University;  
For project design and new feature request:  
* [Qi Zhao](https://github.com/likelet) zhaoqi@sysucc.org.cn, Sun Yat-sen University Cancer Center;  
* [Zhixiang Zuo]() zuozhx@sysucc.org.cn, Sun Yat-sen University Cancer Center;  

> We strongly recommend users open new issues if they have questions or find bugs.


## License
[GPL v3 license](LICENSE)
## Citation 

