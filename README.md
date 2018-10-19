# nf-core/lncpipe
**nf-core LncPipe test**

[![Build Status](https://travis-ci.org/nf-core/lncpipe.svg?branch=master)](https://travis-ci.org/nf-core/lncpipe)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/lncpipe.svg)](https://hub.docker.com/r/nfcore/lncpipe)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable man$


### Documentation
The nf-core/lncpipe pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)


# [LncPipe](https://github.com/likelet/LncPipe) 
[![AUR](https://img.shields.io/aur/license/yaourt.svg)](https://github.com/likelet/LncPipe/blob/master/LICENSE)
 [![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
## Overall
Recently, long noncoding RNA molecules (lncRNA) captured widespread attentions for their critical 
roles in diverse biological process and important implications in variety of human diseases and 
cancers. Identification and profiling of lncRNAs is a fundamental step to advance our knowledge 
on their function and regulatory mechanisms. However, RNA sequencing based lncRNA discovery is 
currently limited due to complicated operations and implementation of the tools involved. Therefore, we present a one-stop multi-tool integrated pipeline called [LncPipe](https://github.com/likelet/LncPipe) focused on characterizing lncRNAs from raw transcriptome sequencing data. 
The pipeline was developed based on a popular workflow framework [Nextflow](https://github.com/nextflow-io/nextflow), composed of four core procedures including reads alignment, assembly, identification and quantification. It contains various unique features such as well-designed lncRNAs annotation strategy, optimized calculating efficiency, diversified classification and interactive analysis report. [LncPipe](https://github.com/likelet/LncPipe) allows users additional control in interuppting the pipeline, resetting parameters from command line, modifying main script directly and resume analysis from previous checkpoint.

## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Schematic diagram](#schematic-diagram)
- [Installation](#installation-and-quick-start)
- [Run Docker](#run-docker)
- [Run with example data](https://github.com/likelet/LncPipeTestData)
- [Interactive reports](#interactive-reports)
- [Parameters](#parameters)
- [FAQ](#faq)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)
- [License](#license)

## Schematic diagram


## Installation
[Nextflow](https://github.com/nextflow-io/nextflow)  
LncPipe is implemented with Nextflow pipeline management system. To run LncPipe. [Nextflow](https://github.com/nextflow-io/nextflow) should be pre-installed at  POSIX compatible system (Linux, Solaris, OS X, etc), It requires BASH and Java 7 or higher to be installed. We do not recommend running the pipes in the Windows since most of bioinformatic tools are not supported.

## Quick start
Here, we show step by step installation of [Nextflow](https://github.com/nextflow-io/nextflow) in a linux system as an example (adopted from [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html)).

* 1. Download the NextFlow executable package by pasting the following command into your terminal window:


        wget -qO- get.nextflow.io | bash  
        
    
> It will create the [Nextflow](https://github.com/nextflow-io/nextflow) main executable file in the current directory.

* 2. Optionally, move the nextflow file to a directory accessible by your `$PATH` variable (only required to avoid typing the full path to this file each time you need to run it). Of course, you can download the lastest binary version of NextFlow by yourself from [here](https://github.com/nextflow-io/nextflow/releases) and add the path to your system environment.All those pipelines were written in [Nextflow](https://github.com/nextflow-io/nextflow) commands. For more details, please see [here](https://www.nextflow.io).

* 3. Download the LncPipe github repository by:
```
git clone https://github.com/likelet/LncPipe.git
```

* 4. Configure the design.file with experimental conditions and replicate info

* 5. Configure your data and reference files in *nextflow.config* or *docker.config* or *singularity.config*

* 6. Run LncPipe nextflow pipeline:

       nextflow -c nextflow.config run LncRNAanalysisPipe.nf

   or docker command

       nextflow -c docker.config run LncRNAanalysisPipe.nf
   
   or singularity command  
        
       # create image 
       singularity build lncPipe.image docker://bioinformatist/lncpipe
       # run command 
       nextflow -c singularity.config run LncRNAanalysisPipe.nf
* __7.Run with test data __ . 

   PlZ go to https://github.com/likelet/LncPipeTestData 
   
### Prepare input files 

#### References, index and annotation files(Mandatory).
* **:blush:Please keep the consistency of your genome sequence,index library and annotation files (Important!): genome version, chromosome format, gtf coordinated e.g. The dependent third-party softwares may stop for any discrepencies in file-formatting.**
*  Genome reference (genome fasta file with suffix `.fa` etc. )

*  Genome Index for alignment (hisat2 or tophat or STAR)

*  GENCODE gene annotation file in GTF format

*  LNCipedia gene annotation file in GTF format.(set null if not available for your species)

*  Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed



#### Species

    >Currently, LncPipe has been tested for detection of lncRNAs in 'humans' only.
    However, LncPipe can be manually configured to run the anlysis for other species as well and requires additional files  "known_protein_coding.gtf" and  "known_lncRNA.gtf" for coding probability calculations. More information on usage for non-human species can be found here.  

* Reference files for humans 

    1. hisat index built from Genome: 
    http://cancerbio.info/pub/hg38_hisat_index.tar.gz 
    
    2. Genome reference: 
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz 
    
    3. GENCODE gene annotation: 
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz 
     
    4. LNCipedia gene annotation: 
     https://lncipedia.org/downloads/lncipedia_5_0_hc_hg38.gtf  
     
    5. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed

* Reference files for mouse 

    1. hisat index built from Genome
    2. Genome reference:   
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz  
    
    3. GENCODE gene annotation: 
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
    
    4. LNCipedia gene annotation: null
    5. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed

## Run Docker 

1. Prepare input files as mentioned earlier.
2. Modify the `docker.config` in `mandatory` section.
3. Install docker and download the latest LncPipe build using:
       ```
       docker pull bioinformatist/lncpipe
       ```
4. Run LncPipe using the following command:


        nextflow -c docker.config run LncRNAanalysisPipe.nf

>The docker image for LncPipe is available on the docker-hub (https://hub.docker.com/r/bioinformatist/lncpipe/tags/). 
> Alternatively, nextflow can automatically pull image from docker.io. `Dockerfile` recorded  that what we have done with the image. For user from local China looking to pull the docker image can use this [mirror site instead](https://github.com/likelet/Blogs_tips/blob/master/README.md#setting-docker-download-mirror-site).

## Dependencies 

TO Install softwares locally on your machine, please see install instructions [here](https://github.com/likelet/LncPipe/blob/master/InstallSoftwareLocally.md)

## Interactive reports

The results of LncPipe are summarized and visualized via interactive plots by our novel R package [LncPipeReporter](https://github.com/bioinformatist/LncPipeReporter). Users can also try LncPipeReporter as stand-alone for visualizing known and novel lncRNAs.

## Configuration
As a nextflow-based analysis pipeline, LncPipe allow users edit configure file `nextflow.config` to set the index files and default file path parameters instead of typing them into the command line.

To configure, please go to `params` line, and set the following information of various file locations and system environment settings

        params {
        /*
            User setting options (mandatory)
             */
        // input file and genome reference
            fastq_ext = '*_{1,2}.fq.gz'
            fasta_ref = '/data/database/hg38/genome.fa'
            design = 'design.file'
            hisat2_index = '/data/database/hg38/hisatIndex/grch38_snp_tran/genome_snp_tran'
            cpatpath='/opt/CPAT-1.2.3'
            //human gtf only
            gencode_annotation_gtf = "/data/database/hg38/Annotation/gencode.v24.annotation.gtf"
            lncipedia_gtf = "/data/database/hg38/Annotation/lncipedia_4_0_hg38.gtf" // set "null" if you are going to perform analysis on other species

        // additional options for non-human species, else leaving them unchanged
            species="human"// mouse , zebrafish, fly
            known_coding_gtf=""
            known_lncRNA_gtf=""
            //for test
            cpatpath = '/home/zhaoqi/software/CPAT/CPAT-1.2.2/'


        /*
            User setting options (optional)
             */
            // tools setting
            star_idex = ''//set if star used
            bowtie2_index = ''//set if tophat used
            aligner = "hisat" // or "star","tophat"
            sam_processor="sambamba"//or "samtools(deprecated)"
            qctools ="fastp"  // or "afterqc","fastp","fastqc"
            detools = "edger"//or "deseq2","noiseq" not supported yet
            quant = "kallisto"// or 'htseq'

            //other setting
            singleEnd = false
            unstrand = false
            skip_combine = false
            lncRep_Output = 'reporter.html'
            lncRep_theme = 'npg'
            lncRep_cdf_percent = 10
            lncRep_max_lnc_len = 10000
            lncRep_min_expressed_sample = 50
            mem=60
            cpu=30
        }

        manifest {
            homePage = 'https//github.com/likelet/LncPipe'
            description = 'LncPipe:a Nextflow-based Long non-coding RNA analysis PIPELINE'
            mainScript = 'LncRNAanalysisPipe.nf'
        }


        timeline {
            enabled: true
            file: "timeline.html"
        }



## Parameters 
> Those parameters would cover the setting from `nextflow.config` file
* Mandatory(plz configure those options in *nextflow.config* or *docker.config* file)

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
|--star_index/--bowtie2_index/--hisat2_index  | -| Path to STAR?bowtie2/hisat2(mutually exclusive) index(required if not set in config file)  |
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
|--qctools |  `fastp` | Tools for assess raw reads quality or filtered by *fastp*, *fastqc*, *afterqc* or *none*(skip qc step)|

* LncPipeReporter options

| Name | Default value | Description |
|-----------|--------------|-------------|
|--lncRep_Output  | `reporter.html` | Specify report file name. |
|--lncRep_theme  | `npg` | Plot theme setting in interactive plot. Values from [ggsci](https://github.com/road2stat/ggsci) |
|--lncRep_min_expressed_sample  | `50` |Minimum expressed gene allowed in each sample, 50 default. Samples not passed were filtered from analysis|

`--fastq_ext`
> Raw fastq files are required for de-novo analysis.This parameters should be set according to your paired or singled reads file names.

For example:

        Sample1_1.fq.gz
        Sample1_2.fq.gz
        Sample2_1.fq.gz
        Sample2_2.fq.gz

Then you can input pattern `*_{1,2}.fq.gz` to make the all paired-end file recognized by [LncPipe](https://github.com/likelet/LncPipe) .

For singled reads file, file pattern should be fed with `--singleEnd` parameter specified


`--star_idex?--bowtie2_index/--hisat2_index`

> This parameter is *required* when not configured in nextflow.config file. It specify the star/tophat/hisat2(mutually exclusive) index folder built before running [LncPipe](https://github.com/likelet/LncPipe) .
If you don't know what it is?You can use `--fasta` to specify the reference sequence data. The index file would be built by [LncPipe](https://github.com/likelet/LncPipe)  automatically.


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
 `Result` folder under current path(default) or output_folder set by user. A typical structure of `Result` is follows:

        Result/
            ├── QC
            │   ├── N1141_1.clean_fastqc.html
            │   ├── N1141_2.clean_fastqc.html
            │   ├── N1177_1.clean_fastqc.html
            │   └── N1177_2.clean_fastqc.html
            ├── Identified_lncRNA
            │   ├── all_lncRNA_for_classifier.gtf
            │   ├── final_all.fa
            │   ├── final_all.gtf
            │   ├── lncRNA.fa
            │   ├── protein_coding.fa
            │   └── protein_coding.final.gtf
            ├── LncReporter
            │   ├── Differential_Expression_analysis.csv
            │   └── Report.html
            ├── Quantification
            │   ├── kallisto.count.txt
            │   └── kallisto.tpm.txt
            └── Star_alignment
                ├── STAR_N1141
                │   ├── N1141Aligned.sortedByCoord.out.bam
                │   ├── N1141Log.final.out
                │   ├── N1141Log.out
                │   ├── N1141Log.progress.out
                │   └── N1141SJ.out.tab
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

Qi Zhao, Yu Sun, Dawei Wang, Hongwan Zhang, Kai Yu, Jian Zheng, Zhixiang Zuo. LncPipe: A Nextflow-based pipeline for identification and analysis of long non-coding RNAs from RNA-Seq data. Journal of Genetics and Genomics. 2018. (_In press_)
