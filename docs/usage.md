# nf-core/lncpipe: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`docker`](#docker)
        * [`awsbatch`](#awsbatch)
        * [`standard`](#standard)
        * [`none`](#none)
    * [`--reads`](#--reads)
    * [`--singleEnd`](#--singleend)
    * [`--unstrand`](#--unstrand)
* [Reference Genomes](#reference-genomes)
    * [`--genome`](#--genome)
    * [`--fasta`](#--fasta)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_emails`](#--plaintext_emails)
    * [`--sampleLevel`](#--sampleLevel)
    * [`--multiqc_config`](#--multiqc_config)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run  LncRNAanalysisPipe.nf  -profile docker,test --fastq_ext '*_R{1,2}.fastq.gz'
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
Result         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Prepare Input file 

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
    
    6. `design` file (optional)

* Reference files for mouse 

    1. hisat index built from Genome  
    
    2. Genome reference:   
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz  
    
    3. GENCODE gene annotation: 
    ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
    
    4. Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed  
    
    5. `design` file (optional)



### Reproducibility
As a nextflow-based analysis pipeline, LncPipe allow users edit configure file `config` to set the index files and default file path parameters instead of typing them into the command line.

To configure, please copy the `conf/temp.config` into a new config file in `conf` folder. For example, we named it to `run.config`. Then open the file, and go to `params` line, and set the following information of various file locations and system environment settings

```groovy
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

```

After the create your own `run.config` file, please go to the `nextflow.config` file to add your profile in the main config.   
In the `profile` chunk, add your `run.config` information to the end as follows:  

```groovy 
profiles {

  standard {
    includeConfig 'conf/base.config'
  }
  docker {
    includeConfig 'conf/base.config'
    includeConfig 'conf/docker.config'
  }
  
  singularity {
    includeConfig 'conf/base.config'
    includeConfig 'conf/sing.config'
  }
  
  test {
    includeConfig 'conf/base.config'
    includeConfig 'conf/test.config'
  }
  
   run {
        includeConfig 'conf/base.config'
        includeConfig 'conf/run.config'
   }

  debug { process.beforeScript = 'echo $HOSTNAME' }
  none {
    // Don't load any config (for use with custom home configs)
  }

}
```

Then you can run the pipeline with your own data  the following command : 
 PlZ go to https://github.com/likelet/LncPipeTestData 

```bash
 nextflow run  LncRNAanalysisPipe.nf -profile standard,docker,run
```


## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub.
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
    * NOTE: in lncPipe, test data should be downloaded separately, plz see details in [here]()
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.


## Reference Genomes

The lncPipe currently support `human` species only, other species will be supported in the near feature.  

We recommended users adopted the reference file from [GENCODE](https://www.gencodegenes.org/)


## Annotation file 

### `--gencode_annotation_gtf`  
An annotation file from GENCODE database for annotating lncRNAs(required if not set in config file). e.g. [gencode.v26.annotation.gtf](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz)   

### `--lncipedia_gtf`    
An annotation file from LNCipedia database for annotating lncRNAs(required if not set in config file) e.g. [lncipedia_4_0_hc_hg38.gtf](http://www.lncipedia.org/downloads/lncipedia_4_0_hc_hg38.gtf) |

### `--star_index/--bowtie2_index/--hisat2_index`

> This parameter is *required* when not configured in nextflow.config file. It specify the star/tophat/hisat2(mutually exclusive) index folder built before running [LncPipe](https://github.com/likelet/LncPipe) .
If you don't know what it is?You can use `--fasta` to specify the reference sequence data. The index file would be built by [LncPipe](https://github.com/likelet/LncPipe)  automatically.


### `--fasta`
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--design`
> Experimental design file matrix for differential expression analysis. Default: `null`
Format:

    WT:Sample1,Sample2,Sample3
    KO:Sample1,Sample2,Sample3

While `KO/WT` represents the two experimental condition, and sample1, sample2, sample3 are replicates which should be comma-delimited in the same line .

For sample names, it should be the sample as the prefix of fastq files which was trimmed by `--fastq_ext`.

For example:

 if fastq file names are `Sample1_1.fq.gz, Sample1_2.fq.gz` that comes from one sample and your `--fastq_ext` is set as `*_{1,2}.fq.gz`, the sample name
should be Sample1.


## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
