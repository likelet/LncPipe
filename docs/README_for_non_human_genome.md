# Example usage for non-human species 

## Introduction 
LncPipe accepts raw reads, annotations and genome reference as input to conduct the whole analysis, 
and is also applicable for selected referenced species. In the first version, we mainly focus on human 
because of well-organized lncRNA annotation files (with .gtf suffixed) from GENCODE and LNCipedia. 
For non-human species, user are required to provide both protein coding annotation file and lncRNA annotation file separately, 
which could be download from GENCODE(human or mouse only) or Ensemble databases. However, not all non-human species are supported 
at present, since one essential tool CPAT included in LncPipe only available for 4 species (human, mouse, fly and zebrafish). 

## Example usage for mouse 
### Step 1. Prepare input files. The following files required by LncPipe
* hisat index: ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz
```shell
    #if database not exsit 
    mkdir ~/database/mouse
    cd ~/database/mouse  
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz
    tar -xzvf grcm38_tran.tar.gz
```
* Genome reference: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz
```shell
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz
    gunzip GRCm38.p5.genome.fa.gz
```
* Long non-coding RNA GTF files:ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.long_noncoding_RNAs.gtf.gz
```shell
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.long_noncoding_RNAs.gtf.gz
    gunzip gencode.vM16.long_noncoding_RNAs.gtf.gz
```
* GTF files including all features :ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
```shell
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
    gunzip gencode.vM16.annotation.gtf.gz
```
* Protein coding GTF file
> This file should be generated from gtf file that contains all features. Users can run following code to get this file:
```shell
    # in the same folder 
    zcat   |grep "protein_coding" > gencode.vM16.proteincoding.gtf
```
* Raw sequence file with \*.fastq.gz / \*.fq.gz suffixed
> For example there are 4 samples with eight .gz suffixed fastq files, like
```shell
sample1_rep1_1.fq.gz,
sample1_rep1_2.fq.gz,
sample1_rep2_1.fq.gz,
sample1_rep2_2.fq.gz,
sample2_rep1_1.fq.gz,
sample2_rep1_2.fq.gz,
sample2_rep2_1.fq.gz,
sample2_rep2_2.fq.gz
```

* `Design.file` are required if you are going to perform differential expression analysis between groups. 

### Step 2. Edit `nextflow.config` file 
Leave the other line unchanged, modified the following sentences like below (According to the location of download files):
> If your are using docker in your server, plz modify `docker.config` instead.  

```shell
    fastq_ext = '*_{1,2}.fq.gz'
    fasta_ref = '~/database/mouse/GRCm38.p5.genome.fa'
    design = 'design.file'
    hisat2_index = '~/database/mouse/grcm38_tran/grcm38_tran'
    cpatpath='/opt/CPAT-1.2.3'
    //human gtf only
    gencode_annotation_gtf = "~/database/mouse/gencode.v24.annotation.gtf"
    lncipedia_gtf = null
    species="mouse"// mouse , zebrafish, fly
    known_coding_gtf="gencode.vM16.proteincoding.gtf"
    known_lncRNA_gtf="gencode.vM16.annotation.gtf"

```
### Step 3. Start your analysis trip with command below   

```shell
    nextflow run -with-trace -with-report report.html -with-timeline timeline.html LncRNAanalysisPipe.nf 
    #or running in a docker image  
    nextflow -c docker.config  LncRNAanalysisPipe.nf 
```
> The default running tools in each step are fastp, hisat, gffcompare, stringtie, cpat, plek, sambamba, kallisto ,edgeR and LncPipeReporter, if you want to change the tool in each step, plz modify `config` file instead.

* Any question, plz open an issue in the issue page, we will reply ASAP :)
