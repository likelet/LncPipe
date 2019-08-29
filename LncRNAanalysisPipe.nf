#!/usr/bin/env nextflow

/*
 * LncPipe was implemented by Dr. Qi Zhao from Sun Yat-sen University Cancer Center, China.
 *
 *
 *   LncPipe is a free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *      See the GNU General Public License for more details.
 *
 *
 */

/*
 * LncPipe: A nextflow-based lncRNA identification and analysis pipeline from RNA sequencing data
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>: design and implement the pipeline.
 * Yu Sun <sun_yu@mail.nankai.edu.cn>: design and implement the analysis report sections.
 * Zhixiang Zuo <zuozhx@sysucc.org.cn>: design the project and perform the testing.
 */

// requirement:
// - fastp/fastqc／AfterQC
// - STAR/tophat2/bowtie2/hisat2/StringTie
// - samtools/sambamba
// - Cufflinks/gffcompare
// - Bedops
// - CPAT
// - PLEK
// - CNCI
// - kallisto [https://pachterlab.github.io/kallisto/starting]

//pre-defined functions for render command
//=======================================================================================
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";


def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }

//Help information
// Nextflow  version
version="v0.2.45"
//=======================================================================================
// Nextflow Version check
if( !nextflow.version.matches('0.30+') ) {
    println print_yellow("This workflow requires Nextflow version 0.26 or greater -- You are running version ")+ print_red(nextflow.version)
// exit 1
}
//help information
params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "LncPipe: a Nextflow-based Long non-coding RNA analysis Pipeline v$version"
    log.info "LncPipe integrates several NGS processing tools to identify novel long non-coding RNAs from"
    log.info "un-processed RNA sequencing data. To run this pipeline, users either need to install required tools manually"
    log.info "or use the docker image for LncPipe that comes with all tools pre-installed. (note: docker needs to be installed on your system). More information on usage can be found at https://github.com/likelet/LncPipe ."
    log.info "Bugs or new feature requests can be reported by opening issues in our github repository."
    log.info print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') +
            print_purple('       Nextflow run LncRNAanalysisPipe.nf \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --input_folder <path>         ') + print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --fastq_ext <*_fq.gz>         ') + print_green('Filename pattern for pairing raw reads, e.g: *_{1,2}.fastq.gz for paired reads\n') +
            print_cyan('      --out_folder <path>           ') + print_green('The output directory where the results will be saved(optional), current path is default\n') +
            print_cyan('      --aligner <hisat>             ') + print_green('Aligner for reads mapping (optional),"hisat"(defalt)/"star"/"tophat"\n') +
            print_cyan('      --qctools <fastp>            ') + print_green('Tools for assess reads quality, fastp(default)/afterqc/fastqc/none(skip QC step)\n') +
            print_cyan('      --detools <edger>             ') + print_green('Tools for differential analysis, edger(default)/deseq/noiseq\n') +
            print_cyan('      --quant <kallisto>            ') + print_green('Tools for estimating abundance of transcript, kallisto(default)/htseq\n') +
            '\n' +
            print_yellow('    Options:                         General options for run this pipeline\n') +
            print_cyan('      --merged_gtf <gtffile>        ') + print_green('Start analysis with assemblies already produced and skip fastqc/alignment step, DEFAOUL NULL\n') +
            print_cyan('      --design <file>               ') + print_green('A flat file stored the experimental design information ( required when perform differential expression analysis)\n') +
            print_cyan('      --singleEnd                   ') + print_green('Reads type, True for single ended \n') +
            print_cyan('      --unstrand                    ') + print_green('RNA library construction strategy, specified for \'unstranded\' library \n') +
            '\n' +
            print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.\n') +
            print_cyan('      --fasta                       ') + print_green('Path to Fasta reference(required)\n') +
            print_cyan('      --gencode_annotation_gtf      ') + print_green('An annotation file from GENCODE database in GTF format (required)\n') +
            print_cyan('      --lncipedia_gtf               ') + print_green('An annotation file from LNCipedia database in GTF format (required)\n') +
            '\n' +
            print_yellow('    LncPipeReporter Options:         LncPipeReporter setting  \n') +
            print_cyan('      --lncRep_Output                ') + print_green('Specify report file name, \"report.html\" default.\n') +
            print_cyan('      --lncRep_theme                 ') + print_green('Plot theme setting in interactive plot, \"npg\" default.\n') +
            print_cyan('      --lncRep_min_expressed_sample  ') + print_green('Minimum expressed gene allowed in each sample, 50 default.\n') +
            '\n' +
            print_yellow('    Other options:                   Specify the email and \n') +
            print_cyan('      --sam_processor                ') + print_green('program to process samfile generated by hisat2 if aligner is hisat2. Default \"sambamba\". \n') +
            print_cyan('      --mail                         ') + print_green('email info for reporting status of your LncPipe execution  \n') +



            log.info '------------------------------------------------------------------------'
    log.info print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    log.info print_yellow('Copyright (c) 2013-2017, Sun Yat-sen University Cancer Center.')
    log.info '------------------------------------------------------------------------'
    exit 0
}

//check parameters
/*
allowed_params = ["input_folder","fastq_ext","out_folder","aligner","qctools","detools","quant",
                   "merged_gtf","design","singleEnd","unstrand",
                   "fasta","gencode_annotation_gtf","lncipedia_gtf",
                   "lncRep_Output", "lncRep_theme","lncRep_min_expressed_sample",
                   "sam_processor","mail"]
params.each { entry ->
    if (! allowed_params.contains(entry.key)) {
        println("The parameter <${entry}.key> is not known");
        System.exit(2);
    }
}
*/

//default values
params.input_folder = './'
params.out_folder = './'
params.multiqc_config = "$baseDir/assets/multiqc_config.yaml" // for generate qc and alignment result
params.merged_gtf = null// dose merged_gtf provided
singleEnd = params.singleEnd ? true : false
skip_combine = params.skip_combine ? true : false
unstrand = params.unstrand ? true : false
//Checking parameters
log.info print_purple("You are running LncPipe with the following parameters:")
log.info print_purple("Checking parameters ...")
log.info print_yellow("=====================================")
log.info print_yellow("Species:                        ") + print_green(params.species)
log.info print_yellow("Fastq file extension:           ") + print_green(params.fastq_ext)
log.info print_yellow("Single end :                    ") + print_green(params.singleEnd)
log.info print_yellow("skip annotation process:        ") + print_green(params.skip_combine)
log.info print_yellow("Input folder:                   ") + print_green(params.input_folder)
log.info print_yellow("Output folder:                  ") + print_green(params.out_folder)
log.info print_yellow("Genome sequence location:       ") + print_green(params.fasta_ref)
log.info print_yellow("STAR index path:                ") + print_green(params.star_index)
log.info print_yellow("HISAT2 index path:               ") + print_green(params.hisat2_index)
log.info print_yellow("bowtie/tophat index path:       ") + print_green(params.bowtie2_index)
if(params.species=="human"){
    log.info print_yellow("GENCODE annotation location:    ") + print_green(params.gencode_annotation_gtf)
    log.info print_yellow("lncipedia annotation location:  ") + print_green(params.lncipedia_gtf)
}else{
    log.info print_yellow("known protein coding annotation location:    ") + print_green(params.known_coding_gtf)
    log.info print_yellow("known lncRNA annotation location:            ") + print_green(params.known_lncRNA_gtf)
}

log.info print_yellow("=====================================")
log.info "\n"

// run information of system file
//automatic set optimize resource for analysis based on current system resources
ava_mem = (double) (Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()
if (params.cpu != null && ava_cpu > params.cpu) {
    ava_cpu = params.cpu
} else if (params.cpu != null && ava_cpu < params.cpu) {
    print print_yellow("Exceeding the max available processors, \n use default parameter to run pipe. ")
}
if (params.mem != null && ava_mem > params.mem) {
    ava_mem = params.mem
} else if (params.mem != null && ava_mem < params.mem) {
    print print_yellow("Exceeding the max available memory, \n use default parameter to run pipe. ")
}
// set individual cpu for fork run
idv_cpu = 40
int fork_number = ava_cpu / idv_cpu
if (fork_number < 1) {
    fork_number = 1
}

// read file
fasta_ref = file(params.fasta_ref)
if (!fasta_ref.exists()) exit 1, "Reference genome not found: ${params.fasta_ref}"
if(params.aligner=='star'){
    star_index = file(params.star_index)
    if (!star_index.exists()) exit 1, "STAR index not found: ${params.star_index}"
}else if(params.aligner =='hisat'){
    hisat2_index = Channel.fromPath("${params.hisat2_index}*")
            .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}else if(params.aligner =='tophat'){
    bowtie2_index = Channel.fromPath("${params.bowtie2_index}*")
            .ifEmpty { exit 1, "bowtie2 index for tophat not found: ${params.bowtie2_index}" }
}

input_folder = file(params.input_folder)
multiqc_config = file(params.multiqc_config)

/*
*Step 1: Prepare Annotations
 */
 

if (params.species=="human") {
    println print_purple("Combining known annotations from GTFs")
    gencode_annotation_gtf = file(params.gencode_annotation_gtf)
    if (!gencode_annotation_gtf.exists()) exit 1, "GENCODE annotation file not found: ${params.gencode_annotation_gtf}"
    lncipedia_gtf = file(params.lncipedia_gtf)
    if (!lncipedia_gtf.exists()) exit 1, "lncipedia annotation file not found: ${params.lncipedia_gtf}"
//Prepare annotations
    annotation_channel = Channel.from(gencode_annotation_gtf, lncipedia_gtf)
    annotation_channel.collectFile { file -> ['lncRNA.gtflist', file.name + '\n'] }
            .set { LncRNA_gtflist }

    process combine_public_annotation {
        storeDir { params.out_folder + "/Combined_annotations" }
        input:
        file lncRNA_gtflistfile from LncRNA_gtflist
        file gencode_annotation_gtf
        file lncipedia_gtf
        output:
        file "gencode_protein_coding.gtf" into proteinCodingGTF, proteinCodingGTF_forClass
        file "known.lncRNA.gtf" into KnownLncRNAgtf
        file "*_mod.gtf" into mod_file_for_rename
        file "gencode_annotation_gtf_mod.gtf" into gencode_annotation_gtf_for_assemble,gencode_annotation_gtf_for_merge,gencode_annotation_gtf_for_filter

        shell:
        cufflinks_threads = ava_cpu- 1


        if(params.aligner=='hisat'){//fix the gtf format required by hisat
            '''
            set -o pipefail
            touch filenames.txt
            
            perl -lpe 's/ ([^"]\\S+) ;/ "$1" ;/g' !{gencode_annotation_gtf} > gencode_annotation_gtf_mod.gtf 
            perl -lpe 's/ ([^"]\\S+) ;/ "$1" ;/g' !{lncipedia_gtf} > lncipedia_mod.gtf 
            
            echo  gencode_annotation_gtf_mod.gtf   >>filenames.txt
            echo lncipedia_mod.gtf   >>filenames.txt
           
            stringtie --merge -o merged_lncRNA.gtf  filenames.txt
            cat gencode_annotation_gtf_mod.gtf   |grep "protein_coding" > gencode_protein_coding.gtf
            gffcompare -r gencode_protein_coding.gtf -p !{cufflinks_threads} merged_lncRNA.gtf
            awk '$3 =="u"||$3=="x"{print $5}' gffcmp.merged_lncRNA.gtf.tmap |sort|uniq|perl !{baseDir}/bin/extract_gtf_by_name.pl merged_lncRNA.gtf - > merged.filter.gtf
            mv  merged.filter.gtf known.lncRNA.gtf
        
            '''
        }
        else {

            '''
            set -o pipefail
            
            cuffmerge -o merged_lncRNA  !{lncRNA_gtflistfile}
            cat !{gencode_annotation_gtf} |grep "protein_coding" > gencode_protein_coding.gtf
            cuffcompare -o merged_lncRNA -r gencode_protein_coding.gtf -p !{cufflinks_threads} merged_lncRNA/merged.gtf
            awk '$3 =="u"||$3=="x"{print $5}' merged_lncRNA/merged_lncRNA.merged.gtf.tmap  |sort|uniq|perl !{baseDir}/bin/extract_gtf_by_name.pl merged_lncRNA/merged.gtf - > merged.filter.gtf
            mv  merged.filter.gtf known.lncRNA.gtf
        
        '''
        }
    }
} else {// for mouse or other species, user should provide known_protein_coding and known_lncRNA GTF file for analysis
    println print_purple("skip combination step due to non-human analysis")
    Known_LncRNAgtf=file(params.known_lncRNA_gtf)
  //  if (!KnownLncRNAgtf.exists()) exit 1, print_red("In non-human mode, known lncRNA GTF annotation file not found: ${params.known_lncRNA_gtf}")
    known_coding_gtf=file(params.known_coding_gtf)
    //if (!known_coding_gtf.exists()) exit 1, print_red("In non-human mode, known protein coding GTF annotation file not found: ${params.known_coding_gtf}")

    process process_non_human_annotation{
        storeDir { params.out_folder + "/Combined_annotations" }
        input:
        file Known_LncRNAgtf
        file known_coding_gtf
        output:
        file "know_lnc.gtf" into KnownLncRNAgtf
        file "known_coding.gtf" into proteinCodingGTF, proteinCodingGTF_forClass
        file "non_human_mod.gtf" into mod_file_for_rename,gencode_annotation_gtf_for_assemble,gencode_annotation_gtf_for_merge,gencode_annotation_gtf_for_filter
        shell:
        '''
        set -o pipefail
        # add quote and remove gene terms avoiding malformed error by stringtie and bedops
         perl -lpe 's/ ([^"]\\S+) ;/ "$1" ;/g'  !{Known_LncRNAgtf} | grep -w gene -v > know_lnc.gtf
         perl -lpe 's/ ([^"]\\S+) ;/ "$1" ;/g' !{known_coding_gtf} | grep -w gene -v > known_coding.gtf 
         cat know_lnc.gtf known_coding.gtf  > non_human_mod.gtf
        
        
        '''
    }





}


// whether the merged gtf have already produced.
if (!params.merged_gtf) {
    /*
     * Step 2: Build read aligner (STAR/tophat/HISAT2) index, if not provided
     */
    //star_index if not exist
    /*if (params.aligner == 'star' && params.star_index == false && fasta_ref) {
        process Make_STARindex {
            tag fasta_ref

            storeDir { params.out_folder + "/STARIndex" }

            input:
            file fasta_ref from fasta_ref
            file gencode_annotation_gtf

            output:
            file "star_index" into star_index

            shell:
            star_threads = ava_cpu- 1
            """
                mkdir star_index
                STAR \
                    --runMode genomeGenerate \
                    --runThreadN ${star_threads} \
                    --sjdbGTFfile $gencode_annotation_gtf \
                    --sjdbOverhang 149 \
                    --genomeDir star_index/ \
                    --genomeFastaFiles $fasta_ref
                """
        }
    } else if (params.aligner == 'star' && params.star_index == false && !fasta_ref) {
        println print_red("No reference fasta sequence loaded! please specify ") + print_red("--fasta_ref") + print_red(" with reference.")

    } else if (params.aligner == 'tophat' && params.bowtie2_index == false && !fasta_ref) {
        process Make_bowtie2_index {
            
            tag fasta_ref
            storeDir { params.out_folder + "/bowtie2Index" }

            input:
            file fasta_ref from fasta_ref

            output:
            file "genome_bt2.*" into bowtie2_index

            shell:
            """
                bowtie2-build !{fasta_ref} genome_bt2
                """
        }
    } else if (params.aligner == 'tophat' && !fasta_ref) {
        println print_red("No reference fasta equence loaded! please specify ") + print_red("--fasta_ref") + print_red(" with reference.")
    } else if (params.aligner == 'hisat' && !fasta_ref) {
        process Make_hisat_index {
            
            tag fasta_ref

            storeDir { params.out_folder + "/hisatIndex" }

            input:
            file fasta_ref from fasta_ref
            file gencode_annotation_gtf

            output:
            file "genome_ht2.*" into hisat2_index

            shell:
            hisat2_index_threads = ava_cpu- 1
            """
                #for human genome it will take more than 160GB memory and take really  long time (6 more hours), thus we recommand to down pre-build genome from hisat website
                extract_splice_sites.py !{gencode_annotation_gtf} >genome_ht2.ss
                extract_exons.py !{gencode_annotation_gtf} > genome_ht2.exon
                hisat2-build -p !{hisat2_index_threads} --ss genome_ht2.ss --exo genome_ht2.exon !{fasta_ref} genome_ht2
                """
        }
    } else if (params.aligner == 'tophat' && params.hisat_index == false && !fasta_ref) {
        println print_red("No reference fasta sequence loaded! please specify ") + print_red("--fasta_ref") + print_red(" with reference.")
    }*/

    println print_purple("Analysis from fastq file")
    //Match the pairs on two channels

    reads = params.input_folder + params.fastq_ext

    /*
    * Step 3: QC (FastQC/AfterQC/Fastp) of raw reads
    */
    println print_purple("Perform quality control of raw fastq files ")
    if (params.qctools == 'fastqc') {
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\n")
        }
        .into { reads_for_fastqc; readPairs_for_discovery;readPairs_for_kallisto}
        process Run_fastQC {
            tag { fastq_tag }
            label 'qc'

            publishDir pattern: "*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "*.html" into fastqc_for_waiting
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            '''
            fastqc -t !{fastq_threads} !{fastq_file[0]} !{fastq_file[1]}
        '''
        }
    }
    else if (params.qctools == 'afterqc'){
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }.set { reads_for_fastqc}
        process Run_afterQC {

            tag { fastq_tag }
            label 'qc'
            publishDir pattern: "QC/*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "QC/*.html" into fastqc_for_waiting
            set val(fastq_tag), file('*.good.fq.gz')  into readPairs_for_discovery,readPairs_for_kallisto
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            if (params.singleEnd) {
                '''
            after.py -z -1 !{fastq_file[0]} -g ./
            '''
            } else {
                '''
            after.py -z -1 !{fastq_file[0]} -2 !{fastq_file[1]} -g ./
            '''
            }
        }
    }
    else if (params.qctools == 'fastp'){
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }
        .set { reads_for_fastqc}
        process Run_FastP {

            tag { fastq_tag }
            label 'qc'

            publishDir pattern: "*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "*.html" into fastqc_for_waiting
            set val(fastq_tag), file('*qc.fq.gz')  into readPairs_for_discovery,readPairs_for_kallisto
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            if (params.singleEnd) {
                '''
            fastp -i !{fastq_file[0]} -o !{samplename}.qc.fq.gz  -h !{samplename}_fastp.html
           
            '''
            } else {
                '''
            fastp -i !{fastq_file[0]}  -I !{fastq_file[1]} -o !{samplename}_1.qc.fq.gz  -O !{samplename}_2.qc.fq.gz -h !{samplename}_fastp.html
            '''
            }
        }
    }else{
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }
                .into{readPairs_for_discovery; readPairs_for_kallisto;fastqc_for_waiting}
    }
    fastqc_for_waiting = fastqc_for_waiting.first()

    /*
    * Step 4: Initialize read alignment (STAR/HISAT2/tophat)
    */
    if (params.aligner == 'star') {
        process fastq_star_alignment_For_discovery {
            
            tag { file_tag }

            publishDir pattern: "",
                    path: { params.out_folder + "/Result/Star_alignment" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(pair) from readPairs_for_discovery
            file tempfiles from fastqc_for_waiting // just for waiting
            file fasta_ref
            file star_index

            output:
            set val(file_tag_new), file("${file_tag_new}Aligned.sortedByCoord.out.bam") into mappedReads,forHtseqMappedReads
            file "${file_tag_new}Log.final.out" into alignment_logs
            shell:
            println print_purple("Start mapping with STAR aligner " + samplename)
            file_tag = samplename
            file_tag_new = file_tag
            star_threads = ava_cpu - 1

            if (params.singleEnd) {
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in single-end mode")
                """
                         STAR --runThreadN !{star_threads} \
                            --twopassMode Basic \
                            --genomeDir !{star_index} \
                            --readFilesIn !{pair} \
                            --readFilesCommand zcat \
                            --outSAMtype BAM SortedByCoordinate \
                            --chimSegmentMin 20 \
                            --outFilterIntronMotifs RemoveNoncanonical \
                            --outFilterMultimapNmax 20 \
                            --alignIntronMin 20 \
                            --alignIntronMax 1000000 \
                            --alignMatesGapMax 1000000 \
                            --outFilterType BySJout \
                            --alignSJoverhangMin 8 \
                            --alignSJDBoverhangMin 1 \
                            --outFileNamePrefix !{file_tag_new} 
                    """
            } else {
                println print_purple("Initial reads mapping of " + samplename + " performed by STAR in paired-end mode")
                '''
                            STAR --runThreadN !{star_threads}  \
                                 --twopassMode Basic --genomeDir !{star_index} \
                                 --readFilesIn !{pair[0]} !{pair[1]} \
                                 --readFilesCommand zcat \
                                 --outSAMtype BAM SortedByCoordinate \
                                 --chimSegmentMin 20 \
                                 --outFilterIntronMotifs RemoveNoncanonical \
                                 --outFilterMultimapNmax 20 \
                                 --alignIntronMin 20 \
                                 --alignIntronMax 1000000 \
                                 --alignMatesGapMax 1000000 \
                                 --outFilterType BySJout \
                                 --alignSJoverhangMin 8 \
                                 --alignSJDBoverhangMin 1 \
                                 --outFileNamePrefix !{file_tag_new} 
                    '''
            }
        }
    }
    else if (params.aligner == 'tophat')
    {
        process fastq_tophat_alignment_For_discovery {
            
            tag { file_tag }

            publishDir pattern: "",
                    path: { params.out_folder + "/Result/tophat_alignment" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(pair) from readPairs_for_discovery
            file tempfiles from fastqc_for_waiting // just for waiting
            file fasta_ref
            file bowtie2_index from bowtie2_index.collect()
            file gtf from gencode_annotation_gtf

            output:
             set val(samplename),file("${file_tag_new}_thout/accepted.bam") into mappedReads,forHtseqMappedReads
            file "${file_tag_new}_thout/Alignment_summary.txt" into alignment_logs
            //align_summary.txt as log file
            shell:
            println print_purple("Start mapping with tophat2 aligner " + samplename)
            file_tag = samplename
            file_tag_new = file_tag
            tophat_threads = ava_cpu- 1
            index_base = bowtie2_index[0].toString() - ~/.\d.bt2/
            strand_str="fr-firststrand"
            if(unstrand){
                strand_str="fr-unstranded"
            }
            if (params.singleEnd) {
                println print_purple("Initial reads mapping of " + samplename + " performed by Tophat in single-end mode")
                '''
                         tophat -p !{tophat_threads} -G !{gtf} -–no-novel-juncs -o !{samplename}_thout --library-type !{strand_str} !{index_base} !{pair} 
                         
                '''
            } else {
                println print_purple("Initial reads mapping of " + samplename + " performed by Tophat in paired-end mode")
                '''
                     tophat -p !{tophat_threads} -G !{gtf} -–no-novel-juncs -o !{samplename}_thout --library-type !{strand_str} !{index_base} !{pair[0]} !{pair[1]} 
                '''
            }
        }
    }
    else if (params.aligner == 'hisat') {
        process fastq_hisat2_alignment_For_discovery {
            
            tag { file_tag }
            label 'para'
            publishDir pattern: "",
                    path: { params.out_folder + "/Result/hisat_alignment" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(pair) from readPairs_for_discovery
            file tempfiles from fastqc_for_waiting // just for waiting
            file fasta_ref
            file hisat2_id from hisat2_index.collect()

            output:
            set val(file_tag_new),file("${file_tag_new}.sort.bam") into hisat_mappedReads,forHtseqMappedReads
            file "${file_tag_new}.hisat2_summary.txt" into alignment_logs
            //align_summary.txt as log file
            shell:
            println print_purple("Start mapping with hisat2 aligner " + samplename)
            file_tag = samplename
            file_tag_new = file_tag
            hisat2_threads = ava_cpu- 2
            index_base = hisat2_id[0].toString() - ~/.\d.ht2/

            if(unstrand){
                if (params.singleEnd) {
                    println print_purple("Initial reads mapping of " + samplename + " performed by hisat2 in single-end mode")
                    '''
                   mkdir tmp
                   hisat2  -p !{hisat2_threads} --dta  -x  !{index_base}  -U !{pair}  -S !{file_tag_new}.sam 2>!{file_tag_new}.hisat2_summary.txt
                  sambamba view -S -f bam -t !{hisat2_threads} !{file_tag_new}.sam -o temp.bam 
                  sambamba sort -o !{file_tag_new}.sort.bam --tmpdir ./tmp -t !{hisat2_threads} temp.bam
                  rm !{file_tag_new}.sam
                  rm temp.bam
                  
                '''
                } else {
                    println print_purple("Initial reads mapping of " + samplename + " performed by hisat2 in paired-end mode")
                    '''
                    mkdir tmp
                  hisat2  -p !{hisat2_threads} --dta  -x  !{index_base}  -1 !{pair[0]}  -2 !{pair[1]}  -S !{file_tag_new}.sam 2> !{file_tag_new}.hisat2_summary.txt
                  sambamba view -S -f bam -t !{hisat2_threads} !{file_tag_new}.sam -o temp.bam
                  sambamba sort -o !{file_tag_new}.sort.bam --tmpdir ./tmp -t !{hisat2_threads} temp.bam
                  rm !{file_tag_new}.sam
                '''
                }
            }else {
                if (params.singleEnd) {
                    println print_purple("Initial reads mapping of " + samplename + " performed by hisat2 in single-end mode")
                    '''
                   mkdir tmp
                   hisat2  -p !{hisat2_threads} --dta --rna-strandness !{params.hisat_strand} -x  !{index_base}  -U !{pair}  -S !{file_tag_new}.sam 2>!{file_tag_new}.hisat2_summary.txt
                  sambamba view -S -f bam -t !{hisat2_threads} !{file_tag_new}.sam -o temp.bam 
                  sambamba sort -o !{file_tag_new}.sort.bam --tmpdir ./tmp -t !{hisat2_threads} temp.bam
                  rm !{file_tag_new}.sam
                  rm temp.bam
                  
                '''
                } else {
                    println print_purple("Initial reads mapping of " + samplename + " performed by hisat2 in paired-end mode")
                    '''
                 mkdir tmp
                  hisat2  -p !{hisat2_threads} --dta --rna-strandness !{params.hisat_strand} -x  !{index_base}  -1 !{pair[0]}  -2 !{pair[1]}  -S !{file_tag_new}.sam 2> !{file_tag_new}.hisat2_summary.txt
                  sambamba view -S -f bam -t !{hisat2_threads} !{file_tag_new}.sam -o temp.bam
                  sambamba sort -o !{file_tag_new}.sort.bam --tmpdir ./tmp -t !{hisat2_threads} temp.bam
                  rm !{file_tag_new}.sam
                '''
                }
            }
        }
    }

    /*
    * Step 5: Transcript assembly using Stringtie
    */
    if(params.aligner == 'hisat'){
        process StringTie_assembly {
            
            tag { file_tag }

            input:
            set val(samplename),file(alignment_bam) from hisat_mappedReads
            file fasta_ref
            file gencode_annotation_gtf from gencode_annotation_gtf_for_assemble

            output:

            file "stringtie_${file_tag_new}_transcripts.gtf" into stringTieoutgtf, StringTieOutGtf_fn

            shell:
            file_tag = samplename
            file_tag_new = file_tag
            stringtie_threads = ava_cpu- 2

            if(unstrand){
                '''
            #run stringtie
            stringtie -p !{stringtie_threads} -G !{gencode_annotation_gtf} -l stringtie_!{file_tag_new} -o stringtie_!{file_tag_new}_transcripts.gtf !{alignment_bam}
            '''
            }else{
                '''
            #run stringtie
            stringtie -p !{stringtie_threads} -G !{gencode_annotation_gtf} --rf -l stringtie_!{file_tag_new} -o stringtie_!{file_tag_new}_transcripts.gtf !{alignment_bam}
            '''
            }

        }
// Create a file 'gtf_filenames' containing the filenames of each post processes cufflinks gtf
        stringTieoutgtf.collectFile { file -> ['gtf_filenames.txt', file.name + '\n'] }
                .set { GTFfilenames }
        /*
        * Step 6: Merged GTFs into one
        */
        process StringTie_merge_assembled_gtf {
            
            tag { file_tag }
            label 'para'
            publishDir pattern: "merged.gtf",
                    path: { params.out_folder + "/Result/Merged_assemblies" }, mode: 'copy', overwrite: true

            input:
            file gtf_filenames from GTFfilenames
            file cufflinksgtf_file from StringTieOutGtf_fn.toList() // not used but just send the file in current running folder
            file fasta_ref


            output:
            file "merged.gtf" into mergeTranscripts_forCompare, mergeTranscripts_forExtract, mergeTranscripts_forCodeingProtential
            shell:

            stringtie_threads = ava_cpu- 1

            '''
            stringtie --merge -p !{stringtie_threads} -o merged.gtf !{gtf_filenames}
            
            
            '''
        }
    }
    else{
        process Cufflinks_assembly {
            
            tag { file_tag }

            input:
            set val(file_tag), file(alignment_bam) from mappedReads
            file fasta_ref
            file gencode_annotation_gtf from gencode_annotation_gtf_for_assemble

            output:

            file "Cufout_${file_tag_new}_transcripts.gtf" into cuflinksoutgtf, cuflinksoutgtf_fn

            shell:
            file_tag_new = file_tag
            cufflinks_threads = ava_cpu- 1
            strand_str="fr-firststrand"
            if(unstrand){
                strand_str="fr-unstranded"
            }
            if (params.aligner == 'tophat') {
                '''
            #run cufflinks
            
            cufflinks -g !{gencode_annotation_gtf} \
                      -b !{fasta_ref} \
                      --library-type !{strand_str}\
                      --max-multiread-fraction 0.25 \
                      --3-overhang-tolerance 2000 \
                      -o Cufout_!{file_tag_new} \
                      -p !{cufflinks_threads} !{alignment_bam}
                      
            mv Cufout_!{file_tag_new}/transcripts.gtf Cufout_!{file_tag_new}_transcripts.gtf
            '''

            } else if (params.aligner == 'star') {
                '''
            #run cufflinks
            
            cufflinks -g !{gencode_annotation_gtf} \
                      -b !{fasta_ref} \
                      --library-type !{strand_str} \
                      --max-multiread-fraction 0.25 \
                      --3-overhang-tolerance 2000 \
                      -o Cufout_!{file_tag_new} \
                      -p !{cufflinks_threads} !{alignment_bam}
                      
            mv Cufout_!{file_tag_new}/transcripts.gtf Cufout_!{file_tag_new}_transcripts.gtf
            '''

            }


        }

// Create a file 'gtf_filenames' containing the filenames of each post processes cufflinks gtf

        cuflinksoutgtf.collectFile { file -> ['gtf_filenames.txt', file.name + '\n'] }
                .set { GTFfilenames }

        /*
        * Step 6: Merged GTFs into one
        */
        process cuffmerge_assembled_gtf {
            
            tag { file_tag }
            label 'para'
            publishDir pattern: "CUFFMERGE/merged.gtf",
                    path: { params.out_folder + "/Result/All_assemblies" }, mode: 'copy', overwrite: true

            input:
            file gtf_filenames from GTFfilenames
            file cufflinksgtf_file from cuflinksoutgtf_fn.toList() // not used but just send the file in current running folder

            file fasta_ref


            output:
            file "CUFFMERGE/merged.gtf" into mergeTranscripts_forCompare, mergeTranscripts_forExtract, mergeTranscripts_forCodeingProtential
            shell:

            cufflinks_threads = ava_cpu- 1

            '''
            mkdir CUFFMERGE
            cuffmerge -o CUFFMERGE \
                      -s !{fasta_ref} \
                      -p !{cufflinks_threads} \
                         !{gtf_filenames}
            
            '''
        }
    }


}
else {
    println print_yellow("Raw reads quality check step was skipped due to provided ") + print_green("--merged_gtf") + print_yellow(" option\n")
    println print_yellow("Reads mapping step was skipped due to provided ") + print_green("--merged_gtf") + print_yellow(" option\n")

    merged_gtf = file(params.merged_gtf)
    Channel.fromPath(merged_gtf)
            .ifEmpty { exit 1, "Cannot find merged gtf : ${merged_gtf}" }
            .into {
        mergeTranscripts_forCompare; mergeTranscripts_forExtract; mergeTranscripts_forCodeingProtential
    }

    // add fastq when do quantification
    reads = params.input_folder + params.fastq_ext
    if (params.qctools == 'fastqc') {
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Fastq file not found, plz check your file path : ${reads}\n")
        }
        .into { reads_for_fastqc; readPairs_for_discovery;readPairs_for_kallisto}
        process Run_fastQC_2 {
            tag { fastq_tag }
            label 'qc'

            publishDir pattern: "*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "*.html" into fastqc_for_waiting
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            '''
            fastqc -t !{fastq_threads} !{fastq_file[0]} !{fastq_file[1]}
        '''
        }
    }
    else if (params.qctools == 'afterqc'){
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Fastq file not found :  ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }
        .set { reads_for_fastqc}
        process Run_afterQC_2 {

            tag { fastq_tag }
            label 'qc'

            publishDir pattern: "QC/*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "QC/*.html" into fastqc_for_waiting
            set val(fastq_tag), file('*.good.fq.gz')  into readPairs_for_discovery,readPairs_for_kallisto
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            if (params.singleEnd) {
                '''
            after.py -z -1 !{fastq_file[0]} -g ./
            '''
            } else {
                '''
            after.py -z -1 !{fastq_file[0]} -2 !{fastq_file[1]} -g ./
            '''
            }
        }
    }
    else if (params.qctools == 'fastp'){
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Fastq file not found :  ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }
        .set { reads_for_fastqc}
        process Run_FastP_2 {

            tag { fastq_tag }
            label 'qc'

            publishDir pattern: "*.html",
                    path: { params.out_folder + "/Result/QC" }, mode: 'copy', overwrite: true

            input:
            set val(samplename), file(fastq_file) from reads_for_fastqc

            output:
            file "*.html" into fastqc_for_waiting
            set val(fastq_tag), file('*qc.fq.gz')  into readPairs_for_discovery,readPairs_for_kallisto
            shell:
            fastq_tag = samplename
            fastq_threads = idv_cpu - 1
            if (params.singleEnd) {
                '''
            fastp -i !{fastq_file[0]} -o !{samplename}.qc.gz -h !{samplename}_fastp.html
           
            '''
            } else {
                '''
            fastp -i !{fastq_file[0]}  -I !{fastq_file[1]} -o !{samplename}_1.qc.fq.gz  -O !{samplename}_2.qc.fq.gz -h !{samplename}_fastp.html
            '''
            }
        }
    }
    else{
        Channel.fromFilePairs(reads, size: params.singleEnd ? 1 : 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: ${reads}\nPlz check your fasta_ref string in nextflow.config file \n")
        }
        .into{readPairs_for_discovery; readPairs_for_kallisto;fastqc_for_waiting}
    }
    fastqc_for_waiting2 = fastqc_for_waiting.first()

}

/*
*Step 7: Compare assembled gtf with known annotations (GENCODE)
*/
    process Merge_assembled_gtf_with_GENCODE {

        tag { file_tag }
        input:
        file mergeGtfFile from mergeTranscripts_forCompare
        file gencode_annotation_gtf from gencode_annotation_gtf_for_merge

        output:
        file "merged_lncRNA.merged.gtf.tmap" into comparedGTF_tmap
        shell:

        gffcompare_threads = ava_cpu- 1
        '''
        #!/bin/sh
        gffcompare -r !{gencode_annotation_gtf} -p !{gffcompare_threads} !{mergeGtfFile} -o merged_lncRNA
        '''
    }



/*
*Step 8: Filter GTFs to distinguish novel lncRNAs
*/
process Identify_novel_lncRNA_with_criterions {

    input:
    file comparedTmap from comparedGTF_tmap
    file fasta_ref
    file mergedGTF from mergeTranscripts_forExtract

    output:
    file "novel.gtf.tmap" into noveltmap
    file "novel.longRNA.fa" into novelLncRnaFasta
    file "novel.longRNA.exoncount.txt" into novelLncRnaExonCount

    shell:
    '''
        # filtering novel lncRNA based on cuffmerged trascripts
        set -o pipefail
        awk '$3 =="x"||$3=="u"||$3=="i"{print $0}' !{comparedTmap} > novel.gtf.tmap
        #   excluding length smaller than 200 nt
        awk '$10 >200{print}' novel.gtf.tmap > novel.longRNA.gtf.tmap
        #   extract gtf
        awk '{print $5}' novel.longRNA.gtf.tmap |perl !{baseDir}/bin/extract_gtf_by_name.pl !{mergedGTF} - >novel.longRNA.gtf
        awk '{if($3=="exon"){print $0}}' novel.longRNA.gtf > novel.longRNA.format.gtf 
        perl !{baseDir}/bin/get_exoncount.pl novel.longRNA.format.gtf  > novel.longRNA.exoncount.txt
        # gtf2gff3
        #check whether required
        # get fasta from gtf
        gffread novel.longRNA.gtf -g !{fasta_ref} -w novel.longRNA.fa -W
     '''
}

/*
*Step 9: Predict coding potential abilities using CPAT and PLEK (CNCI functionality coming soon!)
*/
novelLncRnaFasta.into { novelLncRnaFasta_for_PLEK; novelLncRnaFasta_for_CPAT; }

process Predict_coding_abilities_by_PLEK {
    
    // as PLEK can not return valid exit status even run smoothly, we manually set the exit status into 0 to promote analysis
    validExitStatus 0, 1, 2
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_PLEK
    output:
    file "novel.longRNA.PLEK.out" into novel_longRNA_PLEK_result
    shell:
    plek_threads = ava_cpu- 1
    '''
        PLEK.py -fasta !{novel_lncRNA_fasta} \
                                   -out novel.longRNA.PLEK.out \
                                   -thread !{plek_threads}
	    exit 0
        '''

}
process Predict_coding_abilities_by_CPAT {
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_CPAT
    output:
    file "novel.longRNA.CPAT.out" into novel_longRNA_CPAT_result
    shell:
    if(params.species=="human"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x !{baseDir}/bin/cpat_model/Human_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/Human_logitModel.RData \
                                       -o novel.longRNA.CPAT.out
        '''
    }else if (params.species=="mouse"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x !{baseDir}/bin/cpat_model/Mouse_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/Mouse_logitModel.RData \
                                       -o novel.longRNA.CPAT.out
        '''

    }else if (params.species=="zebrafish"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x !{baseDir}/bin/cpat_model/zebrafish_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/zebrafish_logitModel.RData \
                                       -o novel.longRNA.CPAT.out
        '''
    }else if (params.species=="fly"){
        '''
        cpat.py -g !{novel_lncRNA_fasta} \
                                       -x !{baseDir}/bin/cpat_model/fly_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/fly_logitModel.RData \
                                       -o novel.longRNA.CPAT.out
        '''
    }

}


/*
*Step 9: Merged and filter lncRNAs based on coding potential (CPAT/PLEK)
*/
process Filter_lncRNA_by_coding_potential_result {
    input:
    file novel_longRNA_PLEK_ from novel_longRNA_PLEK_result
    file novel_longRNA_CPAT_ from novel_longRNA_CPAT_result
    file longRNA_novel_exoncount from novelLncRnaExonCount
    file cuffmergegtf from mergeTranscripts_forCodeingProtential
    file gencode_annotation_gtf  from gencode_annotation_gtf_for_filter
    file fasta_ref

    output:
    file "novel.longRNA.stringent.gtf" into Novel_longRNA_stringent_gtf // not used
    file "novel.lncRNA.stringent.gtf" into novel_lncRNA_stringent_gtf
    file "novel.TUCP.stringent.gtf" into novel_TUCP_stringent_gtf // not used

    shell:
    '''
        set -o pipefail
        #merged transcripts
        perl !{baseDir}/bin/integrate_novel_transcripts.pl > novel.longRNA.txt
        awk '$4 >1{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - > novel.longRNA.stringent.gtf
        # retain lncRNA only by coding ability
        awk '$4 >1&&$5=="lncRNA"{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - > novel.lncRNA.stringent.gtf
        awk '$4 >1&&$5=="TUCP"{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - > novel.TUCP.stringent.gtf
        '''
}

/*
*Step 10: Further filtered lncRNAs with known criterion
*/
process Summary_renaming_and_classification {
    publishDir "${params.out_folder}/Result/Identified_lncRNA", mode: 'copy'


    input:
    file knowlncRNAgtf from KnownLncRNAgtf
    file gencode_protein_coding_gtf from proteinCodingGTF
    file novel_lncRNA_stringent_Gtf from novel_lncRNA_stringent_gtf
    file fasta_ref
    file mod_file_for_rename

    output:
//    file "lncRNA.final.v2.gtf" into finalLncRNA_gtf
//    file "lncRNA.final.v2.map" into finalLncRNA_map
    file "protein_coding.final.gtf" into final_protein_coding_gtf
    file "all_lncRNA_for_classifier.gtf" into finalLncRNA_for_class_gtf
    file "final_all.gtf" into finalGTF_for_quantification_gtf, finalGTF_for_annotate_gtf
    file "final_all.fa" into finalFasta_for_quantification_gtf
    file "protein_coding.fa" into final_coding_gene_for_CPAT_fa
    file "lncRNA.fa" into final_lncRNA_for_CPAT_fa
    file "lncRNA_classification.txt" into lncRNA_classification
    file "lncRNA.mapping.file" into rename_mapping_file
    //file "lncRNA.final.CPAT.out" into lncRNA_CPAT_statistic
    //file "protein_coding.final.CPAT.out" into protein_coding_CPAT_statistic

    shell:

    cufflinks_threads = ava_cpu- 1

    if(params.species=="human"){
        '''
        set -o pipefail
        gffcompare -G -o filter \
                    -r !{knowlncRNAgtf} \
                    -p !{cufflinks_threads} !{novel_lncRNA_stringent_Gtf}
        awk '$3 =="u"||$3=="x"{print $5}' filter.novel.lncRNA.stringent.gtf.tmap |sort|uniq| \
                    perl !{baseDir}/bin/extract_gtf_by_name.pl !{novel_lncRNA_stringent_Gtf} - > novel.lncRNA.stringent.filter.gtf
        
        #rename lncRNAs according to neighbouring protein coding genes
        awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} | perl -F'\\t' -lane '$F[8]=~/gene_id "(.*?)";/ && print join qq{\\t},@F[0,3,4],$1,@F[5,6,1,2,7,8,9]' - | \
            sort-bed - > gencode.protein_coding.gene.bed
        gtf2bed < novel.lncRNA.stringent.filter.gtf |sort-bed - > novel.lncRNA.stringent.filter.bed
        gtf2bed < !{knowlncRNAgtf} |sort-bed - > known.lncRNA.bed
        
        perl !{baseDir}/bin/rename_lncRNA_2.pl gencode_annotation_gtf_mod.gtf lncipedia_mod.gtf 
        # mv lncRNA.final.v2.gtf all_lncRNA_for_classifier.gtf
	grep -v "gene_id \"NA-" lncRNA.final.v2.gtf > all_lncRNA_for_classifier.gtf
        perl !{baseDir}/bin/rename_proteincoding.pl !{gencode_protein_coding_gtf} | grep -w exone > protein_coding.final.gtf
        cat all_lncRNA_for_classifier.gtf protein_coding.final.gtf > final_all.gtf
        gffread final_all.gtf -g !{fasta_ref} -w final_all.fa -W
        gffread all_lncRNA_for_classifier.gtf -g !{fasta_ref} -w lncRNA.fa -W
        gffread protein_coding.final.gtf -g !{fasta_ref} -w protein_coding.fa -W
        #classification 
        perl !{baseDir}/bin/lincRNA_classification.pl all_lncRNA_for_classifier.gtf !{gencode_protein_coding_gtf} lncRNA_classification.txt 
        
        
        '''
    }else{
        '''
        set -o pipefail
        gffcompare -G -o filter \
                    -r !{knowlncRNAgtf} \
                    -p !{cufflinks_threads} !{novel_lncRNA_stringent_Gtf}
        awk '$3 =="u"||$3=="x"{print $5}' filter.novel.lncRNA.stringent.gtf.tmap |sort|uniq| \
                    perl !{baseDir}/bin/extract_gtf_by_name.pl !{novel_lncRNA_stringent_Gtf} - > novel.lncRNA.stringent.filter.gtf
        
        #rename lncRNAs according to neighbouring protein coding genes
        awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} | perl -F'\\t' -lane '$F[8]=~/gene_id "(.*?)";/ && print join qq{\\t},@F[0,3,4],$1,@F[5,6,1,2,7,8,9]' - | \
            sort-bed - > gencode.protein_coding.gene.bed
        gtf2bed < novel.lncRNA.stringent.filter.gtf |sort-bed - > novel.lncRNA.stringent.filter.bed
        gtf2bed < !{knowlncRNAgtf} |sort-bed - > known.lncRNA.bed
        perl !{baseDir}/bin/rename_lncRNA_2.pl non_human_mod.gtf
        # mv lncRNA.final.v2.gtf all_lncRNA_for_classifier.gtf
        grep -v "gene_id \"NA-" lncRNA.final.v2.gtf > all_lncRNA_for_classifier.gtf
        perl !{baseDir}/bin/rename_proteincoding.pl !{gencode_protein_coding_gtf}> protein_coding.final.gtf
        cat all_lncRNA_for_classifier.gtf protein_coding.final.gtf > final_all.gtf
        gffread final_all.gtf -g !{fasta_ref} -w final_all.fa -W
        gffread all_lncRNA_for_classifier.gtf -g !{fasta_ref} -w lncRNA.fa -W
        gffread protein_coding.final.gtf -g !{fasta_ref} -w protein_coding.fa -W
        #classification 
        perl !{baseDir}/bin/lincRNA_classification.pl all_lncRNA_for_classifier.gtf !{gencode_protein_coding_gtf} lncRNA_classification.txt 
        
        
        '''
    }

}

/*
*Step 11: Rerun CPAT to evaluate the results
*/
//evaluate lncRNA
process Rerun_CPAT_to_evaluate_lncRNA {
    input:
    file lncRNA_final_cpat_fasta from final_lncRNA_for_CPAT_fa
    output:
    file "lncRNA.final.CPAT.out" into final_lncRNA_CPAT_result
    shell:

    if(params.species=="human"){
        '''
        cpat.py -g !{lncRNA_final_cpat_fasta} \
                                       -x !{baseDir}/bin/cpat_model/Human_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/Human_logitModel.RData \
                                       -o lncRNA.final.CPAT.out
        '''
    }else if (params.species=="mouse"){
        '''
        cpat.py -g !{lncRNA_final_cpat_fasta} \
                                       -x !{baseDir}/bin/cpat_model/Mouse_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/Mouse_logitModel.RData \
                                       -o lncRNA.final.CPAT.out
        '''

    }else if (params.species=="zebrafish"){
        '''
        cpat.py -g !{lncRNA_final_cpat_fasta} \
                                       -x !{baseDir}/bin/cpat_model/zebrafish_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/zebrafish_logitModel.RData \
                                       -o lncRNA.final.CPAT.out
        '''
    }else {
        '''
        cpat.py -g !{lncRNA_final_cpat_fasta} \
                                       -x !{baseDir}/bin/cpat_model/fly_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/fly_logitModel.RData \
                                       -o lncRNA.final.CPAT.out
        '''
   }
}
//evaluate coding
process Rerun_CPAT_to_evaluate_coding {
    input:
    file final_coding_gene_for_CPAT from final_coding_gene_for_CPAT_fa
    output:
    file "protein_coding.final.CPAT.out" into final_coding_gene_CPAT_result
    shell:
    '''
        cpat.py -g !{final_coding_gene_for_CPAT} \
                                       -x !{baseDir}/bin/cpat_model/Human_Hexamer.tsv \
                                       -d !{baseDir}/bin/cpat_model/Human_logitModel.RData \
                                       -o protein_coding.final.CPAT.out
        '''
}
//summary result
process Secondary_basic_statistic {

    input:
    file protein_coding_final_gtf from final_protein_coding_gtf
    file all_lncRNA_for_classifier_gtf from finalLncRNA_for_class_gtf
    file lncRNA_cds from final_lncRNA_CPAT_result
    file coding_gene_cds from final_coding_gene_CPAT_result
    file lncRNA_class from lncRNA_classification
    output:
    file "basic_charac.txt" into statistic_result

    shell:
    '''
        #!/usr/bin/perl -w
        #since CPAT arbitrarily transforms gene names into upper case, we apply 'uc' function to keep the genenames' consistency.  
        use strict;
        open OUT,">basic_charac.txt" or die;
        
        open FH,"all_lncRNA_for_classifier.gtf" or die;
        
        my %class;
        my %g2t;
        my %trans_len;
        my %exon_num;
        while(<FH>){
        chomp;
        my @field=split "\t";
        $_=~/gene_id "(.+?)"/;
        my $gid=$1;
        $_=~/transcript_id "(.+?)"/;
        my $tid=uc($1);
        $class{$tid}=$field[1];
        $g2t{$tid}=$gid;
        my $len=$field[4]-$field[3];
        $trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
        $exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
        }
        open FH,"protein_coding.final.gtf" or die;
        
        while(<FH>){
        chomp;
        my @field=split "\t";
        $_=~/gene_id "(.+?)"/;
        my $gid=uc($1);
        $_=~/transcript_id "(.+?)"/;
        my $tid=$1;
        $class{$tid}="protein_coding";
        $g2t{$tid}=$gid;
        my $len=$field[4]-$field[3];
        $trans_len{$tid}=(exists $trans_len{$tid})?$trans_len{$tid}+$len:$len;
        $exon_num{$tid}=(exists $exon_num{$tid})?$exon_num{$tid}+1:1;
        }
        
        my %lin_class;
        open IN,"lncRNA_classification.txt" or die;                 #change the file name
        while(<IN>){
        chomp;
        my @data = split /\\t/,$_;
        $lin_class{$data[0]} = $data[1];
        }
        open FH,"lncRNA.final.CPAT.out" or die;
        
        <FH>;
        
        while(<FH>){
            chomp;
            my @field=split "\t";
            my $tid=uc($field[0]);
            my $class;
            if (defined($lin_class{$tid})){
                $class = $lin_class{$tid};
            }else{
                $class = 'NA';
            }
            print OUT $g2t{$tid}."\t".$tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\t".$class."\n";
        }
            
        open FH,"protein_coding.final.CPAT.out" or die;
        
        <FH>;
                    
        while(<FH>){
            chomp;
            my @field=split "\t";
            my $tid=uc($field[0]);
            my $class;
            if (defined($lin_class{$tid})){
                $class = $lin_class{$tid};
            }else{
                $class = 'protein_coding';
            }
            print OUT $g2t{$tid}."\t".$tid."\t".$class{$tid}."\t".$field[5]."\t".$trans_len{$tid}."\t".$exon_num{$tid}."\t".$class."\n";
         }

    '''
}



//Keep the channel as constant variable to be used several times in quantification analysis

//The following code is designed for use if the merged_gtf have already been generated previously.
if(!params.merged_gtf){
    /*
*Step 11: Quantification step (Kallisto/Htseq)
*/
    if(params.quant=="htseq"){
        process Run_htseq_for_quantification{
            tag { file_tag }
            input:
            set val(samplename),file(bamfile) from forHtseqMappedReads
            file final_gtf from finalGTF_for_quantification_gtf

            output:
            file "${file_tag_new}.htseq.count " into htseq_tcv_collection

            shell:

            file_tag = samplename
            file_tag_new = file_tag
            if(params.unstrand){
                '''
                sambamba view !{bamfile} > !{samplename}.sam # resolved error caused by bam and htseq version conflicts 
                htseq-count -t exon -i gene_id -s no -r pos -f sam !{samplename}.sam !{final_gtf} > !{samplename}.htseq.count 
                rm !{samplename}.sam
                '''
            }else {
                '''
                sambamba view !{bamfile} > !{samplename}.sam # resolved error caused by bam and htseq version conflicts 
                htseq-count -t exon -i gene_id -r pos -f sam !{samplename}.sam !{final_gtf} > !{samplename}.htseq.count 
                rm !{samplename}.sam
                '''
            }



        }
    }else{
        process Build_kallisto_index_of_GTF_for_quantification {
            
            input:
            file transript_fasta from finalFasta_for_quantification_gtf

            output:
            file "transcripts.idx" into final_kallisto_index

            shell:
            '''
            #index kallisto reference 
            kallisto index -i transcripts.idx !{transript_fasta}
            
            '''
        }
        constant_kallisto_index = final_kallisto_index.first()
        process Run_kallisto_for_quantification {
            

            tag { file_tag }
            label 'para'

            input:
            file kallistoIndex from constant_kallisto_index
            set val(samplename), file(pair) from readPairs_for_kallisto

            output:
            file "${file_tag_new}_abundance.tsv" into kallisto_tcv_collection

            shell:
            file_tag = samplename
            file_tag_new = file_tag
            kallisto_threads = ava_cpu- 1
            if (params.singleEnd) {
                println print_purple("Quantification by kallisto in single end mode")
                '''
                #quantification by kallisto in single end mode
                kallisto quant -i !{kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 --single -l 180 -s 20  !{pair} 
                mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
                '''


            } else {
                println print_purple("quantification by kallisto in paired end mode")
                '''
                #quantification by kallisto 
                kallisto quant -i !{kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 !{pair[0]} !{pair[1]}
                mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
                '''
            }
        }
    }

}else{
    /*
*Step 11: Quantification step (Kallisto/Htseq)
*/
    if(params.quant=="htseq"){
        exit 0, print_red("htseq can not be applicable without mapping step, plz set quant tool using `kallisto`")
    }else {
        process Build_kallisto_index_of_GTF_for_quantification_2 {
            

            input:
            file transript_fasta from finalFasta_for_quantification_gtf

            output:
            file "transcripts.idx" into final_kallisto_index

            shell:
            '''
    #index kallisto reference 
    kallisto index -i transcripts.idx !{transript_fasta}
    
    '''
        }
        constant_kallisto_index = final_kallisto_index.first()
        process Run_kallisto_for_quantification_2 {
            

            tag { file_tag }
            label 'para'

            input:
            file kallistoIndex from constant_kallisto_index
            set val(samplename), file(pair) from readPairs_for_kallisto
            file tempfiles from fastqc_for_waiting2
            output:
            file "${file_tag_new}_abundance.tsv" into kallisto_tcv_collection

            shell:
            file_tag = samplename
            file_tag_new = file_tag
            kallisto_threads = ava_cpu - 1
            if (params.singleEnd) {
                println print_purple("Quantification by kallisto in single end mode")
                '''
                #quantification by kallisto in single end mode
                kallisto quant -i !{kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 --single -l 180 -s 20 !{pair} 
                mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
                
                '''


            } else {
                println print_purple("Quantification by kallisto in paired end mode")
                '''
                #quantification by kallisto 
                kallisto quant -i !{kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 !{pair[0]} !{pair[1]}
                mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
                '''
                    }
        }
    }
}


/*
*Step 12: Generate count matrix for differential expression analysis
*/

if(params.quant=="htseq"){
    process Get_HTseq_matrix {
        tag { file_tag }
        publishDir pattern: "htseq*.txt",
                path: "${params.out_folder}/Result/Quantification/", mode: 'copy'
        input:
        file abundance_tsv_matrix from htseq_tcv_collection.collect()
        file annotated_gtf from finalGTF_for_annotate_gtf
        output:
        file "htseq.count.txt" into expression_matrixfile_count

        shell:
        file_tag = "htseq"
        '''
        perl !{baseDir}/bin/get_map_table.pl  final_all.gtf  > map.file
        R CMD BATCH !{baseDir}/bin/get_htseq_matrix.R
        '''
    }
}else{
    process Get_kallisto_matrix {
        tag { file_tag }
        publishDir pattern: "kallisto*.txt",
                path: "${params.out_folder}/Result/Quantification/", mode: 'copy'
        input:
        file abundance_tsv_matrix from kallisto_tcv_collection.collect()
        file annotated_gtf from finalGTF_for_annotate_gtf
        output:
        file "kallisto.count.txt" into expression_matrixfile_count
        file "kallisto.tpm.txt" into expression_matrixfile_tpm

        shell:
        file_tag = "Kallisto"
        '''
        perl !{baseDir}/bin/get_map_table.pl  --gtf_file=final_all.gtf  > map.file
        R CMD BATCH !{baseDir}/bin/get_kallisto_matrix.R
        '''
    }
}

/*
Step 13: Perform Differential Expression analysis and generate report
 */

// Initialize parameter for lncPipeReporter
lncRep_Output = params.lncRep_Output
lncRep_theme = params.lncRep_theme
lncRep_cdf_percent = params.lncRep_cdf_percent
lncRep_max_lnc_len = params.lncRep_max_lnc_len
lncRep_min_expressed_sample = params.lncRep_min_expressed_sample
detools = params.detools
design=params.design
if(design!=null){
    design = file(params.design)
    if (!design.exists()) exit 1, "Design file not found, plz check your design file path: ${params.design}"

    if(!params.merged_gtf) {
        process Run_LncPipeReporter {
            tag { file_tag }
            publishDir pattern: "*",
                    path: "${params.out_folder}/Result/", mode: 'move'
            input:
            //alignmet log
            file design
            file alignmetlogs from alignment_logs.collect()
            //gtf statistics
            file basic_charac from statistic_result
            //Expression matrix
            file kallisto_count_matrix from expression_matrixfile_count

            output:
            file "*" into final_output
            shell:
            file_tag = "Generating report ..."
            """
            perl ${baseDir}/bin/modifyDesign.pl ${design}  > design.matrix 
            Rscript -e "library(LncPipeReporter);run_reporter(input='.', output = 'reporter.html',output_dir='./LncPipeReports',de.method=\'${detools}\',theme = 'npg',cdf.percent = ${lncRep_cdf_percent},max.lncrna.len = ${lncRep_max_lnc_len},min.expressed.sample = ${lncRep_min_expressed_sample}, ask = FALSE)"
            """
        }
    }else{
        process Run_LncPipeReporter_2 {
            tag { file_tag }
            publishDir pattern: "*",
                    path: "${params.out_folder}/Result/", mode: 'move'
            input:
            //alignment log
            file design
            //gtf statistics
            file basic_charac from statistic_result
            //Expression matrix
            file kallisto_count_matrix from expression_matrixfile_count

            output:
            file "*" into final_output
            shell:
            file_tag = "Generating report ..."
            """
            perl ${baseDir}/bin/modifyDesign.pl ${design}  > design.matrix 
            Rscript -e "library(LncPipeReporter);run_reporter(input='.', output = 'reporter.html',output_dir='./LncPipeReports',de.method=\'${detools}\',theme = 'npg',cdf.percent = ${lncRep_cdf_percent},max.lncrna.len = ${lncRep_max_lnc_len},min.expressed.sample = ${lncRep_min_expressed_sample}, ask = FALSE)"
            """
        }
    }

}else{
    if(!params.merged_gtf) {
        process Run_LncPipeReporter_without_Design {
            tag { file_tag }
            publishDir pattern: "*",
                    path: "${params.out_folder}/Result/", mode: 'move'
            input:
            //alignmet log
            file alignmetlogs from alignment_logs.collect()
            //gtf statistics
            file basic_charac from statistic_result
            //Expression matrix
            file kallisto_count_matrix from expression_matrixfile_count

            output:
            file "*" into final_output
            shell:
            file_tag = "Generating report ..."
            """
             Rscript -e "library(LncPipeReporter);run_reporter(input='.', output = 'reporter.html',output_dir='./LncPipeReports',de.method=\'${detools}\',theme = 'npg',cdf.percent = ${lncRep_cdf_percent},max.lncrna.len = ${lncRep_max_lnc_len},min.expressed.sample = ${lncRep_min_expressed_sample}, ask = FALSE)"
            """
        }
    }else{
        process Run_LncPipeReporter_without_Design_2 {
            tag { file_tag }
            publishDir pattern: "*",
                    path: "${params.out_folder}/Result/", mode: 'move'
            input:
            //alignment log
            //gtf statistics
            file basic_charac from statistic_result
            //Expression matrix
            file kallisto_count_matrix from expression_matrixfile_count

            output:
            file "*" into final_output
            shell:
            file_tag = "Generating report ..."

            """
             Rscript -e "library(LncPipeReporter);run_reporter(input='.', output = 'reporter.html',output_dir='./LncPipeReports',de.method=\'${detools}\',theme = 'npg',cdf.percent = ${lncRep_cdf_percent},max.lncrna.len = ${lncRep_max_lnc_len},min.expressed.sample = ${lncRep_min_expressed_sample}, ask = FALSE)"
            """


        }
    }
}



//pipeline log
if(workflow.success) {
    workflow.onComplete {

        log.info print_green("LncPipe Pipeline Complete!")

        //email information
        if (params.mail) {
            recipient = params.mail
            def subject = 'My LncPipe execution'

            ['mail', '-s', subject, recipient].execute() <<
                    """

            LncPipe execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """
        }


    }
}
workflow.onError {
    println print_yellow("Oops... Pipeline execution stopped with the following message: ")+print_red(workflow.errorMessage)
}

