#!/usr/bin/env nextflow

/*
 * Copyright (c) 2013-2017, Sun Yat-sen University Cancer Center.
 *
 *
 *   LncPipe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *      See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNA-Toy.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * LncPipe: A nextflow-based lncRNA identification and analysis pipeline from RNA sequencing data
 *
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>
 * Yu Sun <sun_yu@mail.nankai.edu.cn>
 * Zhixiang Zuo <zuozhx@sysucc.org.cn>
 */


// requirement:
// - utadapt
// - STAR
// - RSEM
// - Cufflinks
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


def print_red = { String str-> ANSI_RED+str+ANSI_RESET }
def print_black = { String str-> ANSI_BLACK+str+ANSI_RESET }
def print_green = { String str-> ANSI_GREEN+str+ANSI_RESET }
def print_yellow = { String str-> ANSI_YELLOW+str+ANSI_RESET }
def print_blue = { String str-> ANSI_BLUE+str+ANSI_RESET }
def print_cyan= { String str-> ANSI_CYAN+str+ANSI_RESET }
def print_purple = { String str-> ANSI_PURPLE+str+ANSI_RESET }
def print_white = { String str-> ANSI_WHITE+str+ANSI_RESET }

version = '0.0.4'
//Help information
    // Pipeline version



//=======================================================================================
//help information
if (params.help){
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "LncPipe: a Nextflow-based Long non-coding RNA analysis PIPELINE v$version"
    log.info "LncPipe integrates several NGS processing tools to identify novel long non-coding RNAs from"
    log.info "unprocessed RNA sequencing data. Before run this pipeline, users need to install several tools"
    log.info "unprocessed RNA sequencing data. Before run this pipeline, users need to install several tools"
    log.info print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('    The typical command for running the pipeline is as follows:\n') +
            print_purple('       Nextflow lncRNApipe.nf \n') +

            print_yellow('    Mandatory arguments:             Input and output setting\n') +
            print_cyan('      --input_folder                ')+print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --out_folder                  ')+print_green('The output directory where the results will be saved(optional), current path is default\n') +
            print_cyan('      --aligner                     ')+print_green('Aligner for reads mapping (optional), STAR is default\n') +
            '\n'+
            print_yellow('    Options:                         General options for run this pipeline\n') +
            print_cyan('      --singleEnd                   ')+print_green('Specifies that the input is single end reads(optional), paired end mode default \n') +
            print_cyan('      --merged_gtf                  ')+print_green('Start analysis with assemblies already produced and skip fastqc/alignment step, DEFAOUL NULL\n') +
            print_cyan('      --mode                        ')+print_green('Start analysis with fastq or bam mode, can not set with --merged_gtf, values should be \'fastq\' or \'bam\'\n') +
            '\n'+
            print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.\n') +
            print_cyan('      --star_index                  ')+print_green('Path to STAR index(required)\n') +
            print_cyan('      --fasta                       ')+print_green('Path to Fasta reference(required)\n') +
            print_cyan('      --gencode_annotation_gtf      ')+print_green('An annotation file from GENCODE database for annotating lncRNAs(required)\n')+
            print_cyan('      --lncipedia_gtf               ')+print_green('An annotation file from LNCipedia database for annotating lncRNAs(required)\n')+
            print_cyan('      --rRNAmask                    ')+print_green('rRNA GTF for removing rRNA transcript from gtf files(required)\n')+
            '\n'+
            print_yellow('    Skip options:                    Skip the certain step when necessary \n') +
            print_cyan('      --skip_mapping                ')+print_green('Skip mapping step by directly providing bam files\n') +
            print_cyan('      --skip_QC                     ')+print_green('Skip QC step when the reads are clean reads\n') +
            print_cyan('      --skip_DE                     ')+print_green('Skip differential expression analysis step  \n')+
            '\n'+
            print_yellow('    Other options:                   Specify the email and \n') +
            print_cyan('      --email                       ')+print_green('Set  e-mail address to get a summary when the workflow finished\n') +
            print_cyan('      --name                         ')+print_green('Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.\n')
            print_cyan('      --cpu                          ')+print_green('Number of cpu used in analysis. DEFAULT, 16.\n')
            print_cyan('      --mem                          ')+print_green('Memory setting for each analysis run. DEFAULT, 32G.\n')


    log.info '------------------------------------------------------------------------'
    log.info print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    log.info print_yellow('Copyright (c) 2013-2017, Sun Yat-sen University Cancer Center.')
    log.info '------------------------------------------------------------------------'
    exit 0
}






//default values
params.input_folder = './'
params.out_folder = './'

// Reference
// Set reference information here if you don't want to pass them from parameter any longer, we recommand users using the latest reference and the annotation file with the sample genome version.
params.fasta_ref = '/data/database/human/hg38/genome.fa'
params.star_idex = '/data/database/human/hg38/RSEM_STAR_Index'
params.gencode_annotation_gtf = "/data/database/human/hg38/annotation/gencode.v24.annotation.gtf"
params.lncipedia_gtf = "/data/database/human/hg38/annotation/lncipedia_4_0_hg38.gtf"
params.rRNAmask = "/data/database/human/hg38/annotation/hg38_rRNA.gtf";
// software path
params.plekpath = '/home/zhaoqi/software/PLEK.1.2/'
//params.cncipath = '/data/software/CNCI-master'
params.cpatpath = '/home/zhaoqi/software/CPAT/CPAT-1.2.2/'
// other options
params.cpu = null
params.mem = null
params.singleEnd = false


//Checking parameters
log.info print_purple("You are running LncPipe with the following parameters:")
log.info print_purple("Checking parameters ...")
log.info print_yellow("=====================================")
log.info print_yellow("Input folder:                   ")+ print_green(params.input_folder)
log.info print_yellow("Output folder:                  ")+ print_green(params.out_folder)
log.info print_yellow("Genome sequence location:       ")+ print_green(params.fasta_ref)
log.info print_yellow("Star index path:                ")+ print_green(params.star_idex)
log.info print_yellow("GENCODE annotation location:    ")+ print_green(params.gencode_annotation_gtf)
log.info print_yellow("lncipedia ann0tation location:  ")+ print_green(params.lncipedia_gtf)
log.info print_yellow("rRNA annotation location:       ")+ print_green(params.rRNAmask)
log.info print_yellow("=====================================")
log.info "\n"


// fastq file
params.fastq_ext = "fastq.gz"
params.suffix1='_{1,2}'
params.merged_gtf=null
params.mode="fastq"
//aligner
params.aligner="star"


// run information of system file
//automatic set optimize resource for analysis based on current system resources
ava_mem = (double)(Runtime.getRuntime().freeMemory())
ava_cpu = Runtime.getRuntime().availableProcessors()
if(params.cpu!=null && ava_cpu > params.cpu ){
    ava_cpu = params.cpu
}else if(params.cpu!=null && ava_cpu < params.cpu ) {
    print print_red("cpu number set in command is not used for exceeding the max available processors, \n use default parameter to run pipe. ")
}
if(params.mem!=null && ava_mem > params.mem ){
    ava_mem = params.mem
}else if(params.mem!=null && ava_mem < params.mem ) {
    print print_red("memory set in command is not used for exceeding the max available processors, \n use default parameter to run pipe. ")
}
// set individual cpu for fork run
idv_cpu=40
int fork_number=ava_cpu/idv_cpu
if(fork_number<1){
    fork_number=1
}



// read file
fasta_ref = file(params.fasta_ref)
star_idex = file(params.star_idex)
input_folder = file(params.input_folder)
gencode_annotation_gtf = file(params.gencode_annotation_gtf)
lncipedia_gtf = file(params.lncipedia_gtf)
//cncipath = file(params.cncipath)
rRNAmaskfile = file(params.rRNAmask)



//Prepare annotations
annotation_channel = Channel.from(gencode_annotation_gtf, lncipedia_gtf)
annotation_channel.collectFile { file -> ['lncRNA.gtflist', file.name + '\n'] }
        .set { LncRNA_gtflist }


process combine_public_annotation {
    cpus ava_cpu
    input:
    file lncRNA_gtflistfile from LncRNA_gtflist
    file gencode_annotation_gtf
    file lncipedia_gtf

    output:
    file "gencode_protein_coding.gtf" into proteinCodingGTF,proteinCodingGTF_forClass
    file "known.lncRNA.gtf" into KnownLncRNAgtf

    shell:
    cufflinks_threads = ava_cpu.intdiv(2) - 1
    '''
        set -o pipefail
        cuffmerge -o merged_lncRNA \
                    !{lncRNA_gtflistfile}
        cat !{gencode_annotation_gtf} |grep "protein_coding" > gencode_protein_coding.gtf
        cuffcompare -o merged_lncRNA \
                    -r !{gencode_annotation_gtf} \
                    -p !{cufflinks_threads} merged_lncRNA/merged.gtf
        awk '$3 =="u"||$3=="x"{print $5}' merged_lncRNA/merged_lncRNA.merged.gtf.tmap \
        |sort|uniq|perl !{baseDir}/bin/extract_gtf_by_name.pl merged_lncRNA/merged.gtf - > merged.filter.gtf
        mv  merged.filter.gtf known.lncRNA.gtf
        
        '''

}


// whether the merged gtf have already produced.
if (params.merged_gtf==null) {
//Star index
//star_ref = file(params.params.star_idex_ref)


    if (params.mode == 'fastq') {
        println print_purple("Analysis from fastq file")
        println print_purple("Start mapping with STAR aligner")
//

//Match the pairs on two channels
        reads=params.input_folder+'*'+params.suffix1+'.'+params.fastq_ext
        Channel.fromFilePairs(reads,size: params.singleEnd ? 1 : 2)
                .ifEmpty { exit 1, print_red("Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\n" )}
                .into{readPairs_for_discovery;readPairs_for_quatification}

//star_index if not exist
        if(params.aligner == 'star' && params.star_idex==null && fasta_ref){
            process makeSTARindex {
                tag fasta_ref
                cpus ava_cpu

                input:
                file fasta_ref from fasta_ref
                file gtf from gencode_annotation_gtf

                output:
                file "star_index" into star_idex

                script:
                star_threads = ava_cpu.intdiv(2) - 1
                """
                mkdir star_index
                STAR \
                    --runMode genomeGenerate \
                    --runThreadN ${star_threads} \
                    --sjdbGTFfile $gtf \
                    --sjdbOverhang 149 \
                    --genomeDir star_index/ \
                    --genomeFastaFiles $fasta_ref
                """
            }
        }

//Star alignment

            process fastq_star_alignment_For_discovery {
                cpus idv_cpu
                tag { file_tag }
                maxForks fork_number

                input:
                file pair from readPairs_for_discovery
                file fasta_ref
                file star_idex

                output:
                set val(file_tag_new), file("STAR_${file_tag_new}") into STARmappedReads,alignment_logs
                shell:
                file_tag = pair[0].name.replace("${params.suffix1}"+"."+"${params.fastq_ext}", "")
                file_tag_new = file_tag
                star_threads = idv_cpu.intdiv(2) - 1

                if (params.singleEnd) {
                    """
                         STAR --runThreadN !{star_threads} \
                            --twopassMode Basic \
                            --genomeDir !{star_idex} \
                            --readFilesIn !{pair[0]} \
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
                            --outFileNamePrefix !{file_tag_new} > !{file_tag_new}_star_log.txt
                                 
                            mkdir STAR_!{file_tag_new}
                            mv !{file_tag_new}Aligned* STAR_!{file_tag_new}/.
                            mv !{file_tag_new}SJ* STAR_!{file_tag_new}/.
                            mv !{file_tag_new}Log* STAR_!{file_tag_new}/.
                    """
                }else {

                    '''
                            STAR --runThreadN !{star_threads}  \
                                 --twopassMode Basic --genomeDir !{star_idex} \
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
                                 --outFileNamePrefix !{file_tag_new} > !{file_tag_new}_star_log.txt
                                 
                            mkdir STAR_!{file_tag_new}
                            mv !{file_tag_new}Aligned* STAR_!{file_tag_new}/.
                            mv !{file_tag_new}SJ* STAR_!{file_tag_new}/.
                            mv !{file_tag_new}Log* STAR_!{file_tag_new}/.
                    '''
                }
            }



// Cufflinks assembly
        process cufflinks_assembly {
            cpus idv_cpu
            tag { file_tag }
            maxForks fork_number
            input:
            set val(file_tag), file(Star_alignment) from STARmappedReads
            file fasta_ref
            file gencode_annotation_gtf
            file rRNAmaskfile

            output:

            file "Cufout_${file_tag_new}_transcripts.gtf" into cuflinksoutgtf, cuflinksoutgtf_fn

            shell:
            file_tag_new = file_tag
            cufflinks_threads = idv_cpu.intdiv(2) - 1

            '''
            #run cufflinks
            
            cufflinks -g !{gencode_annotation_gtf} \
                      -b !{fasta_ref} \
                      --library-type fr-firststrand \
                      --max-multiread-fraction 0.25 \
                      --3-overhang-tolerance 2000 \
                      -o Cufout_!{file_tag_new} \
                      -p !{cufflinks_threads} !{Star_alignment}/!{file_tag_new}Aligned.sortedByCoord.out.bam
                      
            mv Cufout_!{file_tag_new}/transcripts.gtf Cufout_!{file_tag_new}_transcripts.gtf
            '''

        }

// Create a file 'gtf_filenames' containing the filenames of each post processes cufflinks gtf

        cuflinksoutgtf.collectFile { file -> ['gtf_filenames.txt', file.name + '\n'] }
                .set { GTFfilenames }

// run cuffmerge

        process cuffmerge_assembled_gtf {
            cpus ava_cpu
            tag { file_tag }

            input:
            file gtf_filenames from GTFfilenames
            file cufflinksgtf_file from cuflinksoutgtf_fn.toList() // not used but just send the file in current running folder

            file fasta_ref


            output:
            file "CUFFMERGE/merged.gtf" into cuffmergeTranscripts_forCompare, cuffmergeTranscripts_forExtract, cuffmergeTranscripts_forCodeingProtential
            shell:

            cufflinks_threads = ava_cpu.intdiv(2) - 1

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

if (params.merged_gtf!=null) {

    Channel.fromPath(params.merged_gtf)
            .ifEmpty { exit 1, "Cannot find merged gtf : ${params.merged_gtf}" }
            .into { cuffmergeTranscripts_forCompare;cuffmergeTranscripts_forExtract;cuffmergeTranscripts_forCodeingProtential}
}

// run cuffcompare  merged gtf with gencode annotation
process cuffcompare_GENCODE {
    cpus ava_cpu
    tag { file_tag }
    input:
    file cuffmergefile from cuffmergeTranscripts_forCompare
    file gencode_annotation_gtf

    output:
    file "CUFFCOMPARE_GENCODE" into cuffcomparegencodeDir
    file ""
    shell:

    cufflinks_threads = ava_cpu.intdiv(2) - 1
    '''
        mkdir CUFFCOMPARE_GENCODE
        #!/bin/sh
        cuffcompare -o CUFFCOMPARE_GENCODE/merged_lncRNA \
                    -r !{gencode_annotation_gtf} \
                    -p !{cufflinks_threads} !{cuffmergefile}
        mv merged_lncRNA.merged.gtf.tmap CUFFCOMPARE_GENCODE/
        '''
}

// filter transcript


process ExtractGTF {

    input:
    file cuffcompareDir from cuffcomparegencodeDir
    file fasta_ref
    file mergedGTF from cuffmergeTranscripts_forExtract

    output:
    file "novel.gtf.tmap" into noveltmap
    file "novel.longRNA.fa" into novelLncRnaFasta
    file "novel.longRNA.exoncount.txt" into novelLncRnaExonCount

    shell:
    '''
        # filtering novel lncRNA based on cuffmerged trascripts
        set -o pipefail
        awk '$3 =="x"||$3=="u"||$3=="i"{print $0}' !{cuffcompareDir}/merged_lncRNA.merged.gtf.tmap > novel.gtf.tmap
        #   excluding length smaller than 200 nt
        awk '$11 >200{print}' novel.gtf.tmap > novel.longRNA.gtf.tmap
        #   extract gtf
        awk '{print $5}' novel.longRNA.gtf.tmap |perl !{baseDir}/bin/extract_gtf_by_name.pl !{mergedGTF} - >novel.longRNA.gtf
        perl !{baseDir}/bin/get_exoncount.pl novel.longRNA.gtf > novel.longRNA.exoncount.txt
        # gtf2gff3
        #check whether required
        # get fasta from gtf
        gffread novel.longRNA.gtf -g !{fasta_ref} -w novel.longRNA.fa -W
     '''
}

// predicting coding potential _ parallel
// Divide fasta channel into three for parallel running

novelLncRnaFasta.into { novelLncRnaFasta_for_PLEK; novelLncRnaFasta_for_CPAT; novelLncRnaFasta_for_CNCI }

process run_PLEK {
    cpus ava_cpu
    validExitStatus 0,1,2
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_PLEK
    output:
    file "novel.longRNA.PLEK.out" into novel_longRNA_PLEK_result
    shell:
    plek_threads = ava_cpu.intdiv(2) - 1
    '''
        python !{params.plekpath}/PLEK.py -fasta !{novel_lncRNA_fasta} \
                                   -out novel.longRNA.PLEK.out \
                                   -thread !{plek_threads}
	    exit 0
        '''

}
process run_CPAT {
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_CPAT
    output:
    file "novel.longRNA.CPAT.out" into novel_longRNA_CPAT_result
    shell:
    '''
        python !{params.cpatpath}/bin/cpat.py -g !{novel_lncRNA_fasta} \
                                       -x !{params.cpatpath}/dat/Human_Hexamer.tsv \
                                       -d !{params.cpatpath}/dat/Human_logitModel.RData \
                                       -o novel.longRNA.CPAT.out
        '''
}
//    process run_CNCI{
//        cpus ava_cpu
//        input:
//        file novel_lncRNA_fasta from novelLncRnaFasta_for_CNCI
//        file cncipath
//        output:
//        file lncRNA/CNCI* into novel_longRNA_CNCI_result
//        shell:
//        cnci_threads  = ava_cpu.intdiv(2) - 1
//        '''
//        python !{cncipath}/CNCI.py -f !{novel_lncRNA_fasta} -p !{cnci_threads} -o lncRNA/CNCI -m ve
//        '''
//    }

// merge transcripts and retain lncRNA only by coding ability
process merge_filter_by_coding_potential {
    input:
    file novel_longRNA_PLEK_ from novel_longRNA_PLEK_result
    file novel_longRNA_CPAT_ from novel_longRNA_CPAT_result
    file longRNA_novel_exoncount from novelLncRnaExonCount
    file cuffmergegtf from cuffmergeTranscripts_forCodeingProtential
    file gencode_annotation_gtf
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

//Further  identify novel lncRNA based on annotated database
process Filter_lncRNA_based_annotationbaes {
   // publishDir "${baseDir}/Result/Identified_lncRNA", mode: 'copy'
    cpus ava_cpu

    input:
    file knowlncRNAgtf from KnownLncRNAgtf
    file gencode_protein_coding_gtf from proteinCodingGTF
    file novel_lncRNA_stringent_Gtf from novel_lncRNA_stringent_gtf

    output:
    file "lncRNA.final.v2.gtf" into finalLncRNA_gtf
    file "lncRNA.final.v2.map" into finalLncRNA_map
    file "protein_coding.final.gtf" into final_protein_coding_gtf
    file "all_lncRNA_for_classifier.gtf" into finalLncRNA_for_class_gtf
    file "final_all.gtf" into finalGTF_for_quantification_gtf
    file "final_all.fa" into finalFasta_for_quantification_gtf
    file "lncRNA.final.CPAT.out" into lncRNA_CPAT_statistic
    file "protein_coding.final.CPAT.out" into protein_coding_CPAT_statistic

    shell:

    cufflinks_threads = ava_cpu.intdiv(2) - 1

    '''
        set -o pipefail
        cuffcompare -o filter \
                    -r !{knowlncRNAgtf} \
                    -p !{cufflinks_threads} !{novel_lncRNA_stringent_Gtf}
        awk '$3 =="u"||$3=="x"{print $5}' filter.novel.lncRNA.stringent.gtf.tmap |sort|uniq| \
                    perl !{baseDir}/bin/extract_gtf_by_name.pl !{novel_lncRNA_stringent_Gtf} - > novel.lncRNA.stringent.filter.gtf
        
        #rename lncRNAs according to neighbouring protein coding genes
        #awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} > gencode.protein_coding.gene.gtf
        #gtf2bed < gencode.protein_coding.gene.gtf |sort-bed - > gencode.protein_coding.gene.bed
        awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} | perl -F'\\t' -lane '$F[8]=~/gene_id "(.*?)";/ && print join qq{\\t},@F[0,3,4],$1,@F[5,6,1,2,7,8,9]' - | \
            sort-bed - > gencode.protein_coding.gene.bed
        gtf2bed < novel.lncRNA.stringent.filter.gtf |sort-bed - > novel.lncRNA.stringent.filter.bed
        gtf2bed < !{knowlncRNAgtf} |sort-bed - > known.lncRNA.bed
        perl !{baseDir}/bin/rename_lncRNA_2.pl
        cat novel.lncRNA.stringent.filter.gtf !{knowlncRNAgtf} > all_lncRNA_for_classifier.gtf
        perl rename_proteincoding.pl > protein_coding.final.gtf
        cat all_lncRNA_for_classifier.gtf protein_coding.final.gtf > final_all.gtf
        gffread final_all.gtf -g !{fasta_ref} -w final_all.fa -W
        
        #run statistic 
        perl compare_basic_charac.pl > basic_charac.txt
        '''
}


// star index and quantification by RSEM
// please make sure that RSEM is installed into your $PATH environment

process Reindex_kallisto_index_of_finalGtf_for_quantification{
    input:
    file transript_fasta from finalFasta_for_quantification_gtf

    output:
    file "transcripts.idx" into final_kallisto_index

    shell:
    '''
    #index star reference 
    kallisto index -i transcripts.idx ${transript_fasta}
    
    '''
}

process Run_kallisto_for_quantification{
    cpus idv_cpu
    maxForks fork_number
    tag {file_tag}
    input:
    file kallistoIndex from final_kallisto_index
    file pair from readPairs_for_quatification
    output:
    file("${file_tag_new}_kallisto") into kallisto_tcv_collection

    shell:
    file_tag = pair[0].name.replace("${params.suffix1}"+"."+"${params.fastq_ext}", "")
    file_tag_new = file_tag
    kallisto_threads = ava_cpu.intdiv(2) - 1
    if (params.singleEnd) {
        '''
        #quantification by kallisto in single end mode
        kallisto quant -i ${kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 --single -l 180 -s 20  <(zcat !{pair[0]} ) 
        mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
        '''
    } else {
        '''
        #quantification by kallisto 
        kallisto quant -i ${kallistoIndex} -o !{file_tag_new}_kallisto -t !{kallisto_threads} -b 100 <(zcat !{pair[0]} ) <(zcat !{pair[1]})
        mv !{file_tag_new}_kallisto/abundance.tsv !{file_tag_new}_abundance.tsv
        '''
    }
}

// combine matrix and for further evaluation
process get_kallisto_matrix{
    tag {file_tag}

    input:
    file abundance_tsv_matrix from kallisto_tcv_collection

    output:
    file "kallisto.count.txt" into expression_matrixfile_count
    file "kallisto.tpm.txt" into expression_matrixfile_tpm

    shell:
    '''
    R CMD BATCH !{baseDir}/bin/get_kallisto_matrix.R
    '''
}



// write out command information



//process gathering_Output_and_generate_Report {
//    publishDir "${baseDir}/Result", mode: 'copy'
//    echo true
//
//    input:
//    file ('alignment/*') from alignment_logs.collect().filter( ~/.final.log/ )
//
//    output:
//    file "*multiqc_report.html" into multiqc_report
//    file "*multiqc_data"
//
//    script:
//    """
//    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
//    multiqc -f -n lncRNAanalysis_report . 2>&1
//    """
//}





