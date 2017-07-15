#! /usr/bin/env nextflow

// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"

// requirement:
// - bwa
// - picard
// - samtools/sambamba
// - annovar
// - bcftools


// Pipeline version
version = '0.0.1'

//user options
if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW Exome seq data processing PIPELINE v!{version}'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow ExomeseqPipe.nf '

    exit 1
}




//default values
params.help = null
params.input_folder = './'
params.out_folder = './'

// database file folder genome ref files
params.genome_bwa_idex="/disk/database/human/hg19/bwaIndex"
params.genome_ref_exome="/disk/database/human/hg19"
//params.annovar_db="/disk/soft/Annovar/humandb"


// software
params.annovarpath = ''
//

// fastq file
params.fastq_ext = "fastq.gz"

//params.fastq_ext2 = "fq.gz"
params.suffix1 = "_1"
params.suffix2 = "_2"
// run information of systemfile
params.cpu = 24
params.mem = 64

// read file

genome_bwa_idex = file(params.genome_bwa_idex)
genome_ref_exome = file(params.genome_ref_exome)
annovar_db = file(params.annovar_db)
//annovarpath = file(params.annovarpath)


//specify  weather the merged file was already generated






//Star index
//star_ref = file(params.params.star_idex_ref)

// Check whether fastq file is available
    mode = 'fastq'
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0 ||
            file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext2}/ }.size() > 0) {
        println "fastq files found, proceed with alignment"
    } else {
        if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0) {
            println "BAM files found, proceed with realignment"; mode = 'bam';
            files = Channel.fromPath(params.input_folder + '/*.bam')
        } else {
            println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
        }
    }


    if (mode == 'fastq') {
        println "Analysis from fastq file"
        println "Start mapping with STAR aligner"
//
        keys1 = file(params.input_folder).listFiles().findAll {
            it.name ==~ /.*${params.suffix1}.${params.fastq_ext}/
        }.collect { it.getName() }
                .collect { it.replace("${params.suffix1}.${params.fastq_ext}", '') }
        keys2 = file(params.input_folder).listFiles().findAll {
            it.name ==~ /.*${params.suffix2}.${params.fastq_ext}/
        }.collect { it.getName() }
                .collect { it.replace("${params.suffix2}.${params.fastq_ext}", '') }
        if (!(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1))) {
            println "\n ERROR : There is at least one fastq without its mate, please check your fastq files.";
            System.exit(0)
        }

        println keys1
// parse paired files _1
        reads1 = Channel.fromPath(params.input_folder + '/*' + params.suffix1 + '.' + params.fastq_ext).map { path -> [path.name.replace("${params.suffix1}.${params.fastq_ext}", ""), path] }

// parse paired files _2
        reads2 = Channel.fromPath(params.input_folder + '/*' + params.suffix2 + '.' + params.fastq_ext).map { path -> [path.name.replace("${params.suffix2}.${params.fastq_ext}", ""), path] }

// Match the pairs on two channels
        readPairs = reads1.phase(reads2).map { pair1, pair2 -> [pair1[1], pair2[1]] }


        process bwa_aligment{
            cpus params.cpu
            tag { file_tag }
            input:
            file pair from readPairs
            file genome_bwa_idex


            output:
            set val(file_tag_new), file("${file_tag_new}.bam") into mapped_bam
            file
            shell:
            file_tag = pair[0].name.replace("${params.suffix1}.${params.fastq_ext}", "")
            file_tag_new = file_tag
            bwa_threads = params.cpu.intdiv(2) - 1

            '''
                TEMPDIR=${pwd}/tmp
                mkdir $TEMPDIR
                for every pair in ${pair}
                do
                bwa mem -t !{bwa_threads} !{genome_bwa_idex} !{pair[0]} !{pair[1]}  > !{file_tag_new}.sam
                sambamba  view -S -f bam -t !{bwa_threads} !{file_tag_new}.sam > ${file_tag_new}.bam
                
                rm -rf !{file_tag_new}.sam
                done
                '''
        }

        process bamfile_process{
            cpus params.cpu
            tag { file_tag }
            input:
            set val(file_tag), file(bamfile) from mapped_bam


            output:
            set val(file_tag_new), file("${file_tag_new}_sorted_dedup.bam") into mapped_sorted_dedup_bam_forsamtools
            shell:
            file_tag_new = file_tag
            sambamba_threads = params.cpu.intdiv(2) - 1

            '''
                TEMPDIR=${pwd}/tmp
                mkdir $TEMPDIR
                sambamba sort  -t !{sambamba_threads} -o ${file_tag_new}_sorted.bam ${file_tag_new}.bam
                sambamba markdup -t !{sambamba_threads}  ${file_tag_new}_sorted.bam ${file_tag_new}_sorted_dedup.bam
                '''
        }


        process samtools_bcgtools_process{
            cpus params.cpu
            tag { file_tag }
            //publishDir "${baseDir}/Result", mode: 'copy'
            input:
            set val(file_tag), file(bamfile) from mapped_bam
            file genome_ref_exome

            output:
            set val(file_tag_new), file("${file_tag_new}.vcf") into samtools_vcf_file
            shell:
            file_tag_new = file_tag
            sambamba_threads = params.cpu.intdiv(2) - 1
            '''
            set -o pipefail
                samtools mpileup -C 0 -A -B -d 10000 -v -u -f !{genome_ref_exome} !{bamfile} | \
                !{bcftoolspath}/bcftools call -O v -v -c -n 0.05 -p 1 -A -o !{file_tag_new}.vcf
            '''
        }

    }