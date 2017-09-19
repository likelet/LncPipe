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
params.reads = "*_{1,2}.fastq.gz"
// Match the pairs on two channels
reads=params.input_folder+params.reads
Channel.fromFilePairs(reads,size: params.singleEnd ? 1 : 2)
        .ifEmpty { exit 1, print_red("Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\n" )}
        .into{reads_for_fastqc; readPairs_for_discovery}


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



