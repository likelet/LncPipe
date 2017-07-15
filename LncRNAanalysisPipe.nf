#! /usr/bin/env nextflow

// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"

// requirement:
// - STAR
// - Cufflinks
// - Bedops

// Pipeline version
version = '0.0.4'

//user options
if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW Long non-coding RNA analysis PIPELINE v!{version}'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow lncRNApipe.nf '

    exit 1
}




//default values
params.help = null
params.input_folder = './'
params.out_folder = './'

// HG38 genome ref files
params.fasta_ref = '/data/database/hg38/genome.fa'
params.star_idex = '/data/database/hg38/GENCODE/STARIndex'
params.gencode_annotation_gtf = "/data/database/hg38/GENCODE/gencode.v25.annotation.gtf"
params.lncipedia_gtf = "/data/database/hg38/LNCipedia/lncipedia_4_0_hc_hg38.gtf"
params.rRNAmask = "/data/database/hg38/lncRNAanalysisPipeFile/rRNA_hg38.gtf";
// hg19 genome ref files
//params.fasta_ref = '/data/database/hg19/genome.fa'
//params.star_idex = '/data/database/hg38/GENCODE/STARIndex'
//params.gencode_annotation_gtf = "/data/database/hg38/GENCODE/gencode.v25.annotation.gtf"
//params.lncipedia_gtf = "/data/database/hg38/LNCipedia/lncipedia_4_0_hc_hg38.gtf"
// software
params.plekpath = '/data/software/PLEK.1.2'
params.cncipath = '/data/software/CNCI-master'
params.cpatpath = '/data/software/CPAT-1.2.2/'
//

// fastq file
params.fastq_ext = "fastq.gz"

params.merged_gtf=null
//params.fastq_ext2 = "fq.gz"
params.suffix1 = "_1"
params.suffix2 = "_2"
// run information of systemfile
params.cpu = 16
params.mem = 32

// read file

fasta_ref = file(params.fasta_ref)
star_idex = file(params.star_idex)
input_folder = file(params.input_folder)
gencode_annotation_gtf = file(params.gencode_annotation_gtf)
lncipedia_gtf = file(params.lncipedia_gtf)
plekpath = file(params.plekpath)
cncipath = file(params.cncipath)
cpatpath = file(params.cpatpath)

rRNAmaskfile = file(params.rRNAmask)

//specify  weather the merged file was already generated



annotation_channel = Channel.from(gencode_annotation_gtf, lncipedia_gtf)
annotation_channel.collectFile { file -> ['lncRNA.gtflist', file.name + '\n'] }
        .set { LncRNA_gtflist }
process combine_public_annotation {
    cpus params.cpu

    input:
    file lncRNA_gtflistfile from LncRNA_gtflist
    file gencode_annotation_gtf
    file lncipedia_gtf

    output:
    file "gencode_protein_coding.gtf" into proteinCodingGTF
    file "known.lncRNA.gtf" into KnownLncRNAgtf

    shell:
    cufflinks_threads = params.cpu.intdiv(2) - 1
    '''
        set -o pipefail
        cuffmerge -o merged_lncRNA !{lncRNA_gtflistfile}
        cat !{gencode_annotation_gtf} |grep "protein_coding" > gencode_protein_coding.gtf
        cuffcompare -o merged_lncRNA -r !{gencode_annotation_gtf} -p !{cufflinks_threads} merged_lncRNA/merged.gtf
        awk '$3 =="u"||$3=="x"{print $5}' merged_lncRNA/merged_lncRNA.merged.gtf.tmap \
        |sort|uniq|perl !{baseDir}/bin/extract_gtf_by_name.pl merged_lncRNA/merged.gtf - > merged.filter.gtf
        mv  merged.filter.gtf known.lncRNA.gtf
        
        '''

}

// whether the merged gtf have already produced.
if (params.merged_gtf==null) {
//gene_annotation=file(params.gene_annotation)




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

//prepare annotations
// combine public lnRNA annotation
        //required files
                //#lncRNA.gtflist
                //#gencode.v24.long_noncoding_RNAs.gtf
                //#lncipedia_4_0_hg38.gtf



        //run star iterate over sample or parallel
        paramode = true
// Star alignment
        if (paramode) {
            process fastq_star_alignment_para {
                cpus params.cpu
                tag { file_tag }


                input:
                file pair from readPairs
                file fasta_ref
                file star_idex

                output:
                set val(file_tag_new), file("STAR_${file_tag_new}") into STARmappedReads,alignment_logs
//                file("STAR_${file_tag_new}/*final.out")into alignment_logs

                shell:
                file_tag = pair[0].name.replace("${params.suffix1}.${params.fastq_ext}", "")
                file_tag_new = file_tag
                star_threads = params.cpu.intdiv(2) - 1

                '''
                # create folder for out putfile
                
                STAR --runThreadN !{star_threads}  --twopassMode Basic --genomeDir !{star_idex} \
                --readFilesIn !{pair[0]} !{pair[1]} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
                --chimSegmentMin 20 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 20 \
                --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout \
                --alignSJoverhangMin 8 --alignSJDBoverhangMin 1  --outFileNamePrefix !{file_tag_new} 
                mkdir STAR_!{file_tag_new}
                mv !{file_tag_new}Aligned* STAR_!{file_tag_new}/.
                mv !{file_tag_new}SJ* STAR_!{file_tag_new}/.
                mv !{file_tag_new}Log* STAR_!{file_tag_new}/.
                '''
            }
        } else {
            process fastq_star_alignment_para {
                cpus params.cpu
                tag { file_tag }
/* publishDir params.out_folder, mode: 'move', overwrite: true*/
                input:
                file pair from readPairs
                file fasta_ref
                file star_idex

                output:
                set val(file_tag_new), file("STAR_${file_tag_new}") into STARmappedReads,alignment_logs
//                file "!{file_tag_new}_star_log.txt" into alignment_logs
                shell:
                file_tag = pair[0].name.replace("${params.suffix1}.${params.fastq_ext}", "")
                file_tag_new = file_tag
                star_threads = params.cpu.intdiv(2) - 1

                '''
                for every pair in ${pair}
                do
                    STAR --runThreadN !{star_threads}  --twopassMode Basic --genomeDir !{star_idex} \
                        --readFilesIn !{pair[0]} !{pair[1]} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
                        --chimSegmentMin 20 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 20 \
                        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout \
                        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1  --outFileNamePrefix !{file_tag_new} > !{file_tag_new}_star_log.txt
                    mkdir STAR_!{file_tag_new}
                    mv !{file_tag_new}Aligned* STAR_!{file_tag_new}/.
                    mv !{file_tag_new}SJ* STAR_!{file_tag_new}/.
                    mv !{file_tag_new}Log* STAR_!{file_tag_new}/.
                done
        
                '''
            }

        }

// cufflinks
        process cufflinks_assembly {
            cpus params.cpu
            tag { file_tag }

            input:
            set val(file_tag), file(Star_alignment) from STARmappedReads
            file fasta_ref
            file gencode_annotation_gtf
            file rRNAmaskfile

            output:

            file "Cufout_${file_tag_new}_transcripts.gtf" into cuflinksoutgtf, cuflinksoutgtf_fn

            shell:
            file_tag_new = file_tag
            cufflinks_threads = params.cpu.intdiv(2) - 1

            '''
            #run cufflinks
            
            cufflinks -g !{    gencode_annotation_gtf} -b !{fasta_ref} --library-type fr-firststrand \
            --max-multiread-fraction 0.25 --3-overhang-tolerance 2000 \
            -o Cufout_!{file_tag_new} -p 25 !{Star_alignment}/!{file_tag_new}Aligned.sortedByCoord.out.bam
            mv Cufout_!{file_tag_new}/transcripts.gtf Cufout_!{file_tag_new}_transcripts.gtf
            '''

        }

// Create a file 'gtf_filenames' containing the filenames of each post processes cufflinks gtf

        cuflinksoutgtf.collectFile { file -> ['gtf_filenames.txt', file.name + '\n'] }
                .set { GTFfilenames }

// run cuffmerge

        process cuffmerge_assembled_gtf {
            cpus params.cpu
            tag { file_tag }

            input:
            file gtf_filenames from GTFfilenames
            file cufflinksgtf_file from cuflinksoutgtf_fn.toList() // not used but just send the file in current running folder

            file fasta_ref


            output:
            file "CUFFMERGE/merged.gtf" into cuffmergeTranscripts_forCompare, cuffmergeTranscripts_forExtract, cuffmergeTranscripts_forCodeingProtential
            shell:

            cufflinks_threads = params.cpu.intdiv(2) - 1

            '''
            mkdir CUFFMERGE
            cuffmerge -o CUFFMERGE \
            -s !{fasta_ref} \
            -p !{cufflinks_threads} \
            !{gtf_filenames}
            #    perl !{baseDir}/bin/get_exoncount.pl CUFFMERGE/merge.gtf CUFFMERGE/transcript_exoncount.txt
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
process cuffcompare_GenCODE {
    cpus params.cpu
    tag { file_tag }
//
    input:
    file cuffmergefile from cuffmergeTranscripts_forCompare
    file gencode_annotation_gtf

    output:
// file "gencode.v25.protein_coding.gtf" into knownProteinCoding
    file "CUFFCOMPARE_GENCODE" into cuffcomparegencodeDir
    file ""
    shell:

    cufflinks_threads = params.cpu.intdiv(2) - 1
//
    '''
        mkdir CUFFCOMPARE_GENCODE
        #!/bin/sh
        cuffcompare -o CUFFCOMPARE_GENCODE/merged_lncRNA -r !{gencode_annotation_gtf} -p !{cufflinks_threads} !{cuffmergefile}
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
    //file "transcript_exoncount.txt" into exoncount
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
        #check wether required
        # get fasta frome gtf
        gffread novel.longRNA.gtf -g !{fasta_ref} -w novel.longRNA.fa -W
        '''
}

// predicting coding potential _ parallel
// copy fasta channel into three
novelLncRnaFasta.into { novelLncRnaFasta_for_PLEK; novelLncRnaFasta_for_CPAT; novelLncRnaFasta_for_CNCI }

process run_PLEK {
    cpus params.cpu
    validExitStatus 0,1,2
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_PLEK
    file plekpath
    output:
    file "novel.longRNA.PLEK.out" into novel_longRNA_PLEK_result
    shell:
    plek_threads = params.cpu.intdiv(2) - 1
    '''
        python !{plekpath}/PLEK.py -fasta !{novel_lncRNA_fasta} -out novel.longRNA.PLEK.out -thread !{plek_threads}
	#exit 0
        '''

}
process run_CPAT {
    input:
    file novel_lncRNA_fasta from novelLncRnaFasta_for_CPAT
    file cpatpath
    output:
    file "novel.longRNA.CPAT.out" into novel_longRNA_CPAT_result
    shell:
    '''
        python !{cpatpath}/bin/cpat.py -g !{novel_lncRNA_fasta} -x !{cpatpath}/dat/Human_Hexamer.tsv \
            -d !{cpatpath}/dat/Human_logitModel.RData -o novel.longRNA.CPAT.out
        '''
}

//    process run_CNCI{
//        cpus params.cpu
//        input:
//        file novel_lncRNA_fasta from novelLncRnaFasta_for_CNCI
//        file cncipath
//        output:
//        file lncRNA/CNCI* into novel_longRNA_CNCI_result
//        shell:
//        cnci_threads  = params.cpu.intdiv(2) - 1
//        '''
//        python !{cncipath}/CNCI.py -f !{novel_lncRNA_fasta} -p !{cnci_threads} -o lncRNA/CNCI -m ve
//        '''
//    }

// merge transcripts and retian lncRNA only by coding ability
process merge_filter_by_coding_potential {
    input:
    file novel_longRNA_PLEK_ from novel_longRNA_PLEK_result
    file novel_longRNA_CPAT_ from novel_longRNA_CPAT_result
    file longRNA_novel_exoncount from novelLncRnaExonCount
    file cuffmergegtf from cuffmergeTranscripts_forCodeingProtential
    file gencode_annotation_gtf

    output:
    file "novel.longRNA.stringent.gtf" into Novel_longRNA_stringent_gtf // not used
    file "novel.lncRNA.stringent.gtf" into novel_lncRNA_stringent_gtf
    file "novel.TUCP.stringent.gtf" into novel_TUCP_stringent_gtf // not used


    shell:

    '''
        set -o pipefail
        #merged transcripts
        perl !{baseDir}/bin/integrate_novel_transcripts.pl > novel.longRNA.txt
        awk '$4 >1{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - \
        >novel.longRNA.stringent.gtf
        # retain lncRNA only by coding ability
        awk '$4 >1&&$5=="lncRNA"{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - \
        >novel.lncRNA.stringent.gtf
        awk '$4 >1&&$5=="TUCP"{print $1}' novel.longRNA.txt|perl !{baseDir}/bin/extract_gtf_by_name.pl !{cuffmergegtf} - \
        >novel.TUCP.stringent.gtf
        '''
}

//Further  novel lncRNA based on annotated database
process Filter_lncRNA_based_annotationbaes {
    publishDir "${baseDir}/Result", mode: 'copy'
    cpus params.cpu

    input:
    file knowlncRNAgtf from KnownLncRNAgtf
    file novel_lncRNA_stringent_Gtf from novel_lncRNA_stringent_gtf
    file gencode_protein_coding_gtf from proteinCodingGTF

    output:
    file "lncRNA.final.v2.gtf" into finalLncRNA_gtf
    file "lncRNA.final.v2.map" into finalLncRNA_map
    shell:
    cufflinks_threads = params.cpu.intdiv(2) - 1
    '''
        set -o pipefail
        cuffcompare -o filter -r !{knowlncRNAgtf} -p !{cufflinks_threads} !{novel_lncRNA_stringent_Gtf}
        awk '$3 =="u"||$3=="x"{print $5}' filter.novel.lncRNA.stringent.gtf.tmap |sort|uniq| \
        perl !{baseDir}/bin/extract_gtf_by_name.pl !{novel_lncRNA_stringent_Gtf} - > novel.lncRNA.stringent.filter.gtf
        
        #rename lncRNAs according to neighbouring protein coding genes
        #awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} > gencode.protein_coding.gene.gtf
        #gtf2bed < gencode.protein_coding.gene.gtf |sort-bed - > gencode.protein_coding.gene.bed
        awk '$3 =="gene"{print }' !{gencode_protein_coding_gtf} | perl -F'\\t' -lane '$F[8]=~/gene_id "(.*?)";/ && print join qq{\\t},@F[0,3,4],$1,@F[5,6,1,2,7,8,9]' - | sort-bed - > gencode.protein_coding.gene.bed
        gtf2bed < novel.lncRNA.stringent.filter.gtf |sort-bed - > novel.lncRNA.stringent.filter.bed
        gtf2bed < !{knowlncRNAgtf} |sort-bed - > known.lncRNA.bed
        perl !{baseDir}/bin/rename_lncRNA_2.pl
        #perl rename_proteincoding.pl > protein_coding.final.gtf
        '''

}


//process multiqc {
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





