# nf-core/lncpipe: Troubleshooting

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4.  If you are running Single end data make sure to specify `--singleEnd`

If the pipeline can't find your files then you will get the following error

```
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will.


## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/lncpipe) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).


## Tips & FAQs

## Tips	 

* :blush:Plz keep the consistency of your genome sequence, index library and annotation files: genome version, chromosome format, gtf coordinated e.g. The third-party software may stop for any of the above reasons. 
* :confused:Setting your analysis parameters always in config file, differ project should  corresponding to differ configurations for reproductive analysis. To rerun a project, you can just specify -c `your.config` in your command, which can also help you to record analysis parameters.
* :open_mouth:Run analysis on docker container, no much to say.
* :grimacing:Always use the latest version to be away from the known bugs. 


## FAQ  
In local mode: 
* *1. PLEK throws an error "/data/software/PLEK.1.2/PLEK.py:line12: $'\r': can not find command", how to fix?*
>A: using the follow command as suggested in the installation section. 
    
        perl -CD -pi -e'tr/\x{feff}//d && s/[\r\n]+/\n/' *.py 
    
* *2. IOError: [Errno 2] No such file or directory: '/opt/CPAT-1.2.3/dat/Human_Hexamer.tsv'?*
>A: The cpat command required  the `Human_Hexamer.tsv` to predict lncRNA coding potential, plz check your `cpatpath` parameters. 
* *3. When using htseq to quantify transicript, it throws "Error occured when reading beginning of SAM/BAM file. 'csamtools.AlignedRead' object has no attribute 'reference_start' "*
>A: It's a version conflict caused by htseq and hisat generated bamfile, a possible solution for this is to install the old version of htseq


