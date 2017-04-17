# Workflow
## 1. Description
A collection of self-developed bioinformatics analysis pipeline coding by a DSL language, nextflow.
To run the pipeline files with nf suffixed, plz just type the following command.

```
nextflow <your nf file> -c nextflow.config -with-trace
```

If the pipeline fails at any point and you fix the issue, you can revoke the processes with job avoidance using this command:
```
nextflow <your nf file> -c nextflow.config -with-trace -resume
```

All those pipelines were written in Nextflow commands. For more details, please see the following link:
https://www.nextflow.io/



## 2. Installation of dependencies.
#### Install NextFlow
To run our pipelines. NextFlow should be preinstalled at  POSIX compatible system (Linux, Solaris, OS X, etc), It requires BASH and Java 7 or higher to be installed. We do not recommendate running the pipes in the Windows since most of bioinformatic tools are not supported.
Here, we show the step by step installation of NextFlow in linux system as an example, which adapted from [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html).

1. Download the executable package by copying and pasting the following command in your terminal window: ```wget -qO- get.nextflow.io | bash```. It will create the nextflow main executable file in the current directory.
2. Optionally, move the nextflow file in a directory accessible by your `$PATH` variable (this is only required to avoid to remember and type the Nextflow full path each time you need to run it).
3. Download the lastest binary verion of NextFlow from the https://github.com/nextflow-io/nextflow/releases and add the path into your system environment.
#### Install third-party software and database required by each pipe.
###### Pipe 1
