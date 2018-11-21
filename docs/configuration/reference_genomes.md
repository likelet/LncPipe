# nf-core/lncpipe: Reference Genomes Configuration

We are sorry that the current lncPipe only support few types of organism: human, mouse, fly and zebrafish, as one essential tool CPAT included in LncPipe only available for 4 species

The nf-core/lncpipe pipeline needs a reference genome for alignment and annotation.

These paths can be supplied on the command line at run time (see the [usage docs](../usage.md)),
but for convenience it's often better to save these paths in a nextflow config file.
See below for instructions on how to do this.
Read [Adding your own system](adding_your_own.md) to find out how to set up custom config files.

## Adding paths to a config file
Specifying long paths every time you run the pipeline is a pain.
To make this easier, the pipeline comes configured to understand reference genome keywords which correspond to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline.

Note that this genome key can also be specified in a config file if you always use the same genome.

To use this system, add paths to your config file using the following template:

```nextflow
params {
  genomes {
    'YOUR-ID' {
      fasta  = '<PATH TO FASTA FILE>/genome.fa'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```
