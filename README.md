# FiberGrowth

## Usage

```
Rscript FiberGrowth.R [options] --out=<output_folder> --gff=<gff_file> --proteins=<fasta_file>

Rscript FiberGrowth.R [options] --out=<output_folder> --genome=<fasta_file>

Rscript FiberGrowth.R (-h | --help)

Options:
  --gff=<gff_file>                          Gene locations in gff format.
  --proteins=<fasta_file>                   Amino acid sequences in fasta format.
  --genome=<fasta_file>                     Genome in fasta format.
  --lib=<path_to_library>                   Path to library of PUL models
                                            [default: PUL_models in installation directory].
  --out=<output_folder>                     Output directory (will be created).
  --gffFeatureName=<gffFeatureName>         Name of feature to use in gff file [default: CDS]
  --gffAttributeKey=<gffAttributeKey>       Name of attribute key to use in gff file [default: Name]
  -t --threads=<threads>                    Number of threads to run [default: 1].
  -h, --help                                Show help.

```

## Installation

### Dependencies

  - Hmmer (>3.1b2) (http://hmmer.org)
  - Prodigal (>2.6.3) (https://github.com/hyattpd/Prodigal)
  - R (>4)
  - R packages:
      - data.table
      - DT
      - docopt
      - gggenes
      - ggplot2
      - knitr
      - magrittr
      - rhmmer
      - rmarkdown


## Quick start

Running FiberGrowth on a given genome sequence in fasta format:
```
Rscript FiberGrowthDir/FiberGrowth.R --out=fibergrowth_output --genome=genome.fasta

```

Alternatively, FiberGrowth can be run on amino acid sequences in fasta format and gene locations in gff format:
```
Rscript FiberGrowthDir/FiberGrowth.R --out=fibergrowth_output --gff=genome.gff --proteins=proteins.faa

```
