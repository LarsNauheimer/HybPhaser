# HybPhaser
Detecting and phasing of hybrid accessions in target capture datasets.

HybPhaser was developed to deal with hybrids (and polyploids) in target capture datasets. 

It detects hybrids by measuring heterozygosity in the dataset and phase hybrid accessions by separating reads according to similarity with selected taxa that represent parental clades. 

HybPhaser is built as an extension to the efficient assembly pipeline HybPiper. 

## Installation

HybPhaser Installation

HybPhaser scripts can be downloaded from GitHub

git clone https://github.com/larsnauheimer/HybPhaser.git

Software dependencies

    R (v4.0)
    R-packages:
        ape (v5.4)
        sequinR (v4.2)
        stringR (v1.4)
    BWA
    SAMtools
    Bcftools
    BBSplit (BBMap v.38.87)
    HybPiper (which has additional dependencies)

## Data Preparation

Prior to run HybPhaser, sequence assembly has to be performed using HybPiper (Johnson et al. 2016).

The use of HybPiper is well-explained in the HybPiper-Wiki.

HybPiper requires the sequence reads and a fasta file with the target sequences.

In short, HybPiper pre-selects reads that match to a target gene using BWA or BLAST and then performs a de novo assmebly using Spades of the pre-selected read files. Extronerate is then used to extract and concatenate exon regions to generate gene sequences. Optionally intronerate can be used to recover any available intron regions and concatenate these with the exons to 'supercontigs'.

HybPiper generates one folder for each sample, which contains subfolder for each gene that contains assembly files as well as the assembled contig sequences. HybPhaser uses several of the output files of HybPiper for further analyses. The 'cleanup.py' script of HybPiper can be used before running HybPhaser.
