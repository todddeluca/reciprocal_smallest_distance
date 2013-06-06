# Reciprocal Smallest Distance

Authors: Todd F. DeLuca, Dennis P. Wall  
Organization: Wall Laboratory, Center for Biomedical Informatics, Harvard Medical School, USA, Earth, Sol System, Orion Arm, Milky Way.  
Date: 2011/08/29  


## Introduction

Wall, D.P., Fraser, H.B. and Hirsh, A.E. (2003) Detecting putative orthologs, Bioinformatics, 19, 1710-1711.

The reciprocal smallest distance (RSD) (Wall, et al., 2003.
http://bioinformatics.oxfordjournals.org/content/19/13/1710) algorithm
accurately infers orthologs between pairs of genomes by considering global
sequence alignment and maximum likelihood evolutionary distance between
sequences.  Orthologs inferred with RSD for many species are available at
Roundup (http://roundup.hms.harvard.edu/), which provides multi-species
clusters of orthologous genes, output in formats for other phylogenetics
packages, and sequence metadata such as Gene Ontology terms and database
cross-references.

This package contains source code, scripts for running RSD, and example input
and output files.

- README.md:  the file you are reading now
- bin/rsd_search: a script that runs the reciprocal smallest distance (RSD)
  algorithm to search for orthologs.
- bin/rsd_blast: a script that computes and saves BLAST hits for use in
  multiple runs of RSD.
- bin/rsd_format: a script that turns FASTA-formatted genomes into
  BLAST-formatted indexes.
- rsd/: python package implementing the RSD algorithm.  
- rsd/jones.dat, rsd/codeml.ctl:  used by codeml/paml to compute the
  evolutionary distance between two sequences.
- examples/:  a directory containing examples of inputs and outputs to rsd,
  including fasta-formatted genome protein sequence files, a query sequence id
  file (for --ids), and an orthologs output file.


## Installing RSD

### Prerequisites

RSD depends on Python, NCBI BLAST, PAML, and Kalign.  It has been tested to
work with the versions below.  It might work with other versions too.

Install:

- Python 2.7: http://www.python.org/download/
- NCBI BLAST 2.2.24: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- PAML 4.4: http://abacus.gene.ucl.ac.uk/software/paml.html
- Kalign 2.04: http://msa.sbc.su.se/cgi-bin/msa.cgi

Add the executables for python (version 2.7), makeblastdb, blastp, codeml, and
kalign, to your PATH.


### Install Using Pip

Use pip ( http://www.pip-installer.org/ ) to easily install
reciprocal\_smallest\_distance:

    pip install reciprocal_smallest_distance


### Install From Github

Installing from github is an easy way to make the `examples/` directory
accessible for running the examples in this README file.

Clone the latest version from GitHub:

    cd ~  # Optional
    git clone https://github.com/todddeluca/reciprocal_smallest_distance

Install reciprocal\_smallest\_distance, making sure to use Python 2.7:

    cd reciprocal_smallest_distance
    python setup.py install


## Using RSD to Find Othologs

The following example commands demonstrate the main ways to run `rsd_search`.
Every invocation of `rsd_search` requires specifying the location of a
FASTA-formatted sequence file for two genomes, called the query and subject
genomes.  Their order is arbitrary, but if you use the `--ids` option, the ids
must come from the query genome.  You must also specify a file to write the
results of the orthologs found by the RSD algorithm.  The
[formats](#output_formats) of the output file are described in detail below.
You can optionally specify a file containing ids using the `--ids` option.
Then rsd will only search for orthologs for those ids.  Using `--de` gives you
the option of using different divergence and evalue thresholds from the
defaults.


Get help on how to run `rsd_search`, `rsd_blast`, or `rsd_format`:

    rsd_search -h
    rsd_blast -h
    rsd_format -h


Find orthologs between all the sequences in the query and subject genomes,
using the default divergence and evalue thresholds:

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt


Find orthologs using several divergence and evalue threshold combinations:

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa.several.orthologs.txt \
    --de 0.2 1e-20 --de .5 0.00001 --de 0.8 0.1


It is not necessary to format a FASTA file for BLAST or compute BLAST hits
because `rsd_search` does it for you.  However if you plan on running
`rsd_search` multiple times for the same genomes, especially for large genomes,
you can save time by using `rsd_format` to preformatting the FASTA files and
`rsd_blast` to precomputing the BLAST hits.  When running `rsd_blast`, make
sure to use an --evalue as large as the largest evalue threshold you intend to
give to `rsd_search`.

Here is how to format a pair of FASTA files in place:

    rsd_format -g examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa
    rsd_format -g examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa

And here is how to format the FASTA files, putting the results in another
directory (the current directory in this case):

    rsd_format -g examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa -d .
    rsd_format -g examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa -d .

Here is how to compute forward and reverse blast hits (using the default
evalue):

    rsd_blast -v -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    --forward-hits q_s.hits --reverse-hits s_q.hits

Here is how to compute forward and reverse blast hits for `rsd_search`, using
genomes that have already been formatted for blast and a non-default evalue:

    rsd_blast -v -q Mycoplasma_genitalium.aa \
    --subject-genome=Mycobacterium_leprae.aa \
    --forward-hits q_s.hits --reverse-hits s_q.hits \
    --no-format --evalue 0.1

Find orthologs between all the sequences in the query and subject genomes using
genomes that have already been formatted for blast:

    rsd_search -q Mycoplasma_genitalium.aa \
    --subject-genome=Mycobacterium_leprae.aa \
    -o Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt \
    --no-format


Find orthologs between all the sequences in the query and subject genomes using
hits that have already been computed.  Notice that --no-format is included,
because since the blast hits have already been computed the genomes do not need
to be formatted for blast:

    rsd_search -v --query-genome Mycoplasma_genitalium.aa \
    --subject-genome=Mycobacterium_leprae.aa \
    -o Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa.default.orthologs.txt \
    --forward-hits q_s.hits --reverse-hits s_q.hits --no-format


Find orthologs for specific sequences in the query genome.  For finding
orthologs for only a few sequences, using `--no-blast-cache` can speed up
computation.  YMMV.

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o examples/Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt \
    --ids examples/Mycoplasma_genitalium.aa.ids.txt --no-blast-cache

<a name="output_formats"/>
## Output Formats

Orthologs can be saved in several different formats using the `--outfmt` option
of `rsd_search`.  The default format, `--outfmt -1`, refers to `--outfmt 3`.
Inspired by Uniprot dat files, a set of orthologs starts with a parameters
line, then has 0 or more ortholog lines, then has an end line.  The parametes
are the query genome name, subject genome name, divergence threshold, and
evalue threshold.  Each ortholog is on a single line listing the query sequence
id, the subject sequence id, and the maximum likelihood distance estimate.
This format can represent orthologs for multiple sets of parameters in a single
file as well as sets of parameters with no orthologs.  Therefore it is suitable
for use with `rsd_search` when specifying multiple divergence and evalue
thresholds.

Here is an example containing 2 parameter combinations, one of which has no
orthologs:

    PA\tLACJO\tYEAS7\t0.2\t1e-15
    OR\tQ74IU0\tA6ZM40\t1.7016
    OR\tQ74K17\tA6ZKK5\t0.8215
    //
    PA\tMYCGE\tMYCHP\t0.2\t1e-30
    //

The original format of RSD, `--outfmt 1`, is provided for backward
compatibility.  Each line contains an ortholog, represented as subject sequence
id, query sequence id, and maximum likelihood distance estimate.  It can only
represent a single set of orthologs in a file.

Example:

    A6ZM40\tQ74IU0\t1.7016
    A6ZKK5\tQ74K17\t0.8215

Also provided for backward compatibility is a format used internally by Roundup
(http://roundup.hms.harvard.edu/) which is like the original RSD format, except
the query sequence id column is before the subject sequence id.

Example:

    Q74IU0\tA6ZM40\t1.7016
    Q74K17\tA6ZKK5\t0.8215


