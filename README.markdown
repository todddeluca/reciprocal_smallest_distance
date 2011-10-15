Author: Todd F. DeLuca  
Organization: Wall Laboratory, Center for Biomedical Informatics, Harvard Medical School, USA, Earth, Sol System, Orion Arm, Milky Way.  
Date: 2011/08/29  

This package contains the scripts needed to run the Reciprocal Smallest Distance (RSD) ortholog detection algorithm as well as examples of input and output files.

- README.md:  the file you are reading now
- bin/rsd_search: a script that runs the reciprocal smallest distance (RSD) algorithm to search for orthologs.
- bin/rsd_format: a script that turns FASTA-formatted genomes into BLAST-formatted indexes.
- rsd/: package implementing the RSD algorithm.  
- rsd/jones.dat, rsd/codeml.ctl:  used by codeml/paml to compute the evolutionary distance between two sequences.
- examples/:  a directory containing examples of inputs and outputs to rsd, including fasta-formatted genome protein sequence files,
 a query sequence id file (for --ids), and an orthologs output file.


## Installing RSD

### Prerequisites

RSD depends on Python, NCBI BLAST, PAML, and Kalign.  It has been tested to work with the following versions.  It might work with other versions too.

Install:

- Python 2.7: http://www.python.org/download/
- NCBI BLAST 2.2.24: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- PAML 4.4: http://abacus.gene.ucl.ac.uk/software/paml.html
- Kalign 2.04: http://msa.sbc.su.se/cgi-bin/msa.cgi

Add the executables for python (version 2.7), makeblastdb, blastp, codeml, and kalign, to your PATH.

### Installing from a tarball

Download and untar the latest version from github:

    cd ~
    curl https://github.com/downloads/todddeluca/reciprocal_smallest_distance/reciprocal_smallest_distance-VERSION.tar.gz | tar xvz
    
Install reciprocal\_smallest\_distance, making sure to use Python 2.7:

    cd reciprocal_smallest_distance-VERSION
    python setup.py install


## Using RSD to Find Othologs

The following example commands demonstrate the main ways to run `rsd_search`.
Every invocation of `rsd_search` requires specifying the location of a FASTA-formatted sequence file for two genomes,
called the query and subject genomes.  Their order is arbitrary, but if you use the `--ids` option, the ids must come from the query genome.
You must also specify a file to write the results of the orthologs found by the RSD algorithm.
The format of the output file contains one ortholog per line.  Each line contains the query sequence id, subject sequence id,
and distance (calculated by codeml) between the sequences.
You can optionally specify a file containing ids using the `--ids` option.  Then rsd will only search for orthologs for those ids.
Using `--divergence` and `--evalue`, you have the option of using different thresholds from the defaults.


Get help on how to run rsd_search:
    
    rsd_search -h


Find orthologs between all the sequences in the query and subject genomes.

    rsd_search -d 0.2 -e 1e-20 -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o examples/Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.2_1e-20.orthologs.txt


Find orthologs using non-default divergence and evalue thresholds

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o examples/Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt \
    -d 0.8 -e 1e-5


Find orthologs between all the sequences in the query and subject genomes using genomes that have already been formatted for blast.

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o examples/Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt \
    --no-format


Find orthologs for specific sequences in the query genome.  For finding orthologs for only a few sequences, using `--no-blast-cache` can
speed up computation.  YMMV.

    rsd_search -q examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa \
    --subject-genome=examples/genomes/Mycobacterium_leprae.aa/Mycobacterium_leprae.aa \
    -o examples/Mycoplasma_genitalium.aa_Mycobacterium_leprae.aa_0.8_1e-5.orthologs.txt \
    --ids examples/Mycoplasma_genitalium.aa.ids.txt --no-blast-cache


## Formatting FASTA Files For BLAST

It is not necessary to format a fasta file for BLAST, because rsd_search does it for you.  However should you find yourself running
the program multiple times, especially for large genomes which take some time to format, preformatting the fasta files and 
then running `rsd_search` with the `--no-format` option can save you time.  Here is how to format a FASTA file in place:

    rsd_format -g examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa -v
    
And here is how to format the FASTA file, putting the results in a specific directory:

    rsd_format -g examples/genomes/Mycoplasma_genitalium.aa/Mycoplasma_genitalium.aa -v -d /my/place/for/genomes


## Extras

A user might also want to change the alpha shape parameter for the gamma distribution used in the likelihood calculations of the codeml package of paml.  This can be done by editing the codeml.ctl file included in the distribution.  See the documentation for codeml for more details.
