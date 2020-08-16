<img align="middle" href="https://github.com/biocorecrg/exon_intron_orthology_pipeline" src="https://github.com/biocorecrg/exon_intron_orthology_pipeline/blob/master/docs/logo_s.png?raw=true" />


Summary
-------

ExOrthist is a Nextflow-based pipeline to obtain groups of exon orthologous at all evolutionary timescales.

Requirements
------------

ExOrthist requires the following software:
 * [Nextflow](https://www.nextflow.io/)
 * A linux container engine (either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/cli/singularity_apps.html))
 * R 3.5 or higher, with the [igraph](https://igraph.org/) package installed.
 * Perl 5.10.1 or higher
 * [Mafft](https://mafft.cbrc.jp/alignment/software/) for protein alignment.
 
 Additionally, [liftOver](https://genome-store.ucsc.edu/) and [bedtools](https://bedtools.readthedocs.io/en/latest/) are required to run XXXX-getLif.


Installation
------------

Install Nextflow (version 19.10.0):

~~~~
curl -s https://get.nextflow.io | bash`
~~~~

Clone the ExOrthist repository:
~~~~
git clone --depth 1 https://github.com/biocorecrg/ExOrthist.git
~~~~

Install Docker:

* Docker: https://docs.docker.com/install/ (version 10.03 or later is required)
* Singularity: https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps (version 2.6.1 or 3.2.1 is required)

A Docker image was uploaded to [Docker Hub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/exon_intron_pipe). 


## NextFlow pipeline
Test data are in the folder **test**. Input files are:
* Genomic sequences in fasta files (can be gzipped). The prefix, i.e. the genome identifier, must match the corresponding annotation file:
* Annotation files in GTF format (can be gzipped). The prefix must match the corresponding GTF file:
Example:

```bash
     Ame.fasta.gz Ame.gtf.gz
     Dme.fasta.gz Dme.gtf.gz 
```
* Gene cluster file with this format:
**PLEASE MANU OR YAMILE ADD SOME TEXT ON THIS FILE***

```bash
GF000002.01	Dme	FBgn0000810	fs(1)K10
GF000002.01	Ame	GB41489	fs(1)K10=Hypp632

GF000006.01	Dme	FBgn0031540	CG3238
GF000006.01	Ame	GB54729	CG3238=PIF1

GF000007.01	Dme	FBgn0004828	His3.3B
GF000007.01	Ame	GB54084	His3.3A=H3F3A
GF000007.01	Ame	GB54085	His3.3A=H3F3A
```
The prefix of the file should match the combination of genome identifiers like:
```bash
Ame_Dme.tab.gz
```

The pipeline can be launched in this way:
```
nextflow run main.nf --clusters "test/*.tab.gz" --genomes "test/GENOMES/*_gDNA.fasta.gz" --annotations "$baseDir/test/GTF/*_annot.gtf.gz" --output output -bg > log.txt
```
