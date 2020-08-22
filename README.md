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
 
 Additionally, [liftOver](https://genome-store.ucsc.edu/), [bedtools](https://bedtools.readthedocs.io/en/latest/) as well as specific pairwise [liftOver files](http://hgdownload.soe.ucsc.edu/downloads.html#liftover) are required to run GetLiftOverFile.pl (see [below](#adding-manually-curated-exon-orthology-pairs)).


Installation
------------

Install Nextflow (version 19.10.0):

~~~~
curl -s https://get.nextflow.io | bash
~~~~

Clone the ExOrthist repository:
~~~~
git clone https://github.com/biocorecrg/ExOrthist.git
~~~~

Install Docker:

* Docker: https://docs.docker.com/install/ (version 10.03 or later is required)
* Singularity: https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps (version 2.6.1 or 3.2.1 is required)

A Docker image was uploaded to [Docker Hub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/exon_intron_pipe). 


Running ExOrthist main module
------------

The pipeline can be launched in this way:
```bash
nextflow run main.nf -bg > log.txt
```

A config file "params.config" is need, with the following information:

```
params {
    cluster      = "$baseDir/test/Ame_Cdi_Dme_s.tab.gz"
    genomes      = "$baseDir/test/GENOMES/*_gDNA.fasta.gz"
    annotations  = "$baseDir/test/GTF/*_annot.gtf.gz"
    clusternum   = 100
    extraexons   = ""
    liftover     = ""
    orthofolder  = ""
    vastdb       = ""
    intcons      = 2
    idexons      = 0.2
    maxsize      = 0.65
    output       = "$baseDir/output_small"
    email        = "yourmail@yourdomain"
}

```


Test data are in the folder **test**. Input files are:

* Genomic sequences in fasta files (can be gzipped). The prefix, i.e. the genome identifier, must match the corresponding annotation file:
```bash
     Hsa.fasta.gz
     Mmu.fasta.gz
     Bta.fasta.gz
```

* Annotation files in GTF format (can be gzipped). The prefix must match the corresponding GTF file:
Example:

```bash
     Hsa.gtf.gz
     Mmu.gtf.gz 
     Bta.gtf.gz
```

* Gene cluster file with the following format (tsv): ClusterID Species GeneID
```bash
GF000001	Hsa	ENSG00000151690
GF000001	Mmu	ENSMUSG00000041439
GF000001	Bta	ENSBTAG00000007719
GF000002	Hsa	ENSG00000091127
GF000002	Mmu	ENSMUSG00000057541
GF000002	Bta	ENSBTAG00000007743
GF000003	Hsa	ENSG00000029534
GF000003	Mmu	ENSMUSG00000031543
GF000003	Bta	ENSBTAG00000003275
```

The Species identifier should match that of the other files. Moreover, the prefix of the file should match the combination of genome identifiers like:
```bash
Hsa_Mmu_Bta.tab.gz
```


To provide additional non-annotated exons... The following format (tsv): XXXXX.

```bash
???
```


Adding manually curated exon orthology pairs
------------

To add pairs of manually curated or liftOver-based pairwise orthology associations...

The file format is as follows (tsv; header is expected):  GeneID_Sp1 Exon_Coord_Sp1 GeneID_Sp2 Exon_Coord_Sp2

```bash
ENSG00000171055	chr2:36552056-36552268:-	ENSMUSG00000056121	chr17:78377717-78377890:-
ENSG00000171055	chr2:36578597-36578832:-	ENSMUSG00000056121	chr17:78400630-78400865:-
ENSG00000171055	chr2:36558438-36558513:-	ENSMUSG00000056121	chr17:78384744-78384819:-
```

ExOrthist includes a custom [script](https://github.com/biocorecrg/ExOrthist/blob/master/bin/GetLiftOverFile.pl) to create lists of exon files using liftOver.

```bash
perl GetLiftOverFile.pl -annot_sp1 test/Hsa.exons -annot_sp2 test/Mmu.exons -gene_clusters XXXX -chain_file hg38ToMm10.over.chain
```

Basic cluster stats are provided with Get_stats_exon_cls.pl => DEVELOP FURTHER AND RUN AUTOMATICALLY AT THE END


Basic cluster statistics
------------

XXXXXX



Pairwise species re-clustering
------------

XXXX


Evolutionary comparison of exon lists
------------

XXXXXX


Exon-intron gene plots
------------

XXXXXX

