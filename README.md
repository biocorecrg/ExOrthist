# Please add a name for the pipeline

## License
Please choose a license for this project!

## Docker image
Docker image was uploaded to [Docker Hub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/exon_intron_pipe). 


## NextFlow pipeline
Test data are in the folder **test**. Input files are:
* Genomic sequences in fasta files (can be gzipped). The prefix must match the corresponding annotation file
* Annotation files in GTF format (can be gzipped). The prefix must match the corresponding GTF file
Example:

```bash
     Ame.fasta.gz Ame.gtf.gz
     Dme.fasta.gz Dme.gtf.gz 
```
