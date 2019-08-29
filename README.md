# Please add a name for the pipeline

## License
Plase choose a license for this project!

## Docker image
Docker image was uploaded to [Docker Hub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/exon_intron_pipe). 


## New flow (MI)

0) To add exons from VTS (creating a GTF with fake transcripts)

     a) IntroduceEXON_toGTF_v3.pl 

`Usage: IntroduceEXON_toGTF_v3.pl REFERENCE GTF Sp Dir_vastb_files genome_file`

* deprecates: join_annot_fake_files.pl
* comments: add proper get_options functionality?

     b) cat Fake.gtf >> Sp_annot.gtf


1) get_exint_file_ext.pl 

`Usage: Run_Module_I.pl -GTF path_to_gtfs/ -G path_to_genomes/`

`OPTIONS
     -GTF          Path where GTFs are stored (they should be named Sp1_annot.gtf)
     -G            Path where gDNAs are stores (they should be named Sp1_gDNA.fasta)
     -sp Sp1,Sp2   String of species.
     -h/--help     This help message.`

**LUCA** I think that asking for naming conventions of gtf and fasta file is not a good idea... how many genomes do you expect? why not indicating them like a list of files? -GTF Sp1_annot.gtf,Sp2_annot.gtf.. -G Sp1_gDNA.fasta,Sp2_gDNA.fasta etc? 

* deprecates: get_ref_proteins.pl, get_exint_file.pl, get_ref_prot.pl, get_ref_prot_exint_file.pl

# Nextflow pipeline
2 modules.

* Module 1
     * input:
     * output:
     * options:
     * scripts that have to be called and the command line for their use.
          * script 1
          * script 2
          * etc...
