# Please add a name for the pipeline

## License
Please choose a license for this project!

## Docker image
Docker image was uploaded to [Docker Hub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/exon_intron_pipe). 


## New flow (MI)

0) To add exons from VTS (creating a GTF with fake transcripts)

     a) IntroduceEXON_toGTF_v3.pl 

```
Usage: IntroduceEXON_toGTF_v3.pl REFERENCE GTF Sp Dir_vastb_files genome_file
```

* deprecates: join_annot_fake_files.pl
* comments: add proper get_options functionality?

     b) cat Fake.gtf >> Sp_annot.gtf


1) generate_annotations.pl

```
Usage: generate_annotations.pl -GTF path_to_gtfs/ -G path_to_genomes/ -sp Sp1,Sp2 [-EX_DB path_to_EXONS_DB/ -vastdb REF1,REF2]

Script that creates all annotation files needed for the second module of the pipeline

COMPULSORY
     -GTF              Path where GTFs are stored (they should be named Sp1_annot.gtf)
     -G                Path where gDNAs are stores (they should be named Sp1_gDNA.fasta)
     -sp Sp1,Sp2       String of species.

OPTIONAL
     -vastdb r1,r2     Comma-separated list of VASTDB reference files (must match species list in -sp)
                           If a species is missing the reference file, "NA" should be provided
     -EX_DB            Path to EXONS_DB/ folder (default ./; i.e. working directory)
                           If it does not exit, it will create a EXONS_DB/ folder in the working directory
     -verbose T/F      Verbose (default TRUE) 
     -h/--help         This help message.

```

     

* deprecates: all other scripts from Module I
* Example to run it on test
perl ../../bin/generate_annotations.pl -GTF ../GTF/ -G ../GENOMES/ -sp Dme,Ame -vastdb ../VASTDB_files/REFERENCE-ALL_ANNOT-Dme66.tab,../VASTDB_files/REFERENCE-ALL_ANNOT-Ame101.tab

* It creates all files, but some may have different names and other for sure can be deleted after running.


**LUCA**
```
 I think that asking for naming conventions of gtf and fasta file is not a good idea... how many genomes do you expect? why not indicating them like a list of files? -GTF Sp1_annot.gtf,Sp2_annot.gtf.. -G Sp1_gDNA.fasta,Sp2_gDNA.fasta etc? 
```
**MANU**
```
It could a possibility. It would require changing the code a bit and it's a bit more tedious (and ugly) to run. For instance, Yamile ran it for 16 species.
```

**LUCA**
```
Yes but this then is easy to wrap with NF and to even launch the 16 species in parallel. Moreover we can allow the use of gzipped input files as an option. 
```


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
