# Setting up for Variant calling with GATK

## Environment managememt with module system

Here, we'll again use the module system to load the software we want to use.  
You need to do this every time you log into the server.

```bash
$ source /usr/local/ngseq/etc/lmod_profile
```

To check which modules (software) are available, type the following.  
Use arrow keys to navigate and type `q` to exit.

```bash
$ module avail
```

## It's all relative

We will set up the directories in a specific way during the exercise and I have the commands set up to run based on this directory structure.  
I use relative paths a lot and most things will be run **assuming you are in the VarCall/ directory in YOUR directory on the server.**

The most common reason things don't work are:

- Trying to run things from the wrong directory
- Unpaired quotation marks or brackets
- Misspellings

Make sure you are in the directory you think you are in.  
`pwd`: tell you which directory you currently are in.

## Make project directory

Start from your directory on the server for the course and make a new folder for this exercise:

```bash
$ ssh your_BFabric_account_name@172.23.30.6
$ cd /scratch/bio373_2021/YOUR_USERNAME
$ mkdir VarCall
```

## Set up input by creating symlinks to input files

We use symlinks to help save diskspace on the server.

```bash
$ cd VarCall
$ mkdir 00_input
$ cd 00_input
$ dataDir="/scratch/bio373_2021/data/SNPcalling/inputs"
$ ln -s ${dataDir}/MedtrChr2.fa
$ ln -s ${dataDir}/516950_chr2_R1.fastq.gz
$ ln -s ${dataDir}/516950_chr2_R2.fastq.gz
$ ln -s ${dataDir}/660389_chr2_R1.fastq.gz
$ ln -s ${dataDir}/660389_chr2_R2.fastq.gz
```

### File types

#### FASTA

This is a file with nucleotide sequence(s) in it. For this exercise, it contains the reference sequence from part of chromosome 2 in _M. truncatula_ to which we will map the reads.

The line beginning with `>` is the information about the sequence that follows on the next line.

```bash
$ less MedtrChr2.fa
```

More info: <https://en.wikipedia.org/wiki/FASTA_format>

#### FASTQ

This is the output file given to you after sequencing and contain the reads.

```bash
$ zless 516950_chr2_R1.fastq.gz
```

More info: <https://en.wikipedia.org/wiki/FASTQ_format>

For the quality score encoding: <https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm>
