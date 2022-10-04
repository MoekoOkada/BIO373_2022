# Answers of Excercises: Command line

## Question 1: Create symbolic link of this GFF file to your directory.

```bash
$ cd /scratch/bio373_2022/YOUR_DIRECTORY
$ ln -s /scratch/bio373_2022/data/SNPcalling/reference/Ahal.gff
$ ls
Ahal.gff
$
```

* * *

## Question 2: How many lines are in Ahal.gff? Characters?

```bash
# lines
$ less Ahal.gff | wc -l
712704
$
# characters
$ less Ahal.gff | wc -m
57913032
$
```

* * *

## Question 3: How many CDS’s are there within scaffold_1? Beware! There’s also scaffold_10, scaffold_11, etc

You need to search only scaffold_1, not scaffold_10, scaffold_11, etc. By searching with the tab character, only scaffold_1 can be detected. In general, the tab character is `\t`, but if you want to use `grep` to detect the tab character, write `$'\t'` with `-e` option as describe below.

```bash
$ less Ahal.gff | grep -e "scaffold_1"$'\t'| grep -c "CDS"
6085
$
```

- `-e`: option to use regular expression.

You can use `cat` command instead of `less`.  
You can skip `less` command. `grep Ahal.gff -e "scaffold_1"$'\t' | grep -c "CDS"`.

* * *

## Question 4 (Extra): How many of each type of feature (column 3) are there?

```bash
$ cut -f 3 Ahal.gff | sort -u 
##Ahal.gff
##Split from Akam_synth.gff 2019-07-30
##gff-version 3
CDS
exon
gene
intron
mRNA
start_codon
stop_codon
transcription_end_site
transcription_start_site
$ cut -f 3 Ahal.gff | grep -c "CDS"
177558
$ cut -f 3 Ahal.gff | grep -c "exon"
187838
$ cut -f 3 Ahal.gff | grep -c "gene"
32553
$ cut -f 3 Ahal.gff | grep -c "intron"
143350
$ cut -f 3 Ahal.gff | grep -c "mRNA"
34553
$ cut -f 3 Ahal.gff | grep -c "start_codon"
34273
$ cut -f 3 Ahal.gff | grep -c "stop_codon"
34260
$ cut -f 3 Ahal.gff | grep -c "transcription_end_site"
34192
$ cut -f 3 Ahal.gff | grep -c "transcription_start_site"
34124
$
```

Of course, you can use `wc` command with pipe.

```bash
$ cut -f 3 Ahal.gff | grep "CDS" | wc -l
177558
$ 
```

* * *


## Question 5 (Extra): ipZ the file and modify your commands to answer the same questions!

### 1. zip the file

   Normally, you can compress file with `gzip` command. This command compresses file and overwrite it with similar name.  

```bash
$ ls
file.txt
$ gzip file.txt
$ ls
file.txt.gz
```

   When you use `gzip` command to symlink without any options, `gzip` try to replace the original file itself with the compressed one, so you end up with a looping symbolic link reference.  
   To get the compressed file, we first need to display the compression result on the standard output with the `--stdout` command, and then redirect it with new name.

```bash
$ gzip Ahal.gff --stdout > Ahal.gff.gz
$ ls
Ahal.gff  Ahal.gff.gz
$
```

* * *

### 2. How many lines are in Ahal.gff? Characters?

```bash
# lines
$ zless Ahal.gff.gz | wc -l
712704
$ zless Ahal.gff.gz | wc -m
57913032
$
```

You can use `zcat` command instead of `zless`.

* * *

### 3. How many CDS’s are there within scaffold_1? Beware! There’s also scaffold_10, scaffold_11, etc

```bash
$ zless Ahal.gff | grep -E "scaffold_1"$'\t'| grep -c "CDS"
6085
$
```

You can use `zcat` command instead of `zless`.

* * *

### 4. How many of each type of feature (column 3) are there?

```bash
$ zless Ahal.gff.gz | cut -f 3 | sort -u 
##Ahal.gff
##Split from Akam_synth.gff 2019-07-30
##gff-version 3
CDS
exon
gene
intron
mRNA
start_codon
stop_codon
transcription_end_site
transcription_start_site
$ zless Ahal.gff.gz | cut -f 3 | grep -c "CDS"
177558
$ zless Ahal.gff.gz | cut -f 3 | grep -c "exon"
187838
$ zless Ahal.gff.gz | cut -f 3 | grep -c "gene"
32553
$ zless Ahal.gff.gz | cut -f 3 | grep -c "intron"
143350
$ zless Ahal.gff.gz | cut -f 3 | grep -c "mRNA"
34553
$ zless Ahal.gff.gz | cut -f 3 | grep -c "start_codon"
34273
$ zless Ahal.gff.gz | cut -f 3 | grep -c "stop_codon"
34260
$ zless Ahal.gff.gz | cut -f 3 | grep -c "transcription_end_site"
34192
$ zless Ahal.gff.gz | cut -f 3 | grep -c "transcription_start_site"
34124
$
```

Of course, you can use `wc` command with pipe.

```bash
$ zless Ahal.gff.gz | cut -f 3 | grep "CDS" | wc -l
177558
$ 
```

You can use `zcat` command instead of `zless`.

* * *
