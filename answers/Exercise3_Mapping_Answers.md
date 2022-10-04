# Answers of Excercises: Mapping

## Question 1: 

> Find bitwise flags for a few reads in any bam file and decode them using the link above. On the website, you can check only one box at a time to see what an individual property’s value is.

While you can use `grep` to this task, it would have to read the whole file before returning the results. A more efficient way in this case is to pipe the output from `samtools view` to `less`. The bit flag is in the second column.

```bash
$ samtools view 01_aligned/516950.sorted.bam | less
```

* * *

## Question 2: 

> Using samtools tview, go to chr2:7271 in 516950.sorted.bam by pressing ‘g’ then typing chr2:7271[Enter]. Compare that with chr2:1018541. Why do you think they are different in terms of coverage and mapping quality? Find a few other areas that look interesting to you and take note of their position if you'd like to go back after the SNP calling step and see if a SNP was indeed called.

To open a file, you need specify both the bam file and the reference.

```bash
$ samtools tview 01_aligned/516950.sorted.bam 00_input/MedtrChr2.fa
```

The mapping quality in chr2:1018541 is lower than chr2:7271. There are some possibility.

- Genes can be duplicated in the genome. The second gene is a gene that has multiple copies in the genome, and it is possible that reads that should have been mapped elsewhere were mapped to this region.
- The gene locates on chr2:1018541 is putative cyanidin 3-O-galactoside 2''-O-xylosyltransferase FGGT1. This is related to anthocyanin production. There might be a lot of copies of this genes with variants in the _Medicago trancatula_ genome.
- The first region (chr2:7271) codes putative double-strand break repair protein MRE11. This might be single copy gene in the _Medicago trancatula_ genome. This is important gene. If gene duplication does occur, mutations are unlikely to accumulate.

* * *

## Question 3

> Do these regions look the same in sample 660389?

Yes, both regions look quite similar in 660389.

* * *

## Question 4

> 1. How many unique bitwise flags are there in 516950.sorted.bam file? 

You can combine samtools and unix commands to count unique flags in bam file.

```bash
# 516950.sorted.bam
$ samtools view 01_aligned/516950.sorted.bam | cut -f 2 | sort -u | wc -l
40
$
```

> 2. How many unique bitwise flags are there in the dedupped bam file? 

```bash
# 516950.dedup.bam
$ samtools view 02_dedup/516950.dedup.bam | cut -f 2 | sort -u | wc -l
54
$
```

> 3. How many reads were marked as duplicates?

```bash
$ samtools view -f 1024 02_dedup/516950.dedup.bam | wc -l
2654
$
```

* * *
* * *

## Question 5 (Extra): Practice writing a bash script to run the alignment and dedup steps on both genotypes. 

To encourage organization and reproducibility, make a directory to keep your script(s) in and run from there :)  

Example of bash script

```bash
#!/bin/sh
# HaplotypeCaller.sh
# Project: BIO373
# Date: 28-Sep-2022
# Author: Moeko Okada

## SOFTWARES
source /usr/local/ngseq/etc/lmod_profile
module load Aligner/BWA/0.7.17
module load Tools/samtools/1.11
module load Variants/GATK/4.2.0.0

bio373=/sratch/bio373_2022
workdir=${bio373}/moeko

cd ${workdir}
pwd

# indexing
bwa index 00_input/MedtrChr2.fa
samtools faidx 00_input/MedtrChr2.fa
gatk CreateSequenceDictionary -R 00_input.MedtrChr2.fa

# mapping
fMAP=0
if [ ${fMAP} -eq 1 ]; then
    for acc in "516950" "660389"; do
        bwa mem -M -t 2 -R "@RG\tID:CAV90ANXX.6\tPL:Illumina\tLB:${acc}\tSM:${acc}" \
            00_input/MedtrChr2.fa 00_input/${acc}_chr2_R{1,2}.fastq.gz | samtools sort -m 16G \
            -T /scratch/bio373_2021/moeko -o 01_aligned/${acc}.sorted.bam
    done
fi

# deduplication
fDEDUP=0
if [ ${fDEDUP} -eq 1 ]; then
    for acc in "516950" "660389"; do
        gatk MarkDuplicates -I 01_aligned/${acc}.sorted.bam \
            -O 02_dedup/${acc}.dedup.bam \
            -M 02_dedup/${acc}.metrics
        samtools index 02_dedup/${acc}.dedup.bam
    done
fi
```

#### `for` loop

We frequently use the `for` loop to do the same step for multiple files.  

#### Create switch with `flag` and `if` statement

We can set `flags`, which is switch in shell scripting, combining `flags` and `if` statement.  
Set number of `flags` by writing, for example, `fMAP=0` (flag of mapping equal 0).  
Then when the number in `if [ ${fDEDUP} -eq 1 ]; then` is the same with number you set, the script between the `if` statement (from `then` to `fi`) runs.  
So you can set switch just writing these three lines, and ON/OFF the switch by changing the number either in `flag` or in `if` statement.