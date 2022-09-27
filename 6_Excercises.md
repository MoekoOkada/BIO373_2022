# Exercises

## 1. Command line

GFF files contains information on features of a sequence: genes, introns, etc. Take a look and familiarize yourself with the format.

<http://www.ensembl.org/info/website/upload/gff.html>

1. How many lines are in athal_genome.gtf?
2. How many characters are in athal_genome.gtf?
3. How many CDS’s are there within scaffold_1 of Ahal.gff? Beware! There’s also scaffold_10, scaffold_11, etc...
4. Extra: How many of each type of feature (column 3) in athal_genome.gtf are there?
5. Extra: zip the athal_genome.gtf file in your own directory and modify the commands to answer some of the question above.

## Answers 1

Will be available on the GitHub: [https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise1_CommandLine_Answers.md](https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise1_CommandLine_Answers.md)


* * *

## 2. FASTA and FASTQ

1. How many sequences are there in the reference sequence file (MedtrChr2.fa)?
2. How many reads does each fastq file have (\*\_R1.fastq.gz)? 
3. Does each sample have the same number of R1 and R2 reads? (Caution: Q scores can be + or @)
4. How many bases are in the reference sequence?
5. How many missing bases (N)? Don’t forget ‘\\n’ is considered a character!

## Answers 2

Will be available on the GitHub: [https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise2_QualityControl_Answers.md](https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise2_QualityControl_Answers.md)

* * *

## 3. Mapping

These exercises are really just designed to try to get you to understand what mapping does with the reads, the almost dizzying amount of information encoded in the files, and what a potential variant might look like in a BAM file. These are fairly detailed.

Bitwise flag meaning: <https://broadinstitute.github.io/picard/explain-flags.html>

1. Find bitwise flags for a few reads in any bam file and decode them using the link above. On the website, you can check only one box at a time to see what an individual property’s value is.
2. Using samtools tview, go to chr2:7271 in 516950.deduped.bam by pressing ‘g’ then typing chr2:7271[Enter]. Compare that with chr2:1018541. Why do you think they are different in terms of coverage and mapping quality? Find a few other areas that look interesting to you and take note of their position.
3. Do these regions look the same in sample 660389?
4. 
   1. How many unique bitwise flags are there in 516950.sorted.bam file?
   2. How many unique bitwise flags are there in the dedupped bam file?
   3. How many reads were marked as duplicates? (Hint; the flags are the sum of the value of each individual property assigned to a read; duplicates = 1024)
5. Extra: Practice writing a bash script to run the alignment and dedup steps on both genotypes.

### Answer 3

Will be available on the GitHub: [https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise3_Mapping_Answers.md](https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise3_Mapping_Answers.md)

## 4. Variant call

1. Here, we'll look in the VCF (04_raw_variants.vcf.gz) and take note of the information contained in the file (which is an overwhelming amount!). I like to get to the variants by searching for CHROM (`/CHROM`). You can look at any SNP, but I suggest searching for 7317, then 1018580. Those sites correspond to where we looked at in the BAM file in the mapping exercises. Take note of the variant quality (QD in INFO field). For an individual, take note of the genotype quality (GQ) and depth (AD and DP) as well. If you'd like, view the BAM file again using `samtools tview` and observe how the results we get from GATK compare to what you can see at those positions in a BAM file. Are the genotypes what you would expect just by looking at the BAM file?  

2. How many variants were discovered in this sample set?

3. Count the number of each genotype (0/0, 0/1, etc) for each sample.

### Answer 4

Will be available on the GitHub: [https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise4and5_VariantCalling_Answers.md](https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise4and5_VariantCalling_Answers.md)

## 5. Filter variants

1. Here, we'll look in the filtered VCF (05_variants_filtered.vcf.gz). This time, see if you notice what changed after the filtration step. For example, the FILTER field should now have a value (not just '.').

2. Count the number of variants that failed each of the filters we applied (those applied to stats in the INFO field).

### Answer 5

Will be available on the GitHub: [https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise4and5_VariantCalling_Answers.md](https://github.com/MoekoOkada/BIO373_2022/blob/main/answers/Exercise4and5_VariantCalling_Answers.md)

