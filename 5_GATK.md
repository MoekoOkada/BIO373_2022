# Variant calling using GATK

From here, we use GATK4 software.  
If you logged out of the server, you'll need to load the samtools again.

```bash
$ ssh your_BFabric_account_name@172.23.30.6
$ source /usr/local/ngseq/etc/lmod_profile
$ module load Tools/samtools/1.11
$ module load Variants/GATK/4.2.0.0
```

| :warning: Note that there are huge difference between GATK3 and GATK4. |
| :--------------------------------------------------------------------- |

 GATK4 is the latest version. You cannot use the same commands with GATK3.

* * *

## Prep the reference sequence

GATK needs to use various index and dictionaries of the reference sequence in order to run smoothly.  
You need to use samtools and gatk to generate these.

| :warning: Make sure you're in your VarCall/ directory |
| ----------------------------------------------------- |

```bash
$ samtools faidx 00_input/MedtrChr2.fa
$ gatk CreateSequenceDictionary -R 00_input/MedtrChr2.fa
```

* * *

## Mark duplicates

In order to avoid false-negative/false-positive in variant calling process, you need to remove the reads derived from PCR or sequencing process.

```bash
$ mkdir 02_dedup
$ acc=516950
$ gatk MarkDuplicates -I 01_aligned/${acc}.sorted.bam -O 02_dedup/${acc}.dedup.bam -M 02_dedup/${acc}.metrics    
$ samtools index 02_dedup/${acc}.dedup.bam
```

| :warning: Run for the other sample (660389) ! |
| --------------------------------------------- |

:computer: Forgot which sample number was stored in the acc variable? The echo command will show you!

```bash
$ echo ${acc}
```

* * *

## Find some variants

Now we use several tools in GATK to discover variants, group the samples, and finally filter the dataset to be used in downstream analyses.

First, you need to call variants in each samples. This creates a gVCF file. The step looks at the base that is called in the read that mapped at each site and determines whether there is a variant or not.

```bash
$ mkdir 03_callSNPs
$ acc=516950
$ gatk HaplotypeCaller -R 00_input/MedtrChr2.fa -I 02_dedup/${acc}.dedup.bam -ERC GVCF -O 03_callSNPs/${acc}.g.vcf.gz
```

This takes ~10 minutes for each sample. Maybe add this to your bash script and grab a coffee??

| :warning: Run again for sample 660389 ! |
| --------------------------------------- |

* * *

## Merge GVCFs

Merge GVCFs from multiple samples and create genotype database which will be refered during joint genotyping.

```bash
$ gatk GenomicsDBImport -R 00_input/MedtrChr2.fa -V 03_callSNPs/516950.g.vcf.gz -V 03_callSNPs/660389.g.vcf.gz -L chr2 --genomicsdb-workspace-path Mt_DB 

```

Need to define at least one intervals (where to operate) with `-L` option.

* * *

## Joint genotyping

This steps takes the input of gVCFs from each sample and combines them. This step recalculates some statistics and gives a more confident call of the variant at a particular site.

```bash
$ gatk GenotypeGVCFs -R 00_input/MedtrChr2.fa -V gendb://Mt_DB -O 03_callSNPs/04_raw_variants.vcf.gz
```

### VCF format

All lines beginning with `##` contain information about how and when the VCF was generated and information about the flags included in the file. The single line with `#` tells you what each column of the following lines contains. Every line _without_ a `#` is a variant.

The Genotype field (column 9) is important. Other flags may appear, but these are the minimum that should be included:

| Fields | Mean                                                                                                   |
| ------ | ------------------------------------------------------------------------------------------------------ |
| GT     | genotype. <br>0/0 = homozygous ref; 0/1 = heterozygous; 1/1 = homozygous alt                           |
| AD     | allele depth<br> number of reads containing ref, alt base                                              |
| DP     | total depth<br> total number of reads covering site                                                    |
| GQ     | genotype quality<br> difference between lowest and second lowest PL                                    |
| PL     | genotype likelihood<br> what’s the probability it is NOT the correct genotype. The lowest is always 0. |

:computer: You can view the contents of a VCF using zless. It is a regular text file that has been compressed. To search in the file using `zless`, type `/` when the file is loaded on the screen and type your search. i.e. `/CHROM` will take you to the line which shows the column names and the variants called in each sample.

## Filtering out low quality data

This last step flags the variants that are low confidence. It will put the filterName in the INFO column and the G_filterName in the sample specific column if the site failed the filter.

```bash
$ gatk VariantFiltration -R 00_input/MedtrChr2.fa -O 03_callSNPs/05_variants_filtered.vcf.gz -V 03_callSNPs/04_raw_variants.vcf.gz -filter "QD < 2.0" --filter-name "QD" -filter "MQ < 30.0" --filter-name "MQ" -filter "MQRankSum < -15.0" --filter-name "MQRankSum" -filter "GQ < 20 || DP == 0 " --filter-name "GQ"
```

If at least one sample fails the G_filter, a semicolon-separated list of codes for filters that fail shows up in FORMAT field.  
i.g. chr2:8248 failed the filter named "QD".

![Filtered](https://user-images.githubusercontent.com/58171601/135255461-3cba16e7-2269-4585-9ff9-3eae3068549f.png)
