# Answers of Excercises: Variant call

## Question 1

> Here, we'll look in the VCF (04_raw_variants.vcf.gz) and take note of the information contained in the file (which is an overwhelming amount!). I like to get to the variants by searching for CHROM (`/CHROM`). You can look at any SNP, but I suggest searching for 7317, then 1018580. Those sites correspond to where we looked at in the BAM file in the mapping exercises. Take note of the variant quality (QD in INFO field). For an individual, take note of the genotype quality (GQ) and depth (AD and DP) as well. If you'd like, view the BAM file again using `samtools tview` and observe how the results we get from GATK compare to what you can see at those positions in a BAM file. Are the genotypes what you would expect just by looking at the BAM file?  

```bash
$ zless 03_callSNPs/04_raw_variants.vcf.gz
```

Then type `/7317`.

```text
chr2    7317    .       C       T       1530.13 .       AC=4;AF=1.00;AN=4;DP=40;ExcessHet=3.0103;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=60.00;QD=31.04;SOR=1.121        GT:AD:DP:GQ:PGT:PID:PL:PS       1/1:0,23:23:69:.:.:976,69,0     1|1:0,14:14:42:1|1:7309_A_T:570,42,0:7309
```

There is a C/T SNP just like when you look at bam in tview.

Type `/1018580`.

```text
chr2    1018580 .       T       A       611.50  .       AC=2;AF=0.500;AN=4;BaseQRankSum=0.862;DP=48;ExcessHet=4.7712;FS=11.379;MLEAC=2;MLEAF=0.500;MQ=54.22;MQRankSum=-3.380e+00;QD=13.29;ReadPosRankSum=-5.540e-01;SOR=2.165   GT:AD:DP:GQ:PGT:PID:PL:PS       0|1:11,6:17:99:0|1:1018580_T_A:219,0,444:1018580        0/1:18,11:29:99:.:.:402,0,720
```

There is T/A SNP just like when you look at bam in tview.

* * *

## Question 2: How many variants were discovered in this sample set?

```bash
$ zgrep -c -v "^#" 03_callSNPs/04_raw_variants.vcf.gz
40194
$
```

`-v`: Show non-matching lines

* * *

## Question 3: Count the number of each genotype (0/0, 0/1, etc) for each sample.

The genotype of sample 516950 is in the 10th column.

```bash
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "0/0"
12472
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "0/1"
3104
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "1/1"
8669
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "0|0"
150
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "0|1"
3901
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 10 | grep -c "1|1"
10444
$
```

The genotype of sample 660389 is in the 11th column.

```bash
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "0/0"
11844
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "0/1"
3275
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "1/1"
7823
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "0|0"
75
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "0|1"
4949
$ zless 03_callSNPs/04_raw_variants.vcf.gz | cut -f 11 | grep -c "1|1"
9978
$
```

Another option is to use `sort` and `uniq -c` commands.

```bash
$ zgrep -v "^#" 05_variants_filtered.vcf.gz | cut -f 10 | cut -f 1 -d “:” | sort | uniq -c 
   1452 ./.
     17 .|.
  12472 0/0
   3104 0/1
     10 0/2
    150 0|0
   3672 0|1
      5 0|2
   8669 1/1
    134 1/2
  10442 1|1
     63 1|2
      2 2/2
      2 2|2
$
```

* * *

# Exercise 5

## Question 1:

> Here, we'll look in the filtered VCF (05_variants_filtered.vcf.gz). This time, see if you notice what changed after the filtration step. For example, the FILTER field should now have a value (not just '.'). 

If you search for the name of filter you defined (i.g. `MQ`, `QD`, etc.), you'll see which position failed which filters. GATK4 conduct joint genotyping, so if the position failed some filters, variants can be called.

* * *

## Question 2:

> Count the number of variants that failed each of the filters we applied (those applied to stats in the FILTER field).

```bash
$ zless 03_callSNPs/05_variants_filtered.vcf.gz | cut -f 7 | grep -c -v "PASS"
1301
$
```
