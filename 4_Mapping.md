# Mapping reads to reference genome: DNA edition

Here we map reads to a reference genome.  
These correspond to DNA, whereas Masa did RNA mapping/alignment with you last week.  
They are similar in the types of steps to execute, but the underlying concepts are a bit different.

Before you map, you usually want to check the quality of the reads using something like `fastqc` and then you definitely want to trim the reads. **No use in trying to map poor quality data.**  
These steps have already been done, so we start with the mapping  

We will use the module system to set up our environment:

```bash
$ source /usr/local/ngseq/etc/lmod_profile
```

* * *

## Map reads with BWA

| :warning: Make sure you're back in your `VarCall/` directory |
| ------------------------------------------------------------ |

First, set up the environment with the software we will use: `BWA` and `samtools`.  
In order to run the analysis, we need to index the reference sequence

```bash
$ module load Aligner/BWA/0.7.17
$ module load Tools/samtools/1.11
$ bwa index 00_input/MedtrChr2.fa
```

We have two samples for this exercise, called `516950` and `660389`.  
Run the mapping tool BWA and pipe the output through samtools, which will sort the output.  
We can use it in subsequent steps.

```bash
$ mkdir 01_aligned
$ acc=516950
$ bwa mem -M -t 2 -R "@RG\tID:CAV90ANXX.6\tPL:Illumina\tLB:${acc}\tSM:${acc}" 00_input/MedtrChr2.fa 00_input/${acc}_chr2_R{1,2}.fastq.gz | samtools sort -m 16G -T /scratch/bio373_2021/YOUR_USERNAME -o 01_aligned/${acc}.sorted.bam
$ samtools index 01_aligned/${acc}.sorted.bam
```

Then, change the value in ${acc} variable above for the other sample (660389) and re-run the whole bwa|samtools command.

:computer: Pro tip: if you press the up or down arrow keys, you can scroll through the commands you have already run during your shell session.

You should have 4 files in 01_aligned now: 2 `.bam` files and 2 `.bam.bai` files.

* * *

### BAM files and relative path practice

After mapping, you get a BAM file of where in the reference sequence the read mapped to. Because BAM files are compressed in a special way (not gzipped!) you can only open the BAM with certain tools, here `samtools`.

Below, I've given the basic commands but also a bit of thinking about where you are in the directory structure and how you can modify the input of the command so the computer knows where to look.

#### Trouble shooting when you get errors

- Check where you are: `pwd`
- Move the correct directory: `cd`
- If you don't want to move, Modify the path to the file.

Look at the header of the BAM file:

```bash
$ samtools view -H 516950.sorted.bam 
```

Now look at the actual mapped reads:

```bash
$ samtools view 516950.sorted.bam | less
```

Masa used IGV to view the alignments. You can also do that (remember you have to `scp` the BAM file onto the local computer!) Here's the command line alternative:

```bash
$ samtools tview 516950.sorted.bam MedtrChr2.fa
```

(**Think of where the reference sequence is and where the bam file is RELATIVE to one another**)

While in the samtools tview, use the `?` to open the menu and `q` to exit. Once tview is open, you can type `g` and then within the prompt, type a location you'd like to go to. For example chr2:7271 then press `Enter`. `Esc` will get you out of the location search prompt.
