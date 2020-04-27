# iciHHV6_reconstruction

# download this tool from GitHub
git clone https://github.com/shohei-kojima/iciHHV6_reconstruction

# quick usage guide for impatient

### when you use your BAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-bwaindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
In this case, you need to specify your BAM file with '-b' option. You also need to specify '-alignmentin' option as well.


### when you use your CRAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-c test.cram \
-fa /path/to/reference/genome/of/cram/GRCh38.fa \
-vref /path/to/viral_genomic_seq.fa \
-bwaindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
In this case, you need to specify your CRAM file with '-c' option. You also need to specify '-alignmentin' option and '-fa' option as well.

### when you use only unmapped reads in BAM/CRAM file for mapping to virus genome (alignmentin option + only_unmapped option)
```
python main.py \
-alignmentin \
-only_unmapped \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-bwaindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
By default, this tool use all discordant reads (e.g. reads without sam flag '2') for mapping to viruses when BAM or CRAM file are specified as an input. Those discordant reads includes read-pairs with distant mapped positions and ones with low MAPQ; therefore, the number of discordant reads are far higher than unmapped reads. If you want to use only unmapped reads to reduce computational burden, you can specify '-only_unmapped' option. This option is only available when '-alignmentin' option was specified. This option is primarily intended for quick screening of HHV-6 positive samples. Please use it at your own risk.

### when you use your fastq files as inputs (fastqin option)
```
python main.py \
-fastqin \
-fq1 test_1.fastq \
-fq2 test_2.fastq \
-vref /path/to/viral_genomic_seq.fa \
-bwaindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
You can specify your fastq files as inputs by '-fastqin' option. This tool does NOT support single-end reads. Please specify paired-end reads as separated files ('-fq1' and '-fq2' options).

# output files


# other options

