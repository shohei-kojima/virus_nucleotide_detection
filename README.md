# This is a README file for iciHHV6_reconstruction_pipeline

# 1. download this tool from GitHub
git clone https://github.com/shohei-kojima/iciHHV6_reconstruction

# 2. quick usage guide for impatient

### when you use your BAM file as an input (alignmentin option)
```
python main.py \
-alignmentin \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_bwa_index \
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
-vrefindex /path/to/viral_genomic_seq_bwa_index \
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
-vrefindex /path/to/viral_genomic_seq_bwa_index \
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
-vrefindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
You can specify your fastq files as inputs by '-fastqin' option. This tool does NOT support single-end reads. Please specify paired-end reads as separated files ('-fq1' and '-fq2' options).

### when you perform de novo assembly using metaspades (denovo option)
```
python main.py \
-alignmentin \
-denovo \
-b test.bam \
-vref /path/to/viral_genomic_seq.fa \
-vrefindex /path/to/viral_genomic_seq_bwa_index \
-picard /path/to/picard.jar \
-p 4
```
'-denovo' option performs de novo assembly of HHV-6 sequence, in addition to variant call and consensus sequence generation. 


# 3. output files
### 'virus_detection_summary.txt'
This is a main result file for most users. This contains read coverage information of each virus genome.

- 1st column: RefSeq ID (fasta header without attribution)
- 2nd column: Whether a virus-derived reads are abundantly detected or not
- 3rd column: Coverage information
    - genome_length: Length of the virus genome
    - mapped_length: Length of the virus genome with one or more reads
    - perc_genome_mapped: Percent of virus genome with one or more reads
    - average_depth: Average of mapped read depth (average of whole genome)
    - average_depth_of_mapped_region: Average of mapped read depth of mapped regions (average of only mapped regions)
- 4th column: attribution of fasta header

### 'mapped_to_virus_dedup.bam'
BAM file containing alignment with virus genomes.

### 'mapped_to_virus.bedgraph'
Read depth of virus genomes.

### 'mark_duplicate_metrix.txt'
Summary of picard MarkDuplicates running with 'mapped_to_virus_dedup.bam.'

### 'high_coverage_viruses.pdf'
If there is one or more viruses in your sample, this tool outputs read coverage of those viruses.

### 'hhv6a.vcf.gz', 'hhv6b.vcf.gz'
This file contains variant calls. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B.

### 'hhv6a_reconstructed.fa', 'hhv6b_reconstructed.fa'
This file contains HHV-6 sequence reconstructed with called variants. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B.

### 'hhv6a_metaspades_assembly', 'hhv6b_metaspades_assembly'
This file contains HHV-6 sequence reconstructed by metaspades. This tool outputs this file only when your sample contained either HHV-6A or HHV-6B. You need to specify '-denovo' option to obtain this result.

### 'for_debug.log'
Log file. Stores all information, including input arguments and parameter setting.


# 4. other options
### '-p [integer]'
Number of threads. 3 or more is recommended. Default = 2.

### '-bwa'
When you specify this, this tool will use BWA MEM instead of hisat2 when mapping read to virus genomes. You also need to specify bwa index with '-vrefindex' option.

### '-outdir [out_dir_name]'
You can specify arbitrary name for output directory.

### '-overwrite'
In default, this tool does not overwrite directories and files when output directories and files already exist. When you specify this option, this pipeline overwrite existing directories and files.

### '-keep'
You can keep intermediate files when specifying this option.
