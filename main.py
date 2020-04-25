#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
python main.py -overwrite -c test_data/NA18999.final.cram -fa /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa -vref /home/kooojiii/Documents/NCBI_database/all_viral_genome_nt/200406_1/viral_genomic_200405.fa -bwaindex /home/kooojiii/Documents/NCBI_database/all_viral_genome_nt/200406_1/bwa_index/viral_genomic_200405 -keep -p 4

python main.py -overwrite -unmappedin -unmap1 ./test_data/unmapped_merged_1.fq -unmap2 ./test_data/unmapped_merged_2.fq -fa /home/kooojiii/Documents/genomes/hg38/1kGP/GRCh38_full_analysis_set_plus_decoy_hla.fa -vref /home/kooojiii/Documents/NCBI_database/all_viral_genome_nt/200406_1/viral_genomic_200405.fa -bwaindex /home/kooojiii/Documents/NCBI_database/all_viral_genome_nt/200406_1/bwa_index/viral_genomic_200405 -picard /home/kooojiii/bin/picard.jar -keep -p 4
'''


# version
version='2020/04/24'


# HHV-6 refseq IDs
hhv6a_refid='NC_001664.4'
hhv6b_refid='NC_000898.1'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-alignmentin', help='Optional. Specify if you use BAM/CRAM file for input. You also need to specify either -b or -c.', action='store_true')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')  # , required=True
parser.add_argument('-unmappedin', help='Optional. Specify if you use unmapped reads for input instead of BAM/CRAM file. You also need to specify -unmap1 and -unmap2.', action='store_true')
parser.add_argument('-unmap1', metavar='str', type=str, help='Specify unmapped fastq file, read-1 of read pairs.')
parser.add_argument('-unmap2', metavar='str', type=str, help='Specify unmapped fastq file, read-2 of read pairs.')
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-vref', metavar='str', type=str, help='Required. Specify reference of virus genomes, including HHV-6A and B. Example: viral_genomic_200405.fa')
parser.add_argument('-bwaindex', metavar='str', type=str, help='Required. Specify hisat2 index of virus genomes, including HHV-6A and B. Example: viral_genomic_200405')
parser.add_argument('-picard', metavar='str', type=str, help='Required. Specify full path to picard.jar. Example: /path/to/picard/picard.jar')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
parser.add_argument('-overwrite', help='Optional. Specify if you overwrite previous results.', action='store_true')
parser.add_argument('-keep', help='Optional. Specify if you do not want to delete temporary files.', action='store_true')
parser.add_argument('-p', metavar='int', type=int, help='Optional. Number of threads. 4 is recommended. Default: 1', default=1)
args=parser.parse_args()


# start
import init
init.init(args)


# logging
import log
args.logfilename='for_debug.log'
if os.path.exists(os.path.join(args.outdir, args.logfilename)) is True:
    os.remove(os.path.join(args.outdir, args.logfilename))
log.start_log(args)
log.logger.debug('Logging started.')


# initial check
import initial_check
log.logger.debug('This is version %s' % version)
print()
log.logger.info('Initial check started.')
initial_check.check(args, sys.argv)


# set up
import setup
setup.setup(args, init.base)
params=setup.params


# output file names
import utils
filenames=utils.empclass()

filenames.unmapped_1          =os.path.join(args.outdir, 'unmapped_1.fq')
filenames.unmapped_2          =os.path.join(args.outdir, 'unmapped_2.fq')
filenames.unmapped_3          =os.path.join(args.outdir, 'unmapped_3.fq')
filenames.unmapped_4          =os.path.join(args.outdir, 'unmapped_4.fq')
filenames.unmapped_5          =os.path.join(args.outdir, 'unmapped_5.fq')
filenames.unmapped_6          =os.path.join(args.outdir, 'unmapped_6.fq')
filenames.unmapped_merged_1   =os.path.join(args.outdir, 'unmapped_merged_1.fq')
filenames.unmapped_merged_2   =os.path.join(args.outdir, 'unmapped_merged_2.fq')
filenames.mapped_unsorted_bam =os.path.join(args.outdir, 'mapped_to_virus.bam')
filenames.mapped_to_virus_bam =os.path.join(args.outdir, 'mapped_to_virus_sorted.bam')
filenames.mapped_to_virus_bai =os.path.join(args.outdir, 'mapped_to_virus_sorted.bai')
filenames.bedgraph            =os.path.join(args.outdir, 'mapped_to_virus.bedgraph')
filenames.high_cov_pdf        =os.path.join(args.outdir, 'high_coverage_viruses.pdf')

filenames.tmp_bam             =os.path.join(args.outdir, 'tmp.bam')
filenames.tmp_rg_bam          =os.path.join(args.outdir, 'tmp_rg.bam')
filenames.tmp_fa              =os.path.join(args.outdir, 'tmp.fa')
filenames.tmp_masked_fa       =os.path.join(args.outdir, 'tmp_masked.fa')
filenames.tmp_fa_dict         =os.path.join(args.outdir, 'tmp.dict')
filenames.hhv6a_vcf_gz        =os.path.join(args.outdir, 'hhv6a.vcf.gz')
filenames.hhv6a_norm_vcf_gz   =os.path.join(args.outdir, 'hhv6a_norm.vcf.gz')
filenames.hhv6a_gatk_naive    =os.path.join(args.outdir, 'hhv6a_reconstructed.fa')


# 0. Unmapped read retrieval
if args.alignmentin is True:
    import retrieve_unmapped
    log.logger.info('Unmapped read retrieval started.')
    retrieve_unmapped.retrieve_unmapped_reads(args, params, filenames)
elif args.unmappedin is True:
    log.logger.info('Unmapped read retrieval skipped. Read1=%s, read2=%s.' % (args.unmap1, args.unmap2))
    filenames.unmapped_merged_1=args.unmap1
    filenames.unmapped_merged_2=args.unmap2

# 1. mapping
import mapping
log.logger.info('Mapping of unmapped reads started.')
mapping.map_to_viruses(args, filenames)
if args.alignmentin is True:
    utils.gzip_or_del(args, params, filenames.unmapped_merged_1)
    utils.gzip_or_del(args, params, filenames.unmapped_merged_2)
log.logger.info('BAM to bedgraph conversion started.')
#mapping.bam_to_bedgraph(args, params, filenames)

# 2. identify high coverage viruses
import identify_high_cov
log.logger.info('Identification of high-coverage viruses started.')
identify_high_cov.identify_high_cov_virus_from_bedgraph(args, params, filenames)

# 3. reconstruct HHV-6
import reconstruct_hhv6
if identify_high_cov.hhv6a_highcov is True:
    log.logger.info('HHV-6A sequence reconstruction started.')
    reconstruct_hhv6.pileup(args, params, filenames, hhv6a_refid)

