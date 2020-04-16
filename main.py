#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,datetime,argparse,glob,logging

'''
time python main.py -overwrite -b test_data/NA12878.chr22.bam -fa /home/kooojiii/Documents/genomes/hg38/hg38.fa -p 4
'''


# version
version='2020/04/08'


# args
parser=argparse.ArgumentParser(description='')
parser.add_argument('-b', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end BAM file.')  # , required=True
parser.add_argument('-c', metavar='str', type=str, help='Either -b or -c is Required. Specify input mapped paired-end CRAM file.')  # , required=True
parser.add_argument('-fa', metavar='str', type=str, help='Required. Specify reference genome which are used when input reads were mapped. Example: GRCh38DH.fa')
parser.add_argument('-vref', metavar='str', type=str, help='Required. Specify reference of virus genomes, including HHV-6A and B. Example: viral_genomic_200405.fa')
parser.add_argument('-ht2index', metavar='str', type=str, help='Required. Specify hisat2 index of virus genomes, including HHV-6A and B. Example: viral_genomic_200405')
parser.add_argument('-outdir', metavar='str', type=str, help='Optional. Specify output directory. Default: ./result_out', default='./result_out')
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
filenames.mapped_to_virus_bam =os.path.join(args.outdir, 'mapped_to_virus_sorted.bam')
filenames.mapped_to_virus_bai =os.path.join(args.outdir, 'mapped_to_virus_sorted.bai')
filenames.bedgraph            =os.path.join(args.outdir, 'mapped_to_virus.bedgraph')


# 0. Unmapped read retrieval
import retrieve_unmapped
log.logger.info('Unmapped read retrieval started.')
retrieve_unmapped.retrieve_unmapped_reads(args, params, filenames)

# 1. mapping
import mapping
log.logger.info('Unmapped read retrieval started.')
mapping.map_to_viruses(args, filenames)
utils.gzip_or_del(args, params, filenames.unmapped_merged_1)
utils.gzip_or_del(args, params, filenames.unmapped_merged_2)
mapping.bam_to_bedgraph(args, params, filenames)

# 2. identify high coverage viruses
import identify_high_cov
identify_high_cov.identify_high_cov_virus_from_bedgraph(args, params, filenames)

# 3. reconstruct HHV-6





# 7. search for absent MEs
if do_abs is True:
    import find_absent
    print()
    log.logger.info('Absent ME search started.')
    find_absent.find_abs(args, params, filenames)
    # gzip files
    utils.gzip_file(params, filenames.abs_txt)
    utils.gzip_or_del(args, params, filenames.breakpoint_clean)
    # output comments
    log.logger.info('%d absent ME candidates found.' % find_absent.abs_n)
    log.logger.info('Absent ME search finished!')


# remove unnecessary files
os.remove(filenames.repout_bed)
os.remove(filenames.reshaped_rep)
for f in glob.glob(filenames.repdb +'*'):
    os.remove(f)
log.logger.info('All analysis finished!')
