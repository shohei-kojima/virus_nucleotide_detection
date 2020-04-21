#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,subprocess,pysam
import utils
import log,traceback

def map_to_viruses(args, filenames):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        cmd='hisat2 --mp 2,1 -t -x %s -p %d -1 %s -2 %s --no-spliced-alignment | samtools sort -O BAM - > %s' % (args.ht2index, thread_n, filenames.unmapped_merged_1, filenames.unmapped_merged_2, filenames.mapped_to_virus_bam)
        log.logger.debug('hisat2 command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during hisat2 mapping.')
            exit(1)
        pysam.index('-@', str(thread_n), filenames.mapped_to_virus_bam)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def bam_to_bedgraph(args, params, filenames):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        cmd='bamCoverage --outFileFormat bedgraph -p %d --binSize %d -b %s -o %s' % (thread_n, params.bedgraph_bin, filenames.mapped_to_virus_bam, filenames.bedgraph)
        log.logger.debug('bamCoverage command = "'+ cmd +'"')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bamCoverage.')
            exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

