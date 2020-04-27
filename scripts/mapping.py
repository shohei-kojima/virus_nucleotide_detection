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
        if not args.hisat2 is True:
            cmd='bwa mem -Y -t %d %s %s %s | samtools view -Sbh -o %s -' % (thread_n, args.vrefindex, filenames.unmapped_merged_1, filenames.unmapped_merged_2, filenames.mapped_unsorted_bam)
        else:
            cmd='hisat2 -t -x %s -p %d -1 %s -2 %s --no-spliced-alignment | samtools view -Sbh -o %s -' % (args.vrefindex, thread_n, filenames.unmapped_merged_1, filenames.unmapped_merged_2, filenames.mapped_unsorted_bam)
        log.logger.debug('bwa mem command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bwa mem mapping.')
            exit(1)
        pysam.sort('-@', str(thread_n), '-o', filenames.mapped_to_virus_bam, filenames.mapped_unsorted_bam)
        pysam.index('-@', str(thread_n), filenames.mapped_to_virus_bam)
        if not args.keep is True:
            os.remove(filenames.mapped_unsorted_bam)
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

