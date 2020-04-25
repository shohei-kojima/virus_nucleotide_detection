#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,sys,pysam
import utils
import log,traceback


def retrieve_unmapped_reads(args, params, filenames):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        # read unmapped, mate unmapped
        if not args.b is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '12', '-F', '3328', '-N', '-0', '/dev/null', '-1', filenames.unmapped_1, '-2', filenames.unmapped_2, '-s', '/dev/null', args.b)
        elif not args.c is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '12', '-F', '3328', '-N', '-0', '/dev/null', '-1', filenames.unmapped_1, '-2', filenames.unmapped_2, '-s', '/dev/null', '--reference', args.fa, args.c)
        # read mapped, mate unmapped
        if not args.b is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '8', '-F', '3332', '-N', '-0', '/dev/null', '-1', filenames.unmapped_3, '-2', filenames.unmapped_4, '-s', '/dev/null', args.b)
        elif not args.c is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '8', '-F', '3332', '-N', '-0', '/dev/null', '-1', filenames.unmapped_3, '-2', filenames.unmapped_4, '-s', '/dev/null', '--reference', args.fa, args.c)
        # read unmapped, mate mapped
        if not args.b is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '4', '-F', '3336', '-N', '-0', '/dev/null', '-1', filenames.unmapped_5, '-2', filenames.unmapped_6, '-s', '/dev/null', args.b)
        elif not args.c is None:
            pysam.fastq('-@', '%d' % thread_n, '-f', '4', '-F', '3336', '-N', '-0', '/dev/null', '-1', filenames.unmapped_5, '-2', filenames.unmapped_6, '-s', '/dev/null', '--reference', args.fa, args.c)
        # both mapped, low MAPQ
        if not args.b is None:
            pysam.fastq('-@', '%d' % thread_n, '-U', '-q', '%d' % params.sam_mapq_threshold, '-F', '3340', '-N', '-0', '/dev/null', '-1', filenames.unmapped_7, '-2', filenames.unmapped_8, '-s', '/dev/null', args.b)
        elif not args.c is None:
            pysam.fastq('-@', '%d' % thread_n, '-U', '-q', '%d' % params.sam_mapq_threshold, '-F', '3340', '-N', '-0', '/dev/null', '-1', filenames.unmapped_7, '-2', filenames.unmapped_8, '-s', '/dev/null', '--reference', args.fa, args.c)
        # concatenate fastq
        with open(filenames.unmapped_merged_1, 'w') as outfile:
            for f in [filenames.unmapped_1, filenames.unmapped_3, filenames.unmapped_5, filenames.unmapped_7]:
                with open(f) as infile:
                    for line in infile:
                        outfile.write(line)
                utils.gzip_or_del(args, params, f)
        with open(filenames.unmapped_merged_2, 'w') as outfile:
            for f in [filenames.unmapped_2, filenames.unmapped_4, filenames.unmapped_6, filenames.unmapped_8]:
                with open(f) as infile:
                    for line in infile:
                        outfile.write(line)
                utils.gzip_or_del(args, params, f)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

