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
        # retrieve discordant reads
        if not args.b is None:
            pysam.view('-@', '%d' % thread_n, '-f', '1', '-F', '3842', '-b', '-o', filenames.discordant_bam, args.b, catch_stdout=False)
        elif not args.c is None:
            pysam.view('-@', '%d' % thread_n, '-f', '1', '-F', '3842', '-b', '-o', filenames.discordant_bam, '--reference', args.fa, args.c, catch_stdout=False)
        pysam.sort('-@', '%d' % thread_n, '-n', '-O', 'BAM', '-o', filenames.discordant_sort_bam, filenames.discordant_bam)
        if args.keep is False:
            os.remove(filenames.discordant_bam)
        if args.all_discordant is True:
            pysam.fastq('-@', '%d' % thread_n, '-N', '-0', '/dev/null', '-1', filenames.unmapped_merged_1, '-2', filenames.unmapped_merged_2, '-s', '/dev/null', filenames.discordant_sort_bam)
        else:
            pysam.fastq('-@', '%d' % thread_n, '-f', '12', '-F', '3328', '-N', '-0', '/dev/null', '-1', filenames.unmapped_1, '-2', filenames.unmapped_2, '-s', '/dev/null', filenames.discordant_sort_bam)
            pysam.fastq('-@', '%d' % thread_n, '-f', '8', '-F', '3332', '-N', '-0', '/dev/null', '-1', filenames.unmapped_3, '-2', filenames.unmapped_4, '-s', '/dev/null', filenames.discordant_sort_bam)
            pysam.fastq('-@', '%d' % thread_n, '-f', '4', '-F', '3336', '-N', '-0', '/dev/null', '-1', filenames.unmapped_5, '-2', filenames.unmapped_6, '-s', '/dev/null', filenames.discordant_sort_bam)
            # concatenate fastq
            with open(filenames.unmapped_merged_1, 'w') as outfile:
                for f in [filenames.unmapped_1, filenames.unmapped_3, filenames.unmapped_5]:
                    with open(f) as infile:
                        for line in infile:
                            outfile.write(line)
                    utils.gzip_or_del(args, params, f)
            with open(filenames.unmapped_merged_2, 'w') as outfile:
                for f in [filenames.unmapped_2, filenames.unmapped_4, filenames.unmapped_6]:
                    with open(f) as infile:
                        for line in infile:
                            outfile.write(line)
                    utils.gzip_or_del(args, params, f)
        if args.keep is False:
            os.remove(filenames.discordant_sort_bam)

    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)

