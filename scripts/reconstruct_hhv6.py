#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os
import log,traceback
import pysam
import utils

def pileup(args, params, filenames, refseqid):
    pysam.view(b=filenames.mapped_to_virus_bam, o=filenames.tmp_bam, refseqid)
    header,seq=utils.retrieve_only_one_virus_fasta(args.vref, refseqid)
    with open(filenames.tmp_fa, 'w') as outfile:
        outfile.write('>%s\n%s\n' % (header, seq))
    pysam.faidx(filenames.tmp_fa)
    pysam.mpileup(uf=filenames.tmp_fa, filenames.tmp_bam, o=filenames.hhv6a_pileup_naive)
