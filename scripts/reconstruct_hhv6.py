#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,subprocess
import log,traceback
import pysam
import utils

def pileup(args, params, filenames, refseqid):
    log.logger.debug('started.')
    try:
        pysam.view(filenames.mapped_to_virus_bam, '-h', '-o', filenames.tmp_bam, refseqid, catch_stdout=False)
        header,seq=utils.retrieve_only_one_virus_fasta(args.vref, refseqid)
        with open(filenames.tmp_fa, 'w') as outfile:
            outfile.write('>%s\n%s\n' % (header, seq))
        pysam.faidx(filenames.tmp_fa)
        cmd='bcftools mpileup -f %s %s | bcftools call -mv -Oz -o %s' % (filenames.tmp_fa, filenames.tmp_bam, filenames.hhv6a_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools norm -f %s %s -Oz -o %s' % (filenames.tmp_fa, filenames.hhv6a_vcf_gz, filenames.hhv6a_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools index %s' % filenames.hhv6a_norm_vcf_gz
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
        cmd='bcftools consensus -f %s -o %s %s' % (filenames.tmp_fa, filenames.hhv6a_pileup_naive, filenames.hhv6a_norm_vcf_gz)
        log.logger.debug('bcftools command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during bcftools running.')
            exit(1)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
