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
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        pysam.view(filenames.mapped_to_virus_bam, '-h', '-o', filenames.tmp_bam, refseqid, catch_stdout=False)
        header,seq=utils.retrieve_only_one_virus_fasta(args.vref, refseqid)
        with open(filenames.tmp_fa, 'w') as outfile:
            outfile.write('>%s\n%s\n' % (header, seq))
#        cmd='bcftools mpileup -q -10 -Q -10 -f %s %s | bcftools call -mv -Oz -o %s' % (filenames.tmp_fa, filenames.tmp_bam, filenames.hhv6a_vcf_gz)
        if os.path.exists(filenames.tmp_fa_dict) is True:
            os.remove(filenames.tmp_fa_dict)
        cmd='java -jar %s CreateSequenceDictionary R=%s O=%s' % (args.picard, filenames.tmp_fa, filenames.tmp_fa_dict)
        log.logger.debug('gatk command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        cmd='java -jar %s AddOrReplaceReadGroups I=%s O=%s RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20' % (args.picard, filenames.tmp_bam, filenames.tmp_rg_bam)
        log.logger.debug('gatk command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
            exit(1)
        pysam.index('-@', str(thread_n), filenames.tmp_rg_bam)
        cmd='gatk --java-options "-Xmx4g" HaplotypeCaller -R %s -I %s -O %s' % (filenames.tmp_fa, filenames.tmp_rg_bam, filenames.hhv6a_vcf_gz)
        log.logger.debug('gatk command = `'+ cmd +'`')
        out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
        if not out.returncode == 0:
            log.logger.error('Error occurred during gatk running.')
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
