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


def mask_low_depth(args, params, filenames, orig_seq_file, refseqid):
    log.logger.debug('started.')
    try:
        depth=[]
        with open(filenames.bedgraph) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] == refseqid:
                    val=int(ls[3])
                    for _ in range(int(ls[1]), int(ls[2])):
                        depth.append(val)
        orig_seq=[]
        with open(orig_seq_file) as infile:
            for line in infile:
                if '>' in line:
                    header=line.strip()
                    header += ' masked\n'
                else:
                    orig_seq.append(line.strip())
        orig_seq=''.join(orig_seq)
        if not len(depth) == len(orig_seq):
            log.logger.error('Error occurred during making low depth seq. Length of refseq and depth info is not equal.')
            exit(1)
        masked_seq=[]
        for seq,val in zip(orig_seq, depth):
            if val >= params.reconst_minimum_depth:
                masked_seq.append(seq)
            else:
                masked_seq.append('N')
        masked_seq.append('\n')
        with open(filenames.tmp_masked_fa, 'w') as outfile:
            outfile.write(header)
            outfile.write(''.join(masked_seq))
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)



def reconst_any(args, params, filenames, refseqids):
    log.logger.debug('started.')
    try:
        if args.p <= 2:
            thread_n=args.p
        elif args.p >= 3:
            thread_n=args.p - 1
        if args.alignmentin is True:
            sample_name=os.path.basename(args.b) if not args.b is None else os.path.basename(args.c)
        else:
            sample_name=os.path.basename(args.fq1)
        for refseqid in refseqids:
            # output files
            vcf_gz       =os.path.join(args.outdir, '%s.vcf.gz' % refseqid)
            norm_vcf_gz  =os.path.join(args.outdir, '%s_norm.vcf.gz' % refseqid)
            gatk_naive   =os.path.join(args.outdir, '%s_reconstructed.fa' % refseqid)
            metaspades_o =os.path.join(args.outdir, '%s_metaspades_assembly' % refseqid)
            
            pysam.view(filenames.mapped_to_virus_bam, '-h', '-o', filenames.tmp_bam, refseqid, catch_stdout=False)
            _,seq=utils.retrieve_only_one_virus_fasta(args.vref, refseqid)
            with open(filenames.tmp_fa, 'w') as outfile:
                outfile.write('>%s %s\n%s\n' % (refseqid, sample_name, seq))
            pysam.faidx(filenames.tmp_fa)
            
            # mask low depth regions
            mask_low_depth(args, params, filenames, filenames.tmp_fa, refseqid)
            
            if os.path.exists(filenames.tmp_fa_dict) is True:
                os.remove(filenames.tmp_fa_dict)
            cmd='java -jar %s CreateSequenceDictionary R=%s O=%s' % (args.picard, filenames.tmp_fa, filenames.tmp_fa_dict)
            log.logger.debug('picard command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during gatk running.')
                exit(1)
            cmd='java -jar %s AddOrReplaceReadGroups I=%s O=%s RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20' % (args.picard, filenames.tmp_bam, filenames.tmp_rg_bam)
            log.logger.debug('picard command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during gatk running.')
                exit(1)
    #        pysam.index('-@', str(thread_n), filenames.tmp_rg_bam)
            pysam.index(filenames.tmp_rg_bam)
            cmd='gatk --java-options "-Xmx4g" HaplotypeCaller -R %s -I %s -O %s' % (filenames.tmp_fa, filenames.tmp_rg_bam, vcf_gz)
            log.logger.debug('gatk command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during gatk running.')
                exit(1)
            cmd='bcftools norm -c x -f %s %s -Oz -o %s' % (filenames.tmp_masked_fa, vcf_gz, norm_vcf_gz)
            log.logger.debug('bcftools command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during bcftools running.')
                exit(1)
            cmd='bcftools index %s' % norm_vcf_gz
            log.logger.debug('bcftools command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during bcftools running.')
                exit(1)
            cmd='bcftools consensus -f %s -o %s %s' % (filenames.tmp_masked_fa, gatk_naive, norm_vcf_gz)
            log.logger.debug('bcftools command = `'+ cmd +'`')
            out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
            log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
            if not out.returncode == 0:
                log.logger.error('Error occurred during bcftools running.')
                exit(1)
            # remove unnecessary files
            os.remove(filenames.tmp_rg_bam)
            os.remove(filenames.tmp_rg_bam +'.bai')
            os.remove(filenames.tmp_fa)
            os.remove(filenames.tmp_fa +'.fai')
            os.remove(filenames.tmp_masked_fa)
            os.remove(filenames.tmp_masked_fa +'.fai')
            os.remove(filenames.tmp_fa_dict)
            if args.keep is False:
                os.remove(vcf_gz +'.tbi')
                os.remove(norm_vcf_gz)
                os.remove(norm_vcf_gz +'.csi')
            if args.denovo is True:
                pysam.sort('-@', '%d' % thread_n, '-n', '-O', 'BAM', '-o', filenames.tmp_sorted_bam, filenames.tmp_bam)
                pysam.fastq('-@', '%d' % thread_n, '-N', '-f', '1', '-F', '3852', '-0', '/dev/null', '-1', filenames.tmp_bam_fq1, '-2', filenames.tmp_bam_fq2, '-s', '/dev/null', filenames.tmp_sorted_bam)
                cmd='metaspades.py -1 %s -2 %s -k %s -t %d -m %d -o %s' % (filenames.tmp_bam_fq1, filenames.tmp_bam_fq2, params.metaspades_kmer, thread_n, params.metaspades_memory, metaspades_o)
                log.logger.debug('metaspades command = `'+ cmd +'`')
                out=subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
                log.logger.debug('\n'+ '\n'.join([ l.decode() for l in out.stderr.splitlines() ]))
                if not out.returncode == 0:
                    log.logger.error('Error occurred during metaspades running.')
                    exit(1)
                # remove unnecessary files
                os.remove(filenames.tmp_sorted_bam)
                os.remove(filenames.tmp_bam_fq1)
                os.remove(filenames.tmp_bam_fq2)
            # remove unnecessary files
            os.remove(filenames.tmp_bam)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
