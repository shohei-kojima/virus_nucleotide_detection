#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import os,gzip,shutil
import log,traceback


class empclass:
    pass


def gzip_or_del(args, params, file):
    log.logger.debug('started,file=%s' % file)
    try:
        if args.keep is True:
            with open(file, 'rt') as f_in:
                with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    f_out.flush()
                    os.fdatasync(f_out.fileno())
        os.remove(file)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def gzip_file(params, file):
    log.logger.debug('started,file=%s' % file)
    try:
        with open(file, 'rt') as f_in:
            with gzip.open(file +'.gz', 'wt', compresslevel=params.gzip_compresslevel) as f_out:
                shutil.copyfileobj(f_in, f_out)
                f_out.flush()
                os.fdatasync(f_out.fileno())
        os.remove(file)
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)


def gzip_d(gzfile):
    log.logger.debug('started,infile=%s' % gzfile)
    try:
        outfpath=gzfile[:-3]
        with open(outfpath, 'w') as outfile:
            with gzip.open(gzfile) as infile:
                for line in infile:
                    line=line.decode()
                    outfile.write(line)
            outfile.flush()
            os.fdatasync(outfile.fileno())
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)



def parse_fasta(path_to_file):
    tmp={}
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=seq
                header=line.strip().replace(' ', '_')
                seq=''
            elif '>' in line and not seq:
                header=line.strip().replace(' ', '_')
            else:
                seq += line.strip()
        tmp[header]=seq
    return tmp


def retrieve_only_one_virus_fasta(path_to_file, refseqid):
    seq=''
    keep=False
    with open(path_to_file) as infile:
        for line in infile:
            if refseqid in line:
                header=line.strip().replace('>', '')
                keep=True
            elif keep is True:
                if not '>' in line:
                    seq += line.strip()
                else:
                    break
    return header, seq

