#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import log,traceback

class load:
    def __init__(self, args):
        log.logger.debug('started')
        try:
            # default
            self.gzip_compresslevel=1
            self.sam_mapq_threshold=5
            self.bedgraph_bin=1
            self.read_cov_threshold=3
            self.genome_cov_thresholds=0.05
            self.depth_threshold=1
            
            params_for_debug=[]
            for k,v in self.__dict__.items():
                params_for_debug.append('%s=%s' % (k, str(v)))
            log.logger.debug('parameters:\n'+ '\n'.join(params_for_debug))
        except:
            log.logger.error('\n'+ traceback.format_exc())
            exit(1)
