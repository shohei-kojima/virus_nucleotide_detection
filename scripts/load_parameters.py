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
            self.discordant_reads_clip_len=20
            
            
            params_for_debug=[]
            for k,v in self.__dict__.items():
                params_for_debug.append('%s=%s' % (k, str(v)))
            log.logger.debug('parameters:\n'+ '\n'.join(params_for_debug))
        except:
            log.logger.error('\n'+ traceback.format_exc())
            exit(1)
