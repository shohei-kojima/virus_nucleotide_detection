#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.
'''


import log,traceback
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import text

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


def identify_high_cov_virus_from_bedgraph(args, params, filenames):
    log.logger.debug('started.')
    try:
        prev_id='any'
        global high_cov
        high_cov=[]
        for_plot_d={}
        tmp_retain=[]
        with open(filenames.bedgraph) as infile:
            for line in infile:
                ls=line.split()
                if ls[0] == prev_id:
                    if int(ls[3]) >= params.read_cov_threshold:
                        cov_len += int(ls[2]) - int(ls[1])
                else:
                    if not prev_id == 'any':
                        if (cov_len / total_len) > params.genome_cov_thresholds:
                            high_cov.append(prev_id)
                            for_plot_d[prev_id]=tmp_retain
                        tmp_retain=[]
                    if int(ls[3]) >= params.read_cov_threshold:
                        cov_len= int(ls[2]) - int(ls[1])
                    else:
                        cov_len=0
                    total_len=0
                total_len += int(ls[2]) - int(ls[1])
                prev_id=ls[0]
                tmp_retain.append([ int(i) for i in ls[1:4] ])
        if (cov_len / total_len) > params.genome_cov_thresholds:
            high_cov.append(prev_id)
            for_plot_d[prev_id]=tmp_retain
        if len(high_cov) >= 1:
            log.logger.debug('high_cov_virus=%s' % ';'.join(high_cov))
        global hhv6a_highcov
        global hhv6b_highcov
        high_cov_set=set(high_cov)
        hhv6a_highcov=True if 'NC_001664.4' in high_cov_set else False
        hhv6b_highcov=True if 'NC_000898.1' in high_cov_set else False

        if len(for_plot_d) >= 1:
            # load virus names from virus reference seq file
            virus_names={}
            with open(args.vref) as infile:
                for line in infile:
                    if '>' in line:
                        ls=line.strip().split(' ', 1)
                        virus_names[ls[0].replace('>', '')]=ls[1]
            # plot all high cov viruses
            plt.figure(figsize=(5, len(for_plot_d)+1))
            gs=gridspec.GridSpec(len(for_plot_d), 1, height_ratios=[ 1 for i in range(len(for_plot_d))])
            nums=[ n for n in range(len(for_plot_d)) ]
            labels=[ id +':'+ virus_names[id] for id in for_plot_d ]
            data=[ for_plot_d[id] for id in for_plot_d ]
            for d,n,la in zip(data, nums, labels):
                ax=plt.subplot(gs[n])
                for i in d:
                    if i[2] >= 1:
                        rect=matplotlib.patches.Rectangle((i[0], 0), i[1]-i[0], i[2], color='dodgerblue', ec=None)
                        ax.add_patch(rect)
                ax.set_xlim([0, i[1]])
                ymax=1
                for i in d:
                    if i[2] > ymax:
                        ymax=i[2]
                ax.set_ylim([0, ymax])
                ax.text(0, ymax, la, ha='left', va='top')
            plt.suptitle('Virus(es) with high coverage mapping')
            plt.savefig(filenames.high_cov_pdf)
        
        
    except:
        log.logger.error('\n'+ traceback.format_exc())
        exit(1)
