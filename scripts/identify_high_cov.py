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
        # load virus names from virus reference seq file
        virus_names={}
        with open(args.vref) as infile:
            for line in infile:
                if '>' in line:
                    ls=line.strip().split(' ', 1)
                    virus_names[ls[0].replace('>', '')]=ls[1]
        # identify high cov viruses
        prev_id='any'
        high_cov=[]
        for_plot_d={}
        tmp_retain=[]
        with open(filenames.summary, 'w') as outfile:
            with open(filenames.bedgraph) as infile:
                for line in infile:
                    ls=line.split()
                    if ls[0] == prev_id:
                        if int(ls[3]) >= 1:
                            cov=int(ls[3])
                            for _ in range(int(ls[1]), int(ls[2])):
                                covs.append(cov)
                    else:
                        if not prev_id == 'any':
                            cov_len=len(covs)
                            genome_covered= cov_len / total_len
                            ave_depth= sum(covs) / total_len
                            if cov_len >= 1:
                                ave_depth_norm= sum(covs) / cov_len
                                high_cov_judge='False'
                                if genome_covered >= params.genome_cov_thresholds:
                                    if ave_depth_norm >= params.ave_depth_of_mapped_region_threshold:
                                        high_cov.append([prev_id, genome_covered, ave_depth_norm])
                                        for_plot_d[prev_id]=tmp_retain
                                        high_cov_judge='True'
                            else:
                                ave_depth_norm=0
                                high_cov_judge='False'
                            outfile.write('%s\tvirus_exist=%s\tgenome_length=%d;mapped_length=%d;perc_genome_mapped=%f;average_depth=%f;average_depth_of_mapped_region=%f\tfasta_header=%s\n' % (prev_id, high_cov_judge, total_len, cov_len, 100 * genome_covered, ave_depth, ave_depth_norm, virus_names[prev_id]))
                            tmp_retain=[]
                        total_len=0
                        covs=[]
                        if int(ls[3]) >= 1:
                            cov_len= int(ls[2]) - int(ls[1])
                            for _ in range(int(ls[1]), int(ls[2])):
                                covs.append(cov)
                    total_len += int(ls[2]) - int(ls[1])
                    prev_id=ls[0]
                    tmp_retain.append([ int(i) for i in ls[1:4] ])
            cov_len=len(covs)
            genome_covered= cov_len / total_len
            ave_depth= sum(covs) / total_len
            if cov_len >= 1:
                ave_depth_norm= sum(covs) / cov_len
                high_cov_judge='False'
                if genome_covered >= params.genome_cov_thresholds:
                    if ave_depth_norm >= params.ave_depth_of_mapped_region_threshold:
                        high_cov.append([prev_id, genome_covered, ave_depth_norm])
                        for_plot_d[prev_id]=tmp_retain
                        high_cov_judge='True'
            else:
                ave_depth_norm=0
                high_cov_judge='False'
            outfile.write('%s\tvirus_exist=%s\tgenome_length=%d;mapped_length=%d;perc_genome_mapped=%f;average_depth=%f;average_depth_of_mapped_region=%f\tfasta_header=%s\n' % (prev_id, high_cov_judge, total_len, cov_len, 100 * genome_covered, ave_depth, ave_depth_norm, virus_names[prev_id]))
        if len(high_cov) >= 1:
            log.logger.info('high_cov_virus=%s' % ';'.join([ l[0] for l in high_cov ]))
        global hhv6a_highcov
        global hhv6b_highcov
        high_cov_set=set([ l[0] for l in high_cov ])
        hhv6a_highcov=True if 'NC_001664.4' in high_cov_set else False
        hhv6b_highcov=True if 'NC_000898.1' in high_cov_set else False

        if len(for_plot_d) >= 1:
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
