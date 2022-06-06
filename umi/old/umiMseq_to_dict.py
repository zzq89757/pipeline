#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-15
# Version: 0.2.0.0
# Description: pipline for UMI seq read to pickle type
#########################################################################
import os
import sys
import pickle
from itertools import zip_longest
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def check_file(umiseq_fq1, umiseq_fq2):
    mode='pair'
    if brc.check_file(umiseq_fq2): # umiseq_fq2 != None
        if brc.determainr1r2(umiseq_fq1, umiseq_fq2):
            exit()
        mode='pair'
    else: # umiseq_fq2 == None
        if brc.check_file(umiseq_fq1):
            mode='single'
        else:
            print ("%s not exist"%umiseq_fq1)
            exit()
    return mode

def singleUMI(umiseq_fq1):
    RX_dict={}
    for i, line in enumerate( os.popen('gunzip -c %s'%umiseq_fq1) ):
        if i%4 == 0: #read
            name=line.strip()[1:-2]
            continue
        elif i%4 == 1: #seq
            rx='RX:Z:'+line.strip()
            RX_dict[name]=rx
    return RX_dict

def pairUMI(umiseq_fq1, umiseq_fq2):
    RX_dict={}
    for i, part in enumerate( zip_longest( os.popen('gunzip -c %s'%umiseq_fq1), os.popen('gunzip -c %s'%umiseq_fq2) ) ):
        if i%4 == 0: #read
            name=part[0].strip()[1:-2]
            continue
        elif i%4 == 1: #seq
            rx='RX:Z:'+part[0].strip()+'-'+part[1].strip()
            RX_dict[name]=rx
    return RX_dict

def running(umiseq_fq1, umiseq_fq2, output_dir, output_prefix):
    mode = check_file(umiseq_fq1, umiseq_fq2)
    dirprefix = os.path.join(output_dir, output_prefix)
    brc.createDirectory(output_dir+'/'+'%s_umi_steplog'%output_prefix)
    if mode == 'pair':
        RX_dict = pairUMI(umiseq_fq1, umiseq_fq2)
    elif mode == 'single':
        RX_dict = singleUMI(umiseq_fq1)
    with open('%s_umi_steplog/%s_RX_dict'%(dirprefix, output_prefix), 'wb') as otopen:
        pickle.dump(RX_dict, otopen)

if __name__ == '__main__':
    parser = ArgumentParser('read UMI sequence fastq file to dict Binary pickle type file')
    parser.add_argument('-m1', '--umiseq_fq1', required=True,
                        help='umi fastq file read1')
    parser.add_argument('-m2', '--umiseq_fq2', default=None,
                        help='umi fastq file read2, default=None')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.umiseq_fq1, args.umiseq_fq2, args.output_dir, args.output_prefix)
