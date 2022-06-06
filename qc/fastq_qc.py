#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-11-7 
# Version: 0.2.1.0
# Description: pipeline for raw reads and trimmed reads QC analysis
#########################################################################
import os
import sys
from itertools import zip_longest
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../software/python-site-packages'))
import numpy as np

def GzipOrNot(filename):
    head = os.popen(f"zcat {filename} |head -n 1").readlines()
    if head == []:
        return os.popen(f"cat {filename}")
    elif head != []:
        return os.popen(f"gunzip -c {filename}")

def DealSummary(adapter_list, raw_readnum):
    r1_msgs = adapter_list[0].split(':')[1].strip().split(' ')
    r2_msgs = adapter_list[1].split(':')[1].strip().split(' ')
    r1_num = int(r1_msgs[0].replace(',',''))
    r2_num = int(r2_msgs[0].replace(',',''))
    r1r2_pct = f'R1{r1_msgs[1]};R2{r2_msgs[1]}'
    adapter_ratio = str(round((r1_num+r2_num)/raw_readnum*100, 2))+'%;'+r1r2_pct
    return adapter_ratio

def Answers(fq1, fq2):
    readnum = 0
    basenum = 0
    overq20 = 0
    overq30 = 0
    gccount = 0
    for i,part in enumerate( zip_longest( GzipOrNot(fq1), GzipOrNot(fq2) ) ) : #part[0] is r1 ASCII, part[1] is r2 ASCII
        if (i+1)%4 == 2:
            gc = part[0].count('G')+part[0].count('C')+part[1].count('G')+part[1].count('C')
            gccount += gc
        if (i+1)%4 == 0:
            readnum += 2
            basenum += len(part[0].strip())+len(part[1].strip())
            asciilist = list(part[0].strip())+list(part[1].strip())
            qscore = np.array( list( map( lambda x: ord(x)-33, asciilist ) ) )
            overq20array = qscore[qscore>=20]
            overq20 += len(overq20array)
            overq30array = qscore[qscore>=30]
            overq30 += len(overq30array)
    gcpercent = round(gccount/basenum*100, 2)
    overq20ratio = round(overq20/basenum*100, 2)
    overq30ratio = round(overq30/basenum*100, 2)
    return readnum, basenum, overq20ratio, overq30ratio, gcpercent

def ResultQC(dirprefix, output_prefix, raw_readnum, raw_basenum, raw_overq20ratio, raw_overq30ratio, trim_readnum, gcpercent, adapter_ratio):
    optlist = []
    optlist.append( "\t".join( [ 
        'sample.id', 
        'total.reads', 'total.bases', 'Q20.ratio(%)',
        'Q30.ratio(%)', 'gc.content(%)', 'clean.ratio(%)', 'adapter.ratio(%)'
    ] )+'\n' )
    optlist.append( "\t".join( [ 
        output_prefix, '{:,}'.format(raw_readnum), '{:,}'.format(raw_basenum), 
        str(raw_overq20ratio), str(raw_overq30ratio), str(gcpercent), 
        str(round(trim_readnum/raw_readnum*100, 2)), adapter_ratio
    ] )+'\n' )
    with open(f'{dirprefix}.fastq_QC.xls', 'w') as otopen:
        otopen.writelines(optlist)

def Running(raw_fq1, raw_fq2, trim_fq1, trim_fq2, output_dir, output_prefix, cutadapt_summary):
    dirprefix = os.path.join(output_dir, output_prefix)
    pool = Pool(processes=2)
    raw = pool.apply_async(Answers, args=(raw_fq1, raw_fq2))
    trim = pool.apply_async(Answers, args=(trim_fq1, trim_fq2))
    pool.close()
    pool.join()
    (raw_readnum, raw_basenum, raw_overq20ratio, raw_overq30ratio, raw_gcpercent) = raw.get()
    (trim_readnum, trim_basenum, trim_overq20ratio, trim_overq30ratio, trim_gcpercent) = trim.get()
    if cutadapt_summary is None:
        adapter_ratio = 'NA'
    elif not os.path.exists(cutadapt_summary):
        adapter_ratio = f'NA ({cutadapt_summary} file wrong!)'
    elif os.path.exists(cutadapt_summary):
        adapter_list = os.popen(f"cat {cutadapt_summary} |grep 'with adapter'").readlines()
        if adapter_list == []:
            adapter_ratio = f'NA ({cutadapt_summary} file wrong!)'
        else:
            adapter_ratio = DealSummary(adapter_list, raw_readnum)
    ResultQC(dirprefix, output_prefix, raw_readnum, raw_basenum, raw_overq20ratio, raw_overq30ratio, trim_readnum, raw_gcpercent, adapter_ratio)

if __name__ == "__main__":
    parser = ArgumentParser('trim')
    parser.add_argument('-raw1', '--raw_fq1', required=True,
                        help='fastq file read1')
    parser.add_argument('-raw2', '--raw_fq2', required=True,
                        help='fastq file read2')
    parser.add_argument('-trim1', '--trim_fq1', required=True,
                        help='fastq file read1')
    parser.add_argument('-trim2', '--trim_fq2', required=True,
                        help='fastq file read2')
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-s', '--cutadapt_summary', required=None,
                        help='fastq file cutadapt running result summary file get "with adapter"; default=None;')
    args = parser.parse_args()
    Running(args.raw_fq1, args.raw_fq2, args.trim_fq1, args.trim_fq2, args.output_dir, args.output_prefix, args.cutadapt_summary)
