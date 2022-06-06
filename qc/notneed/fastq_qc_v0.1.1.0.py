#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-13 
# Version: 0.2.0.0
# Description: pipline for raw reads and trimmed reads QC analysis
#########################################################################
import os
import time
import numpy as np
from itertools import zip_longest
from multiprocessing import Pool
from argparse import ArgumentParser

def gzipornot(filename):
    head = os.popen("zcat %s |head -n 1"%filename).readlines()
    if head == []:
        return os.popen("cat %s"%filename)
    elif head != []:
        return os.popen("/usr/bin/gunzip -c %s"%filename)

def answers(fq1, fq2):
    time.sleep(10)
    readnum=0
    basenum=0
    overq30=0
    gccount=0
    for i,part in enumerate( zip_longest( gzipornot(fq1), gzipornot(fq2) ) ) : #part[0] is r1 ASCII, part[1] is r2 ASCII
        if (i+1)%4 == 2:
            gc=part[0].count('G')+part[0].count('C')+part[1].count('G')+part[1].count('C')
            gccount+=gc
        if (i+1)%4 == 0:
            readnum+=2
            basenum+=len(part[0].strip())+len(part[1].strip())
            asciilist=list(part[0].strip())+list(part[1].strip())
            qscore=np.array( list( map( lambda x: ord(x)-33, asciilist ) ) )
            overq30array=qscore[qscore>=30]
            overq30+=len(overq30array)
    gcpercent=round(gccount/basenum*100, 2)
    overq30ratio=round(overq30/basenum*100, 2)
    return readnum, basenum, overq30ratio, gcpercent

def resultQC(dirprefix, output_prefix, raw_readnum, raw_basenum, raw_overq30ratio, trim_readnum, gcpercent):
    optlist=[]
    optlist.append( "\t".join( [ 
        'sample.id', 
        'total.read', 'total.bases', 
        'Q30.ratio(%)', 'gc.content(%)', 'pass.quality.ratio(%)' 
    ] )+'\n' )
    optlist.append( "\t".join( [ 
        output_prefix, '{:,}'.format(raw_readnum), '{:,}'.format(raw_basenum), 
        str(raw_overq30ratio), str(gcpercent), 
        str(round(trim_readnum/raw_readnum*100, 2)) 
    ] )+'\n' )
    with open('%s.fastq_QC.xls'%dirprefix,'w') as otopen:
        otopen.writelines(optlist)

def running(raw_fq1, raw_fq2, trim_fq1, trim_fq2, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    pool = Pool(processes=2)
    raw = pool.apply_async(answers, args=(raw_fq1, raw_fq2))
    trim = pool.apply_async(answers, args=(trim_fq1, trim_fq2))
    pool.close()
    pool.join()
    (raw_readnum, raw_basenum, raw_overq30ratio, gcpercent)=raw.get()
    (trim_readnum, trim_basenum, trim_overq30ratio, trim_gcpercent)=trim.get()
    resultQC(dirprefix, output_prefix, raw_readnum, raw_basenum, raw_overq30ratio, trim_readnum, gcpercent)

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
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.raw_fq1, args.raw_fq2, args.trim_fq1, args.trim_fq2, args.output_dir, args.output_prefix)