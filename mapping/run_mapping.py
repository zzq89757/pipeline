#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15
# Version: 0.2.0.0
# Description: pipeline for trim read mapping
#########################################################################
import os
import sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Running(trimmed_fq1, trimmed_fq2, cfg, output_dir, output_prefix):
    if not bac.CheckFile(trimmed_fq1) or not bac.CheckFile(trimmed_fq2):
        print (f"{trimmed_fq1} or {trimmed_fq2} not exist")
        exit()
    if bac.DetermainR1R2(trimmed_fq1, trimmed_fq2):
        exit()
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    # java=config['language']['java']
    # java_mem=config['language']['java_mem']
    software = config['software']
    bwa = software['bwa']
    samtools = software['samtools']
    # picard=software['picard']
    hg19_fa = config['database']['hg19_fa']
    samtools_threads = config['parameter']['samtools_threads']
    bwa_threads = config['parameter']['bwa_threads']

#Mapping
    cmd_mapping = " ".join([
        bwa, f'mem -t {bwa_threads} -Y -M -R',
        f"'@RG\\tID:{output_prefix}\\tSM:{output_prefix}\\tPL:MGI\\tLB:{output_prefix}\\tPE:100'",
        hg19_fa, trimmed_fq1, trimmed_fq2,
        '|', samtools, 'sort', '--threads', samtools_threads, '-m', '6G',
        '-o', f'{dirprefix}.sort.bam'
    ])
    # '|',java['java_path'],'-jar -Xmx%sG'%java['java_mem'],picard['picard_path'],'AddOrReplaceReadGroups',
    # '--TMP_DIR','%s_tmp'%dirprefix,'--INPUT /dev/stdin',
    # '--OUTPUT','%s.sort.bam'%dirprefix,
    # '--CREATE_INDEX true --SORT_ORDER coordinate','--RGID %s'%output_prefix,
    # '--RGLB Targetseq --RGPL MGI --RGPU MGIseq','--RGSM %s'%output_prefix])
    bac.RunningProcess(cmd_mapping, output_dir, output_prefix, 'Notice:\tbwa mem used, trimmed fastq mapping genome Running messages')
    cmd_index = " ".join([ samtools, 'index', f'{dirprefix}.sort.bam', f'{dirprefix}.sort.bai' ])
    bac.RunningProcess(cmd_index, output_dir, output_prefix, 'Notice:\tsamtools index used, sorted bam samtools index Running messages')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing mapping')
    parser.add_argument('-f1', '--trimmed_fq1', required=True,
                        help='trimmed fastq file read1')
    parser.add_argument('-f2', '--trimmed_fq2', required=True,
                        help='trimmed fastq file read2')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    Running(args.trimmed_fq1, args.trimmed_fq2, args.cfg, args.output_dir, args.output_prefix)
