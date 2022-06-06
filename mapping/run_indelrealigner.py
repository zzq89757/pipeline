#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15 
# Version: 0.2.0.0
# Description: pipeline for rmdup or umi bam gatk3 indelrealigner
#########################################################################
import os
import sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def running(input_bam, bed, cfg, output_dir, output_prefix):
    if not bac.CheckFile(input_bam):
        print (f"{input_bam} not exist")
        exit()
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    gatk3 = config['software']['gatk3']
    hg19_fa = config['database']['hg19_fa']
    gatk3_threads = config['parameter']['gatk3_threads']
    gatk_db = config['database']
#IndelRealigner
    handle = ''
    # if 'markdup' in input_bam:
    #     handle='markdup_'
    # elif 'rmdup' in input_bam:
    #     handle='rmdup_'
    # else:
    #     handle=''
    cmd_indelrealigner = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        gatk3, '-T RealignerTargetCreator',
        '--num_threads', gatk3_threads,
        '-R', hg19_fa,
        '-I', input_bam,
        '-L', bed,
        '--interval_padding 100',
        '-known', gatk_db['1000G'],
        '-known', gatk_db['Mills'],
        '--downsampling_type NONE --allow_potentially_misencoded_quality_scores', # --fix_misencoded_quality_scores两个参数适用其中任何一个都可以
        '-o', f'{dirprefix}.{handle}realigner.intervals',
        '\n',
        java, f'-jar -Xmx{java_mem}G',
        gatk3, '-T IndelRealigner',
        '-R', hg19_fa,
        '-I', input_bam,
        '-known', gatk_db['1000G'],
        '-known', gatk_db['Mills'],
        '-targetIntervals', f'{dirprefix}.{handle}realigner.intervals',
        '--downsampling_type NONE --allow_potentially_misencoded_quality_scores', # --fix_misencoded_quality_scores
        '-o', f'{dirprefix}.{handle}indelrealigner.bam'
    ])
    bac.RunningProcess(cmd_indelrealigner, output_dir, output_prefix, 'handle:\tgatk3 IndelRealigner running messages')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('IndelRealigner')
    parser.add_argument('-b', '--input_bam', required=True,
                        help='input bam file')
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.input_bam, args.bed, args.cfg, args.output_dir, args.output_prefix)
