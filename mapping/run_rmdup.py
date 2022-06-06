#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15 
# Version: 0.2.0.0
# Description: pipeline for sort bam file markduplicates
#########################################################################
import os
import sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Running(input_bam, cfg, output_dir, output_prefix, remove_dup):
    if not bac.CheckFile(input_bam):
        print (f"{input_bam} not exist")
        exit()
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    picard = config['software']['picard']
#MarkDuplicates
    cmd_dup = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        picard, 'MarkDuplicates',
        f'-TMP_DIR {dirprefix}_tmp',
        '-I', input_bam,
        f'-O {dirprefix}.handle.bam',
        '-CREATE_INDEX true',
        '--TAGGING_POLICY All',
        f'-M {dirprefix}.mkdup.metrics',
        '--REMOVE_DUPLICATES NONE'
    ])
    if remove_dup == 'true':
        handle = 'rmdup'
        cmd_dup = cmd_dup.replace('handle','rmdup').replace('NONE','true')
    elif remove_dup == 'false':
        handle = 'markdup'
        cmd_dup = cmd_dup.replace('handle','markdup').replace('NONE','false')
    bac.RunningProcess(cmd_dup, output_dir, output_prefix, f'Notice:\tpicard MarkDuplicates {handle} running messages')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('MarkDuplicates')
    parser.add_argument('-b', '--input_bam', required=True,
                        help='input bam file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-dup', '--remove_dup', choices=['true', 'false'], default='true',
                        help='picard MarkDuplicates removedup; default=true')
    args = parser.parse_args()
    Running(args.input_bam, args.cfg, args.output_dir, args.output_prefix, args.remove_dup)
