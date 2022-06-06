#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-1-17
# Version: 0.1.0.0
# Description: pipeline for HLA 
#########################################################################

import os
import re
import sys
from glob import glob
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def All_Commands(input_bam, output_dir, output_prefix, hla_bed):
    cmd_samtools = " ".join([
        samtools, 'view', '-b -h', input_bam, '--threads', samtools_threads, '-L', hla_bed, '-o /dev/stdout |', 
        java, f'-jar -Xmx{java_mem}G', picard, 'SamToFastq -I /dev/stdin -F', f'{dirprefix}.HLAregion.1.fastq', 
        '-F2', f'{dirprefix}.HLAregion.2.fastq', '--VALIDATION_STRINGENCY SILENT'
    ])
    cmd_seqtk = " ".join([
        seqtk, 'sample -s 100', f'{dirprefix}.HLAregion.1.fastq', '30000', '>', f'{dirprefix}.HLAregion.3w.1.fastq', '\n',
        seqtk, 'sample -s 100', f'{dirprefix}.HLAregion.2.fastq', '30000', '>', f'{dirprefix}.HLAregion.3w.2.fastq'
    ])
    cmd_optitype = " ".join([
        optitype, '-p', output_prefix, '-i', f'{dirprefix}.HLAregion.3w.1.fastq', f'{dirprefix}.HLAregion.3w.2.fastq', '--dna -v --outdir', output_dir
    ])
    return cmd_samtools, cmd_seqtk, cmd_optitype

def Running(cfg, input_bam, output_dir, output_prefix):
    global dirprefix, java, java_mem, samtools, samtools_threads, picard, seqtk, optitype
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    localdir = os.path.dirname(os.path.abspath(__file__))
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    software = config['software']
    samtools = software['samtools']
    samtools_threads = config['parameter']['samtools_threads']
    picard = software['picard']
    seqtk = software['seqtk']
    optitype = software['optitype']
    hla_bed = config['parameter']['hla_bed']
    (cmd_samtools, cmd_seqtk, cmd_optitype) = All_Commands(input_bam, output_dir, output_prefix, hla_bed)
    bac.RunningProcess(cmd_samtools, output_dir, output_prefix, 'Notice:\tHLA analysis samtools picard step')
    bac.RunningProcess(cmd_seqtk, output_dir, output_prefix, 'Notice:\tHLA analysis seqtk step')
    bac.RunningProcess(cmd_optitype, output_dir, output_prefix, 'Notice:\tHLA analysis optitype step')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing HLA')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-b', '--input_bam', required=True,
                        help='input bam file')
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='library id')
    args = parser.parse_args()
    Running(args.cfg, args.input_bam, args.output_dir, args.output_prefix)
