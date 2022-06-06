#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15 
# Version: 0.2.0.0
# Description: pipeline for fastq trim based on "cutadapt" software
#########################################################################
import os
import sys
import gzip
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def PlatForm(raw_read1, platform):
    if platform == 'auto':
        platform = 'MGI'
        try:
            try:
                r1_head = open(raw_read1,'r').readline().split(' ')[1][0] #illumina
                platform = 'illumina'
            except IndexError:
                pass #MGI
        except UnicodeDecodeError:
            try:
                r1_head = bytes.decode(gzip.open(raw_read1,'r').readline()).split(' ')[1][0] #illumina
                platform = 'illumina'
            except IndexError:
                pass #MGI
    return platform

def AdapterSeq(platform, parameter):
    if platform == 'MGI':
        r1_adapt_seq = parameter['MGI_read1_adapter_sequence']
        r2_adapt_seq = parameter['MGI_read2_adapter_sequence']
    elif platform == 'illumina':
        r1_adapt_seq = parameter['illumina_read1_adapter_sequence']
        r2_adapt_seq = parameter['illumina_read2_adapter_sequence']
    return r1_adapt_seq, r2_adapt_seq

def Prepares(raw_read1, raw_read2, cfg, output_dir, output_prefix, platform):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    cutadapt = config['software']['cutadapt']
    parameter = config['parameter']
    # https://github.com/marcelm/cutadapt
    platform = PlatForm(raw_read1, platform)
    (r1_adapt_seq, r2_adapt_seq) = AdapterSeq(platform, parameter)
    cmd_trim = " ".join([
        cutadapt, '--cores 4',
        '-a', r1_adapt_seq,
        '-A', r2_adapt_seq,
        'options',
        '-o', f'{dirprefix}.trim.R1.fq.gz',
        '-p', f'{dirprefix}.trim.R2.fq.gz',
        raw_read1, raw_read2
    ])
    return cmd_trim
# /data/software/bin/cutadapt -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAAX -A AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTGX --error-rate 0.1 --overlap 3 --quality-cutoff 20 --minimum-length 50 --max-n 0 -o .//trim/test.trim.R1.fq.gz -p .//trim/test.trim.R2.fq.gz test.R1.fastq.gz test.R2.fastq.gz
def Running(raw_read1, raw_read2, cfg, output_dir, output_prefix, platform, error_rate, overlap, quality_cutoff, minimum_length):
    if not bac.CheckFile(raw_read1) or not bac.CheckFile(raw_read2):
        print (f"{raw_read1} or {raw_read2} not exist")
        exit()
    if bac.DetermainR1R2(raw_read1, raw_read2):
        exit()
    raw_cmd = Prepares(raw_read1, raw_read2, cfg, output_dir, output_prefix, platform)
    options = " ".join([
        '--error-rate', str(error_rate),
        '--overlap', str(overlap),
        '--quality-cutoff', str(quality_cutoff),
        '--minimum-length', str(minimum_length),
        '--max-n 0'
    ])
    cmd_trim = raw_cmd.replace('options', options)
    bac.RunningProcess(cmd_trim, output_dir, output_prefix, 'Notice:\tcutadapt trim running messages')
# --errors E; --overlap MINLENGTH; --revcomp add; --nextseq-trim(for MGI200) not used; --quality-cutoff; --minimum-length; --max-n 0; 
if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('cutadapt trim')
    parser.add_argument('-fq1', '--raw_read1', required=True,
                        help='raw fastq file read1')
    parser.add_argument('-fq2', '--raw_read2', required=True,
                        help='raw fastq file read2')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', default=True,
                        help='sample id')
    parser.add_argument('-t', '--platform', choices=['auto', 'MGI', 'illumina'], default='auto',
                        help='sequencing platform illumina MGI or auto, auto means Automatic prediction platform based on fastq file. eg illumina is ":1092 1:N:0:" MGI is "0000/1". default=auto')
    parser.add_argument('-e', '--error_rate', type=float, default=0.1,
                        help='Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for full-length adapter match (if E is an integer >= 1). Error rate = no. of errors divided by length of matching region. Default: 0.1 (10%%)')
    parser.add_argument('-l', '--overlap', type=int, default=3,
                        help='Require minlength overlap between read and adapter for an adapter to be found. Default: 3')
    parser.add_argument('-q', '--quality_cutoff', type=int, default=20,
                        help="Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second. Default: 20")
    parser.add_argument('-m', '--minimum_length', type=int, default=50,
                        help='Discard reads shorter than length. Default: 50')
    args = parser.parse_args()
    Running(args.raw_read1, args.raw_read2, args.cfg, args.output_dir, args.output_prefix, args.platform, args.error_rate, args.overlap, args.quality_cutoff, args.minimum_length)
