#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-14 
# Version: 0.2.0.0
# Description: pipeline for downsample bam file
#########################################################################
import os
import sys
import time
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Running(input_bam, bed, cfg, output_dir, output_prefix, downsample_num, CollectHsMetrics_excel, running_type):
    if not bac.CheckFile(input_bam):
        print (f"{input_bam} not exist")
        exit()
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    picard = config['software']['picard']
    samtools = config['software']['samtools']
    hg19 = config['database']

    cmd_bedtolist = " ".join([
        java, '-jar', f'-jar -Xmx{java_mem}G',
        picard, 'BedToIntervalList',
        '-I', bed,
        '-O', f'{dirprefix}.BedToIntervalList.xls',
        '-SD', hg19['hg19_dict']
    ])
    cmd_collecths = " ".join([
        java, '-jar', f'-jar -Xmx{java_mem}G',
        picard, 'CollectHsMetrics',
        '-I', input_bam,
        '-O', f'{dirprefix}.CollectHsMetrics_tmp.xls',
        '-R', hg19['hg19_fa'],
        '-BAIT_INTERVALS', f'{dirprefix}.BedToIntervalList.xls',
        '-TARGET_INTERVALS', f'{dirprefix}.BedToIntervalList.xls'
    ])
    cmd_downsample = " ".join([
        samtools, 'view -s percent', input_bam,
        '-o', input_bam.replace('.bam', f'.downsample_{str(downsample_num)}.bam')
    ])
# running
    hsmetrics_result = os.popen(f'cat {CollectHsMetrics_excel} |head -n 2 |tail -n 1').read()
    if hsmetrics_result.startswith('# CollectHsMetrics'):
        dp = float(os.popen(f"cat {CollectHsMetrics_excel} |head -n 8 |tail -n 1").read().strip().split('\t')[9])
    else:
        if os.path.exists(bed) is False or os.path.getsize(bed) == 0:
            print ("CollectHsMetrics_excel file not exists! and bed file not exitsts!")
            exit()
        if running_type == 'subprocess':
            bac.RunningProcess(cmd_bedtolist, output_dir, output_prefix, 'Notice:\tpicard BedToIntervalList running messages')
            bac.RunningProcess(cmd_collecths, output_dir, output_prefix, 'Notice:\tpicard CollectHsMetrics running messages')
        elif running_type == 'os_system':
            os.system(cmd_bedtolist)
            os.system(cmd_collecths)
        dp = float(os.popen(f"cat {dirprefix}.CollectHsMetrics_tmp.xls |head -n 8 |tail -n 1").read().strip().split('\t')[9])
    if dp > float(downsample_num):
        percent = float(downsample_num)/dp
        cmd_downsample = cmd_downsample.replace('percent', str(percent))
        if running_type == 'subprocess':
            bac.RunningProcess(cmd_downsample, output_dir, output_prefix, f'Notice:\tsort bam file samtools view -s downsample {str(downsample_num)} running messages')
        elif running_type == 'os_system':
            os.system(cmd_downsample)
        otstr = 'PASS at the time: {}\tdownsample running successed; {} bam file downsample to {}x\n'.format( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), input_bam, str(downsample_num) )
    else:
        otstr = f'FAIL\tdownsample running failed; {input_bam} bam file MEAN_TARGET_COVERAGE less than {str(downsample_num)}x\n'
    with open(f'{dirprefix}.downsample_messages.txt','a') as otopen:
        otopen.write(otstr)

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing bam file downsample')
    parser.add_argument('-b', '--input_bam', required=True,
                        help='sorted bam file or umi group bam file')
    parser.add_argument('-d', '--bed', default="",
                        help='panel design bed file; default=""')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-num', '--downsample_num', type=int, default=2000,
                        help='downsample number (x); default=2000')
    parser.add_argument('-hs', '--CollectHsMetrics_excel', default="",
                        help='using picard CollectHsMetrics output file instead of running twice; default=""')
    parser.add_argument('-run', '--running_type', choices=['subprocess', 'os_system'], default='os_system',
                        help='running type "subprocess"==>"run_command_with_return", "os_system"==>"os.system"; default=os.system')
    args = parser.parse_args()
    Running(args.input_bam, args.bed, args.cfg, args.output_dir, args.output_prefix, args.downsample_num, args.CollectHsMetrics_excel, args.running_type)
