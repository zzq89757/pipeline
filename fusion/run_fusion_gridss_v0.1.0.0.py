#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-14
# Version: 0.1.0.0
# Description: pipeline for gridss software fusion detect
#########################################################################
import os
import sys
import json
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def GridssResult(sort_bam, gridss, gridss_jar, hg19_fa, fusion_black_list_bed, output_dir, output_prefix, dirprefix):
    cmd_fusion = " ".join([
        gridss, '--threads 8',
        '--reference', hg19_fa,
        '--jar', gridss_jar,
        '--blacklist', fusion_black_list_bed,
        '--workingdir', output_dir,
        '--output', f"{dirprefix}.fusion_raw.vcf", sort_bam
    ])
    bac.RunningProcess(cmd_fusion, output_dir, output_prefix, 'Notice:\tgridss fusion analysis running messages')

def Running(sort_bam, cfg, output_dir, output_prefix, hotspot):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    gridss = config['software']['gridss']
    gridss_jar = config['software']['gridss_jar']
    hg19_fa = config['database']['hg19_fa']
    fusion_black_list_bed = config['parameter']['fusion_black_list_bed']
    GridssResult(sort_bam, gridss, gridss_jar, hg19_fa, fusion_black_list_bed, output_dir, output_prefix, dirprefix)
    fusion_raw_vcf = f"{dirprefix}.fusion_raw.vcf"
    os.system(f"chmod 755 {output_dir}")

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing fusion')
    parser.add_argument('-b', '--sort_bam', required=True,
                        help='sort bam file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-s', '--hotspot', type=str, default="",
                        help='COSMIC fusion hotspot bed file, default=""')
    # parser.add_argument('-f', '--polish', type=str, default="",
    #                     help='')
    args = parser.parse_args()
    Running(args.sort_bam, args.cfg, args.output_dir, args.output_prefix, args.hotspot)
