#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-23
# Version: 0.1.0.0
# Description: pipeline report result
#########################################################################

import os
import re
import sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def gnomesort(a, b):
    i,j,size = 1,2,len(a)
    while i < size:
        if a[i-1] <= a[i]:
            i,j = j, j+1
        else:
            b[i-1], b[i] = b[i], b[i-1]
            a[i-1], a[i] = a[i], a[i-1]
            i -= 1
            if i == 0:
                i,j = j,j+1
    return b

def HotDict(hotspot):
    hotdict = {}
    for each in open(hotspot):
        col = each.strip().split('\t')
        cell = col[0]
        gene = col[1]
        hgvsp = col[2]
        hotdict.update({f"{hgvsp},{gene}:":cell})
    return hotdict

def SnvIndelHot(multianno_vcf, hotdict, output_prefix):
    # snvindellist = ['CELL\tGENE\tHGVS.p\tCHROM\tPOS\tREF\tALT\tDP\tAD\tAF\tDP4\n']
    snvindellist = []
    aflist = []
    for line in open(multianno_vcf):
        if not line.startswith('#'):
            keys = re.findall(r'p\.\w\d+[A-Za-z0-9_\*]+?,\w+:', line) # r'p\.\w\d+.+?,\w+:'
            if keys != []:
                uniq_keys = list(set(keys))
                for k in uniq_keys:
                    if hotdict.get(k) is not None:
                        l = line.split('\t')
                        msg = l[7].split(';')
                        # 这里去重前和去重后得到hg19_multianno.vcf文件内容不同
                        # snv_indel_sortbam/*_all.hg19_multianno.vcf DP=4562;AF=0.001534;SB=4;DP4=2378,2163,3,6;
                        # annotation/*.hg19_multianno.vcf AF=0.003681;DP=2445;DP4=1239,1134,4,5;
                        if l[7].startswith('DP'):
                            dp = msg[0].replace('DP=', '')
                            af = str( round( float(msg[1].replace('AF=', ''))*100, 4) )
                            dp4 = msg[3].replace('DP4=', '')
                        elif l[7].startswith('AF'):
                            dp = msg[1].replace('DP=', '')
                            af = str( round( float(msg[0].replace('AF=', ''))*100, 4) )
                            dp4 = msg[2].replace('DP4=', '')
                        aflist.append(float(af))
                        ad = str( int(dp4.split(',')[2]) + int(dp4.split(',')[3]) )
                        optstr = "\t".join( [ output_prefix, hotdict.get(k), k.replace(':','').split(',')[1], k.replace(':','').split(',')[0], l[0], l[1], l[3], l[4], dp, ad, af, dp4 ] ) + '\n'
                        snvindellist.append(optstr)
    return aflist, snvindellist

def Running(multianno_vcf, cfg, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    hotspot = config['parameter']['hotspot_snvindel']
    hotdict = HotDict(hotspot)
    aflist, snvindellist = SnvIndelHot(multianno_vcf, hotdict, output_prefix)
    sort_snvindellist = gnomesort(aflist, snvindellist)
    sort_snvindellist.reverse()
    with open(f'{dirprefix}.snvindel.hotspot.result.xls','w') as otopen:
        otopen.writelines([f'{output_prefix}\tCELL\tGENE\tHGVS.p\tCHROM\tPOS\tREF\tALT\tDP\tAD\tAF%\tDP4\n'] + sort_snvindellist)

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing report')
    parser.add_argument('-v', '--multianno_vcf', required=True,
                        help='snv indel annovar annotation *.hg19_multianno.vcf file')
    # parser.add_argument('-d', '--bed', required=True,
    #                     help='panel design bed file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    Running(args.multianno_vcf, args.cfg, args.output_dir, args.output_prefix)
