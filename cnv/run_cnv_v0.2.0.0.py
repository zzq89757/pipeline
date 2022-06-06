#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15
# Version: 0.2.0.0
# Description: pipeline for cnv detection
#########################################################################
import os
import sys
import math
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Standardcnv(cn):
    if cn >= 4.0:
        vartype = 'gain'
    elif cn <= 1.5:
        vartype = 'loss'
    else:
        vartype = 'normal'
    return vartype

def Result(cnvgene, cnvbin, hotspot):
    # 检查panel是否选错
    checkbin = os.popen("cat %s |awk -F \"\t\" 'NR>1{print $4\"\t\"$5}' |sort |uniq" % cnvbin).readlines()[0].strip()
    checkgene = os.popen("cat %s |awk -F \"\t\" 'NR>1{print $2\"\t\"$3}' |sort |uniq" % cnvgene).readlines()[0].strip()
    if checkbin == '' or checkgene == '':
        print ("this bed file not build cnv base line!")
        exit()

    hotdict = {}
    with open(hotspot, 'r') as inopen:
        for i,line in enumerate(inopen):
            if i > 0:
                gene = line.split('\t')[3]
                hotdict[gene] = line.strip()

    # copynumber = [] 不需要cn值从大到小排序
    cnvgene_msg = {}
    cnvbin_msg = {}
    with open(cnvgene, 'r') as inopen:
        for line in inopen:
            col = line.split('\t') # gene name is col[0]
            if hotdict.get(col[0]) is not None:
                cn = float(col[1])
                if int(cn) == 0:
                    cn = 0.0000000001
                vartype = Standardcnv(cn)
                cn_value = str( round(cn, 1) )
                log2_ratio = str( round( math.log2( cn / 2 ), 4 ) )
                cnvgene_msg.update( {col[0]: [log2_ratio, cn_value, vartype]} )
                cnvbin_msg.update( {col[0]: [] } )
                # cnvbin_msg.update( {col[0]: {'binlist':[], 'binnum':[]} } )
                # copynumber.append(cn_value)
    # n = 1
    # tmp = []
    with open(cnvbin, 'r') as inopen:
        for line in inopen:
            col = line.split('\t')
            name = col[1] # gene name is col[1]
            # tmp.append(name)
            # if len(tmp) > 1:
            #     if tmp[0] == tmp[1]:
            #         n += 1
            #     elif tmp[0] != tmp[1]:
            #         n = 1
            #     ss = tmp.pop(0)
            if cnvgene_msg.get(col[1]) is not None:
                gene_vartype = cnvgene_msg.get(col[1])[2]
                cn = float(col[3])
                if Standardcnv(cn) == gene_vartype:
                    cnvbin_msg.get(name).append(col[0])
                    # cnvbin_msg.get(name).get('binlist').append(col[0])
                    # cnvbin_msg.get(name).get('binnum').append(n) # 这里bin的个数收集的不是全部的而是符合要求的第几个

    optlist = ['sample\tchr\tstart\tend\tgene\tstart_exon\tend_exon\tbase_length\texon_cnt\tlog2_ratio\tcn\tgain/loss\tbin_support_CNV_Rate\tbin_support_CNV\ttotal_bin\n']
    for k,v in cnvgene_msg.items():
        col2_9 = hotdict.get(k)
        col10_12 = "\t".join(v)
        bin_supt = len( cnvbin_msg.get(k) )
        total_bin = int( col2_9.split('\t')[-1] )
        bin_rate = round( float(bin_supt)/total_bin, 2)
        optstr = "\t".join ([otprf, col2_9, col10_12, str(bin_rate), str(bin_supt), str(total_bin)] ) + '\n'
        optlist.append(optstr)
    with open(f"{otdir}/{otprf}.cnv_allhot.xls", 'w') as otopen:
        otopen.writelines(optlist)
    os.system(f"cat {otdir}/{otprf}.cnv_allhot.xls |grep -v 'normal' > {otdir}/csv_{otprf}.gene.xls")

def Running(bam, targetcoverage, bed, cfg, output_dir, output_prefix, vcf, hotspot, plotting):
    global dirprefix, otdir, otprf
    otdir = output_dir
    otprf = output_prefix
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    cnv_for_gene = config['scripts']['cnv_for_gene']
    cnv_for_bin = config['scripts']['cnv_for_bin']
    cnv_gene_baseline = config['parameter']['cnv_gene_baseline']
    cnv_bin_baseline = config['parameter']['cnv_bin_baseline']
    python = config['language']['python']
    # bam targetcoverage bed这三个参数都加入，cnv_for_gene及cnv_for_bin会自行判断三者的存在关系自适应分析条件；
    cmd_cnv_bin = " ".join([
        python, cnv_for_bin,
        '--refbaseline', cnv_bin_baseline,
        '--bam', bam,
        '--bed', bed,
        '--targetcoverage', targetcoverage,
        '--vcf', vcf,
        '--output', output_dir
    ])
    if bac.GetSize(targetcoverage) == 0:
        targetcoverage = f"{dirprefix}.targetcov"
    cmd_cnv_gene = " ".join([
        python, cnv_for_gene,
        '--refbaseline', cnv_gene_baseline,
        '--bam', bam,
        '--bed', bed,
        '--targetcoverage', targetcoverage,
        '--vcf', vcf,
        '--output', output_dir
    ])
    bac.RunningProcess(cmd_cnv_bin, output_dir, output_prefix, 'Notice:\tcnv bin Running messages')
    bac.RunningProcess(cmd_cnv_gene, output_dir, output_prefix, 'Notice:\tcnv gene Running messages')
    cnvgene = f"{output_dir}/{output_prefix}_cnvgene.txt"
    cnvbin = f"{output_dir}/{output_prefix}_cnvbin.txt"
    Result(cnvgene, cnvbin, hotspot)

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Copy Number Variant detection')
    parser.add_argument('-b', '--bam', type=str, default="None",
                        help='bam file of tumor sample, default="None"')
    parser.add_argument('-t', '--targetcoverage', type=str, default="None",
                        help='target coverage of tumor sample, if "bam" parameter is False "targetcoverage" file must exists, default="None"')
    parser.add_argument('-d', '--bed', type=str, default="None",
                        help='panel bed file of tumor sample, if "targetcoverage" parameter exists bed file not need, default="None"')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='.',
                        help='output directory, default="."')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id, required=True')
    parser.add_argument('-v', '--vcf', type=str, default="None",
                        help='vcf file of tumor sample, default="None"')
    parser.add_argument('-s', '--hotspot', type=str, default=f"{cfgdir}/../report/hotspot_cnv_gene_table.xls",
                        help='COSMIC cnv hotspot bed file, default="%s/../report/hotspot_cnv_gene_table.xls"'%cfgdir)
    parser.add_argument('-g', '--plotting', type=bool, default=False,
                        help='whether to visualize the distribution of sample AF, default=False') 
    args = parser.parse_args()
    Running(args.bam, args.targetcoverage, args.bed, args.cfg, args.output_dir, args.output_prefix, args.vcf, args.hotspot, args.plotting)
