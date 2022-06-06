#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-18
# Version: 0.2.1.0
# Description: pipeline for gridss software fusion detect
#########################################################################
import os
import re
import sys
import pandas as pd
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def GridssResult(sort_bam, gridss, gridss_jar, hg19_fa, blacklist):
    cmd_fusion = " ".join([
        gridss, '--threads 8',
        '--reference', hg19_fa,
        '--jar', gridss_jar,
        '--blacklist', blacklist,
        '--workingdir', otdir,
        '--output', f"{dirprefix}.fusion_raw.vcf", sort_bam
    ])
    bac.RunningProcess(cmd_fusion, otdir, otprefix, 'Notice:\tgridss fusion analysis running messages')

def TableAnnovar(perl, table_annovar, annovar_humandb, fusion_vcf):
    cmd_tableannovar = " ".join([
        perl, table_annovar, fusion_vcf, annovar_humandb,
        '-buildver hg19', '-out', dirprefix, '-remove',
        '-protocol refGene -operation g --vcfinput'
    ])
    bac.RunningProcess(cmd_tableannovar, otdir, otprefix, 'Notice:\tfusion raw vcf annovar analysis running messages')

def VcfFilter(input_vcf):
    optlist = []
    for line in open(input_vcf, 'r'):
        if line.startswith('#'):
            optlist.append(line)
        else:
            col = line.split('\t')
            optline = col
            chrn = col[0]
            pos = col[1]
            bk2 = col[4]
            try:
                chrn_pos = re.sub(r'\[|\]', ' ', bk2).split(' ')[1]
            except IndexError:
                continue
            optline[0] = chrn_pos.split(':')[0]
            optline[1] = chrn_pos.split(':')[1]
            optline[4] = f'{chrn}:{pos}'
            optlist.append("\t".join(optline))
    with open(f'{dirprefix}.bk2.vcf', 'w') as otopen:
        otopen.writelines(optlist)

def HotTable(hotspot_table):
    hotdict = {}
    for line in open(hotspot_table):
        gene1 = line.strip().split('\t')[1]
        gene2 = line.strip().split('\t')[3]
        if hotdict.get(gene1) is None:
            hotdict[gene1] = [gene2]
        else:
            hotdict.get(gene1).append(gene2)
        if hotdict.get(gene2) is None:
            hotdict[gene2] = [gene1]
        else:
            hotdict.get(gene2).append(gene1)
    return hotdict

# head:
# VF: Count of fragments supporting the variant breakpoint allele and not the reference allele.
# SR: Count of split reads supporting breakpoint
# RP: Count of read pairs supporting breakpoint
# REFPAIR: Count of reference read pairs spanning this breakend supporting the reference allele
# REF: Count of reads mapping across this breakend
# AF: Allele fraction
def BedResult(anno_vcf):
    pdlist = []
    optlist = []
    for line in os.popen(f"cat {anno_vcf} |grep -v '#'"):
        col = line.strip().split('\t')
        bk1 = col[4]
        bk2 = f'{col[0]}:{col[1]}'
        genes = re.findall(r'Gene.refGene=(.*?);GeneDetail', col[7])
        gene1 = genes[0]
        gene2 = genes[1]
        vf = col[9].split(':')[-1]
        sr = col[9].split(':')[-3]
        rp = col[9].split(':')[-5]
        refpair = col[9].split(':')[-6]
        ref = col[9].split(':')[-7]
        af = col[9].split(':')[1]
        pdlist.append( [ bk1, bk2, gene1, gene2, vf, sr, rp, refpair, ref, af ] )
        optstr = "\t".join( [ bk1, bk2, gene1, gene2, vf, sr, rp, refpair, ref, af ] ) + '\n'
        optlist.append(optstr)
    with open(f'{dirprefix}.all_table.txt','w') as otopen:
        otopen.writelines(optlist)
    return pdlist

def Result(pdlist, hotdict):
    hotresult = f'{dirprefix}.sv_hottable.xls'
    column_list = ['bk1','bk2','gene1','gene2','VF','SR','RP','REFPAIR','REF','AF'] # type: str str str str int int int int int float
    dflist = []
    resultlist = []
    for e in pdlist:
        try:
            case1 = e[3] in hotdict.get(e[2])
        except TypeError:
            case1 = False
        try:
            case2 = e[2] in hotdict.get(e[3])
        except TypeError:
            case2 = False
        if case1 or case2:
            # print (e)
            # resultlist.append("\t".join(e)+'\n')
            dflist.append(e)
    if dflist == []:
        with open(hotresult,'w') as otopen:
            otopen.write('\t'+"\t".join(column_list)+'\n')
        return hotresult
    df = pd.DataFrame(dflist)
    df.columns = column_list
    df_intVF = df.astype({'VF':'int'})
    sort_df = df_intVF.sort_values(by=['VF'], ascending=False)
    sort_df.index = [x+1 for x in range(len(sort_df.index))]
    sort_df.to_csv(hotresult, sep='\t', header=True, index=True)
    return hotresult

def ReportVcf(hotresult, hotspot_table):
    head = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=CTX,Description="Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    title = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcsv_{dirprefix}_Tumor\n"
    optlist = []
    optlist.append(head)
    optlist.append(title)
    for line in open(hotresult, 'r'):
        if not line.startswith('\t'):
            col = line.strip().split('\t')
            chr1 = col[1].split(':')[0].replace('chr','')
            pos1 = col[1].split(':')[1]
            # genes = col[3]+'-'+col[4]
            alt = os.popen(f"cat {hotspot_table} |grep {col[3]} |grep {col[4]} ").read().split('\t')[-1].strip()
            chr2 = col[2].split(':')[0].replace('chr','')
            pos2 = col[2].split(':')[1]
            if alt == 'CTX':
                svlen = 0
            else:
                svlen = str( abs( int(pos1) + int(pos2) )+1 )
            dp = str( int(col[-3]) + int(col[-2]) )
            ad = col[5]
            if int(ad) >= 5: # and int(dp) > int(ad)
                tag = 'PASS'
            else:
                tag = 'LQ'
            if int(dp) == 0:
                af = col[-1]
            else:
                af = str( float(ad) / int(dp) )
            optstr = f"{chr1}\t{pos1}\tN\t.\t<{alt}>\t100\t{tag}\tSOMATIC;BND=LR;SVTYPE={alt};CHR2={chr2};END={pos2};SVLEN={svlen}\tGT:AD:AF\t./.:{dp},{ad}:{af}\n"
            optlist.append(optstr)
    with open(f'{otdir}/csv_{otprefix}.sv.vcf', 'w') as otopen:
        otopen.writelines(optlist)

def Running(sort_bam, cfg, output_dir, output_prefix, hotspot):
    global dirprefix, otdir, otprefix
    otdir = output_dir
    otprefix = output_prefix
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    perl = config['language']['perl']
    software = config['software']
    table_annovar = software['table_annovar']
    gridss = software['gridss']
    gridss_jar = software['gridss_jar']
    hg19_fa = config['database']['hg19_fa']
    annovar_humandb = config['database']['annovar_humandb']
    blacklist = config['parameter']['fusion_black_list_bed']
    hotspot = config['parameter']['hotspot_fusion']
    
    GridssResult(sort_bam, gridss, gridss_jar, hg19_fa, blacklist)
    os.system(f"chmod 755 {output_dir}")
    fusion_raw_vcf = f"{dirprefix}.fusion_raw.vcf"
    
    TableAnnovar(perl, table_annovar, annovar_humandb, fusion_raw_vcf)
    VcfFilter(f'{dirprefix}.hg19_multianno.vcf')
    TableAnnovar(perl, table_annovar, annovar_humandb, f'{dirprefix}.bk2.vcf')
    anno_vcf = f'{dirprefix}.hg19_multianno.vcf'

    pdlist = BedResult(anno_vcf)
    hotdict = HotTable(hotspot)
    hotresult = Result(pdlist, hotdict)
    # hotresult = f'{dirprefix}.sv_hottable.xls'
    ReportVcf(hotresult, hotspot)

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
    parser.add_argument('-s', '--hotspot', type=str, default=f"{cfgfile}/../report/hotspot_fusion_gene_table.xls",
                        help='COSMIC fusion hotspot bed file, default="%s/../report/hotspot_fusion_gene_table.xls"'%cfgfile)
    args = parser.parse_args()
    Running(args.sort_bam, args.cfg, args.output_dir, args.output_prefix, args.hotspot)
