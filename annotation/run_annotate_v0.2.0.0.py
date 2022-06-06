#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-10-31 
# Version: 0.2.0.0
# Description: pipeline for snv indel calling result annotation
#########################################################################
import os
import sys
import json
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Prepare(input_bam, bed, cfg, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    software = config['software']
    snpEff = software['snpEff']
    SnpSift = software['SnpSift']
    table_annovar = software['table_annovar']
    Rscript = config['language']['Rscript']
    perl = config['language']['perl']
    hg19_fa = config['database']['hg19_fa']
    annovar_humandb = config['database']['annovar_humandb']
    return dirprefix, java, java_mem, snpEff, SnpSift, table_annovar, perl, Rscript, hg19_fa, annovar_humandb

def SnpEffGetTrans(snpEff, annovar_humandb, passvcf):
    cmd_gettrans = " ".join([
        java, '-jar', f'-jar -Xmx{java_mem}G',
        snpEff, 'ann hg19 -onlyTr', f'{annovar_humandb}/trans.txt',
        passvcf, f'-ud 3 -s {dirprefix}.summary.html',
        '>', f'{dirprefix}.snpeff.vcf'
    ])
    return cmd_gettrans

def SnpSiftFilter(SnpSift):
    cmd_filter = " ".join([
        java, '-jar', f'-jar -Xmx{java_mem}G',
        SnpSift, 'filter',
        "\"(ANN[0].EFFECT != 'intron_variant') & (ANN[0].EFFECT != 'intergenic_region') & (ANN[0].EFFECT != 'intragenic_variant') & (ANN[0].EFFECT != '3_prime_UTR_variant') & (ANN[0].EFFECT != '5_prime_UTR_variant') & (ANN[0].EFFECT != 'splice_region_variant') & (ANN[0].EFFECT != 'synonymous_variant')\"",
        f'{dirprefix}.snpeff.vcf', '>', f'{dirprefix}.rm_intron.vcf'
    ])
    return cmd_filter

def TableAnnovar(perl, table_annovar, annovar_humandb):
    cmd_tableannovar = " ".join([
        perl, table_annovar, f'{dirprefix}.rm_intron.vcf',
        annovar_humandb, '-buildver hg19',
        '-out', dirprefix, '-remove',
        '-protocol refGene,cytoBand,cosmic92_coding,clinvar_20210123,avsnp150,exac03,ALL.sites.2015_08,EAS.sites.2015_08',
        '-operation g,r,f,f,f,f,f,f', '--vcfinput'
    ])
    return cmd_tableannovar

def ExtractFields(SnpSift):
    cmd_extractfields = " ".join([
        java, '-jar', f'-jar -Xmx{java_mem}G',
        SnpSift, 'extractFields', '-s "," -e "."',
        f'{dirprefix}.hg19_multianno.vcf',
        'CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" GEN[0].AF "ANN[*].RANK" GEN[0].DP GEN[0].SB QUAL GEN[0].DP4[0] GEN[0].DP4[1] GEN[0].DP4[2] GEN[0].DP4[3] cosmic92_coding CLNDN CLNSIG avsnp150 ExAC_ALL ExAC_EAS ALL.sites.2015_08 EAS.sites.2015_08',
        '>', f'{dirprefix}.extractFields.xls'
    ])
    return cmd_extractfields

def Running(passvcf, bed, cfg, output_dir, output_prefix):
    global java, java_mem, dirprefix
    (dirprefix, java, java_mem, snpEff, SnpSift, table_annovar, perl, Rscript, hg19_fa, annovar_humandb) = Prepare(passvcf, bed, cfg, output_dir, output_prefix)

    cmd_gettrans = SnpEffGetTrans(snpEff, annovar_humandb, passvcf)
    cmd_filter = SnpSiftFilter(SnpSift)
    cmd_tableannovar = TableAnnovar(perl, table_annovar, annovar_humandb)
    cmd_extractfields = ExtractFields(SnpSift)

    bac.RunningProcess(cmd_gettrans, output_dir, output_prefix, 'Notice:\tsnpEff ann running messages, *.snpeff.vcf produced')
    bac.RunningProcess(cmd_filter, output_dir, output_prefix, 'Notice:\tSnpSift filter running messages, *.rm_intron.vcf produced')

    lines = os.popen(f"cat {dirprefix}.rm_intron.vcf |grep -v '#'").readlines()
    if len(lines) == 0:
        os.system(f"echo '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' > {dirprefix}.hg19_multianno.vcf")
        with open(f"{dirprefix}_final.xls", "w") as otopen:
            otopen.write("chr\tpos\tref\talt\tgene\ttype\ttrans\tHGVS_c\tHGVS_p\tAF\texon\tDP\tSB\tqual\tref_F\tref_R\talt_F\talt_R\tCOSMIC\tCLNDN\tCLNSIG\tavsnp150\tExAC_ALL\tExAC_EAS\t1000G_ALL\t1000G_EAS\tAF_background\tcount_background\ttotal_background\n")
        exit()

    bac.RunningProcess(cmd_tableannovar, output_dir, output_prefix, 'Notice:\ttable_annovar.pl annotation running messages, *.hg19_multianno.vcf produced')
    bac.RunningProcess(cmd_extractfields, output_dir, output_prefix, 'Notice:\tSnpSift extractFields running messages, *.extractFields.xls produced')
    # deal extractFields file
    os.system(f"sed -i 's/\\\\x3d/=/g' {dirprefix}.extractFields.xls")
    os.system(f"sed -i 's/\\\\x3b/:/g' {dirprefix}.extractFields.xls")
    os.system("cat %s.extractFields.xls |awk '{print \"%s\t\"$0}' |sed -e 1d > %s.query.xls"%(dirprefix, output_prefix, dirprefix))
    # query result add annovar background
    config = bac.ResolveConfig(cfg)
    annovar_background = config['database']['annovar_background']
    diff_T_N_v1 = config['scripts']['diff_T_N_v1']
    cmd_diff_T_N_v1 = " ".join([
        Rscript, diff_T_N_v1, f'{dirprefix}.query.xls', '0', annovar_background, f'{dirprefix}_final.xls'
    ])
    bac.RunningProcess(cmd_diff_T_N_v1, output_dir, output_prefix, 'Notice:\tquery result add annovar background;')
    # final background.xls
    # query_list=os.popen(f"cat {dirprefix}.query.xls")
    # os.system("cp %s.query.xls %s"%(dirprefix, annovar_background))
    # with open('%s/background_raw.xls'%annovar_background, 'w') as otopen:
    #     otopen.write("sample\tCHROM\tPOS\tREF\tALT\tGENE\tEFFECT\ttrans\tHGVS_C\tHGVS_P\tAF\tRANK\tDP\tSB\tQUAL\tref_F\tref_R\talt_F\talt_R\tcosmic92_coding\tCLNDN\tCLNSIG\tavsnp150\tEXAC_ALL\tEXAC_EAS\t1000G_ALL\t1000G_EAS")
    #     otopen.writelines(query_list)
    
    # make_background = config['scripts']['make_background']
    # cmd_make_background = " ".join([
    #     Rscript, make_background, '%s/background_raw.xls'%annovar_background, '%s/background.xls'%annovar_background
    # ])
    # 这一步background是将所有分析过的样本汇总并建立一个背景库，但是个人认为每次分析的样本类型、panel、试剂和分析目的均不同所以如果想做背景库再行设计方案整理数据也未尝不可；
    # 有报错信息
# Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
#   line 1 did not have 53 elements
# Calls: read.table -> scan
# Execution halted
    # bac.RunningProcess(cmd_make_background, output_dir, output_prefix, 'Notice:\tbackground.xls;')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing annotation')
    parser.add_argument('-v', '--passvcf', required=True,
                        help='"run_snvindel.py" script running produced *.PASS.vcf file or *.merged.PASS.vcf file')
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    Running(args.passvcf, args.bed, args.cfg, args.output_dir, args.output_prefix)
