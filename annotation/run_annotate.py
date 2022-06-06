#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-31 
# Version: 0.2.0.0
# Description: pipeline for snv indel calling result annotation
#########################################################################
import os
import sys
import json
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def runningprocessing(cmd, output_dir, output_prefix, other_string):
    (stdout, stderr, timer)=brc.run_command_with_return(cmd)
    brc.export_log(cmd, stdout, stderr, timer, output_dir, output_prefix, other_string)

def prepare(input_bam, bed, cfg_json, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = json.load(open(cfg_json,'r'))
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

def snpeffgettrans(java, java_mem, snpEff, annovar_humandb, dirprefix, passvcf):    cmd_gettrans = " ".join([
        java, '-jar -Xmx%sG'%java_mem,
        snpEff, 'ann hg19 -onlyTr', '%s/trans.txt'%annovar_humandb,
        passvcf, '-ud 3 -s %s.summary.html'%dirprefix,
        '>', '%s.snpeff.vcf'%dirprefix
    ])
    return cmd_gettrans

def snpsiftfilter(java, java_mem, SnpSift, dirprefix):
    cmd_filter = " ".join([
        java, '-jar -Xmx%sG'%java_mem,
        SnpSift, 'filter',
        "\"(ANN[0].EFFECT != 'intron_variant') & (ANN[0].EFFECT != 'intergenic_region') & (ANN[0].EFFECT != 'intragenic_variant') & (ANN[0].EFFECT != '3_prime_UTR_variant') & (ANN[0].EFFECT != '5_prime_UTR_variant') & (ANN[0].EFFECT != 'splice_region_variant')\"",
        '%s.snpeff.vcf'%dirprefix, '>', '%s.rm_intron.vcf'%dirprefix
    ])
    return cmd_filter

def tableannovar(perl, table_annovar, dirprefix, annovar_humandb):
    cmd_tableannovar = " ".join([
        perl, table_annovar, '%s.rm_intron.vcf'%dirprefix,
        annovar_humandb, '-buildver hg19',
        '-out %s'%dirprefix, '-remove',
        '-protocol refGene,cytoBand,cosmic92_coding,clinvar_20210123,avsnp150,exac03,ALL.sites.2015_08,EAS.sites.2015_08',
        '-operation g,r,f,f,f,f,f,f', '--vcfinput'
    ])
    return cmd_tableannovar

def extractfields(java, java_mem, SnpSift, dirprefix):
    cmd_extractfields = " ".join([
        java, '-jar -Xmx%sG'%java_mem,
        SnpSift, 'extractFields', '-s "," -e "."',
        '%s.hg19_multianno.vcf'%dirprefix,
        'CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" AF "ANN[*].RANK" DP SB QUAL DP4[0] DP4[1] DP4[2] DP4[3] cosmic92_coding CLNDN CLNSIG avsnp150 ExAC_ALL ExAC_EAS ALL.sites.2015_08 EAS.sites.2015_08',
        '>', '%s.extractFields.xls'%dirprefix
    ])
    return cmd_extractfields

def running(passvcf, bed, cfg_json, output_dir, output_prefix):
    (dirprefix, java, java_mem, snpEff, SnpSift, table_annovar, perl, Rscript, hg19_fa, annovar_humandb) = prepare(passvcf, bed, cfg_json, output_dir, output_prefix)

    cmd_gettrans = snpeffgettrans(java, java_mem, snpEff, annovar_humandb, dirprefix, passvcf)
    cmd_filter = snpsiftfilter(java, java_mem, SnpSift, dirprefix)
    cmd_tableannovar = tableannovar(perl, table_annovar, dirprefix, annovar_humandb)
    cmd_extractfields = extractfields(java, java_mem, SnpSift, dirprefix)

    runningprocessing(cmd_gettrans, output_dir, output_prefix, 'Notice:\tsnpEff ann running messages, *.snpeff.vcf produced')
    runningprocessing(cmd_filter, output_dir, output_prefix, 'Notice:\tSnpSift filter running messages, *.rm_intron.vcf produced')
    runningprocessing(cmd_tableannovar, output_dir, output_prefix, 'Notice:\ttable_annovar.pl annotation running messages, *.hg19_multianno.vcf produced')
    runningprocessing(cmd_extractfields, output_dir, output_prefix, 'Notice:\tSnpSift extractFields running messages, *.extractFields.xls produced')
    # deal extractFields file
    os.system("sed -i 's/\\\\x3d/=/g' %s.extractFields.xls"%dirprefix)
    os.system("sed -i 's/\\\\x3b/:/g' %s.extractFields.xls"%dirprefix)
    os.system("cat %s.extractFields.xls |awk '{print \"%s\t\"$0}' |sed -e 1d > %s.query.xls"%(dirprefix, output_prefix, dirprefix))
    # query result add annovar background
    config = json.load(open(cfg_json,'r'))
    annovar_background = config['database']['annovar_background']
    diff_T_N_v1 = config['scripts']['diff_T_N_v1']
    cmd_diff_T_N_v1 = " ".join([
        Rscript, diff_T_N_v1, '%s.query.xls'%dirprefix, '0', annovar_background, '%s_final.xls'%dirprefix
    ])
    runningprocessing(cmd_diff_T_N_v1, output_dir, output_prefix, 'Notice:\tquery result add annovar background;')
    # final background.xls
    # query_list=os.popen("cat %s.query.xls"%dirprefix)
    # os.system("cp %s.query.xls %s"%(dirprefix, annovar_background))
    # with open('%s/background_raw.xls'%annovar_background, 'w') as otopen:
    #     otopen.write("sample\tCHROM\tPOS\tREF\tALT\tGENE\tEFFECT\ttrans\tHGVS_C\tHGVS_P\tAF\tRANK\tDP\tSB\tQUAL\tref_F\tref_R\talt_F\talt_R\tcosmic92_coding\tCLNDN\tCLNSIG\tavsnp150\tEXAC_ALL\tEXAC_EAS\t1000G_ALL\t1000G_EAS")
    #     otopen.writelines(query_list)
    
    # make_background = config['scripts']['make_background']
    # cmd_make_background = " ".join([
    #     Rscript, make_background, '%s/background_raw.xls'%annovar_background, '%s/background.xls'%annovar_background
    # ])
    # 有报错信息
# Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
#   line 1 did not have 53 elements
# Calls: read.table -> scan
# Execution halted
    # runningprocessing(cmd_make_background, output_dir, output_prefix, 'Notice:\tbackground.xls;')

if __name__ == '__main__':
    cfgfile=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../configure/parameter_configure_file.json')
    parser = ArgumentParser('Doing annotation')
    parser.add_argument('-v', '--passvcf', required=True,
                        help='"run_snvindel.py" script running produced *.PASS.vcf file')
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg_json', type=str, default=cfgfile,
                        help='"configure.file.json" file, default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.passvcf, args.bed, args.cfg_json, args.output_dir, args.output_prefix)
