#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-10-22 
# Version: 0.2.0.0
# Description: pipeline for snv indel calling
#########################################################################
import os
import sys
import json
import subprocess
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def Prepare(cfg, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    config = bac.ResolveConfig(cfg)
    python = config['language']['python']
    perl = config['language']['perl']
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    software = config['software']
    lofreq = software['lofreq']
    samtools = software['samtools']
    gatk4 = software['gatk4']
    whatshap = software['whatshap']
    table_annovar = software['table_annovar']
    filter_to_pass = config['scripts']['filter_to_pass']
    merge_mnp = config['scripts']['merge_mnp']
    hg19_fa = config['database']['hg19_fa']
    annovar_humandb = config['database']['annovar_humandb']
    return dirprefix, python, perl, java, java_mem, lofreq, samtools, gatk4, whatshap, filter_to_pass, merge_mnp, hg19_fa, table_annovar, annovar_humandb

def IndelqualStep(lofreq, input_bam):
    cmd_indelqual = " ".join([
        lofreq, 'indelqual --dindel', '-f', hg19_fa, input_bam, '-o', f'{dirprefix}_indelqual.bam'
    ])
    return cmd_indelqual

def CallParallel(lofreq, input_bam, bed):
    cmd_callparallel = " ".join([
        lofreq, 'call-parallel', '--pp-threads ', '8',
        input_bam, '-f', hg19_fa, '-l', bed,
        '--call-indels', '-Q 25 -q 25 -a 1 -b 1 -A -B --no-default-filter',
        '-o', f'{dirprefix}.raw.vcf.gz'
    ])
    return cmd_callparallel

def TableAnnovar(perl, table_annovar, annovar_humandb):
    cmd_tableannovar = " ".join([
        perl, table_annovar, f'{dirprefix}.raw.vcf.gz',
        annovar_humandb, '-buildver hg19',
        '-out', f'{dirprefix}_all', '-remove',
        '-protocol refGene,cytoBand,cosmic92_coding,clinvar_20210123,avsnp150,exac03,ALL.sites.2015_08,EAS.sites.2015_08',
        '-operation g,r,f,f,f,f,f,f', '--vcfinput'
    ])
    return cmd_tableannovar

def VariantFilter(java, java_mem, gatk4):
    cmd_variantfilter = " ".join([
        java, '-jar -Xmx%sG'%java_mem,
        gatk4, 'VariantFiltration',
        '--variant', f'{dirprefix}.raw.vcf.gz',
        '--output', f'{dirprefix}.filtered.vcf.gz',
        '-filter "(AF>=0.001&&DP>=500&&QUAL>=50&&SB<=60&&DP4[2]>=3&&DP4[3]>=3)&&((vc.hasAttribute(\'INDEL\')&&HRUN<=8))" --filter-name PASS',
        '-filter "AF<0.001" --filter-name low_AF',
        '-filter "DP<500" --filter-name low_DP',
        '-filter "QUAL<50" --filter-name low_QUAL',
        '-filter "DP4[2]<3||DP4[3]<3" --filter-name low_AD',
        '-filter "SB>60" --filter-name strand_bias',
        '-filter "(vc.hasAttribute(\'INDEL\')&&HRUN>8)" --filter-name indel_err'
    ])
    return cmd_variantfilter

def FiltertoReady(filter_to_pass):
    cmd_filtertoready = " ".join([
        python, filter_to_pass, f'{dirprefix}.filtered.vcf.gz'
    ])
    return cmd_filtertoready

def WhatshapPhase(whatshap, indelqual_bam):
    cmd_whatshap = " ".join([
        whatshap, f'--debug phase --reference={hg19_fa}', f'{dirprefix}.readytophase.vcf', indelqual_bam, f'--indels --ignore-read-groups |bgzip --stdout --force --threads 16 > {dirprefix}.PASS.phased.vcf.gz'
    ])
    return cmd_whatshap

def MergeMnpStep(merge_mnp):
    cmd_mergemnp = " ".join([
        python, merge_mnp, f'{dirprefix}.PASS.phased.vcf.gz', hg19_fa, '--max_distance 10', "|awk '{if($1~/^#/){print $0}else if($7==\"PASS\")print $0}' >", f'{dirprefix}.merged.PASS.vcf'
    ])
    return cmd_mergemnp

def Running(input_bam, bed, cfg, output_dir, output_prefix):
    global python, dirprefix, hg19_fa
    # (dirprefix, python, java, java_mem, lofreq, samtools, gatk4, whatshap, filter_to_pass, merge_mnp, hg19_fa) = Prepare(cfg, output_dir, output_prefix)
    dirprefix, python, perl, java, java_mem, lofreq, samtools, gatk4, whatshap, filter_to_pass, merge_mnp, hg19_fa, table_annovar, annovar_humandb = Prepare(cfg, output_dir, output_prefix)

    cmd_indelqual = IndelqualStep(lofreq, input_bam)
    indelqual_bam = f'{dirprefix}_indelqual.bam'
    cmd_callparallel = CallParallel(lofreq, indelqual_bam, bed)
    cmd_tableannovar = TableAnnovar(perl, table_annovar, annovar_humandb)
    cmd_variantfilter = VariantFilter(java, java_mem, gatk4)
    # filtered_vcf_gz = f'{dirprefix}.filtered.vcf.gz'
    cmd_filtertoready = FiltertoReady(filter_to_pass)
    # readytophase_vcf = f'{dirprefix}.readytophase.vcf'
    cmd_whatshap = WhatshapPhase(whatshap, indelqual_bam)
    # pass_phased_vcf_gz = f'{dirprefix}.PASS.phased.vcf.gz'
    cmd_mergemnp = MergeMnpStep(merge_mnp)

    (tmp1, tmp2) = subprocess.Popen(f"rm {dirprefix}.raw.vcf.gz {dirprefix}_indelqual.ba*", shell = True,  stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate() # Cowardly refusing to overwrite already existing VCF output file ./snv_indel/*.raw.vcf.gz
    bac.RunningProcess(cmd_indelqual, output_dir, output_prefix, 'Notice:\tlofreq indelqual running messages, *_indelqual.bam produced')
    os.system(f'{samtools} index {indelqual_bam} {dirprefix}_indelqual.bai')
    bac.RunningProcess(cmd_callparallel, output_dir, output_prefix, 'Notice:\tlofreq call-parallel running messages, *.raw.vcf.gz produced')
    bac.RunningProcess(cmd_tableannovar, output_dir, output_prefix, 'Notice:\traw var vcf table_annovar running messages, *.hg19_multianno.vcf produced')
    bac.RunningProcess(cmd_variantfilter, output_dir, output_prefix, 'Notice:\tgatk4 VariantFiltration running messages, *.filtered.vcf.gz produced')
    bac.RunningProcess(cmd_filtertoready, output_dir, output_prefix, 'Notice:\t*.readytophase.vcf produced')
    bac.RunningProcess(cmd_whatshap, output_dir, output_prefix, 'Notice:\twhatshap running messages, *.PASS.phased.vcf.gz produced')
    os.system(f'tabix -p vcf {dirprefix}.PASS.phased.vcf.gz')
    bac.RunningProcess(cmd_mergemnp, output_dir, output_prefix, 'Notice:\tmerge_mnp.py running messages, *.merged.PASS.vcf produced')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing snv indel')
    parser.add_argument('-b', '--input_bam', required=True,
                        help='input bam file')
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    Running(args.input_bam, args.bed, args.cfg, args.output_dir, args.output_prefix)
