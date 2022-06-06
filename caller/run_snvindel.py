#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-22 
# Version: 0.2.0.0
# Description: pipeline for snv indel calling
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
    lofreq = software['lofreq']
    samtools = software['samtools']
    gatk4 = software['gatk4']
    bcftools = software['bcftools']
    hg19_fa = config['database']['hg19_fa']
    return dirprefix, java, java_mem, lofreq, samtools, gatk4, bcftools, hg19_fa

def indelqual(lofreq, dirprefix, hg19_fa, input_bam):
    cmd_indelqual = " ".join([
        lofreq, 'indelqual --dindel', '-f', hg19_fa, input_bam, '-o', '%s_indelqual.bam'%dirprefix
    ])
    return cmd_indelqual

def callparallel(lofreq, dirprefix, input_bam, hg19_fa, bed):
    cmd_callparallel = " ".join([
        lofreq, 'call-parallel', '--pp-threads ', '8',
        input_bam, '-f', hg19_fa, '-l', bed,
        '--call-indels', '-Q 25 -q 25 -a 1 -b 1 -A -B --no-default-filter',
        '-o', '%s.raw.vcf.gz'%dirprefix
    ])
    return cmd_callparallel

def variantfilter(java, java_mem, gatk4, dirprefix):
    cmd_variantfilter = " ".join([
        java, '-jar -Xmx%sG'%java_mem,
        gatk4, 'VariantFiltration',
        '--variant', '%s.raw.vcf.gz'%dirprefix,
        '--output', '%s.filtered.vcf.gz'%dirprefix,
        '-filter "(AF>=0.001&&DP>=500&&QUAL>=50&&SB<=60&&DP4[2]>=3&&DP4[3]>=3)&&((vc.hasAttribute(\'INDEL\')&&HRUN<=8))" --filter-name PASS',
        '-filter "AF<0.001" --filter-name low_AF',
        '-filter "DP<500" --filter-name low_DP',
        '-filter "QUAL<50" --filter-name low_QUAL',
        '-filter "DP4[2]<3||DP4[3]<3" --filter-name low_AD',
        '-filter "SB>60" --filter-name strand_bias',
        '-filter "(vc.hasAttribute(\'INDEL\')&&HRUN>8)" --filter-name indel_err'
    ])
    return cmd_variantfilter

def bcffilter(bcftools, dirprefix):
    cmd_bcffilter = " ".join([
        bcftools, 'filter',
        '-i \'filter="PASS"\' %s.filtered.vcf.gz'%dirprefix,
        '-o %s.PASS.vcf'%dirprefix
    ])
    return cmd_bcffilter

def running(input_bam, bed, cfg_json, output_dir, output_prefix):
    (dirprefix, java, java_mem, lofreq, samtools, gatk4, bcftools, hg19_fa) = prepare(input_bam, bed, cfg_json, output_dir, output_prefix)

    cmd_indelqual = indelqual(lofreq, dirprefix, hg19_fa, input_bam)
    indelqual_bam = '%s_indelqual.bam'%dirprefix
    cmd_callparallel = callparallel(lofreq, dirprefix, indelqual_bam, hg19_fa, bed)
    cmd_variantfilter = variantfilter(java, java_mem, gatk4, dirprefix)
    cmd_bcffilter = bcffilter(bcftools, dirprefix)

    runningprocessing(cmd_indelqual, output_dir, output_prefix, 'Notice:\tlofreq indelqual running messages, *_indelqual.bam produced')
    os.system('%s index %s %s_indelqual.bai'%(samtools, indelqual_bam, dirprefix))
    runningprocessing(cmd_callparallel, output_dir, output_prefix, 'Notice:\tlofreq call-parallel running messages, *.raw.vcf.gz produced')
    runningprocessing(cmd_variantfilter, output_dir, output_prefix, 'Notice:\tgatk4 VariantFiltration running messages, *.filtered.vcf.gz produced')
    runningprocessing(cmd_bcffilter, output_dir, output_prefix, 'Notice:\tbcfoolts filter running messages, *.PASS.vcf produced')

if __name__ == '__main__':
    cfgfile=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../configure/parameter_configure_file.json')
    parser = ArgumentParser('Doing snv indel')
    parser.add_argument('-b', '--input_bam', required=True,
                        help='input bam file')
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg_json', type=str, default=cfgfile,
                        help='"configure.file.json" file, default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.input_bam, args.bed, args.cfg_json, args.output_dir, args.output_prefix)
