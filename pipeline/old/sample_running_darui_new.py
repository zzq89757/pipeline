#coding:utf8
import os
import sys
import json
import time
import subprocess
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def os_system_run(cmd, step, output_dir, output_prefix, *check_files):
    ls_file = 'ls %s/log/%s_step/%s.%s_step.txt'%(output_dir, output_prefix, output_prefix, step)
    files = subprocess.Popen(ls_file, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()[0] #stdout
    file_size = []
    for each in check_files:
        file_size.append(os.path.getsize(each))
    if files.decode() == '' and 0 not in file_size:
        answers = os.system(cmd)
        if answers == 0:
            os.system("echo '%s step successed!' > %s/log/%s_step/%s.%s_step.txt"%(step, output_dir, output_prefix, output_prefix, step))
        elif answers != 0:
            if step in ['trim', 'mapping', 'umi']:
                print ('the %s step was wrong! now exit!'%step)
                exit()
            else:
                pass
    else:
        print ('the "%s" step allready exists!\tor %s files size is 0!'%(step, check_files))
        # print (check_files)
        # print ('file size is 0!')

def running(fastq_read1, fastq_read2, bed, umi, umi_read1, umi_read2, UMI_type, cfg_json, output_dir, output_prefix):
    # dirprefix=os.path.join(output_dir, output_prefix)
    config_dict = json.load(open(cfg_json,'r'))
    command = config_dict['scripts']
    python = config_dict['language']['python']
    dirlist = ['log', 'trim', 'QC', 'mapping', 'snv_indel', 'annotation', 'fusion', 'cnv', 'msi']
    for each in dirlist:
        (tmp1, tmp2) = subprocess.Popen('mkdir %s/%s'%(output_dir, each), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    (tmp1, tmp2) = subprocess.Popen('mkdir %s/log/%s_step'%(output_dir, output_prefix), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    # trim
    cmd_trim = " ".join([ python, command['trim'], '--raw_read1', fastq_read1, '--raw_read2', fastq_read2, '--cfg_json',cfg_json, '--output_dir', output_dir+'/trim', '--output_prefix', output_prefix ])
    trimfq1 = output_dir+'/trim/'+output_prefix+'.trim.R1.fq.gz'
    trimfq2 = output_dir+'/trim/'+output_prefix+'.trim.R2.fq.gz'
    # fastq QC
    cmd_fastq_qc = " ".join([ python, command['fastqqc'], '--raw_fq1', fastq_read1, '--raw_fq2', fastq_read2, '--trim_fq1', trimfq1, '--trim_fq2', trimfq2, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--cutadapt_summary', output_dir+'/log/'+output_prefix+'.stdout_log.txt' ])
    # mapping
    cmd_mapping = " ".join([ python, command['mapping'], '--trimmed_fq1', trimfq1, '--trimmed_fq2', trimfq2, '--cfg_json', cfg_json, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    sort_bam = output_dir+'/mapping/'+output_prefix+'.sort.bam'
    # rmdup
    cmd_rmdup = " ".join([ python, command['rmdup'], '--input_bam', sort_bam, '--cfg_json', cfg_json, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    rmdup_bam = output_dir+'/mapping/'+output_prefix+'.rmdup.bam'
    # bam QC
    targetcoverage = output_dir+'/QC/'+output_prefix+'.Hs_pertarget_coverage_rmdup.xls'
    # indelrealigner
    cmd_indelrealigner = " ".join([ python, command['indelrealigner'], '--input_bam', 'useful_bam', '--bed', bed, '--cfg_json', cfg_json, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    indelrealigner_bam = output_dir+'/mapping/'+output_prefix+'.indelrealigner.bam'
    # umi
    cmd_umi = " ".join([ python, command['umi'], '--sort_bam', sort_bam, '--umiseq_fq1', umi_read1, '--umiseq_fq2', umi_read2, '--cfg_json', cfg_json, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix, '--UMI_type', UMI_type ])
    group_bam = output_dir+'/mapping/'+output_prefix+'.group.bam'
    consensus_mapped_bam = output_dir+'/mapping/'+output_prefix+'.consensus.mapped.bam'
    # snv_indel
    cmd_snv_indel = " ".join([ python, command['snvindel'], '--input_bam', indelrealigner_bam, '--bed', bed, '--cfg_json', cfg_json, '--output_dir', output_dir+'/snv_indel', '--output_prefix', output_prefix ])
    pass_vcf = output_dir+'/snv_indel/'+output_prefix+'.merged.PASS.vcf'
    # annotation
    cmd_annotation = " ".join([ python, command['annotate'], '--passvcf', pass_vcf, '--bed', bed, '--cfg_json', cfg_json, '--output_dir', output_dir+'/annotation', '--output_prefix', output_prefix ])
    multianno_vcf = output_dir+'/annotation/'+output_prefix+'.hg19_multianno.vcf'
    # msi
    cmd_msi = " ".join([ python, command['msi'], 'detect', '--output_dir', '%s/msi'%output_dir, '--output_prefix', output_prefix, config_dict['parameter']['msi_basline'], indelrealigner_bam, '>', '%s/msi/%s.msi_result.tsv'%(output_dir, output_prefix) ])
    # cnv
    cmd_cnv = " ".join([ python, command['cnv'], '--refbaseline', config_dict['parameter']['cnv_baseline'], '--targetcoverage', targetcoverage, '--output', output_dir+'/cnv', '--vcf', multianno_vcf ])
    # fusion

    # snv_indel add annotation
    cmd_var_add_ann = " ".join( [ cmd_snv_indel, ';', cmd_annotation ] )

    # all running process
    os_system_run(cmd_trim, 'trim', output_dir, output_prefix, fastq_read1, fastq_read2)
    pool=Pool(processes=2)
    pool.apply_async(os_system_run, args=(cmd_mapping, 'mapping', output_dir, output_prefix, trimfq1, trimfq2))
    pool.apply_async(os_system_run, args=(cmd_fastq_qc, 'fastq_qc', output_dir, output_prefix, fastq_read1, fastq_read2, trimfq1, trimfq2))
    pool.close()
    pool.join()
    if umi == 'true':
        os_system_run(cmd_umi, 'umi', output_dir, output_prefix, sort_bam, umi_read1, umi_read2)
        cmd_indelrealigner = cmd_indelrealigner.replace('useful_bam', consensus_mapped_bam)
        os_system_run(cmd_indelrealigner, 'realigner', output_dir, output_prefix, consensus_mapped_bam)
        if os.path.getsize(group_bam) == 0 or os.path.getsize(consensus_mapped_bam) == 0:
            cmd_bam_qc = "echo 'bam_qc need file not exists!' "
        else:
            cmd_bam_qc = " ".join([ python, command['bamqc'], '--bed', bed, '--cfg_json', cfg_json, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--sample_type ctDNA', '--umi true', '--group_bam', group_bam, '--consensus_mapped_bam', consensus_mapped_bam ])
    elif umi == 'false': #这里没有加入downsample的分析内容
        os_system_run(cmd_rmdup, 'rmdup', output_dir, output_prefix, sort_bam)
        cmd_indelrealigner = cmd_indelrealigner.replace('useful_bam', rmdup_bam)
        os_system_run(cmd_indelrealigner, 'realigner', output_dir, output_prefix, rmdup_bam)
        if os.path.getsize(sort_bam) == 0 or os.path.getsize(rmdup_bam) == 0:
            cmd_bam_qc = "echo 'bam_qc need file not exists!' "
        else:
            cmd_bam_qc = " ".join([ python, command['bamqc'], '--sort_bam', sort_bam, '--bed', bed, '--cfg_json', cfg_json, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--sample_type gDNA', '--removedup_bam', rmdup_bam, '--umi false' ])
    else:
        exit()
    pool=Pool(processes=2)
    pool.apply_async(os_system_run, args=(cmd_bam_qc, 'bam_qc', output_dir, output_prefix)) #分析流程不同需要的文件不同，所以这里把判断条件单独提出来在前面判断
    pool.apply_async(os_system_run, args=(cmd_var_add_ann, 'snv_indel_annotation', output_dir, output_prefix, indelrealigner_bam))
    # pool.apply_async(os_system_run, args=(cmd_snv_indel, 'snv_indel', output_dir, output_prefix, indelrealigner_bam))
    # pool.apply_async(os_system_run, args=(cmd_fusion, 'fusion', output_dir, output_prefix))
    pool.close()
    pool.join()
    # os_system_run(cmd_annotation, 'annotation', output_dir, output_prefix, pass_vcf)
    os_system_run(cmd_cnv, 'cnv', output_dir, output_prefix, multianno_vcf, targetcoverage)
    os_system_run(cmd_msi, 'msi', output_dir, output_prefix, indelrealigner_bam)
    # os.system('mkdir %s/result'%output_dir)
    # os.system('cp %s/snv_indel/%s.PASS.vcf %s/msi/%s.msi_result.tsv %s/result'%(output_dir, output_prefix, output_dir, output_prefix, output_dir ))

if __name__ == '__main__':
    localdir=os.path.dirname(os.path.abspath(__file__))
    parser = ArgumentParser('each pair sample pipline pyflow')
    parser.add_argument('-fq1', '--fastq_read1', required=True,
                        help='sample fastq read1 file, required')
    parser.add_argument('-fq2', '--fastq_read2', required=True,
                        help='sample fastq read2 file, required')
    parser.add_argument('-b', '--bed', type=str, default="%s/../bed/YK.P1.SNP.MSI.bed"%localdir,
                        help='panel design bed file, default=%s/../bed/YK.P1.SNP.MSI.bed'%localdir)
    parser.add_argument('-umi', '--umi', type=str, choices=['true', 'false'], default='false',
                        help='running UMI analysis pipline, default=false')
    parser.add_argument('-u1', '--umi_read1', default="",
                        help='umi read1 file, default=""')
    parser.add_argument('-u2', '--umi_read2', default="",
                        help='umi read2 file, default=""')
    parser.add_argument('-y', '--UMI_type', choices=['single_random', 'single_fix', 'duplex_random', 'duplex_fix'], default='duplex_random',
                        help='fixed umi type will correct umi sequence, and random can not; default=duplex_random')
    parser.add_argument('-c', '--cfg_json', type=str, default="%s/../configure/parameter_configure_file.json"%localdir,
                        help='"configure.file.json" file, default=%s/../configure/parameter_configure_file.json'%localdir)
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory, required')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='output prefix, required')

    args = parser.parse_args()
    running(args.fastq_read1, args.fastq_read2, args.bed, args.umi, args.umi_read1, args.umi_read2, args.UMI_type, args.cfg_json, args.output_dir, args.output_prefix)
