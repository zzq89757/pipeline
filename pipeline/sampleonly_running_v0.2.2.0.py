#coding:utf8
import os
import re
import sys
import json
import time
import subprocess
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def GetProcess(umi, running_step):
    allsteps = "trim\tfastqQC\tmapping\trmdup\tumi\treadligner\tbamQC\tsortbamsnvindel\tsnvindel\tannotation\tfusion\tcnv\tmsi\treport\t"
    if umi == 'false':
        allsteps = allsteps.replace('umi\t','')
    elif umi == 'true':
        allsteps = allsteps.replace('rmdup\t','')
    step_str = re.sub(r'\s+\t+|\s+|\t+|,\s+|,+|\.\s+|\.+', ' ', running_step.strip()).strip().upper()
    rm_step = re.findall(r'!(\w+)|(\w+)!', step_str)
    output_steps = allsteps.upper()
    if rm_step != []:
        for each in rm_step:
            each = ''.join(each)
            output_steps = output_steps.replace(each+'\t', '')
        output_steps = output_steps.strip().split('\t')
    else:
        output_steps = step_str.split(' ')
    if 'all'.upper() in step_str.upper():
        allsteps = allsteps.upper().strip().split('\t')
        return allsteps

    return output_steps

def OsSystemRun(cmd, step, output_dir, output_prefix, *check_files):
    full_file = f"ls {output_dir}/log/{output_prefix}_step/{output_prefix}.{step}_fullstep.txt" # only file return stdout or stderr, directory return none
    files = subprocess.Popen(full_file, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()[0] #stdout
    file_size = []
    for each in check_files:
        file_size.append(bac.GetSize(each))
    if cmd.upper() == 'SKIP':
        os.system(f'echo "{step} step skipped! at time: `date`" >> {output_dir}/log/{output_prefix}_step/{output_prefix}.skippedstep.txt')
    else:
        if files.decode() != '':
            print (f'"{step}" step allready exists!')
        elif 0 in file_size:
            print (f"{check_files} files size is 0!")
        elif files.decode() == '' and 0 not in file_size:
            answers = os.system(cmd)
            if answers == 0:
                os.system(f"echo '{step} step accomplished!' > {output_dir}/log/{output_prefix}_step/{output_prefix}.{step}_fullstep.txt")
            elif answers != 0:
                if step in ['trim', 'mapping', 'umi']:
                    print (f'the {step} step was wrong! now exit!')
                    exit()

def BuildCommands(fastq_read1, fastq_read2, bed, umi, umi_read1, umi_read2, UMI_type, cfg, output_dir, output_prefix):
    config = bac.ResolveConfig(cfg)
    command = config['scripts']
    python = config['language']['python']
    # trim
    cmd_trim = " ".join([ python, command['trim'], '--raw_read1', fastq_read1, '--raw_read2', fastq_read2, '--cfg',cfg, '--output_dir', output_dir+'/trim', '--output_prefix', output_prefix ])
    trimfq1 = output_dir+'/trim/'+output_prefix+'.trim.R1.fq.gz'
    trimfq2 = output_dir+'/trim/'+output_prefix+'.trim.R2.fq.gz'
    # fastq QC
    cmd_fastq_qc = " ".join([ python, command['fastqqc'], '--raw_fq1', fastq_read1, '--raw_fq2', fastq_read2, '--trim_fq1', trimfq1, '--trim_fq2', trimfq2, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--cutadapt_summary', output_dir+'/log/'+output_prefix+'.stdout_log.txt' ])
    # mapping
    cmd_mapping = " ".join([ python, command['mapping'], '--trimmed_fq1', trimfq1, '--trimmed_fq2', trimfq2, '--cfg', cfg, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    sort_bam = output_dir+'/mapping/'+output_prefix+'.sort.bam'
    # rmdup
    cmd_rmdup = " ".join([ python, command['rmdup'], '--input_bam', sort_bam, '--cfg', cfg, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    rmdup_bam = output_dir+'/mapping/'+output_prefix+'.rmdup.bam'
    # umi
    cmd_umi = " ".join([ python, command['umi'], '--sort_bam', sort_bam, '--umiseq_fq1', umi_read1, '--umiseq_fq2', umi_read2, '--cfg', cfg, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix, '--UMI_type', UMI_type ])
    group_bam = output_dir+'/mapping/'+output_prefix+'.group.bam'
    consensus_mapped_bam = output_dir+'/mapping/'+output_prefix+'.consensus.mapped.bam'
    # bam QC
    if umi == 'true':
        cmd_bam_qc = " ".join([ python, command['bamqc'], '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--sample_type ctDNA', '--umi true', '--group_bam', group_bam, '--consensus_mapped_bam', consensus_mapped_bam ])
        bamqc_checkfile = (group_bam, consensus_mapped_bam)
    elif umi == 'false':
        cmd_bam_qc = " ".join([ python, command['bamqc'], '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/QC', '--output_prefix', output_prefix, '--sample_type gDNA', '--umi false', '--sort_bam', sort_bam, '--removedup_bam', rmdup_bam ])
        bamqc_checkfile = (sort_bam, rmdup_bam)
    targetcoverage = output_dir+'/QC/'+output_prefix+'.Hs_pertarget_coverage_rmdup.xls'
    # indelrealigner
    cmd_indelrealigner = " ".join([ python, command['indelrealigner'], '--input_bam', bamqc_checkfile[1], '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/mapping', '--output_prefix', output_prefix ])
    indelrealigner_bam = output_dir+'/mapping/'+output_prefix+'.indelrealigner.bam'
    # snv_indel_sortbam
    cmd_snv_indel_sortbam = " ".join([ python, command['snvindel'], '--input_bam', sort_bam, '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/snv_indel_sortbam', '--output_prefix', output_prefix ])
    sortbam_multianno_vcf = output_dir+'/snv_indel_sortbam/'+output_prefix+'_all.hg19_multianno.vcf'
    # snv_indel
    cmd_snv_indel = " ".join([ python, command['snvindel'], '--input_bam', indelrealigner_bam, '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/snv_indel', '--output_prefix', output_prefix ])
    pass_vcf = output_dir+'/snv_indel/'+output_prefix+'.merged.PASS.vcf'
    # annotation
    cmd_annotation = " ".join([ python, command['annotate'], '--passvcf', pass_vcf, '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/annotation', '--output_prefix', output_prefix ])
    multianno_vcf = output_dir+'/annotation/'+output_prefix+'.hg19_multianno.vcf'
    # msi
    cmd_msi = " ".join([ python, command['msi'], 'detect', '--cfg', cfg, '--output_dir', f'{output_dir}/msi', '--output_prefix', output_prefix, config['parameter']['msi_basline'], indelrealigner_bam, '>', f'{output_dir}/msi/{output_prefix}.msi_result.tsv' ])
    # cnv
    localdir=os.path.dirname(os.path.abspath(__file__))
    cnv_baseline={ 'YK.P1.SNP.MSI.bed': f'{localdir}/../cnv/baselines_files/YK.P1.SNP.MSI_cnv_baseline_v3.txt', 'NanOnco_Plus_Panel_v2.0_Covered_hg19.bed': f'{localdir}/../cnv/baselines_files/NanOnco_Plus_Panel_v2.0_Covered_hg19_cnv_baseline_v3.txt' } #后面要把各个panel的cnv基线都建立好，按照命名链接对应的文件
    # cnv这里加入一个判断'--targetcoverage'这个参数不加入也能分析，因为如果不做QC也要能做cnv，这部分内容在command['cnv']脚本中做了修改
    if cnv_baseline.get(bed.split('/')[-1]) is not None:
        cmd_cnv = " ".join([ python, command['cnv'], '--bam', indelrealigner_bam, '--targetcoverage', targetcoverage, '--bed', bed, '--cfg', cfg, '--output_dir', output_dir+'/cnv', '--output_prefix', output_prefix, '--vcf', multianno_vcf ])
    else:
        cmd_cnv = f"echo '{bed} cnv baseline is not found!' "
    # fusion
    cmd_fusion = " ".join([ python, command['fusion'], '--sort_bam', sort_bam, '--cfg', cfg, '--output_dir', output_dir+'/fusion', '--output_prefix', output_prefix ])
    # report
    cmd_report = " ".join([ python, command['report'], '--multianno_vcf', multianno_vcf, '--cfg', cfg, '--output_dir', output_dir+'/report/'+output_prefix, '--output_prefix', output_prefix, '\n', python, command['report'], '--multianno_vcf', sortbam_multianno_vcf, '--cfg', cfg, '--output_dir', output_dir+'/report/'+output_prefix, '--output_prefix', output_prefix+'_sortbam' ])

    running_name = ['TRIM', 'FASTQQC', 'MAPPING', 'RMDUP', 'UMI', 'BAMQC', 'READLIGNER', 'SORTBAMSNVINDEL', 'SNVINDEL', 'ANNOTATION', 'FUSION', 'CNV', 'MSI', 'REPORT']
    running_process = [cmd_trim, cmd_fastq_qc, cmd_mapping, cmd_rmdup, cmd_umi, cmd_bam_qc, cmd_indelrealigner, cmd_snv_indel_sortbam,  cmd_snv_indel, cmd_annotation, cmd_fusion, cmd_cnv, cmd_msi, cmd_report]
    need_files = (trimfq1, trimfq2, sort_bam, bamqc_checkfile, targetcoverage, indelrealigner_bam, sortbam_multianno_vcf, pass_vcf, multianno_vcf)
    return running_name, running_process, need_files

def Running(fastq_read1, fastq_read2, bed, umi, umi_read1, umi_read2, UMI_type, cfg, output_dir, output_prefix, running_step):
# make directory
    dirlist = ['log', 'trim', 'QC', 'mapping', 'snv_indel_sortbam', 'snv_indel', 'annotation', 'fusion', 'cnv', 'msi', 'report']
    for each in dirlist:
        (tmp1, tmp2) = subprocess.Popen(f'mkdir {output_dir}/{each}', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    (tmp1, tmp2) = subprocess.Popen(f'mkdir {output_dir}/log/{output_prefix}_step', shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
# process command
    (running_name, running_process, need_files) = BuildCommands(fastq_read1, fastq_read2, bed, umi, umi_read1, umi_read2, UMI_type, cfg, output_dir, output_prefix)
    (trimfq1, trimfq2, sort_bam, bamqc_checkfile, targetcoverage, indelrealigner_bam, sortbam_multianno_vcf, pass_vcf, multianno_vcf) = need_files
# running process
    steps = GetProcess(umi, running_step)
    for x in range(len(running_name)):
        if running_name[x] not in steps:
            running_process[x] = 'skip'
    (cmd_trim, cmd_fastq_qc, cmd_mapping, cmd_rmdup, cmd_umi, cmd_bam_qc, cmd_indelrealigner, cmd_snv_indel_sortbam, cmd_snv_indel, cmd_annotation, cmd_fusion, cmd_cnv, cmd_msi, cmd_report) = running_process
# process running
    OsSystemRun(cmd_trim, 'trim', output_dir, output_prefix, fastq_read1, fastq_read2)
    pool=Pool(processes=2)
    pool.apply_async(OsSystemRun, args=(cmd_mapping, 'mapping', output_dir, output_prefix, trimfq1, trimfq2))
    pool.apply_async(OsSystemRun, args=(cmd_fastq_qc, 'fastqQC', output_dir, output_prefix, fastq_read1, fastq_read2, trimfq1, trimfq2))
    pool.close()
    pool.join()
    if umi == 'true':
        OsSystemRun(cmd_umi, 'umi', output_dir, output_prefix, sort_bam, umi_read1, umi_read2)
    elif umi == 'false':
        OsSystemRun(cmd_rmdup, 'rmdup', output_dir, output_prefix, sort_bam)
    OsSystemRun(cmd_indelrealigner, 'realinger', output_dir, output_prefix, bamqc_checkfile[1])
    pool=Pool(processes=4)
    pool.apply_async(OsSystemRun, args=(cmd_bam_qc, 'bamQC', output_dir, output_prefix, bamqc_checkfile[0], bamqc_checkfile[1]))
    pool.apply_async(OsSystemRun, args=(cmd_snv_indel, 'snvindel', output_dir, output_prefix, indelrealigner_bam))
    pool.apply_async(OsSystemRun, args=(cmd_fusion, 'fusion', output_dir, output_prefix, sort_bam))
    pool.apply_async(OsSystemRun, args=(cmd_snv_indel_sortbam, 'sortbamsnvindel', output_dir, output_prefix, sort_bam))
    pool.close()
    pool.join()
    OsSystemRun(cmd_annotation, 'annotation', output_dir, output_prefix, pass_vcf)
    if bac.GetSize(bamqc_checkfile[1]) != 0 or bac.GetSize(targetcoverage) != 0:
        # bamqc_checkfile[1] 和 targetcoverage 存在一者即可
        OsSystemRun(cmd_cnv, 'cnv', output_dir, output_prefix, multianno_vcf)
    else:
        print (f"{bamqc_checkfile[1]} and {targetcoverage} file are not exists!")
    OsSystemRun(cmd_msi, 'msi', output_dir, output_prefix, indelrealigner_bam)
# for report
    os.system(f"mkdir {output_dir}/report/{output_prefix}")
    os.system(f"cp {output_dir}/fusion/{output_prefix}.sv_hottable.xls {output_dir}/report/{output_prefix}")
    os.system(f"cp {output_dir}/cnv/{output_prefix}.cnv_allhot.xls {output_dir}/report/{output_prefix}")
    os.system(f"cp {output_dir}/msi/{output_prefix}.msi_result.tsv {output_dir}/report/{output_prefix}")
    OsSystemRun(cmd_report, 'report', output_dir, output_prefix, sortbam_multianno_vcf, multianno_vcf)


if __name__ == '__main__':
    localdir = os.path.dirname(os.path.abspath(__file__))
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
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
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', default=".",
                        help='output directory, default="."')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='output prefix, required')
    parser.add_argument('-s', '--running_step', type=str, default='all',
                        help="All analysis processes are as follows: 'trim' 'fastqQC' 'mapping' 'rmdup' 'umi' 'readligner' 'bamQC' 'sortbamsnvindel' 'snvindel' 'annotation' 'fusion' 'cnv' 'msi' 'report'\nif the name in the string, it means that the analysis process is added, like 'mapping,rmdup,readligner'\nif the name in the string and add '!' before the name, it means these process removed form all of the analysis process, like '!cnv,!fusion', if running_step is 'snvindel,!annotation' the snvindel process is ignore, default='all'")
                        #所有分析步骤都写在上面，如果umi是true则去掉rmdup分析步骤，如果umi是false则去掉umi分析步骤；都不加"-"表示只分析这些步骤，都加"-"表示将这些步骤从总体分析步骤中减去，减去分析步骤的算法不包括rmdup和umi，例如选择umi==true，running_step='rmdup'或='-rmdup'无效，选择umi==false，running_step='umi'或'-umi'无效。
    args = parser.parse_args()
    Running(args.fastq_read1, args.fastq_read2, args.bed, args.umi, args.umi_read1, args.umi_read2, args.UMI_type, args.cfg, args.output_dir, args.output_prefix, args.running_step)
