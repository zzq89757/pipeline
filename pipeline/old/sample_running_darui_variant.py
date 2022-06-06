#coding:utf8
import os
import sys
import json
import time
from glob import glob
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def cmd_running(cmd, output_dir, output_prefix, log_str):
    os.system(cmd+' > %s/log/%s.%s.tmp.txt 2> %s/log/%s.%s.log.txt'%(output_dir, output_prefix, log_str, output_dir, output_prefix, log_str))

def prepare(dirprefix,output_dir,parameter_configure):
    os.system('cp %s %s.parameter_configure'%(parameter_configure, dirprefix))
    brc.resolveConfig('%s.parameter_configure'%dirprefix)
    cfg_json = '%s.parameter_configure.json'%dirprefix
    os.system('mkdir %s/log'%output_dir)
    config_dict = json.load(open(cfg_json,'r'))
    d = os.path.dirname(os.path.abspath(__file__))
    run_trim = d+'/'+'../trim/cutadapt_trim.py'
    run_mapping = d+'/'+'../mapping/run_only_mapping.py'
    downsample = d+'/'+'../mapping/bam_downsample.py'
    bam_add_RX = d+'/'+'../umi/bam_add_RX.pl'
    run_umi = d+'/'+'../umi/run_umi.sh'
    run_rmdup_realigner = d+'/'+'../mapping/run_rmdup_realigner.py'
    fastq_QC = d+'/'+'../qc/qc_fastq.py'
    bam_QC = d+'/'+'../qc/run_QC.py'
    lofreq_v1 = d+'/'+'../caller/lofreq_v1.sh'
    msp = d+'/'+'../msi/msp.py'
    msi_basline =  d+'/'+'../msi/nad_msi_sites.list_baseline'
    run_fusion = d+'/'+'../fusion/run_fusion.py'
    run_cnv = d+'/'+'../CNV/run_CNV.py'
    python = config_dict['python']['python_path']
    perl = config_dict['perl']['perl_path']
    hg19_fa = config_dict['hg19']['hg19_fa']
    return cfg_json, run_trim, run_mapping, downsample, bam_add_RX, run_umi, run_rmdup_realigner, fastq_QC, bam_QC, lofreq_v1, msp, msi_basline, python, perl, hg19_fa, run_fusion, run_cnv
# trim
def trimcommand(python,run_trim,fq1,fq2,cfg_json,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/trim')
    output_dir=os.path.join(output_dir+'/trim/')
    cmd_trim=" ".join([python,run_trim,'--fastq','%s,%s'%(fq1,fq2),'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix])
    trimfq1=os.path.join(output_dir,'%s.trim.R1.fq.gz'%output_prefix)
    trimfq2=os.path.join(output_dir,'%s.trim.R2.fq.gz'%output_prefix)
    return cmd_trim, trimfq1, trimfq2
# fastq file qc
def fastqqccommand(python,fastq_QC,fq1,fq2,trimfq1,trimfq2,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/QC')
    output_dir=os.path.join(output_dir+'/QC/')
    cmd_fastq_qc=" ".join([python,fastq_QC,'--raw_fq1',fq1,'--raw_fq2',fq2,'--trim_fq1',trimfq1,'--trim_fq2',trimfq2,'--output_dir',output_dir,'--output_prefix',output_prefix])
    return cmd_fastq_qc
# mapping
def mappingcommand(python,run_mapping,trimfq1,trimfq2,cfg_json,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/mapping')
    output_dir=os.path.join(output_dir+'/mapping/')
    cmd_mapping=" ".join([python,run_mapping,'--trimfq1',trimfq1,'--trimfq2',trimfq2,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix])
    sort_bam=os.path.join(output_dir,'%s.sort.bam'%output_prefix)
    return cmd_mapping, sort_bam
# downsample
def downsamplecommand(python,bam_downsample,input_bam,bed,cfg_json,output_dir,output_prefix,downsample_num=2000):
    output_dir=os.path.join(output_dir+'/mapping/')
    cmd_downsample=" ".join([python,bam_downsample,input_bam,bed,cfg_json,output_dir,output_prefix,str(downsample_num)])
    return cmd_downsample
# umi
def bamaddRX(perl,bam_add_RX,input_bam,u1,u2,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/umi')
    output_dir=os.path.join(output_dir+'/umi/')
    dirprefix=os.path.join(output_dir, output_prefix)
    cmd_addRX=" ".join([perl,bam_add_RX,'--bam',input_bam,'--UMI1',u1,'--UMI2',u2,'--out','%s.sort.withUMI.bam'%dirprefix])
    return cmd_addRX

def runningumi(run_umi,hg19_fa,output_dir,output_prefix,umi_seq):
    output_dir=os.path.join(output_dir+'/umi/')
    dirprefix=os.path.join(output_dir, output_prefix)
    cmd_umi=" ".join(['/usr/bin/bash',run_umi,'-s paired -t 5 -r',hg19_fa,'-y duplex_ramdom_UMI -d',output_dir,'-p',output_prefix,'-U',umi_seq,'-o',output_dir]) #修改了run_umi.sh的最后一步，将consensus.unmapped.bam拆分出fastq得到consensus.mapped.bam，unmapped.bam中的read Q值没有减去偏移值33，自己修改后减去偏移值，得到r1和r2的fastq同时用picard AddOrReplaceReadGroups重新比对genome，可以进行indelrealigner；结果在/data/home/zhounan/mission/KPI_assess_method_detect/deadline_October_pipline/real_sample_test_pipline/umi这个地址
    return cmd_umi
# 没有umi的rmdup
def removedup(python,run_rmdup_realigner,input_bam,bed,cfg_json,output_dir,output_prefix,dup='true'):
    output_dir=os.path.join(output_dir+'/mapping/')
    cmd_rmdup_realigner=" ".join([python,run_rmdup_realigner,'--input_bam',input_bam,'--bed',bed,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix,'--remove_dup',dup])
    return cmd_rmdup_realigner
# bam文件的qc分析
def bamqccommand(python,bam_QC,input_bam,bed,cfg_json,output_dir,output_prefix,umi=False,group_bam=None,consensus_mapped_bam=None):
    brc.createDirectory(output_dir+'/QC')
    output_dir=os.path.join(output_dir+'/QC/')
    if umi is True and group_bam is not None and consensus_mapped_bam is not None:
        cmd_QC=" ".join([python,bam_QC,'--sort_bam',input_bam,'--bed',bed,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix,'--umi',umi,'--group_bam',group_bam,'--consensus_mapped_bam',consensus_mapped_bam])
    else:
        cmd_QC=" ".join([python,bam_QC,'--sort_bam',input_bam,'--bed',bed,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix])
    return cmd_QC
# snv indel
def snvindelcommand(lofreq_v1,input_bam,output_dir,output_prefix,bed,name):
    brc.createDirectory(output_dir+'/%s'%name)
    output_dir=os.path.join(output_dir+'/%s/'%name)
    cmd_snvindel=" ".join(['/usr/bin/sh',lofreq_v1,output_dir,output_prefix,input_bam,bed])
    return cmd_snvindel
# msi
def msicommand(python,msp,input_bam,msi_basline,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/msi')
    output_dir=os.path.join(output_dir+'/msi/')
    cmd_msi=" ".join([python,msp,'detect',msi_basline,input_bam,'>','%s/%s.msi_result.tsv'%(output_dir,output_prefix)])
    return cmd_msi
# fusion
def fusioncommand(python,run_fusion,input_bam,output_dir,output_prefix):
    brc.createDirectory(output_dir+'/fusion')
    output_dir=os.path.join(output_dir+'/fusion/')
    cmd_fuion=" ".join([python,run_fusion,'--input_bam',input_bam,'--output_dir',output_dir,'--output_prefix',output_prefix])
    return cmd_fuion
# cnv
def cnvcommand():
    pass

def running(fastq_read1, fastq_read2, bed, umi, umi_read1, umi_read2, umi_seq, output_dir, output_prefix):
    dirprefix=os.path.join(output_dir, output_prefix)
    parameter_configure=os.path.dirname(os.path.abspath(__file__))+'/../configure/parameter_configure_file'
    (cfg_json, run_trim, run_mapping, downsample, bam_add_RX, run_umi, run_rmdup_realigner, fastq_QC, bam_QC, lofreq_v1, msp, msi_basline, python, perl, hg19_fa, run_fusion, run_cnv)=prepare(dirprefix, output_dir, parameter_configure)
    (cmd_trim, trimfq1, trimfq2)=trimcommand(python, run_trim, fastq_read1, fastq_read2, cfg_json, output_dir, output_prefix)
    # os.system(cmd_trim) #trim
    cmd_fastq_qc=fastqqccommand(python, fastq_QC, fastq_read1, fastq_read2, trimfq1, trimfq2, output_dir, output_prefix)
    # os.system(cmd_fastq_qc+' &') #fastq qc
    (cmd_mapping, sort_bam)=mappingcommand(python, run_mapping, trimfq1, trimfq2, cfg_json, output_dir, output_prefix)
    # os.system(cmd_mapping) #mapping
    brc.createDirectory(output_dir+'/bam')
    if umi == 'True':
        print ('UMI dedup')
        if umi_read1 is not None and umi_read2 is not None and umi_seq is not None:
            cmd_addRX=bamaddRX(perl, bam_add_RX, sort_bam, umi_read1, umi_read2, output_dir, output_prefix)
            os.system(cmd_addRX) #add RX tag
            cmd_umi=runningumi(run_umi, hg19_fa, output_dir, output_prefix, umi_seq)
            os.system(cmd_umi) #run umi
            consensus_mapped_bam=os.path.join(output_dir+'/umi/', output_prefix)+'.consensus.mapped.bam'
            consensus_mapped_bai=os.path.join(output_dir+'/umi/', output_prefix)+'.consensus.mapped.bai'
            os.system('cp %s %s/bam/%s.RAW.bam'%(consensus_mapped_bam, output_dir, output_prefix))
            os.system('cp %s %s/bam/%s.RAW.bai'%(consensus_mapped_bai, output_dir, output_prefix))
            # result_bam=os.path.join(output_dir+'/bam/', output_prefix)+'.RAW.bam'
            with open('%s/bam/%s.notice.txt'%(output_dir, output_prefix),'w') as otopen:
                otopen.write('%s/umi/%s.consensus.mapped.ba* -> %s/bam/%s.RAW.ba*\n'%(output_dir, output_prefix, output_dir, output_prefix))
            group_bam=os.path.join(output_dir+'/umi/', output_prefix)+'.group.bam'
            cmd_QC=bamqccommand(python, bam_QC, sort_bam, bed, cfg_json, output_dir, output_prefix, umi, group_bam, consensus_mapped_bam)
            os.system(cmd_QC+' &')
        else:
            pass
    elif umi == 'False':
        print ('normal dedup')
        cmd_downsample=downsamplecommand(python, downsample, sort_bam, bed, cfg_json, output_dir, output_prefix)
        # os.system(cmd_downsample) #downsample
        if glob(output_dir+'/mapping/*.sort.downsample_*.bam') != []:
            downsample_bam=glob(output_dir+'/mapping/*.sort.downsample_*.bam')[0]
            running_bam=downsample_bam
        else:
            running_bam=sort_bam
        cmd_rmdup_realigner=removedup(python, run_rmdup_realigner, running_bam, bed, cfg_json, output_dir, output_prefix)
        # os.system(cmd_rmdup_realigner)
        realigner_bam=os.path.join(output_dir+'/mapping/', output_prefix)+'.removedup_IndelRealigner.bam'
        realigner_bai=os.path.join(output_dir+'/mapping/', output_prefix)+'.removedup_IndelRealigner.bai'
        os.system('cp %s %s/bam/%s.RAW.bam'%(realigner_bam, output_dir, output_prefix))
        os.system('cp %s %s/bam/%s.RAW.bai'%(realigner_bai, output_dir, output_prefix))
        with open('%s/bam/%s.notice.txt'%(output_dir, output_prefix),'w') as otopen:
            otopen.write('%s/mapping/%s.removedup_IndelRealigner.ba* -> %s/bam/%s.RAW.ba*\n'%(output_dir, output_prefix, output_dir, output_prefix))
        cmd_QC=bamqccommand(python, bam_QC, sort_bam, bed, cfg_json, output_dir, output_prefix)
        os.system(cmd_QC+' &')
    RAW_bam='%s/bam/%s.RAW.bam'%(output_dir, output_prefix)
    cmd_snvindel_RAW=snvindelcommand(lofreq_v1, RAW_bam, output_dir, output_prefix, bed, 'snv_indel')
    cmd_snvindel_sort=snvindelcommand(lofreq_v1, sort_bam, output_dir, output_prefix, bed, 'snv_indel_sortbam')
    cmd_fusion=fusioncommand(python, run_fusion, RAW_bam, output_dir, output_prefix)
    cmd_msi=msicommand(python, msp, RAW_bam, msi_basline, output_dir, output_prefix)
    pool=Pool(processes=3)
    pool.apply_async(cmd_running, args=(cmd_snvindel_RAW, output_dir, output_prefix, 'snv_indel'))
    # pool.apply_async(cmd_running, args=(cmd_snvindel_sort, output_dir, output_prefix, 'snv_indel_sortbam'))
    pool.apply_async(cmd_running, args=(cmd_fusion, output_dir, output_prefix, 'fusion'))
    pool.close()
    pool.join()
    os.system(cmd_msi)

if __name__ == '__main__':
    localdir=os.path.dirname(os.path.abspath(__file__))
    parser = ArgumentParser('each pair sample pipline pyflow')
    parser.add_argument('-fq1', '--fastq_read1', required=True,
                        help='sample fastq read1 file')
    parser.add_argument('-fq2', '--fastq_read2', required=True,
                        help='sample fastq read2 file')
    parser.add_argument('-b', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-umi', '--umi', type=str, choices=['True', 'False'], default='False',
                        help='running UMI analysis pipline, default=False')
    parser.add_argument('-u1', '--umi_read1', required=None,
                        help='umi read1 file, default=None')
    parser.add_argument('-u2', '--umi_read2', required=None,
                        help='umi read2 file, default=None')
    parser.add_argument('-us', '--umi_seq', type=str, default="%s/../umi/umiseq"%localdir,
                        help='umi sequence file, default=umi dir path')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='output prefix')

    args = parser.parse_args()
    running(args.fastq_read1, args.fastq_read2, args.bed, args.umi, args.umi_read1, args.umi_read2, args.umi_seq, args.output_dir, args.output_prefix)
