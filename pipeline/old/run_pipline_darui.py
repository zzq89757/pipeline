#coding:utf8
import os
import sys
import json
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'version'))


def makedir(output_dir,string):
    try:
        os.mkdir(output_dir+'/'+string)
    except FileExistsError:
        pass

def sample_config_messages():
    pass

run_trim = '/data/home/zhounan/mission/yunkang_analysis/run_pipline/run_trim_yunkang.py'
run_mapping = '/data/home/zhounan/mission/yunkang_analysis/run_pipline/run_mapping_yunkang.py'
run_caller = '/data/home/zhounan/mission/yunkang_analysis/run_pipline/varscan_caller_tumoronly_yunkang.py'
python = '/data/software/bin/python'

def trim_running(run_trim,fq1,fq2,cfg_json,output_dir,output_prefix,log_dir):
    makedir(output_dir,'trim')
    output_dir=os.path.join(output_dir+'/trim/')
    trim_script=" ".join([python,run_trim,'--fq1',fq1,'--fq2',fq2,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix,'>>',os.path.join(log_dir,'%s.script.log'%output_prefix),'2>',os.path.join(log_dir,'%s.trim.log'%output_prefix)])
    os.system(trim_script)
    trimfq1=os.path.join(output_dir,'%s.fastp_trim.R1.fq'%output_prefix)
    trimfq2=os.path.join(output_dir,'%s.fastp_trim.R2.fq'%output_prefix)
    return trimfq1, trimfq2

def mapping_running(run_mapping,trimfq1,trimfq2,bed,cfg_json,output_dir,output_prefix,log_dir):
    makedir(output_dir,'mapping')
    output_dir=os.path.join(output_dir+'/mapping/')
    mapping_script=" ".join([python,run_mapping,'--trimfq1',trimfq1,'--trimfq2',trimfq2,'--bed',bed,'--cfg_json',cfg_json,'-o',output_dir,'-p',output_prefix,'>>',os.path.join(log_dir,'%s.script.log'%output_prefix),'2>',os.path.join(log_dir,'%s.mapping.log'%output_prefix)])
    os.system(mapping_script)
    tumorbam_markdup=os.path.join(output_dir,'%s.markdup_IndelRealigner.bam'%output_prefix)
    tumorbam_removedup=os.path.join(output_dir,'%s.removedup_IndelRealigner.bam'%output_prefix)
    return tumorbam_markdup, tumorbam_removedup

def varscan_running(run_caller,tumorbam,bed,cfg_json,output_dir,output_prefix,log_dir,namestr): #markdup, removedup
    makedir(output_dir,'varscan_%s'%namestr)
    output_dir=os.path.join(output_dir+'/varscan_%s/'%namestr)
    varscan_script=" ".join([python,run_caller,'--tumorbam',tumorbam,'--bed',bed,'--cfg_json',cfg_json,'--output_dir',output_dir,'--output_prefix',output_prefix,'>>',os.path.join(log_dir,'%s.script.log'%output_prefix),'2>',os.path.join(log_dir,'%s.varscan_%s.log'%(output_prefix, namestr) ) ] )
    os.system(varscan_script)
    raw_snp_vcf=os.path.join(output_dir,'%s.single.varscan.snp.vcf'%output_prefix)
    after_snp_vcf=os.path.join(output_dir,'%s.%s.varscan.snp.vcf'%(output_prefix, namestr))
    os.system('mv %s %s'%(raw_snp_vcf, after_snp_vcf))
    raw_indel_vcf=os.path.join(output_dir,'%s.single.varscan.indel.vcf'%output_prefix)
    after_indel_vcf=os.path.join(output_dir,'%s.%s.varscan.indel.vcf'%(output_prefix, namestr))
    os.system('mv %s %s'%(raw_indel_vcf, after_indel_vcf))
    return after_snp_vcf, after_indel_vcf

def fusion_running(tumorbam,bed,output_dir,output_prefix,log_dir,namestr):
    makedir(output_dir,'fusion_%s'%namestr)
    output_dir=os.path.join(output_dir+'/fusion_%s/'%namestr)
    fusion_script=" ".join(['/usr/bin/perl /data/home/zhounan/software/factera_v1.4.4/factera.pl','-o',output_dir,tumorbam,'/data/home/zhounan/database/ucsc/hg19.ncbiRefSeq.all_CDS_longest_transcript_id_reverse_correct_column4.bed','/data/home/zhounan/software/factera_v1.4.4/hg19.2bit',bed])
    os.system(fusion_script)
    os.system('echo "%s" >> %s'%(fusion_script, os.path.join(log_dir,'%s.script.log'%output_prefix)))
    result_txt=os.path.join(output_dir,'%s.%s_IndelRealigner.factera.fusions.txt'%(output_prefix,namestr))
    opt_list=[]
    msgs=os.popen('cat %s'%result_txt).readlines()
    if msgs != []:
        for i,line in enumerate(msgs):
            opt_list.extend(line)
            if i > 0:
                fusion_af=str( round( ( float(line.split('\t')[5]) + float(line.split('\t')[6]) ) / float(line.split('\t')[16])*100, 2) )
                opt_list.append('fusion_AF:\t{}%\n'.format(fusion_af))
        with open(os.path.join(output_dir,'%s.fusions_result.txt'%output_prefix),'w') as otopen:
            otopen.writelines(opt_list)

def running(fq1,fq2,bed,cfg_json,output_dir,output_prefix):
    makedir(output_dir,'log')
    log_dir=os.path.join(output_dir+'/log/')
    (trimfq1, trimfq2) = trim_running(run_trim,fq1,fq2,cfg_json,output_dir,output_prefix,log_dir)
    (tumorbam_markdup, tumorbam_removedup) = mapping_running(run_mapping,trimfq1,trimfq2,bed,cfg_json,output_dir,output_prefix,log_dir)
    (snp_vcf_markdup, indel_vcf_markdup) = varscan_running(run_caller,tumorbam_markdup,bed,cfg_json,output_dir,output_prefix,log_dir,'markdup')
    (snp_vcf_removedup, indel_vcf_removedup) = varscan_running(run_caller,tumorbam_removedup,bed,cfg_json,output_dir,output_prefix,log_dir,'removedup')
    for each in [snp_vcf_markdup, indel_vcf_markdup, snp_vcf_removedup, indel_vcf_removedup]:
        annovar_script="/usr/bin/perl /data/database/annovar/table_annovar.pl %s /data/database/annovar/humandb/ --buildver hg19 --outfile %s --remove --protocol refGene,refGeneWithVer,MT_ensGene,cytoBand,clinvar_20210123,cosmic70,cosmic92_coding,cosmic92_noncoding,gene4denovo201907,icgc28,avsnp150 --operation g,g,g,r,f,f,f,f,f,f,f --vcfinput --checkfile --nastring . > %s/%s.annovar_tmp.log 2> %s/%s.annovar_run.log"%(each, each, log_dir, output_prefix, log_dir, output_prefix)
        os.system(annovar_script)
        os.system('echo "%s" >> %s'%(annovar_script, os.path.join(log_dir,'%s.script.log'%output_prefix)))
    fusion_running(tumorbam_markdup,bed,output_dir,output_prefix,log_dir,'markdup')
    fusion_running(tumorbam_removedup,bed,output_dir,output_prefix,log_dir,'removedup')

if __name__ == '__main__':
    parser = ArgumentParser('running pipline')
    parser.add_argument('-s_cfg', '--running_configure', default=None,
                        help='pairsample running configure file; details like this')
    parser.add_argument('-v', '--version', default=None,
                        help='pipline lateset description')
    # parser.add_argument('-f1', '--fq1', required=True,
    #                     help='fastq file read1')
    # parser.add_argument('-f2', '--fq2', required=True,
    #                     help='fastq file read2')
    # parser.add_argument('-b', '--bed', required=True,
    #                     help='panel design bed file')
    # parser.add_argument('-c', '--cfg_json', required=True,
    #                     help='"configure.file.json" file')
    # parser.add_argument('-o', '--output_dir', required=True,
    #                     help='output directory')
    # parser.add_argument('-p', '--output_prefix',required=True,
    #                     help='sample id')
    args = parser.parse_args()
    if args.version is not None:
        print (darui_pipline_version_description.txt)
        exit()
    elif args.running_configure is not None:
        running(args.running_configure)
    # running(args.fq1, args.fq2, args.bed, args.cfg_json, args.output_dir, args.output_prefix)
