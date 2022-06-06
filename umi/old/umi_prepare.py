#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-15
# Version: 0.2.0.0
# Description: pipline for trim read mapping and add UMI sequence
#########################################################################
import os
import sys
import pickle
import json
import subprocess
from itertools import zip_longest
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def check_file(trimmed_fq1, trimmed_fq2, RX_dict):
    if not brc.check_file(trimmed_fq1) or not brc.check_file(trimmed_fq2) or not brc.check_file(RX_dict):
        print ("%s or %s or %s not exist"%(trimmed_fq1, trimmed_fq2, RX_dict))
        exit()
    if brc.determainr1r2(trimmed_fq1, trimmed_fq2):
        exit()

def prepare(cfg_json, output_dir, output_prefix):
    dirprefix = os.path.join(output_dir, output_prefix)
    config=json.load(open(cfg_json,'r'))
    bwa=config['software']['bwa']
    bwa_threads=config['parameter']['bwa_threads']
    hg19_fa=config['database']['hg19_fa']
    return dirprefix, bwa, bwa_threads, hg19_fa

def printresult(deal_list, umidict):
    if deal_list == []:
        return []
    readname=deal_list[0].split('\t')[0]
    try:
        rx=umidict.pop(readname)
    except KeyError:
        rx='RX:Z:XXX-XXX'
    if len(deal_list) == 2: #[0] is read1, [1] is read2
        mq='MQ:i:%s'%deal_list[1].split('\t')[4]
        r1=deal_list[0]+'\t'+mq+'\t'+rx
        mq='MQ:i:%s'%deal_list[0].split('\t')[4]
        r2=deal_list[1]+'\t'+mq+'\t'+rx
        return [r1, r2]
    elif len(deal_list) > 2:
        print_list=[]
        mqlist=['MQ:i:0','MQ:i:0']
        for line in deal_list:
            flag=line.split('\t')[1]
            mq='MQ:i:%s'%line.split('\t')[4]
            if int(flag)%128 > 64 and int(flag) < 256: #判断为read1
                mqlist[0]=mq
            elif int(flag)%128 < 64 and int(flag) < 256: #判断为read2
                mqlist[1]=mq
        for line in deal_list:
            flag=line.split('\t')[1]
            if int(flag)%128 > 64:
                mq=mqlist[1]
            elif int(flag)%128 < 64:
                mq=mqlist[0]
            print_list.append(line+'\t'+mq+'\t'+rx)
        return print_list

def running(trimmed_fq1, trimmed_fq2, RX_dict, cfg_json, output_dir, output_prefix):
    check_file(trimmed_fq1, trimmed_fq2, RX_dict)
    (dirprefix, bwa, bwa_threads, hg19_fa) = prepare(cfg_json, output_dir, output_prefix)
    # umiseq_fq1和umiseq_fq2当中有些read在trimmed_fq1和trimmed_fq2中没有，是在cutadapt这一步去掉了这些read信息，所以二者不是一一对应的；
    # 所以这里选择的方法是建立一个存在read名称和tag标签的大的哈希表；虽然巨大的哈希表在调用key的时候会出错，但目前使用这个方法；（存在的问题）
    with open(RX_dict, 'rb') as inopen:
        umidict=pickle.load(inopen)
    cmd_mapping = " ".join([
        bwa, 'mem -t %s -Y -M -R'%bwa_threads,
        "'@RG\\tID:%s\\tSM:%s\\tPL:MGI\\tLB:%s\\tPE:100'"%(output_prefix, output_prefix, output_prefix),
        hg19_fa, trimmed_fq1, trimmed_fq2
    ])
    bwa_process=subprocess.Popen(cmd_mapping, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    #MQ:i:score Mapping quality of the mate/next segment.这里MQ给到的是其对应mate read的mapping quality，而不是他自己的mapping quality
    # bwa mem的输出不能直接接picard FixMateInformation，需要sort后的结果才可以
    # samtools view *.bwa_raw.bam |awk '{print $1}' |uniq -c |awk '{if($1==1)print}' #read都是成paired出现的，没有单独read出现两个位置；samtools view --threads 24 *.bwa_raw.bam |awk '{print $1}' |sort |uniq -c |awk '{if($1==1)print}'
    # 这里要同时考虑到多种情况：第一种，read1,read2；第二种多条read;
    readname=''
    deal_list=[]
    for line in bwa_process.stdout:
        line=bytes.decode(line).strip()
        if not line.startswith('@'): #body
            if line.split('\t')[0] != readname:
                result=printresult(deal_list, umidict)
                for each in result:
                    print (each)
                readname=line.split('\t')[0]
                deal_list=[line]
            elif line.split('\t')[0] == readname:
                deal_list.append(line)
        elif line.startswith('@'): #head
            print (line)
    result=printresult(deal_list, umidict)
    for each in result:
        print (each)

if __name__ == '__main__':
    parser = ArgumentParser('Doing mapping and add UMI sequence')
    parser.add_argument('-f1', '--trimmed_fq1', required=True,
                        help='trimmed fastq file read1')
    parser.add_argument('-f2', '--trimmed_fq2', required=True,
                        help='trimmed fastq file read2')
    parser.add_argument('-d', '--RX_dict', required=True,
                        help="UMI RX_dict file; it's file type is pickle")
    parser.add_argument('-c', '--cfg_json', required=True,
                        help='"configure.file.json" file')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    args = parser.parse_args()
    running(args.trimmed_fq1, args.trimmed_fq2, args.RX_dict, args.cfg_json, args.output_dir, args.output_prefix)
