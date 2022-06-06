#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-9-3
# Version: 0.3.2.0
# Description: pipeline for fusion explore and analysis
#########################################################################

import os
import re
import sys
import gzip
import time
import json
import gc
import pysam
from itertools import zip_longest
from multiprocessing import Pool
from argparse import ArgumentParser

#######################################################
def basecontent(softseq, min_single_base_content):
    """计算mapping结果的soft序列的碱基组成，单一类型的碱基数是否小于等于min_single_base_content"""
    result=True
    for x in ['A','T','C','G','a','t','c','g']:
        if softseq.count(x)/float(len(softseq)) >= min_single_base_content:
            return False
    return result

def mismatchregular(readmsg, min_mismatch):
    """计算mapping结果的mismatch值是否小于等于min_mismatch"""
    try:
        mismatch_num = len( re.findall(r'[0-9][A-Z]', re.findall(r'MD:Z:.*?\t|OC:Z:.*?\t', readmsg)[0]) )
    except IndexError:
        mismatch_num = None #如果bam文件中的这一列没有MD这个标签则可能认为这条read没有比对上，则无法计算mismatch数；
    if mismatch_num <= min_mismatch:
        return True
    elif mismatch_num is None:
        return False
    else:
        return False

def insertregular(cigar, min_insert):
    """计算cigar的I数是否小于等于min_insert"""
    ins_list=re.findall(r'(\d+)I', cigar)
    insert_num=sum( list( map(int, ins_list) ) )
    if insert_num <= min_insert:
        return True
    else:
        return False

def deleteregular(cigar, min_delete):
    """计算cigar的D数是否小于等于min_delete"""
    del_list=re.findall(r'(\d+)D', cigar)
    delete_num=sum( list( map(int, del_list) ) )
    if delete_num <= min_delete:
        return True
    else:
        return False

def groups(cigar): #, min_bamsoft_length
    """计算cigar前后S的长度，如果S值小于11记为0，bowtie1软件不支持序列长度小于11的比对
    输出slenlist包含两个值，分别是cigar前后S的长度"""
    # 112M10S, 15S121M, 12S100M25S
    start_S=re.findall(r'^(\d+)S', cigar)
    end_S=re.findall(r'(\d+)S$', cigar)
    start=[int(x) for x in start_S]
    end=[int(x) for x in end_S]
    if start_S == []:
        start = [0]
    if end_S == []:
        end = [0]
    slenlist=start+end
    return slenlist

def softreadfasta(dirprefix,soft_readlines,min_bamsoft_length,min_mismatch,min_insert,min_delete,min_single_base_content):
    """将符合条件的S序列写成fasta格式，包括mismatch,insert,delete三个条件"""
    opt_list=[]
    for line in soft_readlines:
        col=line.split('\t')
        readname=col[0]
        flag=col[1]
        chrn=col[2]
        pos0=col[3]
        if int(pos0) < 500: #如果这个pos值小于500没有办法顺利得到位点前面一定长度的序列，所以这一步还是要去除一下
            continue
        cigar=col[5]
        misrst=mismatchregular(line, min_mismatch)
        insrst=insertregular(cigar, min_insert)
        delrst=deleteregular(cigar, min_delete)
        if not misrst or not insrst or not delrst: #只有三者返回值全是True时才向下进行
            continue
        readlen=len(col[9])
        slenlist=groups(cigar)
        # if slenlist[0] != 0 and slenlist[1] == 0:
        if slenlist[0] >= min_bamsoft_length and slenlist[1] < min_bamsoft_length: #slenlist 15S121M
            fastaid_s='%s:%s=%s_%s_%s_%s_%s_1S'%(chrn,pos0,readname,flag,chrn,pos0,cigar)
            seq_s=col[9][:slenlist[0]]
            if basecontent(seq_s, min_single_base_content):
                opt_list.append('>'+fastaid_s+'\n'+seq_s+'\n')
        # elif slenlist[0] == 0 and slenlist[1] != 0:
        elif slenlist[0] < min_bamsoft_length and slenlist[1] >= min_bamsoft_length: #slenlist 112M10S
            pos = int(pos0) + readlen - slenlist[1]
            fastaid_s='%s:%s=%s_%s_%s_%s_%s_2S'%(chrn,pos,readname,flag,chrn,pos0,cigar)
            seq_s=col[9][-slenlist[1]:]
            if basecontent(seq_s, min_single_base_content):
                opt_list.append('>'+fastaid_s+'\n'+seq_s+'\n')
        # elif slenlist[0] != 0 and slenlist[1] != 0:
        elif slenlist[0] >= min_bamsoft_length and slenlist[1] >= min_bamsoft_length: #slenlist 12S100M25S
            fastaid_s1='%s:%s=%s_%s_%s_%s_%s_1S'%(chrn,pos0,readname,flag,chrn,pos0,cigar)
            seq_s1=col[9][:slenlist[0]]
            if basecontent(seq_s1, min_single_base_content):
                opt_list.append('>'+fastaid_s1+'\n'+seq_s1+'\n')
            pos = int(pos0) + readlen - slenlist[0]
            fastaid_s2='%s:%s=%s_%s_%s_%s_%s_2S'%(chrn,pos,readname,flag,chrn,pos0,cigar)
            seq_s2=col[9][-slenlist[1]:]
            if basecontent(seq_s2, min_single_base_content):
                opt_list.append('>'+fastaid_s2+'\n'+seq_s2+'\n')
    with open('%s.softread.fasta'%dirprefix,'w') as otopen:
        otopen.writelines(opt_list)
    soft_fa='%s.softread.fasta'%dirprefix
    del opt_list
    gc.collect()
    return soft_fa

def bowtie1mapping(bowtie1, samtools, soft_fa, genome, dirprefix):
    """S序列的fasta文件，用bowtie1比对到genome，得到bam文件"""
    bowtie_cmd=" ".join([bowtie1, '--best --strata -v 2 -m 10 --sam --threads 16', genome, '-f', soft_fa, '|', samtools, 'view --threads 24 -b -S -F 4', '|', samtools, 'sort -n -o', '%s.bowtie1_softread.bam'%dirprefix]) #提示：长度小于等于10bp的read在bam文件中不显示
    os.system(bowtie_cmd)
    bowtie_bam='%s.bowtie1_softread.bam'%dirprefix
    return bowtie_bam
# **最终综合结果得到的SV类型存在信息**
# 1S read forward soft forward: {breakpoint1 -> chrn:pos+cigarM(bowtie1); breakpoint2 -> chrn:pos(bwa)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2, bk2+300)}
# 1S read reverse soft forward: {breakpoint1 -> chrn:pos+cigarM(bowtie1); breakpoint2 -> chrn:pos(bwa)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2, bk2+300)}
    # 1s rf sf bk1 chrn:pos+cigarM bk2 bwa
    # 1s rr sf bk1 chrn:pos+cigarM bk2 bwa
# 1S read reverse soft reverse: {breakpoint1 -> chrn:pos(bowtie1); breakpoint2 -> chrn:pos(bwa)} {fusion-genome	upstream -> (bk1, bk1+300).reverse_complement; downstream -> (bk2, bk2+300)}
# 1S read forward soft reverse: {breakpoint1 -> chrn:pos(bowtie1); breakpoint2 -> chrn:pos(bwa)} {fusion-genome	upstream -> (bk1, bk1+300).reverse_complement; downstream -> (bk2, bk2+300)}
    # 1s rr sr bk1 chrn:pos bk2 bwa
    # 1s rf sr bk1 chrn:pos bk2 bwa
# 2S read forward soft forward: {breakpoint1 -> chrn:pos+cigarM(bwa); breakpoint2 -> chrn:pos(bowtie1)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2, bk2+300)}
# 2S read reverse soft forward: {breakpoint1 -> chrn:pos+cigarM(bwa); breakpoint2 -> chrn:pos(bowtie1)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2, bk2+300)}
    # 2s rf sf bk1 bwa bk2 chrn:pos
    # 2s rr sf bk1 bwa bk2 chrn:pos
# 2S read reverse soft reverse: {breakpoint1 -> chrn:pos+cigarM(bwa); breakpoint2 -> chrn:pos+cigarM(bowtie1)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2-300, bk2).reverse_complement}
# 2S read forward soft reverse: {breakpoint1 -> chrn:pos+cigarM(bwa); breakpoint2 -> chrn:pos+cigarM(bowtie1)} {fusion-genome upstream -> (bk1-300, bk1); downstream -> (bk2-300, bk2).reverse_complement}
    # 2s rr sr bk1 bwa bk2 chrn:pos+cigarM
    # 2s rf sr bk1 bwa bk2 chrn:pos+cigarM
def svtypedefine(bowtie1bamline, break_point_stream_length, depth_extend):
    """通过read bwa比对的结果和soft seq的bowtie1比对结果，确定融合断点，融合基因组序列，及SV类型，计算samtools depth的区域；
    bkgenomedp_dict中所包含的信息：bk1上游\t序列形态\tbk2下游\t序列形态\tSV类型\t上游区域dp\t下游区域dp"""
    col=bowtie1bamline.split('\t')
    bwatmp=col[0].split('=')[0] #包括chrn:pos+cigarM和chrn:pos
    bwa_site=(bwatmp.split(':')[0], int(bwatmp.split(':')[1]))
    bowtie1_flag=col[1] # '0'=soft fowrard or '16'=soft reverse
    bowtie1_chrn=col[2]
    bowtie1_pos0=int(col[3])
    bowtie1_addlen=int(col[5][:-1])
    turnbk=re.findall(r'\dS$', col[0])[0] # 1S$ 2S$
    bkresult_dict={ # each condition breakpoint1 and breakpoint2
        ('1S','0'):( (bowtie1_chrn, bowtie1_pos0+bowtie1_addlen), bwa_site), 
        ('1S','16'):( (bowtie1_chrn, bowtie1_pos0), bwa_site), 
        ('2S','0'):(bwa_site, (bowtie1_chrn, bowtie1_pos0) ), 
        ('2S','16'):(bwa_site, (bowtie1_chrn, bowtie1_pos0+bowtie1_addlen) )
    }
    # bkpair=( (bk_raw.split(':')[0], int(bk_raw.split(':')[1]) ), (chrn, pos0) )
    bkpair=bkresult_dict.get((turnbk, bowtie1_flag)) # breakpoint1 breakpoint2
    # 设计思路和判别方法参考：最终综合结果得到的SV类型存在信息
    bkpairsite='%s:%s-%s:%s\t'%(bkpair[0][0], str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]))

    genome_1S_0='%s:%s-%s\trf\t%s:%s-%s\trf\t%s'%( bkpair[0][0], str(bkpair[0][1]-break_point_stream_length+1), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+break_point_stream_length-1), ('deletion' if bkpair[0][0] == bkpair[1][0] else 'translocation') )
    dp_1S_0='\t%s_%s_%s\t%s_%s_%s'%( bkpair[0][0], str(bkpair[0][1]-depth_extend), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+depth_extend) )

    genome_1S_16='%s:%s-%s\trc\t%s:%s-%s\trf\t%s'%( bkpair[0][0], str(bkpair[0][1]), str(bkpair[0][1]+break_point_stream_length-1), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+break_point_stream_length-1), ('inverstion' if bkpair[0][0] == bkpair[1][0] else 'translocation') )
    dp_1S_16='\t%s_%s_%s\t%s_%s_%s'%( bkpair[0][0], str(bkpair[0][1]), str(bkpair[0][1]+depth_extend), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+depth_extend) )

    genome_2S_0='%s:%s-%s\trf\t%s:%s-%s\trf\t%s'%( bkpair[0][0], str(bkpair[0][1]-break_point_stream_length+1), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+break_point_stream_length-1), ('deletion' if bkpair[0][0] == bkpair[1][0] else 'translocation') )
    dp_2S_0='\t%s_%s_%s\t%s_%s_%s'%( bkpair[0][0], str(bkpair[0][1]-depth_extend), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]), str(bkpair[1][1]+depth_extend) )

    genome_2S_16='%s:%s-%s\trf\t%s:%s-%s\trc\t%s'%( bkpair[0][0], str(bkpair[0][1]-break_point_stream_length+1), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]-break_point_stream_length+1), str(bkpair[1][1]), ('inverstion' if bkpair[0][0] == bkpair[1][0] else 'translocation') )
    dp_2S_16='\t%s_%s_%s\t%s_%s_%s'%( bkpair[0][0], str(bkpair[0][1]-depth_extend), str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]-depth_extend), str(bkpair[1][1]) )

    bkgenomedp_dict={
        ('1S','0'):bkpairsite+genome_1S_0+dp_1S_0, 
        ('1S','16'):bkpairsite+genome_1S_16+dp_1S_16, 
        ('2S','0'):bkpairsite+genome_2S_0+dp_2S_0, 
        ('2S','16'):bkpairsite+genome_2S_16+dp_2S_16
    }
    bkgenedp=bkgenomedp_dict.get((turnbk, bowtie1_flag)) #上下游基因组序列位置，breakpoint位点计算深度的depth_extend区域
    return bkpair, bkgenedp

def joinsite(bowtie1bamlist):
    """对于同一个bk位点其softread能够比对多位置，选一个最长的soft序列其他位置不保留；
    这一步可以规划为都保留各个bk组合，都输出各个组合的M及S的长度，之后先按照bk1位点相同排序各个bk2选择soft最长的结果，再按照bk2位点相同排序各个bk1选择soft最长的结果，同时这一步可以输出最长的S的序列，通过找到原始的bam文件找到最长S序列的read的M序列"""
    otsitelist=[]
    tmpbk=[]
    for i in range(len(bowtie1bamlist)):
        current_bk=bowtie1bamlist[i].split('\t')[0].split('=')[0]+'_'+bowtie1bamlist[i].split('\t')[0].split('_')[-1]
        current_mlen=int(bowtie1bamlist[i].split('\t')[5][:-1])
        front_bk=bowtie1bamlist[i-1].split('\t')[0].split('=')[0]+'_'+bowtie1bamlist[i-1].split('\t')[0].split('_')[-1]
        front_mlen=int(bowtie1bamlist[i-1].split('\t')[5][:-1])
        if current_bk == front_bk:
            if current_mlen > tmpmlen:
                tmpbk.append(bowtie1bamlist[i])
                tmpmlen=current_mlen
        else:
            if tmpbk != []: # 第一个重复的先加进来，如果后面发现更大的值则替换，如果后面没有更大的值则tmpbk为空
                otsitelist[-1]=tmpbk[-1]
            otsitelist.append(bowtie1bamlist[i])
            tmpbk=[]
            tmpmlen=current_mlen
    return otsitelist

def numpooling(bkpairpooldict, bk_extend_size):
    newbkpairpooldict={}
    for eachnum in bkpairpooldict.values():
        for bkpair,bkgenedp in eachnum.items():
            numpool='%s:%s-%s:%s'%( bkpair[0][0], str(bkpair[0][1])[:2]+'L'+str(len(str(bkpair[0][1]))), bkpair[1][0], str(bkpair[1][1])[:2]+'L'+str(len(str(bkpair[1][1]))) )
            if newbkpairpooldict.get(numpool) is not None:
                for pairs in list(newbkpairpooldict.get(numpool).keys()):
                    bk1range = (bkpair[0][1] >= pairs[0][1] - bk_extend_size) and (bkpair[0][1] <= pairs[0][1] + bk_extend_size)
                    bk2range = (bkpair[1][1] >= pairs[1][1] - bk_extend_size) and (bkpair[1][1] <= pairs[1][1] + bk_extend_size)
                    if bkpair[0][0] == pairs[0][0] and bk1range and bkpair[1][0] == pairs[1][0] and bk2range:
                        continue
                    else:
                        newbkpairpooldict.get(numpool).update( {bkpair:bkgenedp} )
            else:
                newbkpairpooldict[numpool]={bkpair:bkgenedp}
    return newbkpairpooldict

def bkannotation(pairmsg, bkannodict):
    bk1anno='None'
    bk2anno='None'
    bk1=pairmsg.split('-')[0] #chr1:22103-chr19:41736913
    bk1key1=bk1.split(':')[0]+':'+bk1.split(':')[1][:2]+'L'+str(len(bk1.split(':')[1]))
    if bkannodict.get(bk1key1) is not None:
        for region in list(bkannodict.get(bk1key1).keys()):
            if int(bk1.split(':')[1]) >= region[0] and int(bk1.split(':')[1]) <= region[1]:
                bk1anno=bkannodict.get(bk1key1).get(region)
                break
    bk2=pairmsg.split('-')[1]
    bk2key1=bk2.split(':')[0]+':'+bk2.split(':')[1][:2]+'L'+str(len(bk2.split(':')[1]))
    if bkannodict.get(bk2key1) is not None:
        for region in list(bkannodict.get(bk2key1).keys()):
            if int(bk2.split(':')[1]) >= region[0] and int(bk2.split(':')[1]) <= region[1]:
                bk2anno=bkannodict.get(bk2key1).get(region)
                break
    optresult='%s\t%s-%s\n'%(pairmsg, bk1anno, bk2anno)
    return optresult

def breakpointclassification(bowtie1bamlist, bk_extend_size, dirprefix, break_point_stream_length, depth_extend, bkannodict):
    """通过得到bowtie1比对的bam文件，找出融合的bk1和bk2断点，通过bkpairpooldict对重复或邻近的断点位置进行合并，这里合并并没有讲求其断点位置信息是否是rc的"""
    test_all_site_list=[]
    bkpairpooldict={} #{'chr1:23L7-chr2:29L8':{(('chr1', 2340161),('chr2', 29684468)):'chr1:2340161-chr2:29684468',(('chr1', 2340177),('chr2', 29684488)):'chr1:2340177-chr2:29684488'}}
    for line in bowtie1bamlist:
        col=line.split('\t')
        ids=col[0]
        (bkpair, bkgenedp)=svtypedefine(line, break_point_stream_length, depth_extend)
        #breakpoint1 breakpoint2
        numpool='%s:%s-%s:%s'%( bkpair[0][0], str(bkpair[0][1])[:2]+'L'+str(len(str(bkpair[0][1]))), bkpair[1][0], str(bkpair[1][1])[:2]+'L'+str(len(str(bkpair[1][1]))) )
        # bkstr='%s:%s-%s:%s'%( bkpair[0][0], str(bkpair[0][1]), bkpair[1][0], str(bkpair[1][1]))
        test_all_site_list.append(bkgenedp+'\t'+ids.split('=')[1].split('_')[0]+'\n')
        if bkpairpooldict.get(numpool) is not None:
            for pairs in list(bkpairpooldict.get(numpool).keys()): #(('chr1', 2340161),('chr2', 29684468)):'chr1:2340161-chr2:29684468'
                bk1range = (bkpair[0][1] >= pairs[0][1] - bk_extend_size) and (bkpair[0][1] <= pairs[0][1] + bk_extend_size)
                bk2range = (bkpair[1][1] >= pairs[1][1] - bk_extend_size) and (bkpair[1][1] <= pairs[1][1] + bk_extend_size)
                if bkpair[0][0] == pairs[0][0] and bk1range and bkpair[1][0] == pairs[1][0] and bk2range:
                    continue
                else:
                    bkpairpooldict.get(numpool).update( {bkpair:bkgenedp} )
        else:
            bkpairpooldict[numpool]={bkpair:bkgenedp}

    bkpairpooldict=numpooling(bkpairpooldict, bk_extend_size)
    with open('%s.test_bksite_andread.txt'%dirprefix,'w') as otopen:
        otopen.writelines(test_all_site_list)
    bkresult_list=[]
    addannotation=[]
    for each in bkpairpooldict.values():
        for site in each.values():
            bkresult_list.append(site+'\n')
            #在这里加入annotation注释结果并生成一个新的文件
            pairmsg=site.split('\t')[0]
            anno=bkannotation(pairmsg, bkannodict)
            addannotation.append(anno)
    with open('%s.bkresult.txt'%dirprefix,'w') as otopen:
        otopen.writelines(bkresult_list)
    with open('%s.bk_annotate.txt'%dirprefix,'w') as otopen:
        otopen.writelines(addannotation)
    del bkpairpooldict
    gc.collect()

def dealannobed(break_point_annotation_bed, dirprefix):
    bkannodict={}
    ot={}
    with open(break_point_annotation_bed,'r') as inopen:
        for line in inopen:
            col=line.split('\t')
            key1=col[0]+':'+str(col[1])[:2]+'L'+str(len(col[1]))
            key2=(int(col[1])+1, int(col[2]))
            link='%s\t%s'%(str(int(col[1])+1), col[2])
            value=col[3]
            if bkannodict.get(key1) is None:
                bkannodict.update({key1:{key2:value}})
                ot.update({key1:{link:value}})
            elif bkannodict.get(key1) is not None: 
                bkannodict.get(key1).update({key2:value})
                ot.get(key1).update({link:value})
    with open('%s.bkannodict.json'%dirprefix,'w') as otopen:
        otopen.write(json.dumps(ot))
    return bkannodict

def buildingbreakpoint(input_bam, output_dir, output_prefix, samtools, bed, bowtie1, genome, min_bamsoft_length, min_mismatch, min_insert, min_delete, bk_extend_size, break_point_stream_length, depth_extend, min_single_base_content, break_point_annotation_bed):
    """找到融合断点"""
    chrstring="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    try:
        os.mkdir('%s/01softread'%output_dir)
    except FileExistsError:
        pass
    output_dir01=output_dir+'/01softread'
    dirprefix=os.path.join(output_dir01, output_prefix)
    region_cmd=" ".join( [samtools, "view -F 780 -L bed --threads 24", input_bam, chrstring, "|awk '{if($6~/S/)print}'"] ) # flag 780 means: read unmapped (0x4), mate unmapped (0x8), not primary alignment (0x100), read fails platform/vendor quality checks (0x200)
    if bed is not None and os.path.exists(bed):
        region_cmd=region_cmd.replace('bed',bed)
    else:
        region_cmd=region_cmd.replace('-L bed ','')
    soft_readlines=os.popen(region_cmd).readlines()
    if soft_readlines == []:
        print ("%s not have index file"%input_bam)
        exit()
    soft_fa=softreadfasta(dirprefix, soft_readlines, min_bamsoft_length, min_mismatch, min_insert, min_delete, min_single_base_content)
    del soft_readlines
    gc.collect()
    bowtie_bam=bowtie1mapping(bowtie1, samtools, soft_fa, genome, dirprefix)
    bowtie_bam='%s.bowtie1_softread.bam'%dirprefix
    bowtie1_cmd=" ".join([samtools, "view -F 4 --threads 24", bowtie_bam, "|awk '{if($3!~/chrM/&&$3!~/_/)print}'"]) #chrstring
    bowtie1bamlist=os.popen(bowtie1_cmd).readlines()
    bkannodict=dealannobed(break_point_annotation_bed, dirprefix)
    bowtie1bamlist=joinsite(bowtie1bamlist)
    breakpointclassification(bowtie1bamlist, bk_extend_size, dirprefix, break_point_stream_length, depth_extend, bkannodict)

# 第一部分 寻找潜在的breakpoint分析内容结束：主要包括buildingbreakpoint主函数，用于确定断点信息及SV类型的svtypedefine函数，合并临近位点breakpointclassification函数，提取soft序列用于bowtie1比对的softreadfasta函数；
#######################################################
def allcrossfastq(samtools, output_dir, output_prefix, input_bam, break_point_stream_length):
    """在原始bam文件中找到可能跨过位点的read"""
    chrstring="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    try:
        os.mkdir('%s/01softread'%output_dir)
    except FileExistsError:
        pass
    output_dir01=output_dir+'/01softread'
    dirprefix=os.path.join(output_dir01, output_prefix)
    fusionread_cmd=" ".join( [samtools, "view --threads 24 -h -F 780", input_bam, chrstring, "|awk '{if($1~/^@/||$6~/S/||$9==0||$9>=%d||$9<=-%d)print $0}'"%(break_point_stream_length, break_point_stream_length), "|", samtools, "view --threads 24 -b -S -o", "%s.crossgap.bam"%dirprefix] ) # flag 780 means: read unmapped (0x4), mate unmapped (0x8), not primary alignment (0x100), read fails platform/vendor quality checks (0x200)
    os.system(fusionread_cmd)
    """crossgap.bam当中得到的read1和read2并不是对应的，在后面用bwa mem比对融合基因组后会由于read位置不对应而报错([mem_sam_pe] paired reads have different names:)，所以这里不区分read1和read2而放到一个整体当中"""
    bamtofastq_cmd=" ".join( [samtools, "fastq --threads 24 -c 4 -o", "%s.partone.fq.gz"%dirprefix, "-N %s.crossgap.bam"%dirprefix] )
    os.system(bamtofastq_cmd)
    """注意partone.fq.gz中得到的所有read数并不是crossgap.bam当中的所有read，事实上缺少了一部分flag为256的secondary read，需要通过read名称把这些read找回来"""
    bamread=[x.split('\t')[0] for x in os.popen("%s view --threads 24 %s.crossgap.bam"%(samtools, dirprefix) )]
    partone=map(lambda x:x[1:-3], os.popen("zcat %s.partone.fq.gz |awk '{if(NR%%4==1)print}'"%dirprefix))
    diff_read_name=[x+'\n' for x in list(set(bamread)-set(partone))] #没有发现不同的reads
    if diff_read_name != []:
        with open('%s.diff.readname'%dirprefix,'w') as otopen:
            otopen.writelines(diff_read_name)
        diff_reads=os.popen("%s view --threads 24 -N %s.diff.readname %s.crossgap.bam"%(samtools, dirprefix, dirprefix)).readlines()
        # bam文件read flag值余128得到结果大于64是read1，read flag值余128小于64是read2;
        # bam文件read flag值余32得到结果大于16是reverse，read flag值余32小于16是forward;
        # 这里不应该分出read1和read2，因为这里的read1和read2不是配对的，他们不是对应关系
        # 并且如果一条read(r1或r2)出现forward mapping和reverse mapping则只取forward
        diff_reads_fq=[]
        for line in diff_reads:
            col=line.split('\t')
            readname=col[0]
            flag=col[1]
            seq=col[9]
            qual=col[10]
            if int(flag)%32 > 16 : #如果是forward则不动，如果是reverse则需要序列反向互补，质量值字符倒序
                seq=seq.replace('A','X').replace('C','Y').replace('T','A').replace('G','C').replace('X','T').replace('Y','G').replace('a','x').replace('c','y').replace('t','a').replace('g','c').replace('x','t').replace('y','g')[::-1]
                qual=qual[::-1]
            msg=[readname, seq, '+', qual]
            ids='@'+readname
            msg[0]=ids
            if int(flag)%128 > 64: #判断为read1
                ids='@'+readname+'/1'
                msg[0]=ids
            elif int(flag)%128 < 64: #判断为read2
                ids='@'+readname+'/2'
                msg[0]=ids
            diff_reads_fq.append("\n".join(msg)+'\n')
        otreads=list( set ( diff_reads_fq ) )
        codelines=[x.encode() for x in otreads] #写入gz文件中需要b''形式的字符串
        with gzip.open('%s.diffread.fq.gz'%dirprefix,'wb') as otopen:
            otopen.writelines(codelines)
        os.system("/usr/bin/cat %s.diffread.fq.gz %s.partone.fq.gz > %s.cross.fq.gz"%(dirprefix, dirprefix, dirprefix))
    elif diff_read_name == []:
        os.system("/usr/bin/cp %s.partone.fq.gz %s.cross.fq.gz"%(dirprefix, dirprefix))
    # os.system("/usr/bin/rm %s.diffread.fq.gz %s.partone.fq.gz %s.diff.readname "%(dirprefix, dirprefix, dirprefix))
    del bamread
    gc.collect()

# 第二部分 寻找潜在的能够跨越断点的softread或pairread分析内容结束：仅包括allcrossfastq主函数，寻找跨越断点read信息的逻辑依据是fusionread_cmd方法；
#######################################################
def seqrc(seqstr):
    """用于形成反向互补序列"""
    rc_seq=seqstr.replace('A','X').replace('C','Y').replace('T','A').replace('G','C').replace('X','T').replace('Y','G').replace('a','x').replace('c','y').replace('t','a').replace('g','c').replace('x','t').replace('y','g')[::-1]
    return rc_seq

def unitestreamfa(bksite_list, upstreamfa, downstreamfa, dirprefix, upstrand, downstrand):
    """合并bk1上游及bk2下游的序列，同时注意哪些序列是reverse_complement"""
    otlist=[]
    for i,part in enumerate(zip_longest( open(upstreamfa,'r'), open(downstreamfa,'r') ) ):
        if part[0].startswith('>') and part[1].startswith('>'):
            ids='>'+bksite_list[i]
            otlist.append(ids+'\n')
        elif not part[0].startswith('>') and not part[1].startswith('>'):
            # 这里加入strand判断，sequence序列是否需要rc(reverse_complement)，一共存在下面几种形式
            # 设计思路和判别方法参考：最终综合结果得到的SV类型存在信息
            if upstrand[i] == 'rf' and downstrand[i] == 'rf':
                otlist[-1]=otlist[-1]+part[0].strip()+part[1]
            elif upstrand[i] == 'rc' and downstrand[i] == 'rf':
                rc_seq=seqrc(part[0].strip())
                otlist[-1]=otlist[-1]+rc_seq+part[1]
            elif upstrand[i] == 'rf' and downstrand[i] == 'rc':
                rc_seq=seqrc(part[1].strip())
                otlist[-1]=otlist[-1]+part[0].strip()+rc_seq+'\n'

    otlist=list(set(otlist))
    with open('%s.fusion_allstream.fa'%dirprefix,'w') as otopen:
        otopen.writelines(otlist)
    allstreamfa='%s.fusion_allstream.fa'%dirprefix
    return allstreamfa

def fusionfastaindex(output_dir, output_prefix, samtools, bwa, genome, step1bkresult, break_point_stream_length):
    """通过上下游的信息形成最终的融合基因组"""
    try:
        os.mkdir('%s/01softread'%output_dir)
    except FileExistsError:
        pass
    output_dir02=output_dir+'/01softread'
    dirprefix=os.path.join(output_dir02, output_prefix)
    bksite_list=[]
    upstream_list=[]
    upstrand=[]
    downstream_list=[]
    downstrand=[]
    with open(step1bkresult,'r') as inopen:
        for line in inopen:
            col=line.strip().split('\t')
            # 设计思路和判别方法参考：最终综合结果得到的SV类型存在信息
            bksite_list.append(col[0])
            bksite_list.append(col[0])
            # bk1 upstream
            upstream_list.append(col[1]+'\n')
            upstrand.append(col[2]) # 最终形成上下游的fasta文件，其奇数行为序列id，偶数行为fasta序列
            upstrand.append(col[2]) # upstrand及downstrand不管对应奇数行还是偶数行，都应该判断是正向rf或反向互补rc，对于序列需要反向互补，对于id则判断前后哪一个标记是bk位点；
            # bk2 downstream
            downstream_list.append(col[3]+'\n')
            downstrand.append(col[4])
            downstrand.append(col[4])
    with open('%s.upstream_region_file.txt'%dirprefix,'w') as otopen:
        otopen.writelines(upstream_list)
    with open('%s.downstream_region_file.txt'%dirprefix,'w') as otopen:
        otopen.writelines(downstream_list)
    upstream_fa_cmd=" ".join([samtools, 'faidx', genome, '-r', '%s.upstream_region_file.txt'%dirprefix, '-n', str(break_point_stream_length), '-o', '%s.upstream.fa'%dirprefix])
    downstream_fa_cmd=" ".join([samtools, 'faidx', genome, '-r', '%s.downstream_region_file.txt'%dirprefix, '-n', str(break_point_stream_length), '-o', '%s.downstream.fa'%dirprefix])
    os.system(upstream_fa_cmd)
    os.system(downstream_fa_cmd)
    upfa='%s.upstream.fa'%dirprefix
    downfa='%s.downstream.fa'%dirprefix
    allstreamfa=unitestreamfa(bksite_list, upfa, downfa, dirprefix, upstrand, downstrand)
    # os.system('/usr/bin/rm %s.upstream_region_file.txt %s.downstream_region_file.txt %s.upstream.fa %s.downstream.fa'%(dirprefix,dirprefix,dirprefix,dirprefix))
    index_cmd=" ".join([bwa, 'index', allstreamfa]) #step1allstreamfa='%s.fusion_allstream.fa'%dirprefix
    os.system(index_cmd)

#第三部分 合并断点上下游一定区域的genome序列分析内容结束：包括fusionfastaindex主函数用于得到上游断点的一段序列和下游断点的一段序列，unitestreamfa用于合并两个断点的序列得到融合基因组同时区分正向及反向互补序列；
#######################################################
def samtoolsdp(samtools, input_bam, step1bkresult, depth_extend, output_dir, output_prefix):
    """通过samtools depth计算两个breakpoint断点两侧一定区域深度的平均值"""
    try:
        os.mkdir('%s/02remapping'%output_dir)
    except FileExistsError:
        pass
    output_dir02=output_dir+'/02remapping'
    dirprefix=os.path.join(output_dir02, output_prefix)
    depthregionbed=[]
    bkpairsite=[]
    upsite=[]
    upstrand=[]
    downsite=[]
    downstrand=[]
    with open(step1bkresult,'r') as inopen:
        for line in inopen:
            col=line.strip().split('\t')
            bkpairsite.append(col[0])
            upbed=col[6].replace('_','\t')
            upsite.append(upbed)
            upstrand.append(col[2])
            downbed=col[7].replace('_','\t')
            downsite.append(downbed)
            downstrand.append(col[4])
            depthregionbed.append('%s\n%s\n'%(upbed,downbed))
    with open('%s.region.bed'%dirprefix,'w') as otopen:
        otopen.writelines(depthregionbed)
    depth_cmd=" ".join([samtools, 'depth -b', '%s.region.bed'%dirprefix, '-a -d 0 -s -J -G UNMAP -o', '%s.region.depth'%dirprefix, input_bam])
    os.system(depth_cmd)
    sitedp_dict={}
    with open('%s.region.depth'%dirprefix,'r') as inopen:
        for line in inopen:
            col=line.strip().split('\t')
            sitedp_dict[( col[0], int(col[1]) )]=int(col[2])
    otpairbkdp={}
    for i,key in enumerate(bkpairsite):
        part=( ( key.split('-')[0].split(':')[0], int(key.split('-')[0].split(':')[1]) ), ( key.split('-')[1].split(':')[0], int(key.split('-')[1].split(':')[1]) ) )
        # part[0] is upstream; part[1] is downstream
        # 序列最终id同样也要加入strand判断，rc(reverse_complement)则断点位置信息不同，一共存在下面几种形式
        # 设计思路和判别方法参考：最终综合结果得到的SV类型存在信息
        if upstrand[i] == 'rf' and downstrand[i] == 'rf':
            upsitelist=[(part[0][0], part[0][1] + x) for x in range(-(depth_extend-1),1)]
            downsitelist=[(part[1][0], part[1][1] + x) for x in range(0,depth_extend)]
        elif upstrand[i] == 'rc' and downstrand[i] == 'rf':
            upsitelist=[(part[0][0], part[0][1] + x) for x in range(0,depth_extend)]
            downsitelist=[(part[1][0], part[1][1] + x) for x in range(0,depth_extend)]
        elif upstrand[i] == 'rf' and downstrand[i] == 'rc':
            upsitelist=[(part[0][0], part[0][1] + x) for x in range(-(depth_extend-1),1)]
            downsitelist=[(part[1][0], part[1][1] + x) for x in range(-(depth_extend-1),1)]
        updp=[sitedp_dict.get(x) if sitedp_dict.get(x) is not None else 0 for x in upsitelist]
        downdp=[sitedp_dict.get(x) if sitedp_dict.get(x) is not None else 0 for x in downsitelist]
        otpairbkdp[key]=(round(sum(updp)/len(updp), 2), round(sum(downdp)/len(downdp), 2))
    with open('%s.each_bk.depth.json'%dirprefix,'w') as otopen:
        otopen.write(json.dumps(otpairbkdp))

#第四部分 用于计算断点的深度分析内容结束：仅包括samtoolsdp主函数同时区分正向及反向互补序列；
#######################################################
def fusionseqremapping(output_dir, output_prefix, samtools, bwa, step1allstreamfa, step1crossfq):
    """可能存在支持融合的read比对到融合基因组"""
    try:
        os.mkdir('%s/02remapping'%output_dir)
    except FileExistsError:
        pass
    output_dir02=output_dir+'/02remapping'
    dirprefix=os.path.join(output_dir02, output_prefix)
    bwamem_cmd=" ".join([bwa, 'mem -t 16 -Y -M', step1allstreamfa, step1crossfq, '|', samtools, 'view --threads 24 -b -S -F 12 |', samtools, 'sort --threads 24 -o', '%s.support.fusion.sort.bam'%dirprefix])
    os.system(bwamem_cmd)
    """这里在考虑使用bowtie还是用bwa做分析；最终还是选择bwa，因为允许read存在比对insertion和deletion的存在"""
    index_cmd=" ".join([samtools, 'index', '%s.support.fusion.sort.bam'%dirprefix, '%s.support.fusion.sort.bai'%dirprefix])
    os.system(index_cmd)
    supportfusion_bam='%s.support.fusion.sort.bam'%dirprefix

#第五部分 用于将潜在支持融合的read比对到融合后的基因组用于融合的验证分析内容结束：仅包含fusionseqremapping主函数
#######################################################
def acrossbreakpoint(mstart, mlen, break_point_stream_length, min_corss_breakpoint_length):
    """认定为跨越断点的read的条件"""
    across=False
    if mstart < (break_point_stream_length - min_corss_breakpoint_length) and (mstart + mlen) > break_point_stream_length:
        across=True
    return across

def fusionsupportread(supportfusion_bam, step1bkresult, output_dir, output_prefix, break_point_stream_length, bk_depth_json, min_corss_breakpoint_length):
    """查找融合基因组的read重新比对后是否有支持断点的read
    判断标准：1)read过断点，2)read1在断点的上游或下游read2在其断点的另一边"""
    try:
        os.mkdir('%s/03report'%output_dir)
    except FileExistsError:
        pass
    output_dir03=output_dir+'/03report'
    dirprefix=os.path.join(output_dir03, output_prefix)
    bkcrossreads={} #{bk:[reads]}
    samfile=pysam.AlignmentFile(supportfusion_bam, "rb")
    with open(step1bkresult,'r') as inopen:
        for line in inopen:
            chrn=line.split('\t')[0]
            try:
                msg=samfile.fetch(chrn)
            except ValueError:
                continue # bam文件中没有融合位点信息
            # 对于每一个bk位点都开始记录
            upreads=[]
            downreads=[]
            splitreads=[]
            crossreads=[]
            for each in msg:
            # each.query_name, each.reference_start, each.cigarstring分别表示比对的read名称，比对的起始位置(0base，计算长度时需要从0base开始)，比对的cigar值(用于计算除S外的长度)
                readname=each.query_name # bam file column 1
                mstart=each.reference_start # bam file column 4
                cigarlist=each.cigar #M 0 I 1 D 2 S 4
                mlen=sum([x[1] for x in cigarlist if x[0] != 4])
                cross=acrossbreakpoint(mstart, mlen, break_point_stream_length, min_corss_breakpoint_length)
                if cross:
                    splitreads.append(readname)
                    # crossreads.append(readname)
                else:
                    if (mstart + mlen) <= break_point_stream_length:
                        upreads.append(readname)
                    elif mstart > break_point_stream_length:
                        downreads.append(readname)
                try:
                    lotsofxa=each.get_tag('XA') #这里XA标签可能有多个位置
                except KeyError:
                    continue
                for xa in lotsofxa.split(';')[:-1]:
                    xa_chrn=xa.split(',')[0]
                    xa_mstart=int(xa.split(',')[1][1:])
                    xa_cigar=xa.split(',')[2]
                    xa_mlen=sum([int(x[:-1]) for x in re.findall(r'\d+[A-Z]', xa_cigar) if 'S' not in x])
                    xa_cross=acrossbreakpoint(xa_mstart, xa_mlen, break_point_stream_length, min_corss_breakpoint_length)
                    if xa_cross:
                        if bkcrossreads.get(xa_chrn) is None: #XA标签中只寻找splitreads
                            bkcrossreads.update({xa_chrn:{'splitreads':[readname], 'crossreads':[]}})
                        else:
                            bkcrossreads.get(xa_chrn).get('splitreads').append(readname)
            streampair=list( set(upreads) & set(downreads) ) #read1或read2比对断点上游，同名的read比对断点下游则认为pairread跨越了断点
            crossreads.extend(streampair)
            if bkcrossreads.get(chrn) is None:
                bkcrossreads.update({chrn:{'splitreads':splitreads, 'crossreads':crossreads}})
            else:
                bkcrossreads.get(chrn).get('splitreads').extend(splitreads)
                bkcrossreads.get(chrn).get('crossreads').extend(crossreads)
    with open('%s.bkcrossreads.json'%dirprefix,'w') as otopen:
        otopen.write(json.dumps(bkcrossreads))
    result_list=['fusion.site(breakpoint1-breakpoint2)\tcross.support.splitread.count\tpairs.support.count\tbreakpoint1.depth\tbreakpoint2.depth\tsplitread.list\tpairsread.list\n']
    bkdpdict=json.load(open(bk_depth_json,'r')) #{bkpairsite:(bk1dp, bk2dp)}
    for each in list(bkcrossreads.keys()):
        bk1dp=bkdpdict.get(each)[0]
        bk2dp=bkdpdict.get(each)[1]
        uniqsplit=list(set(bkcrossreads.get(each).get('splitreads')))
        uniqcross=list(set(bkcrossreads.get(each).get('crossreads')))
        splitcount=len(uniqsplit)
        crosscount=len(uniqcross)
        splitlist=",".join(uniqsplit)
        crosslist=",".join(uniqcross)
        result_list.append('%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(each,splitcount,crosscount,str(bk1dp),str(bk2dp),splitlist,crosslist))
    with open('%s.breakpoint.count.txt'%dirprefix,'w') as otopen:
        otopen.writelines(result_list)

#第六部分 计算最终融合基因组的支持read数融合断点深度：包括fusionsupportread主函数分析重新比对的bam文件，对于每个融合genome都通过acrossbreakpoint判断这条read是否支持融合基因组；
#######################################################
def resultreport(bkresult, bkannotate, breakpointcount, output_dir, output_prefix):
    finalresult=['fusion.type\tgenes\tbreakpoint\treadcount\tbk1.dp\tbk2.dp\n']
    msg_dict={}
    rfrclist=[]
    for part in zip_longest( open(bkresult,'r'), open(bkannotate,'r') ):
        rfrclist.append('%s\t%s\t%s\n'%(part[0].split('\t')[0], part[0].split('\t')[2], part[0].split('\t')[4]))
        key=part[0].split('\t')[0]
        types=part[0].split('\t')[5]
        gene=part[1].strip().split('\t')[1]
        msg_dict.update({key:[types, gene]})
    with open('%s/03report/%s.bk_rf_rc.txt'%(output_dir, output_prefix),'w') as otopen:
        otopen.writelines(rfrclist)
    for line in os.popen("cat %s |awk 'NR>1{print $0}' "%breakpointcount):
        bkpair=line.split('\t')[0]
        readcount=str( int(line.split('\t')[1]) + int(line.split('\t')[2]) )
        bk1dp=line.split('\t')[3]
        bk2dp=line.split('\t')[4]
        types=msg_dict.get(bkpair)[0]
        gene=msg_dict.get(bkpair)[1]

        tmpbk1=bkpair.split('-')[0]
        tmpbk2=bkpair.split('-')[1]
        if tmpbk1.split(':')[0] == tmpbk2.split(':')[0] and ( int(tmpbk1.split(':')[1])-int(tmpbk2.split(':')[1])<100 or int(tmpbk1.split(':')[1])-int(tmpbk2.split(':')[1])>-100 ) :
            continue
        finalresult.append("\t".join([types, gene, bkpair, readcount, bk1dp, bk2dp])+'\n')
    with open('%s/03report/%s.breakpoint.result.txt'%(output_dir, output_prefix),'w') as otopen:
        otopen.writelines(finalresult)

#######################################################
def faindexremapdouble(output_dir, output_prefix, samtools, bwa, genome, step1bkresult, step1crossfq, break_point_stream_length):
    fusionfastaindex(output_dir, output_prefix, samtools, bwa, genome, step1bkresult, break_point_stream_length)
    step1allstreamfa = '%s/01softread/%s.fusion_allstream.fa'%(output_dir, output_prefix)
    fusionseqremapping(output_dir, output_prefix, samtools, bwa, step1allstreamfa, step1crossfq)

#######################################################
def running(input_bam,output_dir,output_prefix,genome,samtools,bowtie1,bwa,bed,min_single_base_content,min_mismatch,min_insert,min_delete,min_corss_breakpoint_length,depth_extend,bk_extend_size,min_bamsoft_length,break_point_stream_length,break_point_annotation_bed):
    """分析方法总共分为以下几步：
    1.找到bam文件中存在的SA tag的soft clipped read，通过cigar的M值和S值来确定fusion的breakpoint1和breakpoint2；同时建立融合区域合并池；
    （最终还是选择放弃这一步，主要包括下面的原因：cigar及SA标签中会带有两个S的，例如7S61M32S,40S58M2S或38S62M,51S41M8S或13S60M27S,52M48S等，同时还会存在SA标签中cigar有两个结果的例如33S30M37S SA:Z:chr14,65369458,-,50S50M,8,0;chr16,18724050,+,53S47M,8,2; 对于这种情况考虑问题很复杂，所以不如直接将标记有S的序列提取出来用bowtie1比对到genome）
    2.找到bam文件中除存在SA tag的cigar当中带有S的read信息（实际上也没有去除SA tag信息的read），将S的序列提取出用bowtie1比对到genome，得到breakpoint1和breakpoint2（前面两个分析内容并行运行）；再与融合区域合并池合并和更新；这一步就得到最终所有待验证的融合位点；
    3.samtools depth计算融合位点一定区域内的深度，同时建立融合位点breakpoint1上游一定长度和breakpoint2下游一定长度的genome序列；
    4.将bam文件中所有cigar当中带有S的read提取除与这个融合基因组做bwa mem的比对，确定支持融合的read；（与第3步并行运行）
    \... ...

    要求，分析函数从下向上写，每一步写明具体的含义和技术要求；
    总体分析思路参考文档E:\工作文件夹\工作任务\3.软件开发\开发fusion分析软件方案设计.pptx
    """
    pool=Pool(processes=2)
    pool.apply_async(buildingbreakpoint, args=(input_bam, output_dir, output_prefix, samtools, bed, bowtie1, genome, min_bamsoft_length, min_mismatch, min_insert, min_delete, bk_extend_size, break_point_stream_length, depth_extend, min_single_base_content, break_point_annotation_bed))
    pool.apply_async(allcrossfastq, args=(samtools, output_dir, output_prefix, input_bam, break_point_stream_length))
    pool.close()
    pool.join()
    # buildingbreakpoint(input_bam, output_dir, output_prefix, samtools, bed, bowtie1, genome, min_bamsoft_length, min_mismatch, min_insert, min_delete, bk_extend_size, break_point_stream_length, depth_extend, min_single_base_content, break_point_annotation_bed)
    step1bkresult = '%s/01softread/%s.bkresult.txt'%(output_dir, output_prefix)
    step1bkannotate = '%s/01softread/%s.bk_annotate.txt'%(output_dir, output_prefix)
    # allcrossfastq(samtools, output_dir, output_prefix, input_bam, break_point_stream_length)
    step1crossfq = '%s/01softread/%s.cross.fq.gz'%(output_dir, output_prefix)

    pool=Pool(processes=2)
    pool.apply_async(faindexremapdouble, args=(output_dir, output_prefix, samtools, bwa, genome, step1bkresult, step1crossfq, break_point_stream_length))
    pool.apply_async(samtoolsdp, args=(samtools, input_bam, step1bkresult, depth_extend, output_dir, output_prefix))
    pool.close()
    pool.join()
    # fusionfastaindex(output_dir, output_prefix, samtools, bwa, genome, step1bkresult, break_point_stream_length)
    step1allstreamfa = '%s/01softread/%s.fusion_allstream.fa'%(output_dir, output_prefix)
    # samtoolsdp(samtools, input_bam, step1bkresult, depth_extend, output_dir, output_prefix)
    step2bkdepth = '%s/02remapping/%s.each_bk.depth.json'%(output_dir, output_prefix)
    # fusionseqremapping(output_dir, output_prefix, samtools, bwa, step1allstreamfa, step1crossfq)
    step2remappingbam = '%s/02remapping/%s.support.fusion.sort.bam'%(output_dir, output_prefix)

    fusionsupportread(step2remappingbam, step1bkresult, output_dir, output_prefix, break_point_stream_length, step2bkdepth, min_corss_breakpoint_length)

    step3breakpointcount = '%s/03report/%s.breakpoint.count.txt'%(output_dir, output_prefix)
    resultreport(step1bkresult, step1bkannotate, step3breakpointcount, output_dir, output_prefix)

if __name__ == '__main__':
    localdir=os.path.dirname(os.path.abspath(__file__))
    parser = ArgumentParser('fusion detect methods')
    parser.add_argument('-m', '--input_bam', required=True,
                        help='mapping result "bam" file')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-g', '--genome', type=str, default="%s/../../../database/gatk-bundle-hg19/ucsc.hg19.fasta"%localdir,
                        help='bowtie1 build genome fasta file, default=database path')
    parser.add_argument('-tools', '--samtools', type=str, default="%s/../../../software/samtools-1.13/samtools"%localdir,
                        help='"samtools" software, default=software path')
    parser.add_argument('-bwa', '--bwa', type=str, default="%s/../../../software/bwa-0.7.17/bwa"%localdir,
                        help='"bwa" software, default=software path')
    parser.add_argument('-bowtie1', '--bowtie1', type=str, default="%s/../../../software/bowtie-1.3.0/bowtie"%localdir,
                        help='"bowtie1" software, default=software path')
    parser.add_argument('-b', '--bed', default=None,
                        help='add fusion region bed file')
    parser.add_argument('-minc', '--min_single_base_content', type=float, default=0.5,
                        help='allow break point match reads soft sequence only single base content less than 50%% percent. default=0.5')
    parser.add_argument('-mins', '--min_mismatch', type=int, default=2,
                        help='support break point match reads minimum mismatch base number in each of read. default=2')
    parser.add_argument('-mini', '--min_insert', type=int, default=3,
                        help='support break point match reads minimum insert base number in each of read. default=3')
    parser.add_argument('-mind', '--min_delete', type=int, default=3,
                        help='support break point match reads minimum delete base number in each of read. default=3')
    parser.add_argument('-minr', '--min_corss_breakpoint_length', type=int, default=5,
                        help='cross break point minimum length (bp). default=5')
    parser.add_argument('-dpextend', '--depth_extend', type=int, default=5,
                        help='samtools depth calculate range (bp). default=5')
    parser.add_argument('-breg', '--bk_extend_size', type=int, default=10,
                        help='merge break point upstream and downstream range length (bp). default=10')
    parser.add_argument('-minss', '--min_bamsoft_length', type=int, default=11,
                        help='minimum bam file soft clipped length (bp), over this view the bases will bowtie1 to genome DNA sequence, bowtie1 minimum support length is 11bp. default=11')
    parser.add_argument('-stream', '--break_point_stream_length', type=int, default=300,
                        help='break point in reference upstream or downstream length (bp). default=300')
    parser.add_argument('-anno', '--break_point_annotation_bed', type=str, default="%s/hg19.ncbiRefSeq.all_CDS_longest_transcript_id_reverse_correct_add_intron.bed"%localdir,
                        help="break point annotation bed file. default=local directory")
    # 增加一个continue的功能，能够查找关键步骤的文件是否存在可以继续跑
    args = parser.parse_args()
    running(args.input_bam,args.output_dir,args.output_prefix,args.genome,args.samtools,args.bowtie1,args.bwa,args.bed,args.min_single_base_content,args.min_mismatch,args.min_insert,args.min_delete,args.min_corss_breakpoint_length,args.depth_extend,args.bk_extend_size,args.min_bamsoft_length,args.break_point_stream_length,args.break_point_annotation_bed)

# 1S$
# V350015794L1C002R0331224914     339     chr1    537110  0       61S39M  chr19   7143033 0       CGTCACATTCCCAACATCGCCAAGGGACCTGCGTTTCCGAGATGGCCTGGAACGACAGTAGCCGATAATTGTGTCTTTCCATATACACAAAAGTGAAGTC    EE-D<CEG:FGBF2GF3AGEDF9FECCFGGG@?GFFF@FFFGEFFFG=GEF3CGFDFFGEGGFGFEGDFGEEA5FFFDFEGEGEGEFEFGGEEEBEEFFG    SA:Z:chr19,7143056,-,61M39S,60,0;       XA:Z:chr6,+171030203,39M61S,0;chr8,-36733,61S39M,0;chr5,+180872377,39M61S,0;chr1,+445823,39M61S,0;chr17,+81159946,39M61S,0;     MC:Z:84M16S     MD:Z:39 PG:Z:MarkDuplicates     RG:Z:ZLC210116A NM:i:0  AS:i:39 XS:i:39
# chr1:537110=V350015794L1C002R0331224914_339_chr1_537110_61S39M_1S       0       chr19   7143056 255     61M
# 2S$
# V350015794L1C004R0710231277     97      chr1    566447  0       65M35S  chr20   7513609 0       CCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCAGGATGAATAATCCTTATGTGAAAAGATTTGGACC    FFFFFFFFFGFFFFFFFEFFFBCFEFF5ECFFFFGFEFFFFGFGFFFFFGGFFFFFFGFFFFGGFCFF@FF?FGGFFFFEFE>GFDAFFFEGFFF@F>FG    SA:Z:chr20,7513602,+,65S35M,0,0;        XA:Z:chrM,+5898,65M35S,0;       MC:Z:100M       MD:Z:65 PG:Z:MarkDuplicates     RG:Z:ZLC210116A NM:i:0  AS:i:65 XS:i:65
# chr1:566512=V350015794L1C004R0710231277_97_chr1_566447_65M35S_2S	0	chr20	7513602	255	35M	*	0	0	AGGATGAATAATCCTTATGTGAAAAGATTTGGACC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	XA:i:0	MD:Z:35	NM:i:0	XM:i:2