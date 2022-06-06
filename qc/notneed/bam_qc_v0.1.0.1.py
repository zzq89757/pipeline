#coding:utf8
#########################################################################
# Author: ZZZ 
# Created Time: 2021-10-14 
# Version: 0.2.0.0
# Description: pipline for bam QC
#########################################################################
import os
import sys
import json
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def process_running(cmd_process, output_dir, output_prefix, notice):
    (stdout, stderr, timer)=brc.run_command_with_return(cmd_process)
    brc.export_log(cmd_process, stdout, stderr, timer, output_dir, output_prefix, notice)

def picard_removedup(java, picard, input_bam, dirprefix, output_dir, output_prefix, umi):
    rmdup_bam='%s.MarkDuplicates_tmp.bam'%dirprefix
    cmd_markdup=" ".join([
        java['java'], '-jar', '-Xmx%sG'%java['java_mem'],
        picard, 'MarkDuplicates',
        '-TMP_DIR', '%s_tmp'%dirprefix,
        '-I', input_bam,
        '-O', rmdup_bam,
        '-M', '%s.MarkDuplicates.metrics'%dirprefix,
        '-CREATE_INDEX true', '--REMOVE_DUPLICATES true'
    ])
    if umi == 'true':
        cmd_markdup=cmd_markdup+' --BARCODE_TAG RX'
    process_running(cmd_markdup, output_dir, output_prefix, 'Notice:\t*.sort.bam to *.MarkDuplicates_tmp.bam running messages')
    return rmdup_bam

def picard_insertsize(java, picard, input_bam, dirprefix, output_dir, output_prefix):
    cmd_insertsize=" ".join([
        java['java'], '-jar', '-Xmx%sG'%java['java_mem'],
        picard, 'CollectInsertSizeMetrics',
        '-I', input_bam,
        '-W', '600',
        '-O','%s.CollectInsertSizeMetrics.xls'%dirprefix,
        '-H','%s.CollectInsertSizeMetrics.pdf'%dirprefix
    ])
    process_running(cmd_insertsize, output_dir, output_prefix, 'Notice:\tpicard CollectInsertSizeMetrics running messages')

def picard_bedtolist(java, picard, bed, dirprefix, hg19_dict, output_dir, output_prefix):
    cmd_bedtolist=" ".join([
        java['java'], '-jar', '-Xmx%sG'%java['java_mem'],
        picard, 'BedToIntervalList',
        '-I', bed,
        '-O', '%s.BedToIntervalList.xls'%dirprefix,
        '-SD', hg19_dict
    ])
    process_running(cmd_bedtolist, output_dir, output_prefix, 'Notice:\tpicard BedToIntervalList running messages')

def picard_collecths(java, picard, input_bam, dirprefix, hg19_fa, output_dir, output_prefix, times):
    cmd_collecths=" ".join([
        java['java'], '-jar', '-Xmx%sG'%java['java_mem'],
        picard, 'CollectHsMetrics',
        '-I', input_bam,
        '-O', '%s.CollectHsMetrics_%s.xls'%(dirprefix, times),
        '-R', hg19_fa,
        '-BAIT_INTERVALS', '%s.BedToIntervalList.xls'%dirprefix,
        '-TARGET_INTERVALS', '%s.BedToIntervalList.xls'%dirprefix
    ])
    process_running(cmd_collecths, output_dir, output_prefix, 'Notice:\tpicard CollectHsMetrics running messages')

def qctable(output_prefix, dirprefix, markdup_metrics):
    if markdup_metrics is None:
        dup_rst='%s.MarkDuplicates.metrics'%dirprefix
    else: # if not brc.check_file(dup_rst):
        dup_rst=markdup_metrics
    insert='%s.CollectInsertSizeMetrics.xls'%dirprefix
    nondup_hs='%s.CollectHsMetrics_nondup.xls'%dirprefix
    dedup_hs='%s.CollectHsMetrics_dedup.xls'%dirprefix
    dup_ratio=str(round(float(os.popen('cat %s |head -n 8 |tail -n 1'%dup_rst).read().strip().split('\t')[-2])*100, 2))
    # PERCENT_DUPLICATION
    mean_insert_size=str(round(float(os.popen('cat %s |head -n 8 |tail -n 1'%insert).read().strip().split('\t')[5]), 2))
    # MEAN_INSERT_SIZE
    dedup_mean_depth=str(round(float(os.popen('cat %s |head -n 8 |tail -n 1'%dedup_hs).read().strip().split('\t')[33]), 2))
    # MEAN_TARGET_COVERAGE
    qc_list=os.popen('cat %s |head -n 8 |tail -n 1'%nondup_hs).read().strip().split('\t')
    total_reads='{:,}'.format(int(qc_list[22])) # TOTAL_READS
    mapped_ratio=str(round(float(qc_list[32])*100, 2)) # PCT_PF_UQ_READS_ALIGNED
    target_ratio=str(round(float(qc_list[6])*100, 2)) # PCT_SELECTED_BASES
    mean_depth=str(round(float(qc_list[33]), 2)) # MEAN_TARGET_COVERAGE
    coverage_ratio=str(round(float(qc_list[45])*100, 2)) # PCT_TARGET_BASES_1X
    coverage_1000_ratio=str(round(float(qc_list[55])*100, 2)) # PCT_TARGET_BASES_1000X
    try:
        fold_80=str(round(float(qc_list[44]), 2)) # FOLD_80_BASE_PENALTY (这个值计算会出现?号)
    except ValueError:
        fold_80=qc_list[44]
    opt_list=[]
    opt_list.append("\t".join(['sample.id', 'total.read', 'mapped.ratio(%)', 'mean.insert.size', 'target.ratio(%)', 'duplicated.ratio(%)', 'mean.depth', 'dedup.mean.depth', 'coverage.ratio(%)', 'coverage.1000.ratio(%)', 'fold-80'])+'\n')
    opt_list.append("\t".join([output_prefix, total_reads, mapped_ratio, mean_insert_size, target_ratio, dup_ratio, mean_depth, dedup_mean_depth, coverage_ratio, coverage_1000_ratio, fold_80])+'\n')
    with open('%s.bam_QC.xls'%dirprefix, 'w') as otopen:
        otopen.writelines(opt_list)

def running(raw_bam, bed, cfg_json, output_dir, output_prefix, markdup_metrics, removedup_bam, umi, group_bam, consensus_mapped_bam):
    # this_script_dir=os.path.dirname(os.path.abspath(__file__))
    config=json.load(open(cfg_json,'r'))
    dirprefix=os.path.join(output_dir, output_prefix)
    if not brc.check_file(raw_bam):
        raise IOError('%s must be exist'%(raw_bam))
    java=config['language']
    picard=config['software']['picard']
    samtools=config['software']['samtools']
    hg19_fa=config['database']['hg19_fa']
    hg19_dict=config['database']['hg19_dict']
    ###
    picard_bedtolist(java, picard, bed, dirprefix, hg19_dict, output_dir, output_prefix)
    #当raw_bam, markdup_metrics, removedup_bam同时存在时，removedup_bam insertsize; raw_bam collecths; removedup_bam collecths
    #只有raw_bam则需要先做raw_bam collecths, raw_bam removedup之后再做raw_bam_rmdup insertsize, raw_bam_rmdup collecths
    #当有umi存在时需要先做group_bam sort, group_bam removedup, raw_bam collecths之后再做consensus_mapped_bam insertsize, consensus_mapped_bam collecths
    if umi == 'false':
        if brc.check_file(markdup_metrics) and brc.check_file(removedup_bam):
            pool=Pool(processes=3)
            pool.apply_async(picard_insertsize, args=(java, picard, removedup_bam, dirprefix, output_dir, output_prefix))
            pool.apply_async(picard_collecths, args=(java, picard, raw_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'nondup'))
            pool.apply_async(picard_collecths, args=(java, picard, removedup_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'dedup'))
            pool.close()
            pool.join()
            qctable(output_prefix, dirprefix, markdup_metrics)
            otstr='QC analysis with non of UMI\n%s\n%s\n%s\nused\n'%(raw_bam, markdup_metrics, removedup_bam)
            with open('%s_file_used.txt'%dirprefix,'w') as otopen:
                otopen.write(otstr)
        else:
            pool=Pool(processes=2)
            pool.apply_async(picard_removedup, args=(java, picard, raw_bam, dirprefix, output_dir, output_prefix, umi))
            pool.apply_async(picard_collecths, args=(java, picard, raw_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'nondup'))
            pool.close()
            pool.join()
            rmdup_bam='%s.MarkDuplicates_tmp.bam'%dirprefix

            pool=Pool(processes=2)
            pool.apply_async(picard_collecths, args=(java, picard, rmdup_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'dedup'))
            pool.apply_async(picard_insertsize, args=(java, picard, rmdup_bam, dirprefix, output_dir, output_prefix))
            pool.close()
            pool.join()
            qctable(output_prefix, dirprefix, markdup_metrics)
            otstr='QC analysis with non of UMI\n%s\nused\n'%(raw_bam)
            with open('%s_file_used.txt'%dirprefix,'w') as otopen:
                otopen.write(otstr)
    elif umi == 'true':
        if not brc.check_file(group_bam) or not brc.check_file(consensus_mapped_bam):
            raise IOError('%s and %s file must be exist'%(group_bam, consensus_mapped_bam))
        else:
            cmd_sort=" ".join([
                samtools, 'sort --threads %s'%config['parameter']['samtools_threads'], 
                group_bam, '-o', '%s.group.sort.bam'%dirprefix
            ])
            process_running(cmd_sort, output_dir, output_prefix, 'Notice:\tgroup bam samtools sort running messages')
            group_sort_bam='%s.group.sort.bam'%dirprefix
            pool=Pool(processes=2)
            pool.apply_async(picard_removedup, args=(java, picard, group_sort_bam, dirprefix, output_dir, output_prefix, umi))
            pool.apply_async(picard_collecths, args=(java, picard, raw_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'nondup'))
            pool.close()
            pool.join()

            pool=Pool(processes=2)
            pool.apply_async(picard_collecths, args=(java, picard, consensus_mapped_bam, dirprefix, hg19_fa, output_dir, output_prefix, 'dedup'))
            pool.apply_async(picard_insertsize, args=(java, picard, consensus_mapped_bam, dirprefix, output_dir, output_prefix))
            pool.close()
            pool.join()
            qctable(output_prefix, dirprefix, markdup_metrics)
            otstr='QC analysis with UMI\n%s\n%s\n%s\nused\n'%(raw_bam, group_bam, consensus_mapped_bam)
            with open('%s_file_used.txt'%dirprefix,'w') as otopen:
                otopen.write(otstr)

if __name__ == '__main__':
    parser = ArgumentParser('Quality Control')
    parser.add_argument('-b', '--raw_bam', required=True,
                        help="if methods not with UMI and downsample is PASS, this file is '*.sort.downsample_num.bam';\n and other types this file is '*.sort.bam'")
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg_json', required=True,
                        help='"configure.file.json" file')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-mc', '--markdup_metrics', required=None,
                        help='input "*.mkdup.metrics" file, if exists picard MarkDuplicates will not running, if umi parameter is true this will ignore; default=None')
    parser.add_argument('-rb', '--removedup_bam', required=None,
                        help='if parameter "markdup_metrics" is True; this parameter will input "*.rmdup.bam" file or "*.rmdup_IndelRealigner.bam"; default=None')
    parser.add_argument('-umi', '--umi', type=str, choices=['true', 'false'], default='false',
                        help='count dup ratio with UMI or not, default=false')
    parser.add_argument('-g', '--group_bam', default=None,
                        help='*.group.bam file with UMI, this parameter is used with umi is "True", default=None')
    parser.add_argument('-s', '--consensus_mapped_bam', default=None,
                        help='*.consensus.mapped.bam file with UMI, this parameter is used with umi is "True", default=None')
    args = parser.parse_args()
    running(args.raw_bam, args.bed, args.cfg_json, args.output_dir, args.output_prefix, args.markdup_metrics, args.removedup_bam, args.umi, args.group_bam, args.consensus_mapped_bam)

# QC指标				
# 输出指标	输出文件中选择的参数	picard模块	非umi的bam文件	umi后的bam文件
# total.reads	BAIT_TERRITORY	CollectHsMetrics	sort.bam	sort.bam
# mapped.ratio	PCT_PF_UQ_READS_ALIGNED	CollectHsMetrics	sort.bam	sort.bam
# target.ratio	PCT_SELECTED_BASES	CollectHsMetrics	sort.bam	sort.bam
# mean.insert.size	MEDIAN_INSERT_SIZE	CollectInsertSizeMetrics	sort.bam	sort.bam
# duplicated.ratio	PERCENT_DUPLICATION	MarkDuplicates	sort.bam	group.sort.bam(GroupReadsByUmi之后输出的bam文件)
# mean.depth	MEAN_TARGET_COVERAGE	CollectHsMetrics	sort.bam	sort.bam
# dedup.mean.depth	MEAN_TARGET_COVERAGE	CollectHsMetrics	rmdup.bam	consensus.mapped.bam(最终矫正一致性序列后的bam文件)
# coverage.ratio	PCT_TARGET_BASES_1X	CollectHsMetrics	sort.bam	sort.bam
# coverage.1000.ratio	PCT_TARGET_BASES_1000X	CollectHsMetrics	sort.bam	sort.bam
# fold-80	FOLD_80_BASE_PENALTY	CollectHsMetrics	sort.bam	sort.bam
