#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-10-14 
# Version: 0.2.0.0
# Description: pipeline for bam QC
#########################################################################
import os
import sys
from multiprocessing import Pool
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def RemoveDup(input_bam):
    cmd_rmdup = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        picard, 'MarkDuplicates',
        '-TMP_DIR', f'{dirprefix}_tmp',
        '-I', input_bam,
        '-O', f'{dirprefix}.MarkDuplicates.bam',
        '-M', f'{dirprefix}.MarkDuplicates.metrics',
        '--REMOVE_DUPLICATES true'
    ])
    bac.RunningProcess(cmd_rmdup, otdir, otprefix, 'Notice:\tpicard EstimateLibraryComplexity coculate bed region dup.ratio running messages')

def DuplicateRatio(samtools, samtools_threads, input_bam, bed, umi):
    cmd_dupratio = " ".join([
        samtools, 'view', '--threads', samtools_threads,
        '-h', '-L', bed, input_bam, '|',
        java, f'-jar -Xmx{java_mem}G',
        picard, 'EstimateLibraryComplexity',
        '-I', '/dev/stdin',
        '--MIN_MEAN_QUALITY', '10',
        '-O', f'{dirprefix}.bed.EstimateLibraryComplexity.metrics'
    ])
    # samtools view -h --threads 24 -L bed *.sort.bam |java -jar picard-2.25.0.jar EstimateLibraryComplexity -I /dev/stdin -O tmp
    if umi == 'true':
        cmd_dupratio = cmd_dupratio+' --BARCODE_TAG RX'
    bac.RunningProcess(cmd_dupratio, otdir, otprefix, 'Notice:\tpicard EstimateLibraryComplexity coculate bed region dup.ratio running messages')

def InsertSize(input_bam):
    cmd_insertsize = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        picard, 'CollectInsertSizeMetrics',
        '-I', input_bam,
        '-W', '600',
        '-O', f'{dirprefix}.CollectInsertSizeMetrics.xls',
        '-H', f'{dirprefix}.CollectInsertSizeMetrics.pdf'
    ])
    bac.RunningProcess(cmd_insertsize, otdir, otprefix, 'Notice:\tpicard CollectInsertSizeMetrics running messages')

def BedToList(bed, hg19_dict):
    cmd_bedtolist = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        picard, 'BedToIntervalList',
        '-I', bed,
        '-O', f'{dirprefix}.BedToIntervalList.xls',
        '-SD', hg19_dict
    ])
    bac.RunningProcess(cmd_bedtolist, otdir, otprefix, 'Notice:\tpicard BedToIntervalList running messages')

def CollectHs(input_bam, hg19_fa, times):
    cmd_collecths = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        picard, 'CollectHsMetrics',
        '-I', input_bam,
        '-O', f'{dirprefix}.CollectHsMetrics_{times}.xls',
        '-R', hg19_fa,
        '--BAIT_INTERVALS', f'{dirprefix}.BedToIntervalList.xls', 
        '--TARGET_INTERVALS', f'{dirprefix}.BedToIntervalList.xls', 
        '--PER_BASE_COVERAGE', f'{dirprefix}.Hs_perbase_coverage_{times}.xls', 
        '--PER_TARGET_COVERAGE', f'{dirprefix}.Hs_pertarget_coverage_{times}.xls'
    ])
    bac.RunningProcess(cmd_collecths, otdir, otprefix, 'Notice:\tpicard CollectHsMetrics running messages')

def QCTable(sample_type, umi_result_txt = None):
    covnum = {
        'gDNA':{'raw':'500', 'rmdup':'100'}, 
        'ctDNA':{'raw':'5000', 'rmdup':'1000'}, 
        'WES':{'raw':'100', 'rmdup':'50'}
    }
    dup_rst = f'{dirprefix}.bed.EstimateLibraryComplexity.metrics'
    insert = f'{dirprefix}.CollectInsertSizeMetrics.xls'
    raw_hs = f'{dirprefix}.CollectHsMetrics_raw.xls'
    rmdup_hs = f'{dirprefix}.CollectHsMetrics_rmdup.xls'
    
    # tags
    dup_ratio = str(round(float(os.popen(f'cat {dup_rst} |head -n 8 |tail -n 1').read().strip().split('\t')[-2])*100, 2)) # PERCENT_DUPLICATION
    mean_insert_size = str(round(float(os.popen(f'cat {insert} |head -n 8 |tail -n 1').read().strip().split('\t')[5]), 2)) # MEAN_INSERT_SIZE
    median_insert_size = os.popen(f'cat {insert} |head -n 8 |tail -n 1').read().strip().split('\t')[0] # MEDIAN_INSERT_SIZE
    rmdup_mean_depth = str(round(float(os.popen(f'cat {rmdup_hs} |head -n 8 |tail -n 1').read().strip().split('\t')[9]), 2)) # MEAN_BAIT_COVERAGE
    qc_list = os.popen(f'cat {raw_hs} |head -n 8 |tail -n 1').read().strip().split('\t')
    rmdup_qc_list = os.popen(f'cat {rmdup_hs} |head -n 8 |tail -n 1').read().strip().split('\t')

    total_reads = '{:,}'.format(int(qc_list[22])) # TOTAL_READS
    mapped_ratio = str(round(float(qc_list[32])*100, 2)) # PCT_PF_UQ_READS_ALIGNED
    target_ratio = str(round(float(qc_list[6])*100, 2)) # PCT_SELECTED_BASES
    mean_depth = str(round(float(qc_list[9]), 2)) # MEAN_BAIT_COVERAGE
    raw_coverage_ratio = str(round(float(qc_list[45])*100, 2)) # PCT_TARGET_BASES_1X
    raw_coverage_100_ratio = str(round(float(qc_list[52])*100, 2)) # PCT_TARGET_BASES_100X
    raw_coverage_500_ratio = str(round(float(qc_list[54])*100, 2)) # PCT_TARGET_BASES_500X
    raw_coverage_5000_ratio = str(round(float(qc_list[57])*100, 2)) # PCT_TARGET_BASES_5000X
    # raw_coverage_25000_ratio = str(round(float(qc_list[59])*100, 2)) # PCT_TARGET_BASES_25000X
    # rmdup_coverage_100_ratio = str(round(float(rmdup_qc_list[52])*100, 2)) # PCT_TARGET_BASES_100X
    # rmdup_coverage_1000_ratio = str(round(float(rmdup_qc_list[55])*100, 2)) # PCT_TARGET_BASES_1000X
    # rmdup_coverage_2500_ratio = str(round(float(rmdup_qc_list[56])*100, 2)) # PCT_TARGET_BASES_2500X

    covresult = {
        'gDNA':{'raw':raw_coverage_500_ratio, 'rmdup':raw_coverage_500_ratio}, 
        'ctDNA':{'raw':raw_coverage_5000_ratio, 'rmdup':raw_coverage_5000_ratio}, 
        'WES':{'raw':raw_coverage_100_ratio, 'rmdup':raw_coverage_100_ratio}
    }
    # umi freq
    if bac.CheckFile(umi_result_txt):
        tmp1 = os.popen("cat %s |awk '{print $1}' |sort -n |uniq -c |awk '{if($2==1){equal+=$1*$2}else{none+=$1*$2}}END{print equal/(equal+none)*100}'"%umi_result_txt).read() #这里不是只统计umi出现的次数，应该统计umi出现次数的fragment条数(umi出现次数)，因为只统计了read1和read2其中的一条；
        tmp5 = os.popen("cat %s |awk '{print $1}' |sort -n |uniq -c |awk '{if($2<=5){less5+=$1*$2}else{over5+=$1*$2}}END{print less5/(less5+over5)*100}'"%umi_result_txt).read()
        umi_equal1 = str(round(float(tmp1.strip()), 2))
        umi_less5 = str(round(float(tmp5.strip()), 2))
        # os.system("cat %s |awk '{print $1}' |sort -n |uniq -c |awk 'BEGIN{print \"UMI_freq\tread_count\"}{print $2\"\t\"$2*$1}' > %s.umi_frequence.xls"%(umi_result_txt, dirprefix))
    else:
        umi_equal1 = 'NA'
        umi_less5 = 'NA'
    # fold 80
    try:
        fold_80 = str(round(float(qc_list[44]), 2)) # FOLD_80_BASE_PENALTY (这个值计算会出现?号)
    except ValueError:
        fold_80 = 'NA' # qc_list[44]
    # uniformity
    raw_perbscov = f'{dirprefix}.Hs_perbase_coverage_raw.xls'
    overcount = os.popen("cat %s |awk 'NR>1{if($4>%s*0.2)print $0}' |wc -l"%(raw_perbscov, mean_depth)).read()
    allcount = os.popen("cat %s |awk 'NR>1{print $0}' |wc -l"%raw_perbscov).read()
    uniformity = str( round( float(overcount.strip()) / float(allcount.strip()) * 100, 2) )

    opt_list = []
    opt_list.append("\t".join([
        'sample.id', 'sample.type', 
        'mapped.ratio(%)', 'mean.insert.size', 'median.insert.size', 
        'target.ratio(%)', 'mean.depth', 
        'coverage.ratio(%)', f"{covnum.get(sample_type).get('raw')}x.coverage.ratio(%)", 
        'fold-80', 'uniformity(%)', 
        'bed.duplicated.ratio(%)', 'rmdup.mean.depth', 
        'UMI.freq.equal1(%)', 'UMI.freq.less5(%)' # 'rmdup.{}x.coverage.ratio(%)'.format(covnum.get(sample_type).get('rmdup')), 
    ])+'\n')
    opt_list.append("\t".join([
        otprefix, sample_type, 
        mapped_ratio, mean_insert_size, median_insert_size, 
        target_ratio, mean_depth, 
        raw_coverage_ratio, covresult.get(sample_type).get('raw'), 
        fold_80, uniformity, 
        dup_ratio, rmdup_mean_depth, 
        umi_equal1, umi_less5 # covresult.get(sample_type).get('rmdup'), 
    ])+'\n')
    with open(f'{dirprefix}.bam_QC.xls', 'w') as otopen:
        otopen.writelines(opt_list)

    os.system("paste %s.Hs_pertarget_coverage_raw.xls %s.Hs_pertarget_coverage_rmdup.xls |cut -f -8,13-14,20-22,27-28 |awk 'BEGIN{print \"chrom\\tstart\\tend\\tlength\\tname\\traw_%%gc\\traw_mean_coverage\\traw_normalized_coverage\\traw_pct_0x\\traw_read_count\\trmdup_%%gc\\trmdup_mean_coverage\\trmdup_normalized_coverage\\trmdup_pct_0x\\trmdup_read_count\"}NR>1{print $0}' > %s.each_target_result.xls"%(dirprefix, dirprefix, dirprefix))

    # umi type
    # samtools view --threads 16 -L NanOnco_Plus_Panel_v2.0_Covered_hg19.bed -f 64 ZLC210252A.group.bam |awk '{for(i=1;i<=NF;i++){if($i~/MI/)print $i}}' |awk -F "/" '{print $1}' |uniq -c > uniq.result
    # cat uniq.result |awk '{print $1}' |sort -n |uniq -c |awk '{if($2<=5){less5+=$1}else{over5+=$1}}END{print less5/(less5+over5)*100}'

def Running(sort_bam, bed, cfg, output_dir, output_prefix, sample_type, removedup_bam, umi, group_bam, consensus_mapped_bam, only_QC): #, markdup_metrics
### prepare
    global java, java_mem, picard, dirprefix, otdir, otprefix
    otdir = output_dir
    otprefix = output_prefix
    config = bac.ResolveConfig(cfg)
    dirprefix = os.path.join(output_dir, output_prefix)
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    picard = config['software']['picard']
    samtools = config['software']['samtools']
    samtools_threads = config['parameter']['samtools_threads']
    hg19_fa = config['database']['hg19_fa']
    hg19_dict = config['database']['hg19_dict']
### cases
    # 目前计算duplicate.ratio都计算靶区域的dup率
    # case1: 当存在sort_bam, removedup_bam这三个输入参数时, 做如下分析: sort_bam removedup, removedup_bam insertsize; sort_bam collecths; removedup_bam collecths
    # case2: 当存在sort_bam, group_bam, consensus_mapped_bam, umi=True这四个输入参数时, 做如下分析: 先做group_bam sort, group_sort_bam EstimateLibraryComplexity, group_sort_bam collecths之后再做consensus_mapped_bam insertsize, consensus_mapped_bam collecths
    # case3: 当只有sort_bam输入参数时, 做如下分析: 先做sort_bam collecths, sort_bam removedup, 之后再做sort_bam_rmdup insertsize, sort_bam_rmdup collecths
### strategy of cases
    if bac.CheckFile(sort_bam) and bac.CheckFile(removedup_bam) and umi == 'false': #case 1
        strategy = 'clear'
    elif bac.CheckFile(group_bam) and bac.CheckFile(consensus_mapped_bam) and umi == 'true': #case 2
        strategy = 'umi'
    elif bac.CheckFile(sort_bam) and umi == 'false': #case 3
        strategy = 'all'
    elif only_QC == 'true':
        strategy = 'only'
    else:
        print ('types:\tsort_bam removedup_bam umi=false; group_bam consensus_mapped_bam umi=true; sort_bam umi=false')
        exit()
### running
    BedToList(bed, hg19_dict)
### clear
    if strategy == 'clear':
        pool = Pool(processes=4)
        pool.apply_async(InsertSize, args=(removedup_bam, ))
        pool.apply_async(CollectHs, args=(sort_bam, hg19_fa, 'raw'))
        pool.apply_async(CollectHs, args=(removedup_bam, hg19_fa, 'rmdup'))
        pool.apply_async(DuplicateRatio, args=(samtools, samtools_threads, sort_bam, bed, umi))
        pool.close()
        pool.join()
        QCTable(sample_type)
        otstr = f'QC analysis with non of UMI\n{sort_bam}\n{removedup_bam}\nused\n'
        with open(f'{dirprefix}_file_used.txt','w') as otopen:
            otopen.write(otstr)
### umi
    elif strategy == 'umi':
        cmd_umifreq = " ".join([
            samtools, f'view --threads {samtools_threads}', 
            '-L', bed, '-f 64', group_bam, 
            "|awk '{for(i=1;i<=NF;i++){if($i~/MI/)print $i}}' |awk -F \"/\" '{print $1}' |uniq -c > %s.umi.result.txt"%dirprefix
        ])
        umi_result_txt = f'{dirprefix}.umi.result.txt'
        cmd_sort = " ".join([
            samtools, f'sort --threads {samtools_threads}', 
            group_bam, '-o', f'{dirprefix}.group.sort.bam'
        ])
        bac.RunningProcess(cmd_sort, output_dir, output_prefix, 'Notice:\tgroup bam samtools sort running messages')
        group_sort_bam = f'{dirprefix}.group.sort.bam'
        # 后面如果改进了umi的分析方法会将所有split的group_bam合并成为sort的group_bam，修改group_bam的参数为group_sort_bam
        pool = Pool(processes=3)
        pool.apply_async(DuplicateRatio, args=(samtools, samtools_threads, group_sort_bam, bed, umi))
        pool.apply_async(CollectHs, args=(group_sort_bam, hg19_fa, 'raw'))
        pool.apply_async(bac.RunningProcess, args=(cmd_umifreq, output_dir, output_prefix, 'Notice:\tgroup bam get *.umi.result.txt running messages'))
        pool.close()
        pool.join()
        markdup_metrics = f'{dirprefix}.bed.EstimateLibraryComplexity.metrics'

        pool = Pool(processes=2)
        pool.apply_async(CollectHs, args=(consensus_mapped_bam, hg19_fa, 'rmdup'))
        pool.apply_async(InsertSize, args=(consensus_mapped_bam, ))
        pool.close()
        pool.join()
        QCTable(sample_type, umi_result_txt)
        otstr = f'QC analysis with UMI\n{group_bam}\n{consensus_mapped_bam}\nused\n'
        with open(f'{dirprefix}_file_used.txt','w') as otopen:
            otopen.write(otstr)
### all
    elif strategy == 'all':
        pool = Pool(processes=3)
        pool.apply_async(DuplicateRatio, args=(samtools, samtools_threads, sort_bam, bed, umi))
        pool.apply_async(RemoveDup, args=(sort_bam, ))
        pool.apply_async(CollectHs, args=(sort_bam, hg19_fa, 'raw'))
        pool.close()
        pool.join()
        rmdup_bam = f'{dirprefix}.MarkDuplicates.bam'

        pool = Pool(processes=2)
        pool.apply_async(CollectHs, args=(rmdup_bam, hg19_fa, 'rmdup'))
        pool.apply_async(InsertSize, args=(rmdup_bam, ))
        pool.close()
        pool.join()
        QCTable(sample_type)
        otstr = f'QC analysis with non of UMI\n{sort_bam}\nused\n'
        with open(f'{dirprefix}_file_used.txt','w') as otopen:
            otopen.write(otstr)
### only
    elif strategy == 'only':
        dup_file = bac.CheckFile(f'{dirprefix}.bed.EstimateLibraryComplexity.metrics')
        insert_file = bac.CheckFile(f'{dirprefix}.CollectInsertSizeMetrics.xls')
        hsraw_file = bac.CheckFile(f'{dirprefix}.CollectHsMetrics_raw.xls')
        hsdup_file = bac.CheckFile(f'{dirprefix}.CollectHsMetrics_rmdup.xls')
        umi_result_txt = f'{dirprefix}.umi.result.txt'
        if dup_file and insert_file and hsraw_file and hsdup_file:
            if brc.GetSize(umi_result_txt) != 0:
                QCTable(sample_type, umi_result_txt)
            else:
                QCTable(sample_type)
        else:
            print (f"{dup_file} or {insert_file} or {hsraw_file} or {hsdup_file} file not found!")
            exit()

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Quality Control')
    parser.add_argument('-b', '--sort_bam', required=None,
                        help="if methods not with UMI and downsample is PASS, this file is '*.sort.downsample_num.bam';\n and other types this file is '*.sort.bam'; default=None")
    parser.add_argument('-d', '--bed', required=True,
                        help='panel design bed file')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-tp', '--sample_type', type=str, choices=['gDNA', 'ctDNA', 'WES'], default='gDNA',
                        help='sample type; gDNA(FFPE/WBC) ctDNA or WES, default=gDNA')
    parser.add_argument('-rb', '--removedup_bam', required=None,
                        help='this parameter will input "*.rmdup.bam" file or "*.rmdup_IndelRealigner.bam"; default=None')
    parser.add_argument('-umi', '--umi', type=str, choices=['true', 'false'], default='false',
                        help='count dup ratio with UMI or not, default=false')
    parser.add_argument('-gb', '--group_bam', default=None,
                        help='*.group.bam file with UMI, this parameter is used with umi is "True", default=None')
    parser.add_argument('-sb', '--consensus_mapped_bam', default=None,
                        help='*.consensus.mapped.bam file with UMI, this parameter is used with umi is "True", default=None')
    parser.add_argument('-only', '--only_QC', type=str, choices=['true', 'false'], default='false',
                        help='do not running picard or samtools running, based on *.xls file getting *.bam_QC.xls, default=false')
    args = parser.parse_args()
    Running(args.sort_bam, args.bed, args.cfg, args.output_dir, args.output_prefix, args.sample_type, args.removedup_bam, args.umi, args.group_bam, args.consensus_mapped_bam, args.only_QC)

# QC指标				
# 输出指标	输出文件中选择的参数	picard模块	非umi的bam文件	umi后的bam文件
# total.reads	BAIT_TERRITORY	CollectHsMetrics	sort.bam	sort.bam
# mapped.ratio	PCT_PF_UQ_READS_ALIGNED	CollectHsMetrics	sort.bam	sort.bam
# target.ratio	PCT_SELECTED_BASES	CollectHsMetrics	sort.bam	sort.bam
# mean.insert.size	MEDIAN_INSERT_SIZE	CollectInsertSizeMetrics	sort.bam	sort.bam
# duplicated.ratio	PERCENT_DUPLICATION	MarkDuplicates/EstimateLibraryComplexity	sort.bam	group.sort.bam(GroupReadsByUmi之后输出的bam文件)
# mean.depth	MEAN_BAIT_COVERAGE	CollectHsMetrics	sort.bam	sort.bam
# rmdup.mean.depth	MEAN_BAIT_COVERAGE	CollectHsMetrics	rmdup.bam	consensus.mapped.bam(最终矫正一致性序列后的bam文件)
# coverage.ratio	PCT_TARGET_BASES_1X	CollectHsMetrics	sort.bam	sort.bam
# coverage.1000.ratio	PCT_TARGET_BASES_1000X	CollectHsMetrics	sort.bam	sort.bam
# fold-80	FOLD_80_BASE_PENALTY	CollectHsMetrics	sort.bam	sort.bam
