#!/data/software/bin/python
# coding=utf-8

'''
    文件夹格式：
    cnv_detection
        copy_number_detection.py
        input_CNV
            baselines
                build_baseline.py
                input_baseline
                    baseline_samples
                    target.bed
                    ucsc.hg19.fasta
                output_baseline
                    baseline_targetcoverage
                    baseline.txt
            bam
            vcf
            coverage
        output_CNV
            res_CNV
            visual_AF
            visual_CNV
'''

import os, time, vcf
import numpy as np
import pandas as pd
import argparse
import scipy.stats as st
from distutils.spawn import find_executable

# 设置pandas显示所有列
pd.set_option('display.max_columns', None)
# disable chained assignments
pd.options.mode.chained_assignment = None

# 设置参数
def get_parser():
    # 设置脚本参数
    parser = argparse.ArgumentParser(description='copy number variant detection')

    working_dir = os.path.split(os.path.realpath(__file__))[0]

    # 默认的参数变量名是-或者--后面的值，也可以用dest参数指定，在调用的时候用args.dest调用
    parser.add_argument('-r', '--refbaseline', dest="baseline",type=str,required=True,
                        help='copy number baseline file')
    parser.add_argument('-b', '--bam', dest="bam",type=str,required=True,
                        help='bam file of tumor sample')                                        
    parser.add_argument('-o', '--output', dest="output",type=str,required=True,
                        help='file path of output results')
    parser.add_argument('-v', '--vcf', dest="vcf",type=str,required=True, 
                        help='vcf file of tumor sample')                    
    parser.add_argument('-p', '--plotting',dest="visualization",type=bool,required=False, default=False,
                        help='whether to visualize the distribution of sample AF') 

    return parser

# 读取vcf文件, 返回一个字典，保存各个样本过滤后的AF信息
def vcf2dict(sample_vcf):
    # 创建一个Dict，保存结果
    dict_vcf = {}
    list_vcf = []

    # 读取vcf
    vcf_reader = vcf.Reader(filename=f"{sample_vcf}")

    # 根据VCF中的信息过滤
    for record in vcf_reader:
        var_type = record.INFO['ANN'][0].split("|")[1]
        var_AF = record.INFO['AF'][0]

        # 只保留错义突变，设置突变频率区间
        if 'missense_variant' == var_type and var_AF >= 0.011 and var_AF <= 0.329:
            list_vcf.append(var_AF)

        # 保存结果到字典中
        dict_vcf.update({os.path.split(sample_vcf)[1].split(".vcf")[0]: list_vcf})

    return dict_vcf

# 输入一个AF的数组，估计TP
# 得到两个结果，一个均值一个中位数
# 如果需要修改肿瘤细胞浓度计算，修改这个函数
def af2tp(aflist):
    tp = list()

    # 取百分位数
    percentile = np.percentile(aflist, (80,95), interpolation='midpoint')
    aflist = [x for x in aflist if x >= percentile[0]]

    #tp1 = round(np.mean(aflist)*1.2, 4)
    tp2 = round(np.median(aflist)*1.2, 4)
    #tp.append(tp1)
    tp.append(tp2)
    
    return tp2

# 获取picard执行程序
def picard_exe():
    path = os.path.join(os.path.split(os.path.realpath(__file__))[0], '../../../software/')
    picard_exe = str(find_executable("picard-2.25.0.jar", path = path))

    return picard_exe

# 计算肿瘤样本中的coverage
def calc_tumorcov(reference, targetbed, tumor_sample, output):
    picard = picard_exe()

    # 计算tumor sample的Hsmetrics
    path = os.path.split(tumor_sample)[0]
    prefix = os.path.split(tumor_sample)[1].split(".")[0]

    os.system(f"java -jar {picard} CollectHsMetrics I={tumor_sample} \
        O={output}/out_hsmetrics{prefix}.txt \
        R={reference} PER_TARGET_COVERAGE={output}/{prefix}.targetcov \
        BAIT_INTERVALS= {os.path.splitext(targetbed)[0]}.interval_list \
        TARGET_INTERVALS= {os.path.splitext(targetbed)[0]}.interval_list")

def extract_tumorcov(tumorcoverage):
    results = dict()
    bins = dict()

    cov = pd.read_csv(tumorcoverage, sep = "\t", usecols=[0,1,2,3,4,6,7])

    mean_cov = cov['mean_coverage'].tolist()
    norm_cov = cov['normalized_coverage'].tolist()
    prefix = os.path.split(tumorcoverage)[-1].split(".")[0]
    results[f"{prefix}_meancov"] = mean_cov
    results[f"{prefix}_normcov"] = norm_cov

    bins["binsize"] = cov["length"].values
    index_bins = cov['chrom'].str.cat(cov['start'].astype(str),sep=":").str.cat(cov['end'].astype(str),sep="-").str.cat(cov['name'],sep=":").tolist()

    res = pd.DataFrame(results, index=index_bins)
    df_bins = pd.DataFrame(bins, index=index_bins)
    
    return res, df_bins, prefix

def filter(sample):
    zscore = sample.copy()
    
    # 得到zscore矩阵
    zscores = (zscore-zscore.mean())/zscore.std()
    
    # 得到矩阵中的正常值
    condition = zscores.abs()<=2

    # 如果全部都是正常值，返回结果
    if condition.all():
        # print("已找到所有异常值！")
        # 返回布尔Series，指明哪些index对应的是正常值
        return condition
    else:
        # 如果存在异常值，保留正常值再重新过滤
        # print("过滤中...")
        return filter(sample[condition])

# 对肿瘤样本的normalized coverage重新标准化
def renorm(df_coverage, df_normalized_coverage, binsizes, samplename):
    # 对picard的结果重新标准化，得到收敛后的结果
    binsizes = binsizes.iloc[:,0]

    # 修改df_coverage和df_normalized_coverage的列名为样本名称
    df_coverage.columns = [f"{samplename}"]
    df_normalized_coverage.columns = [f"{samplename}"]

    for sample in df_normalized_coverage.columns:
        condition = filter(df_normalized_coverage[sample])
        
        total = df_coverage.loc[condition.index,sample] * binsizes[condition.index]
        mean_rd = total.sum()/binsizes[condition.index].sum()

        df_renorm = df_coverage.loc[condition.index,sample]/mean_rd

        df_normalized_coverage.loc[:,sample][condition.index] = df_renorm[df_renorm.index]


    df_tumor = pd.DataFrame(data=df_normalized_coverage, index = df_norm.index)
    df_tumor = df_tumor.fillna(value=0)

    return df_tumor

def calc_copynumbers(baseline, tumor, estimateTP, output):
    #print(tumor.columns.values[0])
    df_baseline = pd.read_csv(baseline, sep="\t", header=0, index_col=0)
    df_merge = pd.concat([df_baseline, tumor], axis = 1)
    
    zscore = (df_merge[tumor.columns.values[0]]-df_merge["baseline_mean"])/df_merge["baseline_std"]
    # 计算p值
    df_merge = df_merge.assign(pvalue = st.norm.sf(abs(zscore)))
    # 计算cn值
    df_merge = df_merge.assign(copy_numbers = (2*(df_merge[tumor.columns.values[0]]/df_merge["normalized_coverage"]) - 2*(1-float(estimateTP)))/float(estimateTP))

    # 调整列顺序
    df_merge = df_merge[['binsize', 'copy_numbers','pvalue',f'{tumor.columns.values[0]}', 'normalized_coverage', 'mean_coverage', 'baseline_mean','baseline_std','baseline_mad']]
    # 修改列名
    df_merge.columns = ['binsize', 'copy_numbers','pvalue',f'{tumor.columns.values[0]}_normcov','baseline_normcov', 'baseline_meancov', 'baseline_mean','baseline_std','baseline_mad']

    # 拷贝数为负数的修改为0
    df_merge.loc[df_merge["copy_numbers"]<0, "copy_numbers"] = 0
    # p值大于1的修改为1
    df_merge.loc[df_merge["pvalue"]>1, "pvalue"] = 1

    # 计算基因的CN值
    df_geneCN = df_merge

    # 过滤p值大于0.01的bin
    df_geneCN = df_geneCN.loc[df_geneCN['pvalue'] < 0.01]

    df_index = df_geneCN.index
    gene_name = [x.split(":")[-1] for x in df_index]

    # 加上基因名称用于分组计算
    df_geneCN.insert(0,"gene", gene_name)

    # 计算每个bin加权后的copy number
    df_geneCN.insert(2,"weights", df_geneCN['copy_numbers']*df_geneCN.groupby('gene')['binsize'].transform(lambda x:(x/x.sum())))

    # 计算基因的拷贝数
    gene_copynumbers = df_geneCN.groupby('gene')['weights'].sum()

    # 保存target bin的检测结果和gene copy number的检测结果
    df_merge.to_csv(f"{output}/resbin_{tumor.columns.values[0]}.txt", sep="\t", index = True, index_label="bin_info")
    
    df_genes = pd.DataFrame(gene_copynumbers)
    df_genes.to_csv(f"{output}/resgene_{tumor.columns.values[0]}.txt", sep="\t", index = True, header=False)

    

if __name__ == "__main__":
    start = time.time()

    # 获取参数
    parser = get_parser()
    args = parser.parse_args()

    # 获取当前脚本所在路径
    working_dir = os.path.split(os.path.realpath(__file__))[0]

    # 估算肿瘤细胞浓度
    res_vcf = vcf2dict(args.vcf)
    # dictTP的值就是估算的tumor purity, 分别是均值估算和中位数估算
    dict_TP = {}
    
    for key, value in res_vcf.items():
        dict_TP.update({key:af2tp(value)})
    
    prefix = os.path.split(args.bam)[1].split(".")[0]

    # 计算肿瘤样本的bin coverage\
    if os.path.exists(f"{args.output}/{prefix}.targetcov"):
        pass
    else:
        calc_tumorcov(f"{working_dir}/../../../database/gatk-bundle-hg19/ucsc.hg19.fasta", \
        f"{working_dir}/../bed/NanOnco_Plus_Panel_v2.0_Covered_hg19.bed", \
        f"{args.bam}", args.output)
    
    res, binsizes, samplename = extract_tumorcov(f"{args.output}/{prefix}.targetcov")
    
    df_norm = res.iloc[:,[i%2==1 for i in range(len(res.columns))]]
    df_cov = res.iloc[:,[i%2==0 for i in range(len(res.columns))]]

    # 得到重新标准化的肿瘤Normalized coverage
    df_tumor = renorm(df_cov, df_norm, binsizes, samplename)

    # 计算copy numbers 
    purity = list(dict_TP.values())[0]
    print(f"Estimate Tumor Purity: {purity*100}%")
    calc_copynumbers(args.baseline, df_tumor, purity, args.output)

    end = time.time() 
    print("Finish!")
    print(f"Time consumption: {end-start}s")
