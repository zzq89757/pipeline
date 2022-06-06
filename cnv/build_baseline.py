#!/data/software/bin/python
# coding=utf-8

import os, time
import numpy as np
import pandas as pd
import argparse
from distutils.spawn import find_executable

'''
    建基线
'''

# 忽略计算时分母为0的警告信息
np.seterr(invalid='ignore')

def get_parser():
    # 设置参数
    parser = argparse.ArgumentParser(description='baseline building workflow')

    # 获取当前脚本的工作路径
    working_dir = os.path.split(os.path.realpath(__file__))[0]

    parser.add_argument('-s', '--baselineSamples', dest="samples",type=str,required=False, default=f"{working_dir}/input_baseline/baseline_samples",
                        help='file path of samples for building baseline')
    parser.add_argument('-r', '--reference', dest="reference",type=str,required=False, default=f"{working_dir}/input_baseline/ucsc.hg19.fasta",
                        help='genome reference fasta file')
    parser.add_argument('-b', '--targetbed', dest="bed",type=str,required=False, default=f"{working_dir}/input_baseline/NanOnco_Plus_Panel_v2.0_Covered_hg19.bed",
                        help='target bed file')


    return parser

# 获取picard执行程序
def picard_exe():
    path = "/data/software/src/"
    picard_exe = str(find_executable("picard.jar", path = path))

    return picard_exe

def getres_picard(reference, targetbed, samples):
    # 用picard计算hsmetrics
    picard = picard_exe()

    # 根据fasta文件创建picard dict
    if os.path.exists(f"{os.path.splitext(reference)[0]}.dict"):
        pass
    else:
        os.system(f"java -jar {picard} CreateSequenceDictionary R={reference} O={os.path.splitext(reference)[0]}.dict")

        # 把targetbed转换为picard的IntervalList
        os.system(f"java -jar {picard} BedToIntervalList I={targetbed} O={os.path.splitext(targetbed)[0]}.interval_list \
            SD={os.path.splitext(reference)[0]}.dict")


    
    # 计算HsMetrics
    baseline_samples = os.listdir(samples)
    baseline_samples = [file for file in baseline_samples if file.endswith("removedup_IndelRealigner.bam")]
    #print(baseline_samples)

    
    for sample in baseline_samples:
        print(f"正在计算{sample}的coverage信息。。。")

        os.system(f"java -jar {picard} CollectHsMetrics I={samples}/{sample} O=out_hsmetrics{sample}.txt \
        R={reference} PER_TARGET_COVERAGE={samples}/../../output_baseline/baseline_targetcoverage/{sample}.targetcov \
        BAIT_INTERVALS= {os.path.splitext(targetbed)[0]}.interval_list \
        TARGET_INTERVALS= {os.path.splitext(targetbed)[0]}.interval_list")
    
def extract_cov(coverages):
    # 从picard结果中提取mean_cov和norm_cov等信息
    
    results = dict()
    bins = dict()
    sampleNames = []
    
    for (dirpath, dirnames, filenames) in os.walk(coverages):
        for filename in filenames:
            cov = pd.read_csv(os.path.join(dirpath, filename), sep = "\t", usecols=[0,1,2,3,4,6,7])
            
            mean_cov = cov['mean_coverage'].tolist()
            norm_cov = cov['normalized_coverage'].tolist()
            prefix = filename.split(".")[0]
            results[f"{prefix}_meancov"] = mean_cov
            results[f"{prefix}_normcov"] = norm_cov
            sampleNames.append(prefix)
            
            bins["binsize"] = cov["length"].values
            index_bins = cov['chrom'].str.cat(cov['start'].astype(str),sep=":").str.cat(cov['end'].astype(str),sep="-").str.cat(cov['name'],sep=":").tolist()
            
    res = pd.DataFrame(results, index=index_bins)
    df_bins = pd.DataFrame(bins, index=index_bins)

    return res, df_bins, sampleNames

def filter_norm(sample):
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
        return filter_norm(sample[condition])

def renorm(df_coverage, df_normalized_coverage, binsizes, samplenames):
    # 对picard的结果重新标准化，得到收敛后的结果
    binsizes = binsizes.iloc[:,0]
    print(binsizes)

    # 修改df_coverage和df_normalized_coverage的列名为样本名称
    df_coverage.columns = samplenames
    df_normalized_coverage.columns = samplenames


    for sample in df_normalized_coverage.columns:
        condition = filter_norm(df_normalized_coverage[sample])
        
        total = df_coverage.loc[condition.index,sample] * binsizes[condition.index]
        mean_rd = total.sum()/binsizes[condition.index].sum()

        df_renorm = df_coverage.loc[condition.index,sample]/mean_rd

        df_normalized_coverage.loc[:,sample][condition.index] = df_renorm[df_renorm.index]
    
    mean_coverages = []
    medians = []
    meanvalues = []
    standarddeviations = []
    madvalues = []
    
    # 过滤bin中的异常，取剩余值的median
    for bin in df_normalized_coverage.index:
        condition = filter_norm(df_normalized_coverage.loc[bin])

        median = np.median(df_normalized_coverage.loc[bin, condition.index])
        std = np.std(df_normalized_coverage.loc[bin, condition.index])
        mean_norm = np.mean(df_normalized_coverage.loc[bin, condition.index])
        mean_cov = np.mean(df_coverage.loc[bin, condition.index])
        mad = np.median(abs(df_normalized_coverage.loc[bin, condition.index] - median))
        mean_coverages.append(mean_cov)
        medians.append(median)
        meanvalues.append(mean_norm)
        standarddeviations.append(std)
        madvalues.append(mad)
        
    # 参考基线的内容
    baseline = {
        "binsize": binsizes, 
        "mean_coverage" : mean_coverages, 
        "normalized_coverage" : medians,
        "baseline_mean": meanvalues,
        "baseline_std" : standarddeviations,
        "baseline_mad" : madvalues
    }

    df_baseline = pd.DataFrame(data=baseline, index = df_norm.index)
    df_baseline = df_baseline.fillna(value=0)

    return df_baseline

if __name__ == "__main__":
    start = time.time()

    # 获取当前脚本的工作路径
    working_dir = os.path.split(os.path.realpath(__file__))[0]
    
    # 获取输入参数
    parser = get_parser()
    args = parser.parse_args()

    # 第一步：picard计算hsmetrics
    # 测试通过
    # getres_picard(f"{args.reference}", f"{args.bed}", f"{args.samples}")

    # 第二步：从coverage结果中提取相应的信息并重新标准化得到baseline结果
    res, bins, samplenames = extract_cov(f"{working_dir}/output_baseline/baseline_targetcoverage")

    # 读取DF中的标准化列和均值列，分别保存
    df_norm = res.iloc[:,[i%2==1 for i in range(len(res.columns))]]
    df_cov = res.iloc[:,[i%2==0 for i in range(len(res.columns))]]

    # 对picard的结果重新进行标准化
    df_baseline = renorm(df_cov, df_norm, bins, samplenames)

    df_baseline.to_csv(f"{working_dir}/output_baseline/df_baseline.txt", sep="\t", index=True, index_label="bin_info", float_format='%.4f')


    end = time.time()
    print("Finish!")
    print(f"Time consumption: {end-start}s")