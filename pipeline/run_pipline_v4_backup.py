#coding:utf8
#########################################################################
# Author: ZhouNan 
# Created Time: 2021-11-23
# Version: 0.2.0.0
# Description: pipeline for all of the analysis process
#########################################################################
import os
import sys
import json
import pandas as pd #reading excel file need xlrd==1.2.0; pip install xlrd==1.2.0
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import brc

def versions(directory, keystring='parameter_configure_file'): 写入brc
    """用于展示最新版本及判断历史版本信息"""
    allv_list=[]
    version_dict={}
    files=next(os.walk(directory))[2]
    for each in files:
        if keystring in each:
            try:
                pointv=re.findall(r'(\d+)\.(\d+)\.(\d+)\.(\d+)', each)[0]
                allv_list.append('v'+".".join(pointv))
                allnum=0
                for x in list(range(len(pointv)))[::-1]:
                    allnum+=int(pointv[-x-1])*100**x
                version_dict.update({allnum:each})
            except IndexError:
                version_dict.update({0:each})
    compare_num=0
    lateset=keystring
    if version_dict == {}:
        raise KeyError('%s keystring not found in %s directory all files'%(keystring, direcotry))
    else:
        for k,v in version_dict.items():
            if k >= compare_num:
                lateset=v
    return lateset, allv_list


def ConfigureHandle():
    """用于处理配置文件的模块"""
    pass

def LibraryResult(match_table):
    """用于处理生信对接表的模块"""
    try:
        df = pd.read_excel(match_table)
    except ValueError: # ValueError: File is not a recognized excel file
        df = pd.read_table(match_table)
    dropna = df[['Barcode','Slide','Lane']].dropna(axis=0, how='any').index
    df_dpna = df.loc[dropna] # remove the explain message
    # df_dpna = df.dropna(axis=0, how='all').dropna(axis=1, how='all')
    # deal library structure
    ids = df_dpna[(df_dpna['Library_structure']=='/') & (df_dpna['Type']=='PE100')].index
    df_dpna.loc[ids, 'Library_structure'] = 'R100R100I10I10'
    df_dpna.loc[df_dpna[df_dpna['UMI_type']=='/'].index, 'umi_option'] = 'false'
    df_dpna.loc[df_dpna[df_dpna['UMI_type']!='/'].index, 'umi_option'] = 'true'
    # deal slide and lane
    df_sl = df_dpna[['Slide','Lane']].drop_duplicates()
    for x in df_sl.index:
        slide = df_sl.loc[x]['Slide']
        lane = df_sl.loc[x]['Lane']
        tmp = df_dpna[(df_dpna['Slide']==slide) & (df_dpna['Lane']==lane)]
        library = tmp[['Library','Barcode_1','Barcode_2','Library_structure']]
        library.loc[library.index, 'Panel'] = 'N'
        # library.rename(columns={'Libaray':'Sample'}, inplace=True)
        library.to_csv(f'{slide}_{lane}_split_fq_library.txt', index=False, sep='\t')
        sampletable = tmp[['Library','Panel','Pipeline','umi_option']]
        sampletable.to_csv(f'{slide}_{lane}_sample_message.txt', index=False, sep='\t')

def RunningProcess():
    """用于调度各个分析流程模块"""
    pass

def CheckProcess():
    """用于检查各个分析流程模块的运行情况"""
    pass


if __name__ == '__main__':
    cfgdir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../configure/')
    parser = ArgumentParser('all process pipeline')
    parser.add_argument('-match', '--match_table', required=True,
                        help='matching table file; excel file *.xlsx format')
    parser.add_argument('-process', '--analysis_process', choices=['all', 'isolate', 'trim', 'mapping', 'umi', 'qc', 'caller', 'fusion', 'cnv', 'msi'], default='all',
                        help='the processing of analysis')
    parser.add_argument('-cfg', '--configure_file', default=None,
                        help='running configure file; default=latest version')
    parser.add_argument('-v', '--version', default=None,
                        help='print the lateset version number; and other ')
    目前需要的几个参数：
        读取生信对接表match_table (match_table)
        分析功能选择all, isolate, trim, mapping, umi, QC, caller, fusion, cnv, msi等，选择对应选项后执行相应的分析
        读取配置文件的txt格式（存在一个默认的路径，默认选择最新版本，列举出所有版本信息）
        能够输出各个版本号及版本的信息内容
            加入单个样本从rawdata开始的分析过程的功能
        选择continue分析还是从头分析
        选择local模式还是sge模式
    args = parser.parse_args()
    if args.version is not None:
        print (darui_pipeline_version_description.txt)
        exit()
    elif args.running_configure is not None:
        running(args.running_configure)
    # running(args.fq1, args.fq2, args.bed, args.cfg_json, args.output_dir, args.output_prefix)
