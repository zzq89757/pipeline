#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-20
# Version: 0.2.0.0
# Description: pipeline for all of the analysis process
#########################################################################
import os
import sys
import json
import pandas as pd #reading excel file need xlrd==1.2.0; pip install xlrd==1.2.0
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac
# 能够读取cfg文件，为了适合结果找出bed文件的位置
def LibraryResult(match_table):
    """用于处理生信对接表的模块:
    1.单barcode非umi，link原始fastq文件到Rawdata下面
    2.单barcode umi，cutadapt按照Library_structure拆分出M序列和R序列到Rawdata下面
    3.混barcode非umi，构建MGI_demultiplex分析library表格后对undecode的fastq文件拆分得到R1和R2文件，再将原始机器拆分的read fastq文件与拆分后的R1和R2文件cat到一起；
    4.混barcode umi，构建MGI_demultiplex分析library表格后对undecode的fastq文件拆分得到R1和R2文件和M1和M2文件，再将原始机器拆分的read fastq文件利用cutadapt按照Library_structure拆分出M序列和R序列，之后read fastq与拆分后的R1和R2文件cat到一起，cutadapt拆分出的umi fastq与M1和M2文件cat到一起；
    """
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
    for x in df_dpna.index:
        lib = df_dpna.loc[x]['Library']
        b1 = df_dpna.loc[x]['Barcode_1']
        b2 = df_dpna.loc[x]['Barcode_2']
        umi_option = df_dpna.loc[x]['umi_option']
        if ',' in b1 and ',' in b2:
            pass # barcode combine
        elif ',' not in b1 and ',' not in b2:
            if umi_option is 'false':
                # single barcode non umi

        else:
            with open(f'{lib}_waring_msg.txt','a') as otopen:
                otopen.write('barcode type unknown, combined barcode or one barcode only?\n')

    # df_sl = df_dpna[['Slide','Lane']].drop_duplicates()
    # for x in df_sl.index:
    #     slide = df_sl.loc[x]['Slide']
    #     lane = df_sl.loc[x]['Lane']
    #     tmp = df_dpna[(df_dpna['Slide']==slide) & (df_dpna['Lane']==lane)]
    #     library = tmp[['Library','Barcode_1','Barcode_2','Library_structure']]
    #     library.loc[library.index, 'Panel'] = 'N'
    #     # library.rename(columns={'Libaray':'Sample'}, inplace=True)
    #     library.to_csv(f'{slide}_{lane}_split_fq_library.txt', index=False, sep='\t')
    #     sampletable = tmp[['Library','Panel','Pipeline','umi_option']]
    #     sampletable.to_csv(f'{slide}_{lane}_sample_message.txt', index=False, sep='\t')
