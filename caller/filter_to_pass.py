import os
import sys
# python直接写入bgzip文件并建立*.PASS.phased.vcf.gz.tbi
# bgzip -c -f ${path}/filter/${sample}.PASS.phased.vcf > ${path}/filter/${sample}.PASS.phased.vcf.gz
# tabix -p vcf ${path}/filter/${sample}.PASS.phased.vcf.gz
# 得到*.readytophase.vcf文件后面用whatshap做分析
if len(sys.argv) != 2:
    print ("python %s *.filtered.vcf.gz"%sys.argv[0])
    exit()

vcf_file=sys.argv[1]
dirprefix=vcf_file.replace('.filtered.vcf.gz','')
name=dirprefix.split('/')[-1]
header=os.popen("gunzip -c %s |grep '##' "%vcf_file).readlines()
header.append("##contig=<ID=chr1>\n##contig=<ID=chr2>\n##contig=<ID=chr3>\n##contig=<ID=chr4>\n##contig=<ID=chr5>\n##contig=<ID=chr6>\n##contig=<ID=chr7>\n##contig=<ID=chr8>\n##contig=<ID=chr9>\n##contig=<ID=chr10>\n##contig=<ID=chr11>\n##contig=<ID=chr12>\n##contig=<ID=chr13>\n##contig=<ID=chr14>\n##contig=<ID=chr15>\n##contig=<ID=chr16>\n##contig=<ID=chr17>\n##contig=<ID=chr18>\n##contig=<ID=chr19>\n##contig=<ID=chr20>\n##contig=<ID=chr21>\n##contig=<ID=chr22>\n##contig=<ID=chrX>\n##contig=<ID=chrY>\n##contig=<ID=chr6_cox_hap2>\n##contig=<ID=chr6_mann_hap4>\n##contig=<ID=chr6_mcf_hap5>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n")
title=['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n'%name]
body=[]
for line in os.popen("gunzip -c %s |grep -v '#' "%vcf_file):
    line=line.strip()+'.0'
    # line=line.replace('0;DP=',';DP=')
    col=line.split('\t')
    if col[6] == 'PASS':
        msg=col[7].split(';')
        gt='0/1'
        af=msg[0].replace('AF=','') #AF 0.1500 -> 0.15
        # if af[-1] == '0':
        #     af=af[:-1]
        dp=msg[1].replace('DP=','')
        dp4=msg[2].replace('DP4=','')
        sb=msg[-1].replace('SB=','')
        ad=str( int(dp4.split(',')[0]) + int(dp4.split(',')[1]) )+','+str( int(dp4.split(',')[2]) + int(dp4.split(',')[3]) )
        format_str=":".join([gt,ad,af,dp,dp4,sb])
        body.append("\t".join([line, 'GT:AD:AF:DP:DP4:SB', format_str])+'\n')

opt_list=header+title+body
with open('%s.readytophase.vcf'%dirprefix,'w') as otopen:
    otopen.writelines(opt_list)
