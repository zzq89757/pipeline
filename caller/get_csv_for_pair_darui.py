#coding:utf8
import os
import sys

global header
header = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">
##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Germline event">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP or INDEL">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Average base quality for reads supporting alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="RMS mapping quality">
##FORMAT=<ID=SP,Number=1,Type=String,Description="Start point number of reads supporting mutant allele">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of ref-forward, ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=SQ,Number=1,Type=Integer,Description="Snp quality">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=M,length=16571>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
"""

def PassToCsvStep1(pass_vcf, sampleclass, output_dir, output_prefix):
    snv_list=[header]
    indel_list=[header]
    snv_list.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{output_prefix}\n")
    indel_list.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{output_prefix}\n")
    for line in open(pass_vcf):
        if line.startswith('#'):
            continue
        else:
            col=line.strip().split('\t')
            chrm=col[0].replace('chr','')
            pos=col[1]
            ids=col[2]
            ref=col[3]
            alt=col[4]
            qual=col[5]
            filters=col[6]
            formats=col[7].split(';')
            try:
                testaf=float(formats[0].replace('AF=',''))
            except ValueError:
                continue
            ad=formats[2].replace('DP4=','').split(',')[2]+','+formats[2].replace('DP4=','').split(',')[3]
            bq='30.0'
            dp=formats[1].replace('DP=','')
            af=formats[0].replace('AF=','')
            if sampleclass == 'germline' and float(af) > 0.5:
                gt='1/1'
            elif sampleclass == 'germline' and float(af) <= 0.5:
                gt='0/1'
            elif sampleclass == 'somatic':
                gt='0/0'
            gq='99'
            mq='60'
            sp='0'
            dp4=formats[2].replace('DP4=','')
            sq='0'
            opt_formats_str=":".join([ad,bq,dp,af,gt,gq,mq,sp,dp4,sq])
            opt_str="\t".join([chrm,pos,ids,ref,alt,qual,filters,'info','AD:BQ:DP:AF:GT:GQ:MQ:SP:DP4:SQ',f'{opt_formats_str}'])+'\n'
            if sampleclass == 'germline':
                sclass='GERMLINE'
                otfilename=f'csv_{output_prefix}.tmp.germline.vcf'
            elif sampleclass == 'somatic':
                sclass='SOMATIC'
                otfilename=f'csv_{output_prefix}.all_snv.tmp.vcf'
            if 'INDEL' in line:
                opt_str=opt_str.replace('info',f'{sclass};VT=INDEL')
                indel_list.append(opt_str)
            else:
                opt_str=opt_str.replace('info',f'{sclass};VT=SNV')
                snv_list.append(opt_str)

    snp_filename=otfilename.replace('tmp','snp')
    indel_filename=otfilename.replace('tmp','indel')
    with open(f'{output_dir}/{snp_filename}','w') as otopen:
        otopen.writelines(snv_list)
    with open(f'{output_dir}/{indel_filename}','w') as otopen:
        otopen.writelines(indel_list)

def PassToCsvStep2(all_snv_var, germline_var, output_dir, output_preifx, types):
    germline_var_dict={}
    optlist=[]
    for line in os.popen(f'cat {germline_var} |grep -v "#"'):
        col=line.split('\t')
        var_site=':'.join([col[0],col[1],col[3],col[4]])
        germline_var_dict[var_site]=[]
    for line in open(all_snv_var):
        if line.startswith('#'):
            optlist.append(line)
        else:
            col=line.split('\t')
            var_site=':'.join([col[0],col[1],col[3],col[4]])
            if germline_var_dict.get(var_site) == None:
                optlist.append(line)
    with open(f'{output_dir}/csv_{output_preifx}.merge_snv.{types}.vcf', 'w') as otopen:
        otopen.writelines(optlist)

if len(sys.argv) != 6:
    print("python this_script.py somatic.PASS.vcf germline.PASS.vcf output_dir somatic_prefix germline_prefix")
    exit()
somatic_pass_vcf = sys.argv[1]
germline_pass_vcf = sys.argv[2]
output_dir = sys.argv[3]
somatic_prefix = sys.argv[4]
germline_prefix = sys.argv[5]
PassToCsvStep1(somatic_pass_vcf, 'somatic', output_dir, somatic_prefix)
PassToCsvStep1(germline_pass_vcf, 'germline', output_dir, germline_prefix)

all_snv_snp = f"{output_dir}/csv_{somatic_prefix}.all_snv.snp.vcf"
all_snv_indel = f"{output_dir}/csv_{somatic_prefix}.all_snv.indel.vcf"
germline_snp = f"{output_dir}/csv_{germline_prefix}.snp.germline.vcf"
germline_indel = f"{output_dir}/csv_{germline_prefix}.indel.germline.vcf"
PassToCsvStep2(all_snv_snp, germline_snp, output_dir, somatic_prefix, 'snp')
PassToCsvStep2(all_snv_indel, germline_indel, output_dir, somatic_prefix, 'indel')

os.system(f"rm {all_snv_snp} {all_snv_indel}")
# python get_csv_darui.py ../L-KY-022P/snv_indel/L-KY-022P.PASS.vcf ../L-KY-022W/snv_indel/L-KY-022W.PASS.vcf . L-KY-022P L-KY-022W
