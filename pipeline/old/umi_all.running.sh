#!/bin/bash
name=$1
dir=$2
bed=$3
r1=$4
r2=$5
m1=$6
m2=$7
if  [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] || [ -z "$6" ] || [ -z "$7" ]
then
	echo "sh this_script.sh name output_dir bed read1 read2 umiseq1 umiseq2"
	exit
fi
script_dir=/data/home/zhounan/mission/KPI_assess_method_detect/deadline_October_pipeline/command/v0.2.0.0
# trim
mkdir log trim
python ${script_dir}/trim/run_trim.py --raw_read1 ${r1} --raw_read2 ${r2} --output_dir ${dir}/trim --output_prefix ${name}
trimr1=${dir}/trim/${name}.trim.R1.fq.gz
trimr2=${dir}/trim/${name}.trim.R2.fq.gz
# fastq QC
mkdir QC
python ${script_dir}/qc/fastq_qc.py --raw_fq1 ${r1} --raw_fq2 ${r2} --trim_fq1 ${trimr1} --trim_fq2 ${trimr2} --output_dir ${dir}/QC/ --output_prefix ${name} --cutadapt_summary ${dir}/log/${name}.stdout_log.txt &
# mapping
mkdir mapping
python ${script_dir}/mapping/run_mapping.py --trimmed_fq1 ${trimr1} --trimmed_fq2 ${trimr2} --output_dir ${dir}/mapping/ --output_prefix ${name}
sort_bam=${dir}/mapping/${name}.sort.bam
# umi
mkdir umi
python ${script_dir}/umi/run_umi.py --sort_bam ${sort_bam} --umiseq_fq1 ${m1} --umiseq_fq2 ${m2} --output_dir ${dir}/umi --output_prefix ${name} --UMI_type duplex_random
group_bam=${dir}/umi/${name}.group.bam
consensus_mapped_bam=${dir}/umi/${name}.consensus.mapped.bam
# indelrealigner
python ${script_dir}/mapping/run_indelrealigner.py --input_bam ${consensus_mapped_bam} --bed ${bed} --output_dir ${dir}/mapping --output_prefix ${name}
indelrealigner_bam=${dir}/mapping/${name}.IndelRealigner.bam
# bam QC
python ${script_dir}/qc/bam_qc.py --bed ${bed} --output_dir ${dir}/QC --output_prefix ${name} --sample_type ctDNA --group_bam ${group_bam} --consensus_mapped_bam ${consensus_mapped_bam} --umi true &
# snv_indel
mkdir snv_indel
python ${script_dir}/caller/run_snvindel.py --input_bam ${indelrealigner_bam} --bed ${bed} --output_dir ${dir}/snv_indel/ --output_prefix ${name}
# annotation
mkdir annotation
python ${script_dir}/annotation/run_annotate.py --passvcf ${dir}/snv_indel/${name}.PASS.vcf --bed ${bed} --output_dir ${dir}/annotation --output_prefix ${name}
