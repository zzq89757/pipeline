#!/bin/bash
name=$1
dir=$2
bed=$3
r1=$4
r2=$5
if  [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]
then
	echo "sh this_script.sh name output_dir bed read1 read2"
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
# remove duplicate
python ${script_dir}/mapping/run_rmdup.py --input_bam ${sort_bam} --output_dir ${dir}/mapping/ --output_prefix ${name}
rmdup_bam=${dir}/mapping/${name}.rmdup.bam
# bam QC
python ${script_dir}/qc/bam_qc.py --sort_bam ${sort_bam} --bed ${bed} --output_dir ${dir}/QC/ --output_prefix ${name} --sample_type gDNA --removedup_bam ${rmdup_bam} --umi false
meandp=`cat ${dir}/QC/${name}.CollectHsMetrics_raw.xls |head -n 8 |tail -n 1 |awk -F "\t" '{print $10}'`

downsample_number=2000
if [[ ${meandp} > ${downsample_number} ]]
then
    python ${script_dir}/mapping/bam_downsample.py --sort_bam ${sort_bam} --bed ${bed} --output_dir ${dir}/mapping/ --output_prefix ${name} --downsample_num ${downsample_number}
    downsample_bam=${dir}/mapping/${name}.sort.downsample_${downsample_number}.bam
    # downsample remove duplicate
    python ${script_dir}/mapping/run_rmdup.py --input_bam ${downsample_bam} --output_dir ${dir}/mapping/ --output_prefix ${name}.downsample_${downsample_number}
    rmdup_downsample_bam=${dir}/mapping/${name}.downsample_${downsample_number}.rmdup.bam
    # downsample bam QC
    python ${script_dir}/qc/bam_qc.py --sort_bam ${downsample_bam} --bed ${bed} --output_dir ${dir}/QC/ --output_prefix ${name} --sample_type gDNA --removedup_bam ${rmdup_downsample_bam} --umi false &
    # downsample indelrealigner
    python ${script_dir}/mapping/run_indelrealigner.py --input_bam ${rmdup_downsample_bam} --bed ${bed} --output_dir ${dir}/mapping/ --output_prefix ${name}.downsample_${downsample_number}
else
    # indelrealigner
    python ${script_dir}/mapping/run_indelrealigner.py --input_bam ${rmdup_bam} --bed ${bed} --output_dir ${dir}/mapping/ --output_prefix ${name}
fi
