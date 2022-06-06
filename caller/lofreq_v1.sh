path=$1
sample=$2
input_bam=$3
bed=$4

localdir=$(cd "$(dirname "$0")"; pwd)
ref=${localdir}/../../../database/gatk-bundle-hg19/ucsc.hg19.fasta
# bed=/data/home/caiqinglong/TCGA/data_preparing/panels/NanOnco_Plus_Panel_v2.0_Covered_hg19.bed
# bed=/data/home/caiqinglong/Database/panel/LungCancer_Panel_v1.0_Probe_Covered_hg19.bed
annovar_db=${localdir}/../../../database/annovar_humandb_hg19
background=${localdir}/../../../database/annovar_background
diff_T_N_v1=${localdir}/bin/diff_T_N_v1.R
make_background=${localdir}/bin/make_background.R

mkdir -p ${path}/raw_vcf
mkdir -p ${path}/indelqualbam
mkdir -p ${path}/filter
mkdir -p ${path}/annotation
mkdir -p ${path}/extractFields
mkdir -p ${path}/result

#${path}/bam/${sample}.sort.bam
echo "start indelqual for ${sample} at" `date` >> ${path}/time_cost.txt
${localdir}/../../../software/lofreq-2.1.5/bin/lofreq indelqual \
	--dindel \
	-f ${ref} \
	${input_bam} \
	-o ${path}/indelqualbam/${sample}_indelqual.bam
echo "end indelqual for ${sample} at" `date` >> ${path}/time_cost.txt
samtools index ${path}/indelqualbam/${sample}_indelqual.bam

echo "start call-parallel for ${sample} at" `date` >> ${path}/time_cost.txt
${localdir}/../../../software/lofreq-2.1.5/bin/lofreq call-parallel \
       --pp-threads 8 \
       ${path}/indelqualbam/${sample}_indelqual.bam \
       -f ${ref} \
       -l ${bed} \
       --call-indels \
       -Q 25 \
       -q 25 \
       -a 1 \
       -b 1 \
       -A \
       -B \
       --no-default-filter \
       -o ${path}/raw_vcf/${sample}.raw.vcf.gz
echo "end call-parallel for ${sample} at" `date` >> ${path}/time_cost.txt

java -jar -Xmx64G ${localdir}/../../../software/gatk-package-4.2.0.0-local.jar VariantFiltration \
        --variant ${path}/raw_vcf/${sample}.raw.vcf.gz \
        --output ${path}/filter/${sample}.filtered.vcf.gz \
	-filter "(AF>=0.001&&DP>=500&&QUAL>=50&&SB<=60&&DP4[2]>=3&&DP4[3]>=3)&&((vc.hasAttribute('INDEL')&&HRUN<=8))" \
        --filter-name PASS \
        -filter "AF<0.001" \
        --filter-name low_AF \
        -filter "DP<500" \
        --filter-name low_DP \
        -filter "QUAL<50" \
        --filter-name low_QUAL \
        -filter "DP4[2]<3||DP4[3]<3" \
        --filter-name low_AD \
        -filter "SB>60" \
        --filter-name strand_bias \
        -filter "(vc.hasAttribute('INDEL')&&HRUN>8)" \
        --filter-name indel_err

${localdir}/../../../software/bcftools-1.13/bin/bcftools filter \
        -i 'filter="PASS"' ${path}/filter/${sample}.filtered.vcf.gz \
        -o ${path}/filter/${sample}.PASS.vcf
:<<!
#merge variants
vcf-genotype-annotator \
        ${path}/filter/${sample}.PASS.vcf \
        ${sample} \
        0/1 \
        -o ${path}/filter/${sample}.PASS.GT.vcf

whatshap --debug phase \
        -o ${path}/filter/${sample}.PASS.phased.vcf \
        --reference=/data/database/gatk-bundle-hg19/ucsc.hg19.fasta \
        ${path}/filter/${sample}.PASS.GT.vcf \
        ${path}/bam/${sample}.sort.bam \
        --indels \
        --ignore-read-groups
!
#variants annotation
java -jar -Xmx64G ${localdir}/../../../software/snpEff-5.0e/snpEff.jar ann \
        hg19 \
	-onlyTr ${annovar_db}/trans.txt \
        ${path}/filter/${sample}.PASS.vcf \
        -ud 3 \
        -s ${path}/annotation/${sample}.summary.html \
        > ${path}/annotation/${sample}.snpeff.vcf

java -jar -Xmx64G ${localdir}/../../../software/snpEff-5.0e/SnpSift.jar filter \
	"(ANN[0].EFFECT != 'intron_variant') & (ANN[0].EFFECT != 'intergenic_region') & (ANN[0].EFFECT != 'intragenic_variant') & (ANN[0].EFFECT != '3_prime_UTR_variant') & (ANN[0].EFFECT != '5_prime_UTR_variant') & (ANN[0].EFFECT != 'splice_region_variant')" \
        ${path}/annotation/${sample}.snpeff.vcf \
        > ${path}/annotation/${sample}.rm_intron.vcf
#如果${sample}.rm_intron.vcf文件为空则会在table_annovar这一步报错
perl ${localdir}/../../../software/annovar-2020.06.08/table_annovar.pl \
	${path}/annotation/${sample}.rm_intron.vcf \
	${annovar_db} \
	-buildver hg19 \
	-out ${path}/annotation/${sample} \
	-remove \
	-protocol refGene,cytoBand,cosmic92_coding,clinvar_20210123,avsnp150,exac03,ALL.sites.2015_08,EAS.sites.2015_08 \
	-operation g,r,f,f,f,f,f,f \
	--vcfinput
# /data/software/src/snpEff
java -jar -Xmx64G ${localdir}/../../../software/snpEff-5.0e/SnpSift.jar extractFields \
        -s "," \
        -e "." \
        ${path}/annotation/${sample}.hg19_multianno.vcf \
        CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" AF "ANN[*].RANK" DP SB QUAL DP4[0] DP4[1] DP4[2] DP4[3] cosmic92_coding CLNDN CLNSIG avsnp150 ExAC_ALL ExAC_EAS ALL.sites.2015_08 EAS.sites.2015_08 \
        > ${path}/extractFields/${sample}.extractFields.xls
/usr/bin/sed -i 's+\\x3d+=+g' ${path}/extractFields/${sample}.extractFields.xls
/usr/bin/sed -i 's+\\x3b+:+g' ${path}/extractFields/${sample}.extractFields.xls
/usr/bin/awk '{print "'${sample}'\t"$0}' ${path}/extractFields/${sample}.extractFields.xls| /usr/bin/sed -e 1d > ${path}/result/${sample}.query.xls

# /data/home/caiqinglong/scripts/pipeline/diff_T_N_v1.R
Rscript ${diff_T_N_v1} ${sample} ${path} 0 ${background}
#update background
cp ${path}/result/${sample}.query.xls ${background}
echo -e "sample\tCHROM\tPOS\tREF\tALT\tGENE\tEFFECT\ttrans\tHGVS_C\tHGVS_P\tAF\tRANK\tDP\tSB\tQUAL\tref_F\tref_R\talt_F\talt_R\tcosmic92_coding\tCLNDN\tCLNSIG\tavsnp150\tEXAC_ALL\tEXAC_EAS\t1000G_ALL\t1000G_EAS" > ${background}/background_raw.xls
cat ${background}/*.query.xls >> ${background}/background_raw.xls
# /data/home/caiqinglong/scripts/pipeline/make_background.R
Rscript ${make_background} ${background}/background_raw.xls ${background}/background.xls
