#!/bin/env bash
#author:LBB
#Date:2021-09-18
#Version:1.0
# Description: pipline for UMI


####soft###########
localdir=$(cd "$(dirname "$0")"; pwd)
picard=${localdir}/../../software/picard_2.25.0.jar
fgbio=${localdir}/../../software/fgbio.jar
bwa=${localdir}/../../software/bwa-0.7.17/bwa
samtools=${localdir}/../../software/samtools-1.13/samtools

helpFunction()
{
  echo ""
  echo "Usage: $0  -s strategy -t thread  -r reference  -y umi_type -d fastq_pwd -p bam_prefix -U umifile  -o outdir"
  echo -e "\t-s strategy=Adjacency(singleUMI) or Paired(duplexUMI)"
  echo -e "\t-t thread"
  echo -e "\t-r reference the input original fasta file "
  echo -e "\t-y type single_ramdom_UMI | duplex_ramdom_UMI | duplex_fix_UMI "
  echo -e "\t-d bam_prefix"
  echo -e "\t-p bam_prefix"
  echo -e "\t-U  UMI file "
  echo -e "\t-o outdir"
  exit -1 # exit script after pringting help information
}

######### opts##################
while getopts 'h:U:s:t:r:y:p:d:o:' OPT
do 
	case "$OPT" in
    	h) helpFunction ;; # print helpFuntion in case parameter is non-existent
	s) strategy="$OPTARG" ;;   ##singleUMI:adjacency,duplexUMI:paired
	t) thread="$OPTARG" ;;
	r) reference="$OPTARG" ;;
	y) umi_type="$OPTARG" ;;
	p) bam_prefix="$OPTARG" ;;
	U) umifile="$OPTARG" ;;
	d) dir="$OPTARG" ;;
	o) out="$OPTARG" ;;
esac
done
echo $strategy
echo $thread
echo $reference
echo $umi_type
echo $dir
echo $bam_prefix
echo $umifile
echo $out

# Print helpFunction if parameters are empty



if  [ -z "$strategy" ] || [ -z "$thread" ] || [ -z "$reference" ] || [ -z "$umi_type" ] || [ -z "$bam_prefix" ]
then
                echo "Some or all of the parameters are empty";
		helpFunction
fi




########umi deal##########
function sigle_ramdomUMI {
	
	#first step:add UMI 
	#second step:add MQ
	#step:group callconsensusreads realign
	 java -jar -Xmx64G $picard FixMateInformation -I $dir/$bam_prefix.sort.withUMI.bam -O $dir/$bam_prefix.MQ.UMI.bam -MC true -SO coordinate
	 java -jar -Xmx64G $fgbio GroupReadsByUmi -i $dir/$bam_prefix.MQ.UMI.bam \
		  -s ${strategy} -e 0 -m 20  -o $dir/$bam_prefix.group.bam
	#  $samtools sort $dir/$bam_prefix.group.bam -o $dir/$bam_prefix.group.sort.bam
	  java -jar -Xmx64G $fgbio CallMolecularConsensusReads -i $dir/$bam_prefix.group.bam \
		-1 45 -2 30 \
                -N 40 -m 30 \
		-M 1 -o $dir/$bam_prefix.consensus.unmapped.bam 
	 
	 java -jar -Xmx64G $picard SamToFastq -I $dir/$bam_prefix.consensus.unmapped.bam  -F /dev/stdout -INTER true | \
          $bwa mem -p -t ${thread} $reference /dev/stdin  | \
         java -jar -Xmx64G $picard MergeBamAlignment \
                -UNMAPPED $dir/$bam_prefix.consensus.unmapped.bam -O $out/$bam_prefix.consensus.mapped.bam -ALIGNED /dev/stdin -R $reference \
                -SO coordinate  -MAX_GAPS -1 
	${samtools} index $out/$bam_prefix.consensus.mapped.bam $out/$bam_prefix.consensus.mapped.bai
}	    		

function duplex_ramdomUMI {
	
	java -jar -Xmx64G $picard FixMateInformation -I $dir/$bam_prefix.sort.withUMI.bam -O $dir/$bam_prefix.MQ.UMI.bam -MC true -SO coordinate	
	java -jar -Xmx64G $fgbio GroupReadsByUmi -i $dir/$bam_prefix.MQ.UMI.bam \
                 -s $strategy -e 0 -m 20 -o $dir/$bam_prefix.group.bam
	# $samtools sort $dir/$bam_prefix.group.bam -o $dir/$bam_prefix.group.sort.bam
	 java -jar -Xmx64G $fgbio CallDuplexConsensusReads -i $dir/$bam_prefix.group.bam  -o $dir/$bam_prefix.consensus.unmapped.bam \
                -1 45 -2 30 \
                -m 30
	 java -jar -Xmx64G $picard SamToFastq -I $dir/$bam_prefix.consensus.unmapped.bam  -F /dev/stdout -INTER true | \
         $bwa mem -p -t ${thread} $reference /dev/stdin  | \
          java -jar -Xmx64G $picard MergeBamAlignment \
                -UNMAPPED $dir/$bam_prefix.consensus.unmapped.bam -O $out/$bam_prefix.consensus.mapped.bam -ALIGNED /dev/stdin -R $reference \
                -SO coordinate  -MAX_GAPS -1 
        ${samtools} index $out/$bam_prefix.consensus.mapped.bam $out/$bam_prefix.consensus.mapped.bai
}

function duplex_fixUMI {
	
	java -jar -Xmx64G $picard FixMateInformation -I $dir/$bam_prefix.sort.withUMI.bam -O $dir/$bam_prefix.MQ.UMI.bam -MC true -SO coordinate
	java -jar -Xmx64G $fgbio CorrectUmis -i $dir/$bam_prefix.MQ.UMI.bam \
                 -o $dir/$bam_prefix.mapped.fixedumi.bam \
                 -m 3 -d 1 \
                 -M $dir/$bam_prefix.metrics.txt \
                 -r $dir/$bam_prefix.rejected.bam \
                 -t RX -U $umifile
	java -jar -Xmx64G $fgbio GroupReadsByUmi -i $dir/$bam_prefix.mapped.fixedumi.bam \
                -s $strategy -e 0 -m 20  -o $dir/$bam_prefix.group.bam 
	# $samtools sort $dir/$bam_prefix.group.bam -o $dir/$bam_prefix.group.sort.bam
        java -jar -Xmx64G $fgbio CallDuplexConsensusReads -i $dir/$bam_prefix.group.bam -o $dir/$bam_prefix.consensus.unmapped.bam \
                -1 45 -2 30 \
                -m 30
        java -jar -Xmx64G $picard SamToFastq -I $dir/$bam_prefix.consensus.unmapped.bam  -F /dev/stdout -INTER true | \
         $bwa mem -p -t ${thread} $reference /dev/stdin  | \
          java -jar -Xmx64G $picard MergeBamAlignment \
                -UNMAPPED $dir/$bam_prefix.consensus.unmapped.bam -O $out/$bam_prefix.consensus.mapped.bam -ALIGNED /dev/stdin -R $reference \
               -SO coordinate  -MAX_GAPS -1  
        ${samtools} index $out/$bam_prefix.consensus.mapped.bam $out/$bam_prefix.consensus.mapped.bai
}


echo -e "---------------[ `date`: BEGIN ]---------------"  
if [ $umi_type == "single_ramdom_UMI" ];then
        sigle_ramdomUMI
        echo -e "---------------[ `date`: END ]---------------" 

elif [ $umi_type == "duplex_ramdom_UMI" ];then
        duplex_ramdomUMI
        echo -e "---------------[ `date`: END ]---------------" 

elif  [ $umi_type == "duplex_fix_UMI" ];then
        duplex_fixUMI
        echo -e "---------------[ `date`: END ]---------------" 
fi

