#coding:utf8
#########################################################################
# Author:  
# Created Time: 2021-12-15
# Version: 0.2.0.0
# Description: pipline for UMI analysis
#########################################################################
import os
import sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../core'))
import bac

def RunningProcessUpdate(cmd, output_dir, output_prefix, other_string, dirprefix, keyword, existfile):
    # check command running or not, not add step tag
    filesize = bac.GetSize(existfile)
    instepkeyward = os.popen(f'grep {keyword} {dirprefix}.umistep.txt').read()
    if instepkeyward != '' and filesize != 0:
        pass
    # if file exists and keyword in umistep.txt file, command will not run, otherwise command will run
    else:
        bac.RunningProcess(cmd, output_dir, output_prefix, other_string) # this command run and add step tag
        outputfie = cmd.strip().split(' ')[-1] # 切记：所有要生产的文件都放在命令行的最后
        if bac.GetSize(outputfie) != 0:
            os.system(f'echo "{keyword} command is accomplished! at the time: `date`" >> {dirprefix}.umistep.txt')

def CheckFiles(sort_bam, umiseq_fq1, umiseq_fq2, UMI_type):
    if not bac.CheckFile(sort_bam) or not bac.CheckFile(umiseq_fq1):
        print (f"{sort_bam} or {umiseq_fq1} not exist")
        exit()
    if 'duplex' in UMI_type:
        if not bac.CheckFile(umiseq_fq2):
            print (f"UMI_type is {UMI_type} but {umiseq_fq2} not exist")
            exit()
    if bac.DetermainR1R2(umiseq_fq1, umiseq_fq2):
        exit()

def Prepare(sort_bam, umiseq_fq1, umiseq_fq2, cfg, output_dir, output_prefix, UMI_type):
    CheckFiles(sort_bam, umiseq_fq1, umiseq_fq2, UMI_type)
    dirprefix = os.path.join(output_dir, output_prefix)
    if not bac.CheckFile(f'{dirprefix}.umistep.txt'):
        os.system(f'touch {dirprefix}.umistep.txt')
    config = bac.ResolveConfig(cfg)
    python = config['language']['python']
    perl = config['language']['perl']
    java = config['language']['java']
    java_mem = config['language']['java_mem']
    software = config['software']
    bwa = software['bwa']
    samtools = software['samtools']
    picard = software['picard']
    fgbio = software['fgbio']
    hg19_fa = config['database']['hg19_fa']
    bwa_threads = config['parameter']['bwa_threads']
    samtools_threads = config['parameter']['samtools_threads']
    # script_dir = os.path.dirname(os.path.abspath(__file__))
    addRX = config['scripts']['addRX']
    # bam2fq = config['scripts']['bam2fq']
    duplex_fix = config['parameter']['duplex_fixUMI']
    # addRX=os.path.join(script_dir, 'bam_add_RX.pl')
    # bam2fq=os.path.join(script_dir, 'unmapped_bam2fq.py')
    # duplex_fix=os.path.join(script_dir, 'duplex_fixUMI')
    return perl, python, java, java_mem, bwa, samtools, picard, fgbio, hg19_fa, bwa_threads, samtools_threads, dirprefix, addRX, duplex_fix

def CommandList(sort_bam, umiseq_fq1, umiseq_fq2, cfg, output_dir, output_prefix, UMI_type):
    # parpare each software, database and scripts
    (perl, python, java, java_mem, bwa, samtools, picard, fgbio, hg19_fa, bwa_threads, samtools_threads, dirprefix, addRX, duplex_fix) = Prepare(sort_bam, umiseq_fq1, umiseq_fq2, cfg, output_dir, output_prefix, UMI_type)
    # add umi sequence
    cmd_addrx = " ".join([
        perl, addRX, '--bam', sort_bam, 
        '--UMI1', umiseq_fq1, '--UMI2', umiseq_fq2, 
        '--out', f'{dirprefix}.sort.withUMI.bam'
    ])
    # add mate MQ; MQ:i:num
    cmd_mq = " ".join([
        java, f'-jar -Xmx{java_mem}G', 
        picard, 'FixMateInformation', 
        '--INPUT', f'{dirprefix}.sort.withUMI.bam', 
        '--TMP_DIR', f'{dirprefix}_FixMateInformation_tmp', 
        '--ADD_MATE_CIGAR true --SORT_ORDER coordinate',
        '--OUTPUT', f'{dirprefix}.MQ.UMI.bam' 
    ])
    # duplex double umi running fgbio CorrectUmis module
    cmd_correctumis = " ".join([
        java, f'-jar -Xmx{java_mem}G', 
        fgbio, f'--tmp-dir {dirprefix}_CorrectUmis_tmp', 
        'CorrectUmis',
        '--input', f'{dirprefix}.MQ.UMI.bam',
        '--max-mismatches 3 --min-distance 1',
        f'--metrics {dirprefix}.metrics.txt',
        # '-r %s.rejected.bam', #Reject BAM file to save unassigned reads. [Optional].
        '--umi-tag RX --umi-files', duplex_fix,
        '--output', f'{dirprefix}.fixedumi.bam'
    ])
    # fgbio GroupReadsByUmi module
    cmd_groupby = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        fgbio, f'--tmp-dir {dirprefix}_GroupReadsByUmi_tmp', 
        'GroupReadsByUmi',
        '--input', 'input_bam',
        '--strategy', 'strategy_used', '--edits 1 --min-map-q 20',
        '--output', f'{dirprefix}.group.bam'
    ])
    # single UMI fgbio CallMolecularConsensusReads module
    cmd_single_consensus = " ".join([
        java, f'-jar -Xmx{java_mem}G',
        fgbio, f'--tmp-dir {dirprefix}_CallMolecularConsensusReads_tmp',
        'CallMolecularConsensusReads',
        '--input', f'{dirprefix}.group.bam',
        '--error-rate-pre-umi 45 --error-rate-post-umi 30',
        '--min-input-base-quality 20 --min-consensus-base-quality 40',
        '--min-reads 1 1 0',
        '--output', f'{dirprefix}.consensus.unmapped.bam'
    ])
    # paired UMI fgbio CallDuplexConsensusReads module
    cmd_duplex_consensus = cmd_single_consensus.replace(' --min-consensus-base-quality 40 --min-reads 1 1 0',' --min-reads 1 1 0 --threads 24').replace('CallMolecularConsensusReads','CallDuplexConsensusReads')
    cmd_mappend_consensus = " ".join([ #自行转换umi后生成consensus.mapped.bam的碱基Q值会影响picard CollectHsMetrics --PER_TARGET_COVERAGE计算每个target的min_coverage值，所以fgbio建议不做碱基Q值的转换，碱基质量值越高说明Consensus后得到的结果越准确
        # python, bam2fq, samtools, '%s.consensus.unmapped.bam', '|',
        samtools, 'fastq -N --threads %s'%samtools_threads, f'{dirprefix}.consensus.unmapped.bam', '|', 
        bwa, 'mem -p -t %s'%bwa_threads, '-Y -M', '-R',
        "'@RG\\tID:%s\\tSM:%s\\tPL:MGI\\tLB:%s\\tPE:100'"%(output_prefix, output_prefix, output_prefix), hg19_fa, '/dev/stdin', '|',
        samtools, 'sort --threads %s'%samtools_threads, '-o', f'{dirprefix}.consensus.mapped.bam'
    ])
    # samtools index mapped.bam
    cmd_index = " ".join([samtools, 'index', f'{dirprefix}.consensus.mapped.bam', f'{dirprefix}.consensus.mapped.bai'])

    return dirprefix, cmd_addrx, cmd_mq, cmd_correctumis, cmd_groupby, cmd_single_consensus, cmd_duplex_consensus, cmd_mappend_consensus, cmd_index

def Running(sort_bam, umiseq_fq1, umiseq_fq2, cfg, output_dir, output_prefix, UMI_type):
    (dirprefix, cmd_addrx, cmd_mq, cmd_correctumis, cmd_groupby, cmd_single_consensus, cmd_duplex_consensus, cmd_mappend_consensus, cmd_index) = CommandList(sort_bam, umiseq_fq1, umiseq_fq2, cfg, output_dir, output_prefix, UMI_type)
    
    strategy = {'single_random':'Adjacency', 'single_fix':'Adjacency', 'duplex_random':'Paired', 'duplex_fix':'Paired'}.get(UMI_type)

    RunningProcessUpdate(cmd_addrx, output_dir, output_prefix, 'Notice:\tperl bam_add_RX.pl, getting *.sort.withUMI.bam file running messages', dirprefix, 'addRX', f'{dirprefix}.sort.withUMI.bam')
    RunningProcessUpdate(cmd_mq, output_dir, output_prefix, 'Notice:\tpicard FixMateInformation, getting *.MQ.UMI.bam file running messages', dirprefix, 'addMQ', f'{dirprefix}.MQ.UMI.bam')
    
    if 'fix' in UMI_type:
        RunningProcessUpdate(cmd_correctumis, output_dir, output_prefix, 'Notice:\tfgbio CorrectUmis, getting *.fixedumi.bam file running messages', dirprefix, 'CorrectUMI', f'{dirprefix}.fixedumi.bam')
        cmd_groupby = cmd_groupby.replace('input_bam', f'{dirprefix}.fixedumi.bam').replace('strategy_used', strategy)
    elif 'random' in UMI_type:
        cmd_groupby = cmd_groupby.replace('input_bam', f'{dirprefix}.MQ.UMI.bam').replace('strategy_used', strategy)
    
    RunningProcessUpdate(cmd_groupby, output_dir, output_prefix, 'Notice:\tfgbio GroupReadsByUmi, getting *.group.bam file running messages', dirprefix, 'GroupBy', f'{dirprefix}.group.bam')
    
    if 'single' in UMI_type:
        RunningProcessUpdate(cmd_single_consensus, output_dir, output_prefix, 'Notice:\tfgbio CallMolecularConsensusReads, getting *.consensus.unmapped.bam file running messages', dirprefix, 'ConsensusSingleUnmap', f'{dirprefix}.consensus.unmapped.bam')
    elif 'duplex' in UMI_type:
        RunningProcessUpdate(cmd_duplex_consensus, output_dir, output_prefix, 'Notice:\tfgbio CallDuplexConsensusReads, getting *.consensus.unmapped.bam file running messages', dirprefix, 'ConsensusDuplexUnmap', f'{dirprefix}.consensus.unmapped.bam')
    
    RunningProcessUpdate(cmd_mappend_consensus, output_dir, output_prefix, 'Notice:\tbaw mem mapped to genome, getting *.consensus.mapped.bam file running messages', dirprefix, 'ReMap', f'{dirprefix}.consensus.mapped.bam')
    RunningProcessUpdate(cmd_index, output_dir, output_prefix, 'Notice:\tsamtools index *.consensus.unmapped.bam, getting *.consensus.mapped.bai file running messages', dirprefix, 'Index', f'{dirprefix}.consensus.mapped.bai')

if __name__ == '__main__':
    cfgdir = os.path.dirname(os.path.abspath(__file__)) + '/../configure'
    cfgfile = cfgdir +'/b'+ bac.Versions(cfgdir)[0]
    parser = ArgumentParser('Doing UMI analysis')
    parser.add_argument('-b', '--sort_bam', required=True,
                        help='sorted bam file')
    parser.add_argument('-m1', '--umiseq_fq1', required=True,
                        help='umi fastq file read1')
    parser.add_argument('-m2', '--umiseq_fq2', required=None,
                        help='umi fastq file read2; default=None')
    parser.add_argument('-c', '--cfg', type=str, default=cfgfile,
                        help='cfg_version file default=%s'%cfgfile)
    parser.add_argument('-o', '--output_dir', type=str, default='./',
                        help='output directory, default="./"')
    parser.add_argument('-p', '--output_prefix', required=True,
                        help='sample id')
    parser.add_argument('-y', '--UMI_type', choices=['single_random', 'single_fix', 'duplex_random', 'duplex_fix'], default='duplex_random',
                        help='fixed umi type will correct umi sequence, and random can not; default=duplex_random')
    args = parser.parse_args()
    Running(args.sort_bam, args.umiseq_fq1, args.umiseq_fq2, args.cfg, args.output_dir, args.output_prefix, args.UMI_type)

# each UMI_type 'cmd_addrx' and 'cmd_mq' is first and 'cmd_mappend_consensus' is end; if single UMI consensus type is CallMolecularConsensusReads duplex UMI consensus type is CallDuplexConsensusReads; if fix UMI need CorrectUmis and random not need
# UMI_type is single_random, GroupReadsByUmi strategy is Adjacency
    # cmd_groupby (input bamfile is cmd_mqumi output bamfile)
    # cmd_single_consensus
# UMI_type is single_fix, GroupReadsByUmi strategy is Adjacency
    # cmd_correctumis
    # cmd_groupby (input bamfile is cmd_correctumis output bamfile)
    # cmd_single_consensus
# UMI_type is duplex_random, GroupReadsByUmi strategy is Paired
    # cmd_groupby (input bamfile is cmd_mqumi output bamfile)
    # cmd_duplex_consensus
# UMI_type is duplex_fix, GroupReadsByUmi strategy is Paired
    # cmd_correctumis
    # cmd_groupby (input bamfile is cmd_correctumis output bamfile)
    # cmd_duplex_consensus