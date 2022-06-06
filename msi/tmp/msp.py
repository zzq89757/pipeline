# 使用msisensor-pro检测MSI
import os
import io
import sys
import math
import random
import pickle
import argparse
import subprocess
from multiprocessing import Pool
from collections import OrderedDict
from copy import deepcopy
from string import Template
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

def get_args():
    '''
    获取命令行参数
    '''
    parser = argparse.ArgumentParser(description='使用msisensor-pro程序检测MSI')
    parser.add_argument(
        '-c', '--coverage',
        default = 20,
        help = "最小覆盖深度"
    )
    parser.add_argument(
        '-t', '--threads',
        default = 8,
        help = "并行线程数"
    )

    subparsers = parser.add_subparsers()

    parser_baseline = subparsers.add_parser('baseline', help = "建立基线")
    parser_baseline.set_defaults(func = baseline)
    parser_baseline.add_argument(
        '-s', '--scanlist',
        help = "msisensor-pro 生成的基因组scan文件"
    )
    parser_baseline.add_argument(
        '-g', '--genome',
        help = "参考基因组fasta文件，如果没有指定 scanlist 将由基因组文件生成"
    )
    parser_baseline.add_argument(
        '-b', '--bedfile',
        help = "MSI位点的BED文件"
    )
    parser_baseline.add_argument(
        '-p', '--plot',
        help = "画基线图"
    )
    parser_baseline.add_argument(
        '-o', '--output',
        help = "基线输出文件"
    )
    parser_baseline.add_argument(
        'bamfiles',
        nargs = '+',
        help = '用于建立基线的BAM文件'
    )

    parser_detect = subparsers.add_parser('detect', help = "检测样本MSI状态")
    parser_detect.set_defaults(func = detect)
    parser_detect.add_argument(
        '-s', '--sample',
        help = "样本名称"
    )
    parser_detect.add_argument(
        '--html',
        help = '生成html格式的报告'
    )
    parser_detect.add_argument(
        'baseline',
        help = "基线文件"
    )
    parser_detect.add_argument(
        'bamfile',
        help = '样本的BAM文件'
    )

    args = parser.parse_args()
    return args

def baseline(args):
    '''建立基线'''
    if not args.scanlist and not args.genome:
        raise ValueError("scanlist 和 genome 必须指定一个")

    if not args.scanlist and args.genome: # 运行msisensor-pro scan程序
        genome = os.path.basename(args.genome)
        if genome[-3:] == '.fa':
            scanlist = genome[:-3] + '.list'
        elif genome[-6:] == '.fasta':
            scanlist = genome[:-6] + '.list'
        else:
            scanlist = genome + '.list'
        cmdlist = ['msisensor-pro', 'scan', '-d', args.genome, '-o', scanlist]
        proc = subprocess.run(cmdlist, capture_output=True)
        if proc.returncode != 0:
            raise

    if args.scanlist:
        scanlist = args.scanlist

    if not os.path.isfile(scanlist):
        raise ValueError(f'文件不存在: {scanlit}')

    if args.bedfile:
        sites = []
        with open(args.bedfile) as f:
            for line in f:
                lines = line.split()
                sites.append(f'{lines[0]}:{lines[1]}')

        sitelist = os.path.basename(args.bedfile)[:-4] + '.list'
        with open(scanlist) as i, open(sitelist, 'w') as o:
            for line in i:
                lines = line.split()
                if lines[0] == 'chromosome': # header
                    o.write(line)
                    continue
                if f'{lines[0]}:{lines[1]}' in sites:
                    o.write(line)
    else:
        sitelist = scanlist

    # msisensor-pro baseline
    samplelist = 'samples.txt'
    with open(samplelist, 'w') as f:
        i = 1
        for bamfile in args.bamfiles:
            f.write(f'case{i}\t{bamfile}\n')
            i += 1

    cmdlist = ['msisensor-pro', 'baseline', '-d', sitelist, '-i', samplelist, '-o', 'baseline',
               '-c', str(args.coverage), '-b', str(args.threads)]
    proc = subprocess.run(cmdlist, capture_output=True)
    if proc.returncode != 0:
        raise

    prefix = '.'.join(sitelist.split('.')[:-1])
    if args.output:
        baseline = args.output
    else:
        baseline = prefix + '.list_baseline'
    subprocess.run(['cp', f'baseline/{prefix}.list_baseline', baseline], capture_output=True)

    bednames = get_table_cols(args.bedfile, slice(0, 4))
    names = bednames
    osites = plot_baseline(baseline, args.bamfiles, names = names, figout = args.plot, threads = args.threads)
    if args.plot:
        plt.savefig(args.plot)
    with open(baseline+'.pkl', 'wb') as f:
        pickle.dump(osites, f)

    return baseline

def get_table_cols(bedfile, col, header = 0):
    values = []
    with open(bedfile) as f:
        for line in f:
            if header != 0:
                header -= 1
                continue
            lines = line.strip()
            lines = lines.split()
            try:
                values.append(lines[col])
            except IndexError:
                values.append(None)
    return values

def run_mpp(baseline, bamfile):
    '''run msisensor-pro pro'''
    prefix = os.path.basename(bamfile)
    prefix = prefix[:-4]
    cmdlist = ['msisensor-pro', 'pro', '-d', baseline, '-t', bamfile, '-o', prefix]
    proc = subprocess.run(cmdlist, capture_output=True)
    if proc.returncode == 0:
        statfile = prefix
        disfile = prefix + '_dis'
        allfile = prefix + '_all'
        unsfile = prefix + '_unstable'
        return(statfile, disfile, allfile, unsfile)

def extract_dis(disfile, freq = True, returndict = False):
    '''重复数分布'''
    ms_sites = []
    with open(disfile) as f:
        for line in f:
            line = line.strip()
            line = line.split()
            if line[0][:3] == 'chr':
                ms_sites.append(['{}-{}-{}'.format(line[0], line[1], line[3])])
            else:
                ms_sites[-1].append([ int(depth) for depth in line[1:] ])

    ods = []
    odd = {}
    for site in ms_sites:
        od = OrderedDict()
        tx = 1
        total = 0
        for num in site[1]:
            if num != 0:
                od[tx] = num
            total += num
            tx += 1
        if freq:
            for tx in od:
                od[tx] = od[tx] / total
        ods.append(od)
        chrom, pos, _ = site[0].split('-')
        site_position = chrom + ':' + pos
        odd[site_position] = od
    if returndict:
        return(odd)
    else:
        return(ods)

def get_msp_dis(baseline, bamfile):
    result = run_mpp(baseline, bamfile)
    dis = extract_dis(result[1])
    return dis


def plot_baseline(baseline, bamfiles, names = None, figout = 'figout', threads = 8):
    '''画基线图'''
    with Pool(threads) as p:
        args = ( (baseline, bamfile) for bamfile in bamfiles )
        samples_dis = p.starmap(get_msp_dis, args)
    site_number = len(samples_dis[0])
    max_number = 20 # 最多画20个点
    if site_number <= max_number:
        max_number = site_number

    sites = {}
    width = 16
    height = 12
    plt.figure(figsize=(width, height))
    rows = int(math.log(max_number*width/height, 2))
    cols = math.ceil(max_number/rows)
    for sample in samples_dis:
        i = 1
        for site in sample:
            plt.subplot(rows, cols, i)
            x = list(site.keys())
            y = list(site.values())
            plt.plot(x, y)
            if names:
                name = names[i-1][3]
                location = names[i-1][0] + ":" + names[i-1][1]
                plt.title(name, loc = 'left')
                if location not in sites:
                    sites[location] = {}
                for pos in site:
                    if pos in sites[location]:
                        sites[location][pos].append(site[pos])
                    else:
                        sites[location][pos] = [site[pos]]
            i += 1
    osites = OrderedDict()
    for location in sites:
        osites[location] = OrderedDict()
        for pos in sorted(sites[location]):
            osites[location][pos] = sum(sites[location][pos]) / len(sites[location][pos])
    #plt.savefig(figout, bbox_inches='tight')
    return osites

def canon_dis(canon_sites, disfiles):
    rt = {}
    for disfile in disfiles:
        dis = extract_dis(disfile, freq=True, returndict=True)
        for site in canon_sites:
            if site not in rt:
                rt[site] = dis[site]
            else:
                for i in dis[site]:
                    if i not in rt[site]:
                        rt[site][i] = dis[site][i]
                    else:
                        rt[site][i] = (rt[site][i] + dis[site][i]) / 2

    return rt

def detect(args):
    '''检测MSI'''
    if not os.path.isfile(args.baseline):
        sys.exit(f'File No Found: {args.baseline}')
    if not os.path.isfile(args.bamfile):
        sys.exit(f'File No Found: {args.bamfile}')
    if args.sample:
        sample = args.sample
    else:
        sample = os.path.basename(args.bamfile).split('.')[0]


    canon_sites = {
        'chr2:47641559': 'BAT-26',
        'chr4:55598211': 'BAT-25',
        'chr2:39564893': 'Mono-27',
        'chr2:95849361': 'NR-24',
        'chr11:102193508': 'NR-27',
        'chr11:125490765': 'NR-22',
        'chr14:23652346': 'NR-21',
        'chr2:51288511': 'D2S123',
        'chr1:120053340': 'BAT-40',
        'chr5:112213678': 'D5S346',
        'chr17:10269157': 'D17S520',
        'chr17:15359877': 'D17S261',
        'chr18:39248243': 'D18S34',
        'chr3:30691871': 'BAT-RII',
    }

    result = run_mpp(args.baseline, args.bamfile)
    allfile = result[2]
    site_details = get_table_dict(allfile, ['chromosome', 'location'])

    # 判断MSI-H
    canon_high = 0
    other_high = 0
    for site in site_details:
        instable = (site_details[site]['pro_p'] > site_details[site]['threshold'])
        site_details[site]['instable'] = instable
        if instable:
            if site in list(canon_sites)[:6]:  # 这个地方注意字典的顺序, bug!
                canon_high += 1
            else:
                other_high += 1
    msi_high = False
    if canon_high > 3:
        msi_high = True
    elif canon_high == 3 and other_high > 0:
        msi_high = True
    elif canon_high == 2 and other_high > 1:
        msi_high = True
    else:
        if other_high / len(site_details) > 0.5:
            msi_high = True

    site_details2 = deepcopy(site_details) # for html report

    # 生成简单的报告
    header = ['sample']
    sites = [sample]
    for site in canon_sites:
        if site in site_details:
            header.append(canon_sites[site])
            site_statu = site_details.pop(site, None)
            sites.append(site_statu['instable'])
    i = 1
    for site in site_details:
        header.append(f"MSI-{i}")
        sites.append(site_details[site]['instable'])
        i += 1
    header.append('MSI-H')
    sites.append(msi_high)

    # 输出
    print("\t".join(header))
    print("\t".join([ str(i) for i in sites]))
    canon_names = []
    canon_site_pos = []
    disfile = result[1]
    tumor_dis = extract_dis(disfile, returndict = True)
    baseline_dis_file = args.baseline + '.pkl'
    if os.path.exists(baseline_dis_file):
        with open(baseline_dis_file, 'rb') as f:
            baseline_dis = pickle.load(f)
    else:
        baseline_dis = {}


    i = 1
    sites_plot = ['BAT-26', 'BAT-25', 'Mono-27', 'NR-21', 'NR-24']
    for site in site_details2:
        if canon_sites[site] in sites_plot:
            ax = plt.subplot(5, 1, i)
            i += 1

            ax.plot(baseline_dis[site].keys(), [ i*100 for i in baseline_dis[site].values()], 'b-', label='normal')
            x = list(tumor_dis[site].keys())
            y = list(tumor_dis[site].values())

            comm = set(baseline_dis[site].keys()) & set(x)
            #t_x = [ baseline_dis[site][i] for i in comm ]
            #t_y = [ tumor_dis[site][i] for i in comm ]
            t_x = []
            t_y = []
            depth = 1000
            random.seed(100)
            for j in baseline_dis[site].keys():
                t_x.extend([j] * int(depth * baseline_dis[site][j]))
            for j in x:
                t_y.extend([j] * int(depth * tumor_dis[site][j]))
            w, p = stats.ranksums(random.sample(t_x, 50), random.sample(t_y, 50))

            y = [0 for i in range(50) if i < min(x)] + y
            x = [min(x)-i for i in range(50) if i < min(x)] + x
            y.extend([0 for i in range(50) if max(x)+i<60])
            x.extend([max(x)+i for i in range(50) if max(x)+i<60])
            y = [ i*100 for i in y ]
            ax.plot(x, y, 'r-', label='tumor')
            if p < 0.0001:
                ax.text(0.6, 0.4, canon_sites[site]+"\n"+"p<0.0001", transform = ax.transAxes)
            else:
                ax.text(0.6, 0.4, canon_sites[site]+"\n"+"p={:.4f}".format(p), transform = ax.transAxes)
            if i == 2:
                ax.legend()
            if i != 6:
                plt.tick_params(axis='x', labelcolor='none', which='both')
            if i == 4:
                plt.ylabel("Support Reads Numner(normalization)")


    plt.xlabel("Microsates Repeat Times")
    plt.savefig(sample)

    if args.html:
        for site in site_details2:
            site_details2[site]['sitepos'] = site
            site_details2[site]['score'] = "%.6f" % (float(site_details2[site]['threshold']) - float(site_details2[site]['pro_p']))
            if site in list(canon_sites):
                site_details2[site]['name'] = canon_sites[site]
            else:
                site_details2[site]['name'] = 'NA'
            if site in list(canon_sites)[:6]:
                canon_names.append(canon_sites[site])
                canon_site_pos.append(site)
                site_details2[site]['figure'] = plotly_plot(tumor_dis[site], baseline_dis.get(site, None))
            else:
                site_details2[site]['figure'] = ''
        metadata = {
            'sample': sample,
            'statu': 'MSI-H' if msi_high else 'MSS',
            'baseline': args.baseline,
            'bamfile': args.bamfile,
            'total_site': len(site_details2),
            'canon_site': ', '.join(canon_names),
            'canon_site_pos': canon_site_pos,
            'mincov': args.coverage
        }
        html_report(site_details2, metadata = metadata, outfile = args.html)

def plotly_plot(site_dis, baseline_dis = None):
    '''使用plotly画图'''
    import plotly.express as px
    x, y, c= [], [], []
    if baseline_dis:
        x.extend(list(baseline_dis.keys()))
        y.extend(list(baseline_dis.values()))
        c.extend(['normal'] * len(baseline_dis))
    x.extend(list(site_dis.keys()))
    y.extend(list(site_dis.values()))
    c.extend(['tumor'] * len(site_dis))
    df = pd.DataFrame(dict(repeat_times = x, frequency = y, color = c))
    fig =px.line(df, x='repeat_times', y='frequency', color = 'color')
    plot = io.StringIO()
    fig.write_html(plot, full_html=False)
    return(plot.getvalue())



def html_report(sites, metadata, outfile = 'out.html'):
    '''HTML检测报告'''
    html_tmpl = '''
<!doctype html>
<html>
<head>
<meta charset='UTF-8'><meta name='viewport' content='width=device-width initial-scale=1'>
<title>$sample</title>
<style>
    body{text-align: center;}
    #myid{margin:0 auto;border:1px solid #000;width:800px;text-align:left;padding:0 10px;}
    table{width: 100%;border:1px solid #000; border-spacing: 0; border-collapse: collapse;}
    th{border-bottom:1px solid #000; text-align: left; min-width:100px; max-width:600px}
    td{border-bottom:1px solid #000; text-align: left; min-width:100px; max-width:600px}
</style>
</head>
<body>
<div id="myid">
<h4 id='sample-id'>样本编号： $sample</h4>
<h4 id='msi-statu'>MSI状态：$statu</h4>
<h4 id='parameters'>检测参数：</h4>
<figure><table id="table1">
<tr><td>输入文件  :</td><td>$baseline</td></tr>
<tr><td>使用基线  :</td><td>$bamfile</td></tr>
<tr><td>总位点数  :</td><td>$total_site</td></tr>
<tr><td>经典位点  :</td><td>$canon_site</td></tr>
<tr><td>最小覆盖度:</td><td>$mincov</td></tr>
</table></figure>
<p>&nbsp;</p>
<h4 id='detail'>位点详情：</h4>
$details
</div>
</body>
</html>
    '''
    site_tmpl = '''
<figure><table>
<thead>
<tr><th>位置</th><th>名称</th><th>重复单元</th><th>参考重复数</th><th>覆盖度</th><th>分值</th></tr></thead>
<tbody><tr><td>$sitepos</td><td>$name</td><td>$repeat_unit_bases</td><td>$repeat_times</td><td>$CovReads</td><td>$score</td></tr><tr><td colspan=6>$figure</td></tr></tbody>
</table></figure>
'''
    details = ''
    tpl = Template(site_tmpl)
    for site in metadata['canon_site_pos']:
        details += tpl.substitute(sites[site])
    other_site = set(sites.keys()) - set(metadata['canon_site_pos'])
    details += '''<figure><table><thead><tr><th>位置</th><th>名称</th><th>重复单元</th><th>参考重复数</th><th>覆盖度</th><th>分值</th></tr></thead><tbody>'''
    for s in other_site:
        site = sites[s]
        details += "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(site['sitepos'], site['name'], site['repeat_unit_bases'], site['repeat_times'], site['CovReads'], site['score'])
    details += '</tbody></table></figure>'
    metadata['details'] = details
    tpl = Template(html_tmpl)
    html = tpl.substitute(metadata)

    with open(outfile, 'w') as f:
        f.write(html)

def get_table_dict(filename, keys, sep = ':'):
    header = []
    values = {}
    try:
        _ = ( e for e in keys )
    except TypeError:
        raise
    with open(filename) as f:
        for line in f:
            lines = line.strip()
            lines = line.split()
            if not header:
                header = lines
                continue
            item = dict(zip(header, lines))
            key = sep.join([ item.pop(i, None) for i in keys])
            values[key] = item

    return values

def main():
    args = get_args()
    args.func(args)
    #print(args)
    #plot_baseline('nad_msi_sites.list_baseline', args.bamfiles, threads = args.threads)
    #print(args)

    #bamfile = '/data/home/chenshulin/project/MSI/NAD/analysis/L01_UDB-1/2.align/L01_UDB-1.sorted.markdup.bam'
    #baseline = '/data/home/chenshulin/project/MSI/test/baseline/nad.list_baseline'
    #print(get_msp_dis(baseline, bamfile))


if __name__ == '__main__':
    main()
