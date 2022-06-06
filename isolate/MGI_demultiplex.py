#!/bin/env python3
import re
import json
import gzip
import os
import fcntl
#import threading
import queue
import multiprocessing
#import mgzip
from xopen import xopen
#import dnaio
import sys
import textwrap
import datetime
import argparse

#version
version = 'Alpha 0.1a'
date = '2021-10-18'

###### Usage information
usage = '''

    Version: %s
    Author: Fang Peng (fun_peng@foxmail.com)
    Update: %s

    Usage: %s --read1 <undecoded read 1 file> --read2 <undecoded read 1 file> --library <library file> --outdir <outDir>
    
''' % (version, date, os.path.basename(sys.argv[0]))


def parse_library_type(library_type, total_cycle = 220):
    pattern = re.compile(r'(\D\d+)')
    structure = pattern.findall(library_type)
    res = {}
    start = 1
    pattern2 = re.compile(r'(\D)(\d+)')
    for i in structure:
        code = pattern2.findall(i)[0]
        j = 1
        if(not (code[0] in res)):
            res[code[0]] = {code[0]+str(j):[start, start + int(code[1]) - 1]}
            start += int(code[1])
        else:
            while code[0]+str(j) in res[code[0]]:
                j +=1
            res[code[0]][code[0]+str(j)] = [start, start + int(code[1]) - 1]
            start += int(code[1])
    if start -1 > total_cycle:
        raise Exception('文库的组成长度为：{}，程序的指定长度为：{}。请调整参数!'.format(start-1, total_cycle))
    return(res)

def replace_char(string, char, index):
    string = list(string)
    string[index] = char
    return(''.join(string))

def reverse_complement(seq):
    seq = seq.upper()
    seq = seq[::-1]
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('C', 'g')
    seq = seq.replace('G', 'c')
    return(seq.upper())

def build_barcode_pool(barcode1, barcode2, mismatch = 1, reverse = True, version = "ECR4", mode = "PE"):  ## version为ECR4:：SE和PE均先测barcode1再测barcode2；ECR3：SE先测barcode2再测barcode1，PE先测barcode1再测barcode2
    barcode1 = barcode1.upper()
    barcode2 = barcode2.upper()
    
    if barcode1 == "" or barcode2 == "":  # 单barcode
        barcode = ""
        if barcode1 == "":
            barcode = barcode2
        else:
            barcode = barcode1
        if reverse == True:
            barcode = reverse_complement(barcode)
            
        barcode_pools = []
        if mismatch >= 1:
            for i in range(0,len(barcode)):
                    for z in ['A','T', 'C', 'G']:
                        barcode_mis = replace_char(string= barcode, char = z, index = i)
                        if mismatch >1:
                            barcode_pools += build_barcode_pool(barcode1 = barcode, barcode2 = "", mismatch = mismatch -1, reverse = False, version = version, mode = mode)
                        else:
                            if version == "ECR3":
                                barcode_pools.append(barcode_mis)
                            elif version == "ECR4":
                                barcode_pools.append(barcode_mis)
        elif mismatch == 0:
            if version == "ECR3":
                barcode_pools.append(barcode)
            elif version == "ECR4":
                barcode_pools.append(barcode)
    else:  # 双barcode
        if reverse == True:
            barcode1 = reverse_complement(barcode1)
            barcode2 = reverse_complement(barcode2)
            
        barcode_pools = []
        if mismatch >= 1:
            for i in range(0,len(barcode1)):
                for j in range(0,len(barcode2)):
                    for z in ['A','T', 'C', 'G']:
                        barcode1_mis = replace_char(string= barcode1, char = z, index = i)
                        for x in ['A','T', 'C', 'G']:
                            barcode2_mis = replace_char(string=barcode2, char = x, index = j)
                            if mismatch >1:
                                barcode_pools += build_barcode_pool(barcode1 = barcode1_mis, barcode2 = barcode2_mis, mismatch = mismatch -1, reverse = False, version = version, mode = mode)
                            else:
                                if version == "ECR3":
                                    if mode == "PE":
                                        barcode_pools.append(barcode1_mis + barcode2_mis)
                                    elif mode == "SE":
                                        barcode_pools.append(barcode2_mis + barcode1_mis)
                                elif version == "ECR4":
                                    barcode_pools.append(barcode1_mis + barcode2_mis)
        elif mismatch == 0:
            if version == "ECR3":
                if mode == "PE":
                    barcode_pools.append(barcode1 + barcode2)
                elif mode == "SE":
                    barcode_pools.append(barcode2 + barcode1)
            elif version == "ECR4":
                barcode_pools.append(barcode1 + barcode2)
                
    return(list(set(barcode_pools)))

def build_mix_barcode_pool(barcodes1, barcodes2, mismatch = 1, reverse = True, version = "ECR4", mode = "PE"):
    mix_barcode_pool = []
    for i in barcodes1:
        for j in barcodes2:
            mix_barcode_pool += build_barcode_pool(barcode1 = i, barcode2 = j, mismatch = mismatch, reverse = reverse, version = version)
    return(list(set(mix_barcode_pool)))

def read_fq(files, line = 100000):
    fastq = {}
    eof = False
    for i in range(0, line):
        header = files[0].readline()
        if header == "":
            eof = True
            break
        if "R1" in fastq:
            fastq['R1'].append([header.strip('\n'), files[0].readline().strip('\n'), files[0].readline().strip('\n'), files[0].readline().strip('\n')])
            if len(files) >1:
                fastq['R2'].append([files[1].readline().strip('\n'), files[1].readline().strip('\n'), files[1].readline().strip('\n'), files[1].readline().strip('\n')])
        else: 
            fastq["R1"] = [[header.strip('\n'), files[0].readline().strip('\n'), files[0].readline().strip('\n'), files[0].readline().strip('\n')]]
            if len(files) > 1:
                fastq["R2"] = [[files[1].readline().strip('\n'), files[1].readline().strip('\n'), files[1].readline().strip('\n'), files[1].readline().strip('\n')]]
    if len(fastq) != 0:
        fastq['EOF'] = eof
    return(fastq)

def get_barcode(fq1, fq2, start, end): # barcode在R2上
    barcode = [fq2[0], fq2[1][(start - len(fq1[1]) - 1):(end - len(fq1[1]))], fq2[2], fq2[3][(start - len(fq1[1]) - 1):(end - len(fq1[1]))]]
    return(barcode)

def get_barcode_se(fq1, start, end): # barcode在R2上
    barcode = [fq1[0], fq1[1][(start - 1):(end)], fq1[2], fq1[3][(start - 1):(end)]]
    return(barcode)

def format_fq_se(fq1, library):
    fastq = {}
    for i in library:
        for j in library[i]:
            start = int(library[i][j][0])
            end = int(library[i][j][1])
            fastq[j] = [fq1[0], fq1[1][(start - 1):end], fq1[2], fq1[3][(start - 1):end]]
    return(fastq)

def format_fq(fq1, fq2, library):
    fastq = {}
    for i in library:
        for j in library[i]:
            start = int(library[i][j][0])
            end = int(library[i][j][1])
            if start < len(fq1[1]):
                fastq[j] = [fq1[0], fq1[1][(start - 1):end], fq1[2], fq1[3][(start - 1):end]]
            else:
                fastq[j] = [fq2[0], fq2[1][(start - len(fq1[1]) - 1):(end - len(fq1[1]))], fq2[2], fq2[3][(start - len(fq1[1]) - 1):(end - len(fq1[1]))]]
    return(fastq)

def main_split(fq, library, mode = 'PE'): # 仅支持两个位置的barcode，超过两个不行   ## barcode总长度长的优先 
    fastq = {}
    for i in range(len(fq['R1'])):
        ambiguous = 0
        fq_tmp = {}
        
        for j in sorted([int(x) for x in list(library.keys())], reverse=True):
            j = str(j)
            for loc in library[j]:
                if loc.find(',') >= 0:  # 仅支持两个位置的barcode，超过两个不行，即只支持双barcode
                    locs = loc.split(',')
                    locs = locs[0].split('-') + locs[1].split('-')
                    if mode == 'PE':
                        barcode = get_barcode(fq1 = fq['R1'][i], fq2 = fq['R2'][i], start = int(locs[0]), end = int(locs[1]))[1] + get_barcode(fq1 = fq['R1'][i], fq2 = fq['R2'][i], start = int(locs[2]), end = int(locs[3]))[1]
                    elif mode == "SE":
                        barcode = get_barcode_se(fq1 = fq['R1'][i], start = int(locs[0]), end = int(locs[1]))[1] + get_barcode_se(fq1 = fq['R1'][i], start = int(locs[2]), end = int(locs[3]))[1]
                else:
                    locs = loc.split('-')
                    if mode == "PE":
                        barcode = get_barcode(fq1 = fq['R1'][i], fq2 = fq['R2'][i], start = int(locs[0]), end = int(locs[1]))[1]
                    elif mode =="SE":
                        barcode = get_barcode_se(fq1 = fq['R1'][i], start = int(locs[0]), end = int(locs[1]))[1]
                for sample in library[j][loc]:
                    if barcode in library[j][loc][sample][2]:
                        if len(fq_tmp) > 0:
                            ambiguous = 1
                            break
                        if sample in fq_tmp:
                            if mode == "PE":
                                fq_tmp[sample].append(format_fq(fq1 = fq['R1'][i], fq2 = fq['R2'][i], library = library[j][loc][sample][3]))
                            elif mode == "SE":
                                fq_tmp[sample].append(format_fq_se(fq1 = fq['R1'][i], library = library[j][loc][sample][3]))
                        else:
                            if mode == "PE":
                                fq_tmp[sample] = [format_fq(fq1 = fq['R1'][i], fq2 = fq['R2'][i], library = library[j][loc][sample][3])]
                            elif mode == "SE":
                                fq_tmp[sample] = [format_fq_se(fq1 = fq['R1'][i], library = library[j][loc][sample][3])]
#                        break #可加速分析，但可能barcode之间设计的时候就较近
                if ambiguous == 1:
                    break
            if ambiguous == 1:
                break
        if ambiguous == 1:
            if 'ambiguous' in fastq:
                if 'R1' in fastq['ambiguous']:
                    fastq['ambiguous']['R1'].append(fq['R1'][i])
                    if mode == "PE":
                        fastq['ambiguous']['R2'].append(fq['R2'][i])
                else:
                    fastq['ambiguous']['R1'] = [fq['R1'][i]]
                    if mode == "PE":
                        fastq['ambiguous']['R2'] = [fq['R2'][i]]
            else:
                if mode == "PE":
                    fastq['ambiguous']= {'R1':[fq['R1'][i]], 'R2':[fq['R2'][i]]}
                elif mode == "SE":
                    fastq['ambiguous']= {'R1':[fq['R1'][i]]}
        elif len(fq_tmp) == 0:
            if 'undecoded' in fastq:
                if 'R1' in fastq['undecoded']:
                    fastq['undecoded']['R1'].append(fq['R1'][i])
                    if mode == "PE":
                        fastq['undecoded']['R2'].append(fq['R2'][i])
                else:
                    fastq['undecoded']['R1'] = [fq['R1'][i]]
                    if mode == "PE":
                        fastq['undecoded']['R2'] = [fq['R2'][i]]
            else:
                if mode == "PE":
                    fastq['undecoded']= {'R1':[fq['R1'][i]], 'R2':[fq['R2'][i]]}
                elif mode == "SE":
                    fastq['undecoded']= {'R1':[fq['R1'][i]]}
        else:
            for sample in fq_tmp:
                for x in fq_tmp[sample]:
                    for tag in x.keys():
                        if sample in fastq:
                            if tag in fastq[sample]:
                                fastq[sample][tag].append(x[tag])
                            else:
                                fastq[sample][tag] = [x[tag]]
                        else:
                            fastq[sample]= {tag:[x[tag]]}
    return(fastq)

def library_build(library_file, seq_mode = "PE100", mismatch = 1, total_cycle = 220, version = "ECR4"):  
    lib = {}
    flag = {}
    library_info = {}
    if seq_mode.upper().startswith('PE'):
        mode = "PE"
    else:
        mode = "SE"
    if seq_mode == "SE_Proton":
        reverse = False
    else:
        reverse = True
        
    with open(library_file, mode='r', encoding='utf-8') as f:
        for line in f.readlines():
            if line.upper().startswith('SAMPLE'):
                continue
            if line.upper().startswith('LIBRARY'):
                continue
            if line.strip('\n') == "":
                continue
            d = line.strip('\n').split('\t')
            if 'sample' in library_info:
                library_info['sample'].append(d[0].replace(' ', ''))
            else:
                library_info['sample'] = [d[0].replace(' ', '')]
            if 'barcode1' in library_info:
                library_info['barcode1'].append(d[1].replace(' ', ''))
            else:
                library_info['barcode1'] = [d[1].replace(' ', '')]
            if 'barcode2' in library_info:
                library_info['barcode2'].append(d[2].replace(' ', ''))
            else:
                library_info['barcode2'] = [d[2].replace(' ', '')]
            if 'library_type' in library_info:
                library_info['library_type'].append(d[3].replace(' ', ''))
            else:
                library_info['library_type'] = [d[3].replace(' ', '')]
            if 'bed' in library_info:
                library_info['bed'].append(d[4].replace(' ', ''))
            else:
                library_info['bed'] = [d[4].replace(' ', '')]
                
    for i in range(len(library_info['sample'])):
        library = parse_library_type(library_type = library_info['library_type'][i], total_cycle = total_cycle)
        tag = []
        for j in library.keys():
            tag += list(library[j].keys())
            tag = list(set(tag))
        barcode_len = 0
        if library_info['barcode1'][i].find(',') >= 0 or library_info['barcode2'][i].find(',') >= 0:
            barcode = build_mix_barcode_pool(barcodes1 = library_info['barcode1'][i].split(','), barcodes2 = library_info['barcode2'][i].split(','), mismatch = mismatch, reverse = reverse, version = version, mode = mode)
            barcode_len = len(library_info['barcode1'][i].split(',')[0]) + len(library_info['barcode2'][i].split(',')[0])
        else:
            barcode = build_barcode_pool(barcode1 = library_info['barcode1'][i], barcode2 = library_info['barcode2'][i], mismatch = mismatch, reverse = reverse, version = version, mode = mode)
            barcode_len = len(library_info['barcode1'][i]) + len(library_info['barcode2'][i])
        
        index = ''
        if len(library['I'].keys()) >1:
            for n in range(1, len(library['I'].keys()) +1):
                if index == '':
                    index = str(library['I']['I'+str(n)][0])+'-'+str(library['I']['I'+str(n)][1])
                else:
                    index += ','+str(library['I']['I'+str(n)][0])+'-'+str(library['I']['I'+str(n)][1])
        else:
            index = str(library['I']['I1'][0])+'-'+str(library['I']['I1'][1])
        
        if str(barcode_len) in lib:
            if index in lib[str(barcode_len)]:
                lib[str(barcode_len)][index][library_info['sample'][i]] = [library_info['sample'][i], list(library["I"].values()), barcode, library, library_info['bed'][i]]
            else:
                lib[str(barcode_len)][index] = {library_info['sample'][i]:[library_info['sample'][i], list(library["I"].values()), barcode, library, library_info['bed'][i]]}
        else:
            lib[str(barcode_len)] = {index:{library_info['sample'][i]:[library_info['sample'][i], list(library["I"].values()), barcode, library, library_info['bed'][i]]}}
        
        flag[library_info['sample'][i]] = tag
    
    if seq_mode.upper().startswith('PE'):
        flag['undecoded'] = ['R1', 'R2']
        flag['ambiguous'] = ['R1', 'R2']
    else:
        flag['undecoded'] = ['R1']
        flag['ambiguous'] = ['R1']
    return([lib, flag])

#### 多进程
def estimate_compression_threads(cores: int):
    return max(0, min(cores, 4))

def multi_split(r1_file, r2_file, lib, outdir, flag, line = 1000000, cpu = 20, mode = "PE"):
    if r2_file == "" and mode == "SE":
        if r1_file.endswith('.gz'):
            r1 = gzip.open(r1_file, 'rt')
        else:
            r1 = open(r1_file, 'rt')
    else:
        if r1_file.endswith('.gz'):
            r1 = gzip.open(r1_file, 'rt')
            r2 = gzip.open(r2_file, 'rt')
        else:
            r1 = open(r1_file, 'rt')
            r2 = open(r2_file, 'rt')
    
    file_opener = FileOpener(compression_level=6, threads=estimate_compression_threads(cores = cpu))
    readQueue = multiprocessing.Queue(maxsize=5000)
    if mode == "SE":
        readthread = readThread('read', [r1], line, readQueue)
    elif mode == "PE":
        readthread = readThread('read', [r1, r2], line, readQueue)
    readthread.start()
        
    workQueue = {}
    writethreads ={}
    writehandle = {}
    write = 0
    for s in flag.keys():
        for f in flag[s]:
            if s in writehandle:
                writehandle[s][f] = file_opener.xopen(os.path.join(outdir, s+'_'+f+'.fastq.gz'), "w")
            else:
                writehandle[s] = {f:file_opener.xopen(os.path.join(outdir, s+'_'+f+'.fastq.gz'), "w")}
        workQueue[s] = multiprocessing.Queue()
        writethreads[s] = writeThread('write'+str(write), workQueue[s], writehandle[s])
        writethreads[s].start()
        write += 1
        
            # if s in workQueue:
                # workQueue[s][f] = multiprocessing.Queue()
            # else:
                # workQueue[s] = {f:multiprocessing.Queue()}
            # if s in writethreads:
                # writethreads[s][f] = writeThread('write'+str(write), workQueue[s][f], writehandle[s+f])
            # else:
                # writethreads[s] = {f:writeThread('write'+str(write), workQueue[s][f], writehandle[s+f])}
            # writethreads[s][f].start()
            # write += 1
        
    threadID = 1
    threads = []
    Lock = multiprocessing.Lock()
    # Lock = multiprocessing.RLock()
    for i in range(cpu):
        thread = workThread(threadID, readQueue, lib, workQueue, Lock, mode)
        thread.start()
        threads.append(thread)
        threadID += 1
        
    while True:
        if not readthread.is_alive():
            if readQueue.empty():
                readQueue.close()
                del readthread
                del readQueue
                r1.close()
                del r1
                if mode == "PE":
                    r2.close()
                    del r2
                break
        
    while len(threads) > 0:
        for t in range(len(threads)):
            if threads[t].is_alive() == False:
                del(threads[t])
                break
        
    tag = 0
    num_write = len(workQueue.keys())
    while len(threads) == 0 and tag != num_write:
        for s in workQueue.keys():
            # for f in workQueue[s].keys():
            if workQueue[s].empty():
                workQueue[s].put(['EOF'], block=True, timeout=None)
                tag += 1
    
    while len(writethreads) > 0:
        for s in writethreads.keys():
            if not writethreads[s].is_alive():
                del  writethreads[s]
                # for f in writehandle[s].keys():
                    # writehandle[s][f].close()
                break
    while len(workQueue) > 0:
        for s in workQueue.keys():
            if workQueue[s].empty():
                workQueue[s].close()
                del workQueue[s]
            break

class workThread (multiprocessing.Process):
    def __init__(self, threadID, readQueue, lib, workQueue, lock, mode):
        multiprocessing.Process.__init__(self)
        self.threadID = threadID
        self.fastq = readQueue
        self.lib = lib
        self.q = workQueue
        self.lock = lock
        self.mode = mode
        self._running = True
    # def terminate(self):
        # self._running = False
    def run(self):
        while self._running:
            self.lock.acquire()
            try:
                fq = self.fastq.get(block=True, timeout=5) # 5秒
            # except multiprocessing.TimeoutError :
                # self._running = False
            # except ValueError:
                # self._running = False
            # except queue.Empty:
                # self._running = False
            # except multiprocessing.TimeoutError:
                # self.lock.release()
                # continue
                # # self._running = False
            except:
                self._running = False
            self.lock.release()
            if self._running:
                if len(fq) != 0:
                    # self.lock.acquire(block=True, timeout=None)
                    fastq = main_split(fq = fq, library = self.lib, mode = self.mode)
                    self.lock.acquire()
                    for s in fastq.keys():
                        # for f in fastq[s].keys():
                        self.q[s].put([s, fastq[s]], block=True, timeout=None)
                    if fq['EOF']:
                        self._running = False
                    self.lock.release()
                else:
                    self._running = False

class readThread (multiprocessing.Process):
    def __init__(self, threadID, file, line, readQueue):
        multiprocessing.Process.__init__(self)
        self.threadID = threadID
        self.file = file
        self.line = line
        self.q = readQueue
        self._running = True
    # def terminate(self):
        # self._running = False
    def run(self):
        while self._running:
            fastq = read_fq(files = self.file, line = self.line)
            if len(fastq) != 0:
#                queueLock.acquire()
                self.q.put(fastq, block=True, timeout=None)
                if fastq['EOF']:
                    self._running = False
#                queueLock.release()
            else:
                self.q.put({'EOF':True}, block=True, timeout=None)
                self._running = False

class writeThread (multiprocessing.Process):
    def __init__(self, threadID, workQueue, file):
        multiprocessing.Process.__init__(self)
        self.threadID = threadID
        self.q = workQueue
        self.file = file
        self._running = True
    def stop(self):
        self._running = False
    def run(self):
        while self._running:
            # fcntl.flock(self.file.outfile, fcntl.LOCK_EX)
        # while self._running:
            try:
                fq = self.q.get(block=True, timeout=None)
                if len(fq) == 0:
                    self._running = False
                elif fq[0] == "EOF":
                    for tag in self.file:
                        self.file[tag].close()
                    self._running = False
                else:
                    for tag in fq[1].keys(): 
                        for j in fq[1][tag]:
                            # fcntl.flock(self.self.file[tag].outfile, fcntl.LOCK_EX)
                            n = self.file[tag].write("\n".join(j)+"\n")
                            # fcntl.flock(self.self.file[tag].outfile, fcntl.LOCK_UN)
            except:
                self._running = False

class FileOpener:
    def __init__(self, compression_level: int = 6, threads: int = None):
        """
        threads -- no. of external compression threads.
            0: write in-process
            None: min(cpu_count(0, 4)
        """
        self.compression_level = compression_level
        self.threads = threads
    def xopen(self, path, mode):
        threads = self.threads if "w" in mode else 0
        f = xopen(
            path, mode, compresslevel=self.compression_level, threads=self.threads
        )
        return f

def main():
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage, formatter_class = argparse.RawTextHelpFormatter)
    ArgParser.add_argument("--version", action="version", version=version)
    ArgParser.add_argument("--read1", action="store", type=str, default=None, required=True, help="fastq of undecoded read 1. required")
    ArgParser.add_argument("--read2", action="store", type=str, default=None, required=False,  help="fastq of undecoded read 2. required if PE mode")
    ArgParser.add_argument("--library", action="store", type=str, default=None, required=True, 
                        help=textwrap.dedent('''\
                        Barcode information for samples. 
                        format:
                        Sample\tBarcode_1\tBarcode_2\tLibrary_structure\tPanel
                        ZLC210244A\tCGTCGTGTTA,GAAGAGAACA,TTGTGAACTC,AGGACATCAG\tGTCATCGACG,CGCCTTATAT,GAGCACGATC,CCGAAGACTG\tM8S1R91S1M8S1R91S1I10I10\tNad
                        '''))
    ArgParser.add_argument("--mode", action="store", type=str, default='PE100', metavar=['PE100', 'PE150', 'SE100', 'SE50', 'SE_Proton'], help="Seq mode. [%(default)s]")
    ArgParser.add_argument("--outdir", action="store", type=str, default='./fastq', help="Output dir. [%(default)s]")
    ArgParser.add_argument("--cpus", action="store", type=int, default=40, help="Number of cpus. [%(default)s]")
    ArgParser.add_argument("--mismatch", action="store", type=int, default=1,  help="Mismatch of each barcode. [%(default)s]")
    ArgParser.add_argument("--total_cycle", action="store", type=int, default=0,  help="Total cycle, 0 for auto (PE100: 220, PE150: 320, SE100: 120, SE50: 70). [%(default)s]")
    ArgParser.add_argument("--mgi_version", action="store", type=str, default='ECR3', metavar=['ECR3', 'ECR4'], help="MGISEQ Instrument control software version. [%(default)s]")
    
    args = ArgParser.parse_args()
    if args.total_cycle == 0:
        if args.mode == "PE100":
            args.total_cycle = 220
        elif args.mode == "PE150":
            args.total_cycle = 320
        elif args.mode == "SE100":
            args.total_cycle = 120
        elif args.mode == "SE50":
            args.total_cycle = 70
        else:
            raise Exception('--mode:{}, 需要指定--total_cycle,不能为: {}。请调整参数!'.format(args.mode, args.total_cycle))
    
    if(os.path.isfile(args.read1) and os.path.isfile(args.library)):
        if args.mode.upper().startswith('PE'):
            if not os.path.isfile(args.read2):
                ArgParser.error("file not exists in options --read2 in %s mode" % (args.mode))
        startTime = datetime.datetime.now()
        print ('Start Time: %s\n' % startTime.strftime("%Y-%m-%d %H:%M:%S"))
        
        lib, flag = library_build(library_file = args.library, seq_mode = args.mode, mismatch = args.mismatch, total_cycle = args.total_cycle, version = args.mgi_version)
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        if args.mode.upper().startswith('PE'):
            multi_split(r1_file = args.read1, r2_file = args.read2, line = 1000, cpu = args.cpus, lib = lib, outdir = args.outdir, flag = flag, mode = "PE")
        else:
            multi_split(r1_file = args.read1, r2_file = "", line = 1000, cpu = args.cpus, lib = lib, outdir = args.outdir, flag = flag, mode = "SE")
        endTime = datetime.datetime.now()
        print ('End Time: %s\n' % endTime.strftime("%Y-%m-%d %H:%M:%S"))
        print ('Totally spend time %s\n' % (endTime - startTime))
    else:
        ArgParser.error("file not exists in options --read1 or --library")

if __name__ == "__main__":
    main()

