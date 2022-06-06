#coding:utf8
import subprocess
import sys
import os
import re
import gzip
import glob
import json
import time
from time import sleep
#说明：CreateDirectory, CheckFile, DetermainR1R2, RunCommandWithReturn, ExportLog, ResolveConfig, Versions, RunningProcess, GetSize, 验证是正确的

# Description: create directory or not.
# input      : all path direcotry
# return     : None
def CreateDirectory(direcotry):
    if not os.path.exists(direcotry):
        os.mkdir(direcotry)
    else:
        print ('%s directory all exists!'%direcotry)
# Description: check file existed or not.
# input      : file
# return     : True or False
def CheckFile(filename):
    if filename is None:
        return False
    elif not os.path.exists(filename):
        return False
    else:
        return True
# Description: determain fastq file read1 is real or not, read2 is real or not
# input      : read1 and read2
# return     : read1 and read2 are the real file return True otherwize return False
def DetermainR1R2(fastq_read1,fastq_read2):
    r1 = False
    r2 = False
    try:
        try:
            r1_head = open(fastq_read1,'r').readline().split(' ')[1][0] #illumina
            r2_head = open(fastq_read2,'r').readline().split(' ')[1][0]
        except IndexError:
            r1_head = open(fastq_read1,'r').readline().strip().split('/')[1] #MGI
            r2_head = open(fastq_read2,'r').readline().strip().split('/')[1]
    except UnicodeDecodeError:
        try:
            r1_head = bytes.decode(gzip.open(fastq_read1,'r').readline()).split(' ')[1][0] #illumina
            r2_head = bytes.decode(gzip.open(fastq_read2,'r').readline()).split(' ')[1][0]
        except IndexError:
            r1_head = bytes.decode(gzip.open(fastq_read1,'r').readline()).strip().split('/')[1] #MGI
            r2_head = bytes.decode(gzip.open(fastq_read2,'r').readline()).strip().split('/')[1]
    if r1_head != '1':
        r1 = True
        print ('read1 file is wrong; it is %s file'%fastq_read1)
    if r2_head != '2':
        r2 = True
        print ('read2 file is wrong; it is %s file'%fastq_read2)
    result = r1 or r2
    return result
# Description: run command and return stdout and stderr string
# input      : command line
# return     : stdout and stderr string
def RunCommandWithReturn(cmd, shell = True):
    start_time = 'Start time:\t'+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    proc = subprocess.Popen(cmd,
                shell = shell,
                universal_newlines = True,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    if return_value != 0:
        sys.stderr.write('\n\tError Message: ' + stderr + '\n')
        sys.stderr.write('\n\tStdout Message: ' + stdout + '\n')
        raise OSError('Error during: "{}"; at the time {}'.format(cmd, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime() ) ) )
    # else:
    #     stderr = stderr if stderr else ''
    #     sys.stdout.write('\n\tStderror of {}:'.format(cmd) + '\n\t\t' +  stderr + '\n')
    #     stdout = stdout if stdout else ''
    #     sys.stdout.write('\n\tStdout of {}:'.format(cmd) + '\n\t\t' + stdout + '\n')
    end_time = 'End time:\t'+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    timer = [start_time, end_time]
    return stdout, stderr, timer
# Description: command and stdout and stderr string export log directory or current directory
# input      : command and stdout and stderr messages
# return     : strings export to log file
def ExportLog(command, stdout, stderr, timer, output_dir, output_prefix, other_string = ''):
    outputstring = other_string+'\n'+timer[0]+'\n'+'temporary'+'\n'+timer[1]+'\n'
    dirprefix = os.path.join(output_dir, output_prefix)
    log_dir = output_dir+'/'+'../log'
    logdirprefix = os.path.join(log_dir, output_prefix)
    if os.path.exists(log_dir):
        resultprefix = logdirprefix
    else:
        resultprefix = dirprefix
    open("%s.command_log.txt"%resultprefix,'a').write(outputstring.replace('temporary',command))
    open("%s.stdout_log.txt"%resultprefix,'a').write(outputstring.replace('temporary',stdout))
    open("%s.stderr_log.txt"%resultprefix,'a').write(outputstring.replace('temporary',stderr))
# Description: add the RunCommandWithReturn function and ExportLog function
# input      : input runcommand outputdirectory outputprefix and logkeystring
# return     : None only execute Exportlog
def RunningProcess(cmd, output_dir, output_prefix, other_string):
    (stdout, stderr, timer) = RunCommandWithReturn(cmd)
    ExportLog(cmd, stdout, stderr, timer, output_dir, output_prefix, other_string)
# Description: parse configure file
# input      : input configure file
# return     : config dict
def ResolveConfig(configure_file):
    config_dict = {}
    with open(configure_file,'r') as inopen:
        for line in inopen:
            line = line.strip()
            if line.startswith("//"):
                items = line[2:]
                config_dict.update({items:{}})
            elif line != "" and not line.startswith('#'):
                (parameter, value) = re.split(r"\s*=+\s*", line)
                # parameter=line.split(" = ")[0]
                value = value.replace('None','')
                try:
                    basic = re.findall(r'\$\{(.+)\}', value)[0]
                    value = value.replace('${%s}'%basic, config_dict.get('basic').get(basic))
                except IndexError:
                    pass
                config_dict.get(items).update({parameter:value})
    # with open("%s.json"%configure_file,"w") as otopen:
    #     otopen.write(json.dumps(config_dict))
    return config_dict
# Description: checking the directory all file; return the latest version script and all of the version number
# input      : direcrtory and keystring which used to check the file
# return     : latest version script and all of the version number
def Versions(directory, keystring = 'cfg'):
    allv_list = []
    version_dict = {}
    files = next(os.walk(directory))[2]
    for each in files:
        if keystring in each:
            try:
                pointv = re.findall(r'(\d+)\.(\d+)\.(\d+)\.(\d+)', each)[0]
                allv_list.append('v'+".".join(pointv))
                allnum = 0
                for x in list(range(len(pointv)))[::-1]:
                    allnum += int(pointv[-x-1])*100**x
                version_dict.update({allnum:each})
            except IndexError:
                version_dict.update({0:each})
    compare_num = 0
    lateset = keystring
    if version_dict == {}:
        raise KeyError('%s keystring not found in %s directory all files'%(keystring, direcotry))
    else:
        for k,v in version_dict.items():
            if k >= compare_num:
                lateset = v
    return lateset, allv_list
# Description: get filesize
# input      : cmd check file
# return     : size number int type
def GetSize(infile):
    try:
        size = os.path.getsize(infile)
    except Exception:
        size = 0
    return size
###############################################################
# Description: run command
# input      : cmd string, error tip string
# return     : None
    # def run_command(run_cmd, error_log, doprint = True, error_exit = True):
    #     try:
    #         for cmd in run_cmd.split("\n"):
    #             if cmd.strip() == "":
    #                 continue
    #             if doprint:
    #                 print("[Command]:\t" + cmd)
    #             subprocess.check_call(cmd, shell = True)
    #     except:
    #         if error_exit:
    #             sys.exit(error_log)
# Description: check path is not existed. create directory
# input      : path
# return     : None
# def createDirectory(dirctory):
#     if not os.path.exists(dirctory):
#         run_command('mkdir -p %s' % dirctory, ">>> Create dirctory %s failed." % dirctory, False)

# Description: write README for each module
# input      : output file, version, note list
# return     : return a string
# def write_readme(output_file, module_version, note_list):
#     output_line = "# version: %s\n" % module_version
#     for idx, each_note in enumerate(note_list):
#         output_line += "# %d. %s\n" %(idx+1, each_note)
#     if output_file != "":
#         open(output_file, 'w').write(output_line)
#     return output_line

    # def rm_files(path):
    #     flist = glob.glob(path)
    #     if len(flist)>0:
    #         os.system("rm -rf " + path)

# def check_files(files):
#     error_flag = 0
#     for efile in files:
#         if not os.path.isfile(efile):
#             print("ERROR: Could not find file %s." % efile)
#             error_flag = 1
#         elif os.path.getsize(efile) == 0:
#             print("WARNING: File %s size is zero!" % efile)
#     if error_flag:
#         exit(error_flag)

    # def check_dirs(dirs):
    #     error_flag = 0
    #     for edir in dirs:
    #         if not os.path.isdir(edir):
    #             print("ERROR: Could not find directory %s." % edir)
    #             error_flag = 1
    #     if error_flag:
    #         exit(error_flag)

# Description: check files exists or not
    # def command_running_completed(files): #, tests=6
    #     result=False
    #     if not os.path.exists(files):
    #         sleep(30)
    #     else:
    #         md5code_raw='0'
    #         times=0
    #         for x in range(1,361):
    #             md5code = subprocess.getoutput('md5sum %s'%files).split(' ')[0]
    #             if md5code_raw != md5code:
    #                 md5code_raw = md5code
    #                 sleep(10)
    #             elif md5code_raw == md5code:
    #                 times+=1
    #                 sleep(10)
    #                 if times == 6: #tests
    #                     result=True
    #                     break
    #     return result

if __name__ == '__main__':
    pass
