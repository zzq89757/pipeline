import os
import sys

if len(sys.argv) != 3:
    print ('python %s samtools unmapped_bam'%sys.argv[0])
    exit()

samtools=sys.argv[1]
bam=sys.argv[2]
cmd_view=" ".join([samtools, 'view --threads 16', bam])
for line in os.popen(cmd_view):
    # col=bytes.decode(line).split('\t')
    col=line.split('\t')
    qual_list=[ord(x)-33-33 if ord(x)-33-33 >= 0 else 0 for x in list(col[10])]
    result_q="".join([chr(x+33) for x in qual_list])
    if col[1] == '77': #read1
        print ('@%s/1'%col[0])
    elif col[1] == '141': #read1
        print ('@%s/2'%col[0])
    print (col[9])
    print ('+')
    print (result_q)
