#!/usr/bin/env python
### submit cosmomc jobs
#### by  Zhiqi Huang (zqhuang@cita.utoronto.ca)
import re
import os
import sys
import glob
import string

def search_value(fname, pattern):
    fp = open(fname, 'r')
    file_content = fp.read()
    fp.close()
    m = re.search(pattern, file_content, flags = re.M + re.I)
    if m:
        return m.group(1)
    else:
        return ""

if(len(sys.argv) < 2):
    print "Syntax: "
    print "python submit.py [R-1 threshold]"
    sys.exit()
inifile = sys.argv[1]
if( not os.path.isfile(inifile)):
    print inifile + " does not exist"
    sys.exit()
threshold = 0.03
if(len(sys.argv) >=  3):
    try:
        threshold = float(sys.argv[2])
    except:
        print "R-1 threshold format is incorrect"
        sys.exit()
fileroot = search_value(inifile, r'file\_root\s*\=\s*(\S+)')
rootdir = search_value(inifile, r'root\_dir\s*\=\s*(\S+)')
covmat = search_value(inifile, r'propose\_matrix\s*\=\s*(\S+)')
if(fileroot == ''):
    print "ini file does not contain key file_root"
    sys.exit()
fileroot = rootdir+fileroot    
print "file_root = " + fileroot
if(covmat != ''):
    if(os.path.isfile(covmat)):
        print "propose_matrix = " + covmat        
    else:
        print "propose_matrix " + covmat + " does not exist"
        sys.exit()
else:
    print "no propose_matrix"
if(os.path.isfile(fileroot + r'.converge_stat')):
    fp = open(fileroot + r'.converge_stat', 'r')
    conv = fp.read()
    fp.close()    
    try:
        rm = float(conv)
    except:
        rm = 1000.
    if(rm < threshold):
        print "chains are already converged, not submitting the job."
        sys.exit()

print "submitting " + inifile

current_path = os.getcwd()
patterns = [r'.*\/', r'scan\_', r'fixrp(\d\d)\d', r'dpp', r'dpl', r'qcdm\_1param', r'qcdm\_3param', 'lowTEB_', 'plikTTTEEE', 'plikTT', 'plikTEEE', 'plikEE', 'BAO_JLA_HSTlow', 'BAO', 'lens', 'RSD_WL', 'RSD', 'WL', r'\.ini']
repls = ['', '', r'r\1', 'S', 'L', 'q1', 'q3', '', 'A', 'T', 'TE', 'E', 'pr', 'B', 'l', 'RW', 'R','W', '']
jobname = inifile
for i in range(len(patterns)):
    jobname = re.sub(patterns[i], repls[i], jobname)
    
fp = open(r'scripts/' + jobname + r'.jb', 'w')
fp.write(r'#!/bin/csh -f' + "\n" + r'#PBS -N '+jobname + "\n" + r'#PBS -l nodes=8:ppn=16' + "\n" + r'#PBS -q hpq' + "\n" + r'#PBS -l walltime=48:00:00' + "\n" + r'#PBS -r n'  + "\n" + r'#PBS -j oe'+ "\n" + r'cd ' + current_path + "\n" + 'mpirun -pernode ./cosmomc ' + inifile + ' > ./scripts/'+jobname+r'.log' + "\n")

fp.close()

os.chdir('scripts')
os.system('qsub ' + jobname + r'.jb')
os.chdir(current_path)


    


