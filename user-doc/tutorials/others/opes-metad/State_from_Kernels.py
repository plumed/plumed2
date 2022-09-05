#! /usr/bin/env python3

### Generate an OPES STATE file from a KERNELS file ###
# For postprocessing only, do not use for restarting a simulation
# (the idea is to fake a restart with the plumed driver and dump the OPES state)

import sys
import argparse
import subprocess
plumed_exe='plumed'

error='--- ERROR: %s '
#parser
parser = argparse.ArgumentParser(description='Generate an OPES STATE file from a KERNELS file, so that it can be used with FES_from_State.py. DO NOT use the obtained STATE file for restart')
parser.add_argument('--kernels','-f',dest='filename',type=str,default='KERNELS',help='the kernels file name, with the deposited kernels')
parser.add_argument('--outfile','-o',dest='outfile',type=str,default='STATE',help='name of the output file')
parser.add_argument('--tmpname',dest='tmpname',type=str,default='tmp-plumed_driver.dat',help='name of the temporary plumed file')
args = parser.parse_args()
#parsing
filename=args.filename
outfile=args.outfile
tmp_plumed_file=args.tmpname

#get info
f=open(filename,'r')
line=f.readline() #fields
if line.split()[1]!='FIELDS':
  sys.exit(error%(' no FIELDS found in file "'+filename+'"'))
if len(line.split())<7:
  sys.exit(error%(' not enough FIELDS found in file "'+filename+'"'))
if (len(line.split())-5)%2!=0:
  sys.exit(error%(' wrong number of FIELDS found in file "'+filename+'"'))
ncv=int((len(line.split())-5)/2)
cvname=[]
for i in range(ncv):
  cvname.append(line.split()[3+i])
  if line.split()[3+ncv+i]!='sigma_'+cvname[i]:
    sys.exit(error%(' expected "sigma_%s" instead of "%s"'%(cvname[i],line.split()[4+i])))
line=f.readline() #action
if line.split()[3]=='OPES_METAD_kernels':
  action='OPES_METAD'
elif line.split()[3]=='OPES_METAD_EXPLORE_kernels':
  action='OPES_METAD_EXPLORE'
else:
  sys.exit(error%(' this script only works with OPES_METAD or OPES_METAD_EXPLORE KERNELS files'))
line=f.readline() #biasfactor
if line.split()[2]!='biasfactor':
  sys.exit(error%(' biasfactor not found!'))
biasfactor=line.split()[3]
line=f.readline() #epsilon
if line.split()[2]!='epsilon':
  sys.exit(error%(' epsilon not found!'))
epsilon=line.split()[3]
line=f.readline() #kernel_cutoff
if line.split()[2]!='kernel_cutoff':
  sys.exit(error%(' kernel_cutoff not found!'))
kernel_cutoff=line.split()[3]
line=f.readline() #compression_threshold
if line.split()[2]!='compression_threshold':
  sys.exit(error%(' compression_threshold not found!'))
compression_threshold=line.split()[3]
periodic=['NO']*ncv
line=f.readline()
while line.split()[0]=='#!':
  for i in range(ncv):
    if line.split()[2]=='min_'+cvname[i]:
      periodic[i]=line.split()[3]+','
      line=f.readline()
      if line.split()[2]!='max_'+cvname[i]:
        sys.exit(error%(' periodic CVs should have both min and max value!'))
      periodic[i]+=line.split()[3]
  line=f.readline()
f.close()

#unfourtunately plumed does not allow for CVs to be named with a dot
#for now we only support .x .y .z in no periodic CVs
suffix=['']*ncv
for i in range(ncv):
  if cvname[i].find('.')!=-1:
    suffix[i]=cvname[i][cvname[i].find('.'):]
    if ((suffix[i]=='.x' or suffix[i]=='.y' or suffix[i]=='.z') and periodic[i]=='NO'):
      cvname[i]=cvname[i][:cvname[i].find('.')]
    else:
      sys.exit(error%(' %s: you must modify the KERNELS file and remove any "." from CVs names'%cvname[i]))

#create temporary plumed file
plumed_input='# Temporary file used to convert an opes KERNELS file into a STATE file\n'
plumed_input+='# vim:ft=plumed\nUNITS NATURAL\n'
plumed_input+='RESTART\n'
plumed_input+='f: FIXEDATOM AT=0,0,0\n' #fake atom
plumed_input+='d: DISTANCE ATOMS=f,f\n' #unfourtunately the FAKE colvar has issues with PERIODIC
plumed_input+='COMMITTOR ARG=d BASIN_LL1=-1 BASIN_UL1=1\n' #this will kill the driver
for i in range(ncv):
  if suffix[i]=='':
    plumed_input+=cvname[i]+': COMBINE ARG=d PERIODIC='+periodic[i]+'\n' #recreate CVs label doesn't work if they have components!
  else:
    if not cvname[i] in cvname[i+1:]: #avoid duplicates
      plumed_input+=cvname[i]+': DISTANCE ATOMS=f,f COMPONENTS\n'
    cvname[i]+=suffix[i]
plumed_input+='opes: '+action
plumed_input+=' ARG='+cvname[0]
for ii in range(1,ncv):
  plumed_input+=','+cvname[ii]
plumed_input+=' FILE='+filename
plumed_input+=' STATE_WFILE='+outfile
plumed_input+=' BIASFACTOR='+biasfactor
plumed_input+=' EPSILON='+epsilon
plumed_input+=' KERNEL_CUTOFF='+kernel_cutoff
plumed_input+=' COMPRESSION_THRESHOLD='+compression_threshold
plumed_input+=' TEMP=1 PACE=1 BARRIER=1 SIGMA=1'
for ii in range(1,ncv):
  plumed_input+=',1'

f=open(tmp_plumed_file,'w')
f.write(plumed_input)
f.close()

#run driver
cmd_string=plumed_exe+' driver --noatoms --plumed '+tmp_plumed_file
cmd=subprocess.Popen(cmd_string,shell=True)
cmd.wait()
print(' +++ IMPORTANT: do not use the obtained STATE file for restarting a simulation +++')
