#! /usr/bin/env python3

### Get the FES estimate used by OPES, from a dumped state file (STATE_WFILE). 1D or 2D only ###
# usage is similar to plumed sum_hills

import sys
import argparse
import numpy as np
import pandas as pd #much faster reading from file
do_bck=False #requires the bck.meup.sh script
if do_bck:
  bck_script='bck.meup.sh' #e.g. place the script in your ~/bin
  import subprocess

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the FES estimate used by OPES, from a dumped state file (STATE_WFILE). 1D or 2D only')
# files
parser.add_argument('--state','-f',dest='filename',type=str,default='STATE',help='the state file name, with the compressed kernels')
parser.add_argument('--outfile','-o',dest='outfile',type=str,default='fes.dat',help='name of the output file')
# compulsory
kbt_group = parser.add_mutually_exclusive_group(required=True)
kbt_group.add_argument('--kt',dest='kbt',type=float,help='the temperature in energy units')
kbt_group.add_argument('--temp',dest='temp',type=float,help='the temperature in Kelvin. Energy units is Kj/mol')
# grid related
parser.add_argument('--min',dest='grid_min',type=str,required=False,help='lower bounds for the grid')
parser.add_argument('--max',dest='grid_max',type=str,required=False,help='upper bounds for the grid')
parser.add_argument('--bin',dest='grid_bin',type=str,default="100,100",help='number of bins for the grid')
# other options
parser.add_argument('--fmt',dest='fmt',type=str,default='% 12.6f',help='specify the output format')
parser.add_argument('--deltaFat',dest='deltaFat',type=float,required=False,help='calculate the free energy difference between left and right of given c1 value')
parser.add_argument('--all_stored',dest='all_stored',action='store_true',default=False,help='print all the FES stored instead of only the last one')
parser.add_argument('--nomintozero',dest='nomintozero',action='store_true',default=False,help='do not shift the minimum to zero')
parser.add_argument('--der',dest='der',action='store_true',default=False,help='calculate also FES derivatives')
# some easy parsing
args=parser.parse_args()
filename=args.filename
outfile=args.outfile
if args.kbt is not None:
  kbt=args.kbt
else:
  kbt=args.temp*0.0083144621
fmt=args.fmt
calc_deltaF=False
if args.deltaFat is not None:
  calc_deltaF=True
ts=args.deltaFat
mintozero=(not args.nomintozero)
calc_der=args.der
all_stored=args.all_stored
if all_stored:
  if outfile.rfind('/')==-1:
    prefix=''
    outfile_n=outfile
  else:
    prefix=outfile[:outfile.rfind('/')]
    if prefix+'/'==outfile:
      outfile+='fes_rew.dat'
    outfile_n=outfile[outfile.rfind('/'):]
  if outfile_n.rfind('.')==-1:
    suffix=''
  else:
    suffix=outfile_n[outfile_n.rfind('.'):]
    outfile_n=outfile_n[:outfile_n.rfind('.')]
  outfile_n=prefix+outfile_n+'_%d'+suffix
explore='unset'

### Get data ###
# get data and check number of stored states
data=pd.read_table(filename,sep='\s+',header=None)
fields_pos=[]
tot_lines=len(data.iloc[:,1])
for i in range(tot_lines):
  if data.iloc[i,1]=='FIELDS':
    fields_pos.append(i)
if len(fields_pos)==0:
  sys.exit(' no FIELDS found in file "'+filename+'"')
if len(fields_pos)>1:
  print(' a total of %d stored states where found'%len(fields_pos))
  if all_stored:
    print('  -> all will be printed')
  else:
    print('  -> only the last one will be printed. use --all_stored to instead print them all')
    fields_pos=[fields_pos[-1]]
fields_pos.append(tot_lines)

for n in range(len(fields_pos)-1):
  print('   working...   0% of {:.0%}'.format(n/(len(fields_pos)-1)),end='\r')
  l=fields_pos[n]
  dim2=False
  if len(data.iloc[l,:])==6:
    name_cv_x=data.iloc[l,3]
  elif len(data.iloc[l,:])==8:
    dim2=True
    name_cv_x=data.iloc[l,3]
    name_cv_y=data.iloc[l,4]
  else:
    sys.exit(' wrong number of FIELDS in file "'+filename+'": only 1 or 2 dimensional bias are supported')
  action=data.iloc[l+1,3]
  if action=="OPES_METAD_state":
    if explore!='no':
      explore='no'
      print(' building free energy from OPES_METAD')
  elif action=="OPES_METAD_EXPLORE_state":
    if explore!='yes':
      explore='yes'
      print(' building free energy from OPES_METAD_EXPLORE')
  else:
    sys.exit(' This script works only with OPES_METAD_state and OPES_METAD_EXPLORE_state')
  if data.iloc[l+2,2]!='biasfactor':
    sys.exit(' biasfactor not found!')
  sf=1 #scaling factor for explore mode
  if explore=='yes':
    sf=float(data.iloc[l+2,3])
  if data.iloc[l+3,2]!='epsilon':
    sys.exit(' epsilon not found!')
  epsilon=float(data.iloc[l+3,3])
  if data.iloc[l+4,2]!='kernel_cutoff':
    sys.exit(' kernel_cutoff not found!')
  cutoff=float(data.iloc[l+4,3])
  val_at_cutoff=np.exp(-0.5*cutoff**2)
  if data.iloc[l+6,2]!='zed':
    sys.exit(' zed not found!')
  Zed=float(data.iloc[l+6,3])
  if explore=='no':
    if data.iloc[l+7,2]!='sum_weights':
      sys.exit(' sum_weights not found!')
    Zed*=float(data.iloc[l+7,3])
  if explore=='yes':
    if data.iloc[l+9,2]!='counter':
      sys.exit(' counter not found!')
    Zed*=float(data.iloc[l+9,3])
  l+=10 #there are always at least 10 header lines
# get periodicity
  period_x=0
  period_y=0
  while data.iloc[l,0]=='#!':
    if data.iloc[l,2]=='min_'+name_cv_x:
      if data.iloc[l,3]=='-pi':
        grid_min_x=-np.pi
      else:
        grid_min_x=float(data.iloc[l,3])
      l+=1
      if data.iloc[l,2]!='max_'+name_cv_x:
        sys.exit(' min_%s was found, but not max_%s !'%(name_cv_x,name_cv_x))
      if data.iloc[l,3]=='pi':
        grid_max_x=np.pi
      else:
        grid_max_x=float(data.iloc[l,3])
      period_x=grid_max_x-grid_min_x
      if calc_der:
        sys.exit(' derivatives not supported with periodic CVs, remove --der option')
    if dim2 and data.iloc[l,2]=='min_'+name_cv_y:
      if data.iloc[l,3]=='-pi':
        grid_min_y=-np.pi
      else:
        grid_min_y=float(data.iloc[l,3])
      l+=1
      if data.iloc[l,2]!='max_'+name_cv_y:
        sys.exit(' min_%s was found, but not max_%s !'%(name_cv_y,name_cv_y))
      if data.iloc[l,3]=='pi':
        grid_max_y=np.pi
      else:
        grid_max_y=float(data.iloc[l,3])
      period_y=grid_max_y-grid_min_y
      if calc_der:
        sys.exit(' derivatives not supported with periodic CVs, remove --der option')
    l+=1
  if l==fields_pos[-1]:
    sys.exit(' missing data!')
# get kernels
  center_x=np.array(data.iloc[l:fields_pos[n+1],1],dtype=float)
  if dim2:
    center_y=np.array(data.iloc[l:fields_pos[n+1],2],dtype=float)
    sigma_x=np.array(data.iloc[l:fields_pos[n+1],3],dtype=float)
    sigma_y=np.array(data.iloc[l:fields_pos[n+1],4],dtype=float)
    height=np.array(data.iloc[l:fields_pos[n+1],5],dtype=float)
  else:
    sigma_x=np.array(data.iloc[l:fields_pos[n+1],2],dtype=float)
    height=np.array(data.iloc[l:fields_pos[n+1],3],dtype=float)

### Prepare the grid ###
  grid_bin_x=int(args.grid_bin.split(',')[0])
  if period_x==0:
    grid_bin_x+=1 #same as plumed sum_hills
  if args.grid_min is None:
    if period_x==0: #otherwise is already set
      grid_min_x=min(center_x)
  else:
    if args.grid_min.split(',')[0]=='-pi':
      grid_min_x=-np.pi
    else:
      grid_min_x=float(args.grid_min.split(',')[0])
  if args.grid_max is None:
    if period_x==0: #otherwise is already set
      grid_max_x=max(center_x)
  else:
    if args.grid_max.split(',')[0]=='pi':
      grid_max_x=np.pi
    else:
      grid_max_x=float(args.grid_max.split(',')[0])
  grid_cv_x=np.linspace(grid_min_x,grid_max_x,grid_bin_x)
  if dim2:
    if len(args.grid_bin.split(','))!=2:
      sys.exit('two comma separated integers expected after --bin')
    grid_bin_y=int(args.grid_bin.split(',')[1])
    if period_y==0:
      grid_bin_y+=1 #same as plumed sum_hills
    if args.grid_min is None:
      if period_y==0: #otherwise is already set
        grid_min_y=min(center_y)
    else:
      if len(args.grid_min.split(','))!=2:
        sys.exit('two comma separated floats expected after --min')
      if args.grid_min.split(',')[1]=='-pi':
        grid_min_y=-np.pi
      else:
        grid_min_y=float(args.grid_min.split(',')[1])
    if args.grid_max is None:
      if period_y==0: #otherwise is already set
        grid_max_y=max(center_y)
    else:
      if len(args.grid_max.split(','))!=2:
        sys.exit('two comma separated floats expected after --max')
      if args.grid_max.split(',')[1]=='pi':
        grid_max_y=np.pi
      else:
        grid_max_y=float(args.grid_max.split(',')[1])
    grid_cv_y=np.linspace(grid_min_y,grid_max_y,grid_bin_y)
    x,y=np.meshgrid(grid_cv_x,grid_cv_y)
  if calc_deltaF and (ts<=grid_min_x or ts>=grid_max_x):
    print(' +++ WARNING: the provided --deltaFat is out of the CV grid +++')
    calc_deltaF=False

### Calculate FES ###
  max_prob=0
  if not dim2:
    prob=np.zeros(grid_bin_x)
    if calc_der:
      der_prob_x=np.zeros(grid_bin_x)
    for i in range(grid_bin_x):
      print('   working...  {:.0%} of '.format(i/grid_bin_x),end='\r')
      if period_x==0:
        dist_x=(grid_cv_x[i]-center_x)/sigma_x
      else:
        dx=np.absolute(grid_cv_x[i]-center_x)
        dist_x=np.minimum(dx,period_x-dx)/sigma_x
      kernels_i=height*(np.maximum(np.exp(-0.5*dist_x*dist_x)-val_at_cutoff,0))
      prob[i]=np.sum(kernels_i)/Zed+epsilon
      if calc_der:
        der_prob_x[i]=np.sum(-dist_x/sigma_x*kernels_i)/Zed
      if mintozero and prob[i]>max_prob:
        max_prob=prob[i]
  else:
    prob=np.zeros((grid_bin_y,grid_bin_x))
    if calc_der:
      der_prob_x=np.zeros((grid_bin_y,grid_bin_x))
      der_prob_y=np.zeros((grid_bin_y,grid_bin_x))
    for i in range(grid_bin_y):
      print('   working...  {:.0%} of '.format(i/grid_bin_y),end='\r')
      for j in range(grid_bin_x):
        if period_x==0:
          dist_x=(x[i,j]-center_x)/sigma_x
        else:
          dx=np.absolute(x[i,j]-center_x)
          dist_x=np.minimum(dx,period_x-dx)/sigma_x
        if period_y==0:
          dist_y=(y[i,j]-center_y)/sigma_y
        else:
          dy=np.absolute(y[i,j]-center_y)
          dist_y=np.minimum(dy,period_y-dy)/sigma_y
        kernels_ij=height*(np.maximum(np.exp(-0.5*(dist_x**2+dist_y**2))-val_at_cutoff,0))
        prob[i,j]=np.sum(kernels_ij)/Zed+epsilon
        if calc_der:
          der_prob_x[i,j]=np.sum(-dist_x/sigma_x*kernels_ij)/Zed
          der_prob_y[i,j]=np.sum(-dist_y/sigma_y*kernels_ij)/Zed
        if mintozero and prob[i,j]>max_prob:
          max_prob=prob[i,j]
  if not mintozero:
    max_prob=1
  fes=-kbt*sf*np.log(prob/max_prob)
  if calc_der:
    der_fes_x=-kbt*sf/prob*der_prob_x
    if dim2:
      der_fes_y=-kbt*sf/prob*der_prob_y
# calculate deltaF  
# NB: summing is as accurate as trapz, and logaddexp avoids overflows
  if calc_deltaF:
    if not dim2:
      fesA=-kbt*np.logaddexp.reduce(-1/kbt*fes[grid_cv_x<ts])
      fesB=-kbt*np.logaddexp.reduce(-1/kbt*fes[grid_cv_x>ts])
    else:
      fesA=-kbt*np.logaddexp.reduce(-1/kbt*fes[x<ts])
      fesB=-kbt*np.logaddexp.reduce(-1/kbt*fes[x>ts])
    deltaF=fesB-fesA

### Print to file ###
# prepare file
  if all_stored:
    outfile=outfile_n%(n+1)
  if do_bck:
    cmd=subprocess.Popen(bck_script+' -i '+outfile,shell=True)
    cmd.wait()
# actual print
  f=open(outfile,'w')
  fields='#! FIELDS '+name_cv_x
  if dim2:
    fields+=' '+name_cv_y
  fields+=' file.free'
  if calc_der:
    fields+=' der_'+name_cv_x
    if dim2:
      fields+=' der_'+name_cv_y
  f.write(fields+'\n')
  if calc_deltaF:
    f.write('#! SET DeltaF %g\n'%(deltaF))
  f.write('#! SET min_'+name_cv_x+' %g\n'%(grid_min_x))
  f.write('#! SET max_'+name_cv_x+' %g\n'%(grid_max_x))
  f.write('#! SET nbins_'+name_cv_x+' %g\n'%(grid_bin_x))
  if period_x==0:
    f.write('#! SET periodic_'+name_cv_x+' false\n')
  else:
    f.write('#! SET periodic_'+name_cv_x+' true\n')
  if not dim2:
    for i in range(grid_bin_x):
      line=(fmt+'  '+fmt)%(grid_cv_x[i],fes[i])
      if calc_der:
        line+=(' '+fmt)%(der_fes_x[i])
      f.write(line+'\n')
  else:
    f.write('#! SET min_'+name_cv_y+' %g\n'%(grid_min_y))
    f.write('#! SET max_'+name_cv_y+' %g\n'%(grid_max_y))
    f.write('#! SET nbins_'+name_cv_y+' %g\n'%(grid_bin_y))
    if period_y==0:
      f.write('#! SET periodic_'+name_cv_y+' false\n')
    else:
      f.write('#! SET periodic_'+name_cv_y+' true\n')
    for i in range(grid_bin_y):
      for j in range(grid_bin_x):
        line=(fmt+' '+fmt+'  '+fmt)%(x[i,j],y[i,j],fes[i,j])
        if calc_der:
          line+=(' '+fmt+' '+fmt)%(der_fes_x[i,j],der_fes_y[i,j])
        f.write(line+'\n')
      f.write('\n')
  f.close()
