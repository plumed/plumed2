#! /usr/bin/env python3

### Get the FES estimate used by OPES_WT from a dumped state file (STATE_WFILE), 2D only ###
# slightly similar to plumed sum_hills

import sys
import numpy as np
import pandas as pd
import argparse

#parser
parser = argparse.ArgumentParser(description='get the FES estimate used by OPES from a dumped STATE_WFILE, 2D only')
parser.add_argument('--state',dest='filename',type=str,default='State.data',help='the state file name, with the compressed kernels')
parser.add_argument('--kt',dest='kbt',type=float,required=True,help='the temperature in energy units')
parser.add_argument('--angle1',dest='angle1',action='store_true',default=False,help='the cv1 is an angle in the range [-pi,pi]')
parser.add_argument('--angle2',dest='angle2',action='store_true',default=False,help='the cv2 is an angle in the range [-pi,pi]')
parser.add_argument('--min',dest='grid_min',type=str,required=False,help='lower bounds for the grid')
parser.add_argument('--max',dest='grid_max',type=str,required=False,help='upper bounds for the grid')
parser.add_argument('--bin',dest='grid_bin',type=str,default="100,100",help='number of bins for the grid')
parser.add_argument('--mintozero',dest='mintozero',action='store_true',default=False,help='shift the minimum to zero')
parser.add_argument('--no_der',dest='no_der',action='store_true',default=False,help='skip derivatives to run faster')
parser.add_argument('--outfile',dest='outfile',type=str,default='fes.dat',help='name of the output file')
args = parser.parse_args()
#parsing
filename=args.filename
kbt=args.kbt
mintozero=args.mintozero
calc_der=(not args.no_der)

#get kernels
f=open(filename)
line=f.readline()
if len(line.split())!=8:
  sys.exit('  something is wrong with file'+filename)
cv1name=line.split()[3]
cv2name=line.split()[4]
line=f.readline() #action
if line.split()[-1]=="OPES_WT":
  explore=False
  print(' building free energy from OPES_WT')
elif line.split()[-1]=="OPES_WT_EXPLORE":
  explore=True
  print(' building free energy from OPES_WT_EXPLORE')
else:
  sys.exit("This script works onyl with OPES_WT and OPER_WT_EXPLORE")
line=f.readline() #biasfactor
if explore:
  kbt*=float(line.split()[-1])
line=f.readline() #epsilon
epsilon=float(line.split()[-1])
line=f.readline() #kernel_cutoff
cutoff=float(line.split()[-1])
val_at_cutoff=np.exp(-0.5*cutoff**2)
line=f.readline() #compression_threshold
line=f.readline() #zed
Zed=float(line.split()[-1])
f.close()
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3,4,5])
center_x=np.array(data.iloc[:,0])
center_y=np.array(data.iloc[:,1])
sigma_x=np.array(data.iloc[:,2])
sigma_y=np.array(data.iloc[:,3])
height=np.array(data.iloc[:,4])
del data

#set grid
period_x=0
period_y=0
if len(args.grid_bin.split(','))!=2:
  sys.exit('two comma separated integers expected after --bin')
grid_bin_x=int(args.grid_bin.split(',')[0])
grid_bin_y=int(args.grid_bin.split(',')[1])
if args.grid_min is None:
  grid_min_x=min(center_x)
  grid_min_y=min(center_y)
else:
  if len(args.grid_min.split(','))!=2:
    sys.exit('two comma separated floats expected after --min')
  grid_min_x=float(args.grid_min.split(',')[0])
  grid_min_y=float(args.grid_min.split(',')[1])
if args.grid_max is None:
  grid_max_x=max(center_x)
  grid_max_y=max(center_y)
else:
  if len(args.grid_max.split(','))!=2:
    sys.exit('two comma separated floats expected after --max')
  grid_max_x=float(args.grid_max.split(',')[0])
  grid_max_y=float(args.grid_max.split(',')[1])
if args.angle1:
  if calc_der:
    print(' +++ WARNING: derivatives are not supported for periodic CVs +++')
    calc_der=False
  grid_min_x=-np.pi
  grid_max_x=np.pi
  period_x=2*np.pi
  grid_bin_x-=1
if args.angle2:
  if calc_der:
    print(' +++ WARNING: derivatives are not supported for periodic CVs +++')
    calc_der=False
  grid_min_y=-np.pi
  grid_max_y=np.pi
  period_y=2*np.pi
  grid_bin_y-=1
cv_grid_x=np.linspace(grid_min_x,grid_max_x,grid_bin_x)
cv_grid_y=np.linspace(grid_min_y,grid_max_y,grid_bin_y)
x,y=np.meshgrid(cv_grid_x,cv_grid_y)

#calculate
max_prob=0
prob=np.zeros((grid_bin_y,grid_bin_x))
if calc_der:
  der_prob_x=np.zeros((grid_bin_y,grid_bin_x))
  der_prob_y=np.zeros((grid_bin_y,grid_bin_x))
for i in range(grid_bin_y):
  print('    working... {:.0%}'.format(i/grid_bin_y),end='\r')
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
    arg2=dist_x**2+dist_y**2
    prob[i,j]=np.sum(height*(np.maximum(np.exp(-0.5*arg2)-val_at_cutoff,0)))
    prob[i,j]=prob[i,j]/Zed+epsilon
    if calc_der:
      der_prob_x[i,j]=np.sum(dist_x/sigma_x*height*(np.maximum(np.exp(-0.5*arg2)-val_at_cutoff,0)))/Zed
      der_prob_y[i,j]=np.sum(dist_y/sigma_y*height*(np.maximum(np.exp(-0.5*arg2)-val_at_cutoff,0)))/Zed
    if mintozero and prob[i,j]>max_prob:
      max_prob=prob[i,j]
if not mintozero:
  max_prob=1

#print out
output=open(args.outfile,'w')
print('#! FIELDS '+cv1name+' '+cv2name+' file.free der_'+cv1name+' der_'+cv2name,file=output)
print('#! SET min_'+cv1name+' '+str(grid_min_x),file=output)
print('#! SET max_'+cv1name+' '+str(grid_max_x),file=output)
print('#! SET nbins_'+cv1name+' '+str(grid_bin_x),file=output)
if period_x==0:
  print('#! SET periodic_'+cv1name+' false',file=output)
else:
  print('#! SET periodic_'+cv1name+' true',file=output)
print('#! SET min_'+cv2name+' '+str(grid_min_y),file=output)
print('#! SET max_'+cv2name+' '+str(grid_max_y),file=output)
print('#! SET nbins_'+cv2name+' '+str(grid_bin_y),file=output)
if period_y==0:
  print('#! SET periodic_'+cv2name+' false',file=output)
else:
  print('#! SET periodic_'+cv2name+' true',file=output)
for i in range(grid_bin_y):
  for j in range(grid_bin_x):
    if calc_der:
      print(x[i,j],y[i,j],0-kbt*np.log(prob[i,j]/max_prob),0-kbt/prob[i,j]*der_prob_x[i,j],0-kbt/prob[i,j]*der_prob_y[i,j],file=output)
    else:
      print(x[i,j],y[i,j],0-kbt*np.log(prob[i,j]/max_prob),file=output)
  print('',file=output)
output.close()
