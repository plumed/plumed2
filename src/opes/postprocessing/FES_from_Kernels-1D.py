#! /usr/bin/env python3

### Get the running FES estimate used by OPES_WT, 1D only ###
# similar to plumed sum_hills

import sys
import numpy as np
import pandas as pd
import subprocess
import argparse

# Manually change if using RECURSTIVE_MERGE_OFF
recursive_merge=True

#parser
parser = argparse.ArgumentParser(description='get the running FES estimate used by OPES, 1D only')
parser.add_argument('--kernels',dest='filename',type=str,default='KERNELS',help='the kernels file name')
parser.add_argument('--kt',dest='kbt',type=float,required=True,help='the temperature in energy units')
parser.add_argument('--angle',dest='angle',action='store_true',default=False,help='the cv is an angle in the range [-pi,pi]')
parser.add_argument('--min',dest='grid_min',type=float,required=False,help='lower bound for the grid')
parser.add_argument('--max',dest='grid_max',type=float,required=False,help='upper bound for the grid')
parser.add_argument('--bin',dest='grid_bin',type=int,default=100,help='number of bins for the grid')
parser.add_argument('--stride',dest='stride',type=int,default=0,help='how often to print running fes')
parser.add_argument('--mintozero',dest='mintozero',action='store_true',default=False,help='shift the minimum to zero')
parser.add_argument('--no_der',dest='no_der',action='store_true',default=False,help='skip derivative to run faster')
parser.add_argument('--faster',dest='faster',action='store_true',default=False,help='run faster, ingoring Zed and the epsilon cutoff')
parser.add_argument('--outfile',dest='outfile',type=str,default='fes.dat',help='name of the output file')
args = parser.parse_args()
#parsing
filename=args.filename
kbt=args.kbt
mintozero=args.mintozero
faster=args.faster
calc_der=(not args.no_der)

#get kernels
f=open(filename)
line=f.readline() #header
if len(line.split())!=7:
  sys.exit('  something is wrong with file '+filename)
cvname=line.split()[3]
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
line=f.readline() #cutoff
cutoff=float(line.split()[-1])
line=f.readline() #threshold
threshold=float(line.split()[-1])
for n in range(3):
  line=f.readline() #skip CV info if present
line2=f.readline() #second line
pace_to_time=(float(line2.split()[0])-float(line.split()[0]))
f.close()
cutoff2=cutoff**2
val_at_cutoff=np.exp(-0.5*cutoff2)
data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,usecols=[1,2,3])
center=np.array(data.iloc[:,0])
sigma=np.array(data.iloc[:,1])
height=np.array(data.iloc[:,2])
del data
print('  all data loaded')

#set grid
if args.grid_min is None:
  grid_min=min(center)
else:
  grid_min=args.grid_min
if args.grid_max is None:
  grid_max=max(center)
else:
  grid_max=args.grid_max
grid_bin=args.grid_bin+1
period=0
if args.angle:
  if (args.grid_min is not None or args.grid_max is not None):
    sys.exit('do not set min and max if variable is an angle')
  grid_min=-np.pi
  grid_max=np.pi
  period=2*np.pi
  grid_bin-=1
  if calc_der:
    print(' +++ WARNING: derivative is not supported for periodic CV +++')
    calc_der=False
cv_grid=np.linspace(grid_min,grid_max,grid_bin)

#output files
head='#! FIELDS '+cvname+' file.free der_'+cvname
head+='\n#! SET min_'+cvname+' '+str(grid_min)
head+='\n#! SET max_'+cvname+' '+str(grid_max)
head+='\n#! SET nbins_'+cvname+' '+str(grid_bin)
if period==0:
  head+='\n#! SET periodic_'+cvname+' false'
else:
  head+='\n#! SET periodic_'+cvname+' true'
print_stride=args.stride
if print_stride==0:
  print_stride=len(center)+1
outfile=args.outfile
if print_stride<=len(center):
  file_ext=outfile.split('.')[-1]
  if len(file_ext)>1:
    file_ext='.'+file_ext
  current_fes_running='fes_running/'+outfile[:-len(file_ext)]+'.t-%d'+file_ext
  cmd=subprocess.Popen('mkdir fes_running',shell=True)
  cmd.wait()

# useful functions
z_center=[center[0]]
z_sigma=[sigma[0]]
z_height=[height[0]]

def get_merge_candidate(c,self):
  min_dist=threshold
  min_j=-1
  for j in range(len(z_center)):
    if j==self:
      continue
    dist=abs(c-z_center[j])/z_sigma[j] #merging through periodicity is not implemented in the code
    if dist<min_dist:
      min_j=j
  return min_j

def merge(j,m_height,m_center,m_sigma):
    h=z_height[j]+m_height
    c=(z_height[j]*z_center[j]+m_height*m_center)/h
    s2=(z_height[j]*(z_sigma[j]**2+z_center[j]**2)+m_height*(m_sigma**2+m_center**2))/h-c**2
    z_height[j]=h
    z_center[j]=c
    z_sigma[j]=np.sqrt(s2)

def build_fes(c,s,h):
  nker=len(c)
  prob=np.zeros(grid_bin)
  der_prob=np.zeros(grid_bin)
  for x in range(len(cv_grid)):
    if period==0:
      dist=(cv_grid[x]-c)/s
    else:
      dx=np.abs(cv_grid[x]-c)
      dist=np.minimum(dx,period-dx)/s
    prob[x]=np.sum(h*np.maximum(np.exp(-0.5*dist**2)-val_at_cutoff,0))
    if calc_der:
      der_prob[x]=np.sum(dist/s*h*np.maximum(np.exp(-0.5*dist**2)-val_at_cutoff,0))
  if not faster:
    Zed=0
    for j in range(nker):
      if period==0:
        dist=(c[j]-c)/s
      else:
        dx=np.abs(c[j]-c)
        dist=np.minimum(dx,period-dx)/s
      Zed+=np.sum(h*np.maximum(np.exp(-0.5*dist**2)-val_at_cutoff,0))
    Zed/=nker
    prob=prob/Zed+epsilon
    der_prob/=Zed
  norm=1
  if mintozero:
    norm=max(prob)
  return 0-kbt*np.log(prob/norm),0-kbt/prob*der_prob #adding zero avoids printing -0

# compression
for i in range(1,len(center)):
  print('    working... {:.0%}'.format(i/len(center)),end='\r')
  j=get_merge_candidate(center[i],-1)
  if j>=0:
    merge(j,height[i],center[i],sigma[i])
    if recursive_merge:
      jj=get_merge_candidate(z_center[j],j)
      while jj>=0:
        merge(jj,z_height[j],z_center[j],z_sigma[j])
        z_height.pop(j)
        z_center.pop(j)
        z_sigma.pop(j)
        if j<jj:
          j=jj-1
        else:
          j=jj
        jj=get_merge_candidate(z_center[j],j)
  else:
    z_center.append(center[i])
    z_sigma.append(sigma[i])
    z_height.append(height[i])
  if (i+1)%print_stride==0:
    fes,der_fes=build_fes(z_center,z_sigma,z_height)
    np.savetxt(current_fes_running%((i+1)*pace_to_time),np.c_[cv_grid,fes,der_fes],header=head,comments='',fmt='%14.9f')
print(' total kernels read from file: %d'%len(center))
print(' total kernels in compressed FES: %d'%len(z_center))

#build and print final
fes,der_fes=build_fes(z_center,z_sigma,z_height)
np.savetxt(outfile,np.c_[cv_grid,fes,der_fes],header=head,comments='',fmt='%14.9f')

