#!/usr/bin/python
from math import *
from itertools import ifilterfalse

Temperature=1700   # Temperature in Kelvins
dtmd=0.1           # timestep in COLVAR file (ps)
vesbiascol=5       # Column of VES bias (index starting from 1)
metabiascol=7      # Column of metadynamics bias (index starting from 1)
fname='COLVAR'     # name of COLVAR file

vesbiascol=vesbiascol-1    # python index starts from 0
metabiascol=metabiascol-1  # python index starts from 0

ftimes=open('fpt.dat','w')

def iscomment(s):
    return s.startswith('#')

kT=(2.479/298.0/96.485)*Temperature
beta=1.0/kT

for i in range(30):
   string=fname+'-'+str(i+1)
   f=open(string,'r')

   timesum=0.0
   for line in ifilterfalse(iscomment,f):
      line=line.strip()
      columns=line.split()
      bias=float(columns[vesbiascol])+float(columns[metabiascol])
      timesum=timesum+exp(beta*bias)*dtmd

   conversion=(1e-12)
   ftimes.write(str(timesum*conversion)+'\n')
   print i
   f.close()

ftimes.close()
