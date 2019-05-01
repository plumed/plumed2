#!/usr/bin/python
from math import *
from itertools import ifilterfalse

Temperature=1700   # Temperature in Kelvins
dtmd=0.1           # timestep in COLVAR file (ps)
CVcol=4            # Column of CV (index starting from 1)
vesbiascol=5       # Column of VES bias (index starting from 1)
metabiascol=7      # Column of metadynamics bias (index starting from 1)
fname='COLVAR'     # name of COLVAR file
fout='COLVAR-RW'   # name of output COLVAR

CVcol=CVcol-1              # python index starts from 0
vesbiascol=vesbiascol-1    # python index starts from 0
metabiascol=metabiascol-1  # python index starts from 0

# for skipping headers
def iscomment(s):
    return s.startswith('#')

# Conversion from energy in eV to kBT
kT=(2.479/298.0/96.485)*Temperature
beta=1.0/kT

conversion=1e-12   # conversion factor ps ---> s

f=open(fname,'r')
fout=open(fout,'w')

timesum=0.0
for line in ifilterfalse(iscomment,f):
   line=line.strip()
   columns=line.split()
   cv=float(columns[CVcol])
   bias=float(columns[vesbiascol])+float(columns[metabiascol])
   timesum=timesum+exp(beta*bias)*dtmd
   fout.write(str(timesum*conversion)+' '+str(cv)+'\n')


print 'Unbiased first passage time (s)',timesum*conversion

f.close()
fout.close()
