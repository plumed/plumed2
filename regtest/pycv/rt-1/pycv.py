
import numpy as np
log=open("pycv.log","w")

print("Imported.",file=log)

def cv1(x):
    print(x,file=log)
    d = x[0,:]-x[1,:]
    d2 = np.dot(d,d)
    return np.sqrt(d2)
