import numpy as np

# In reality one should not log stuff here...
log=open("log.txt","w",1)

print("At import",file=log)



# Define the distance function
def dist(x):
    print("In dist()",file=log)
    print(x,file=log)
    r = x[0,:]-x[1,:]
    d2 = np.dot(r,r)
    return np.sqrt(d2)

def grad_dist(x):
    print("In grad_dist()",file=log)
    d = dist(x)
    r = x[0,:]-x[1,:]
    g = r/d
    return np.array([g,-g])

# The CV function actually called
def cv(x):
    return dist(x), grad_dist(x)

