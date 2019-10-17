import numpy as np

# In reality one should not log stuff here...
log=open("log.txt","w",1)

print("At import",file=log)


# Define the distance function
def dist(x):
    r = x[0,:]-x[1,:]
    d2 = np.dot(r,r)
    return np.sqrt(d2)

def grad_dist(x):
    d = dist(x)
    r = x[0,:]-x[1,:]
    g = r/d
    return np.array([g,-g])


# The CV function actually called
def cv(X):
    print("At CV",file=log)
    row12 = np.array([0,1])
    row13 = np.array([0,2])
    
    X12 = X[row12,:]
    X13 = X[row13,:]

    rf = {
        'd12': dist(X12),
        'd13': dist(X13)
        }

    gd12 = grad_dist(X12)
    gd13 = grad_dist(X13)

    rg = {
        'd12': np.array([ gd12[0,:], gd12[1,:], [0.,0.,0.] ]),
        'd13': np.array([ gd13[0,:], [0.,0.,0.], gd13[1,:] ])
        }
    
    return rf, rg

