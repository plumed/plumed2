import numpy as np
import plumedCommunications

# In reality one should not log stuff here...
log = open("log.txt", "w", 1)

print("At import", file=log)

plumedInit={"Value": plumedCommunications.defaults.COMPONENT_NODEV,}

# The CV function actually called
def cv(X):
    #0.36 is arbitrary
    return 0.36*(X.getStep()+1)
