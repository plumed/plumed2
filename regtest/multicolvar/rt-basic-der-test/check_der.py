import numpy as np

d1c_raw = np.loadtxt("deriv_d1c.reference")[:,2]
c1 = np.loadtxt("deriv_c1.reference")[:,2]

dbase, cbase = 0, 0
for i in range(5) :
    # Atom 1
    for j in range(3) : 
        totder = d1c_raw[dbase+j] + d1c_raw[dbase+6+j] + d1c_raw[dbase+12+j] + d1c_raw[dbase+18+j]
        if not np.isclose(totder,c1[cbase+j],rtol=1e-3) : print("Found difference in derivatives on atom 1", totder, c1[cbase+j])
    # Atom 2
    twopos = dbase + 6*4  
    for j in range(3) : 
        totder = d1c_raw[twopos+j] + d1c_raw[twopos+6+j] + d1c_raw[twopos+12+j] + d1c_raw[twopos+18+j]
        if not np.isclose(totder,c1[cbase+3+j],rtol=1e-3) : print("Found difference in derivatives on atom 2") 
    # Other atoms
    for k in range(4) : 
        astart=dbase+k*6
        for j in range(3) :
            totder = d1c_raw[astart+3+j] + d1c_raw[astart+4*6+3+j]
            if not np.isclose(totder,c1[cbase+6+3*k+j],atol=1e-3 ) : print("Found difference in derivatives on atom", 3+k, totder,c1[cbase+6+3*k+j] )   

    # Check the virial
    for j in range(9) : 
        if not np.isclose(d1c_raw[dbase+8*6+j],c1[cbase+6*3+j]) : print("Found difference in virial")
    dbase, cbase = dbase+57, cbase+27
    
