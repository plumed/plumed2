import numpy as np

cmap_data = np.loadtxt("contact_matrix.reference")
cmap2_data = np.loadtxt("square_matrix.reference")
eig = np.loadtxt("colvar.reference")

for row in range(cmap_data.shape[0]) : 
    matrix = cmap_data[row,1:].reshape([7,7])
    matrix2 = cmap2_data[row,1:].reshape([7,7])
    tmat2 = np.linalg.matrix_power( matrix, 2 )
    diffmat = np.fabs( tmat2 - matrix2 ) 
    if not (diffmat < 1E-3 ).all() : raise Exception("Differences in matrices") 
    eigval, eigvec = np.linalg.eig( matrix2 )
    print( max(eigval), eig[row] )
