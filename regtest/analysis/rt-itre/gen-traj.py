import numpy.random as rnd

for i in range(50) :
    vec = rnd.multivariate_normal( [0.95527742, -0.0073749 ], [[0.01089105, 0.00065146],[0.00065146, 0.0261473 ]] )
    print( 1 )
    print("100. 100. 100.")
    print( "X", vec[0], vec[1], 0.0 )

for i in range(50) :
    vec = rnd.multivariate_normal( [ 0.00674389, -0.95262677], [[0.02834656, 0.00011025],[0.00011025, 0.01101571]]  )
    print( 1 )
    print("100. 100. 100.")
    print( "X", vec[0], vec[1], 0.0 )

for i in range(50) :
    vec = rnd.multivariate_normal( [ 0.03124637, -0.02893218], [[0.02697116, 0.01789088],[0.01789088, 0.02588994]]  )
    print( 1 )
    print("100. 100. 100.")
    print( "X", vec[0], vec[1], 0.0 )
