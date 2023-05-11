#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "analysis")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
AnalysisBase.cpp
Average.cpp
AverageVessel.cpp
Committor.cpp
DataCollectionObject.cpp
EuclideanDissimilarityMatrix.cpp
FarthestPointSampling.cpp
Histogram.cpp
LandmarkSelectionBase.cpp
LandmarkStaged.cpp
OutputColvarFile.cpp
OutputPDBFile.cpp
PrintDissimilarityMatrix.cpp
ReadAnalysisFrames.cpp
ReadDissimilarityMatrix.cpp
ReselectLandmarks.cpp
SelectRandomFrames.cpp
SelectWithStride.cpp
WhamHistogram.cpp
WhamWeights.cpp
)
ADDMODULENEEDS(${module_name}
	core tools reference vesselbase gridtools multicolvar bias
)
ADDMODULEDEPENDENCIES(${module_name}
	reference #core
	#TODO: resolve THIS circular dependency:
	#vesselbase gridtools multicolvar# it does actually depend on them, but it is better to avid cyclic dependencies
	bias #tools core
)
