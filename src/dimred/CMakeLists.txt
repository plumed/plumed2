#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "dimred")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
ClassicalMultiDimensionalScaling.cpp
DimensionalityReductionBase.cpp
OutputPCAProjections.cpp
PCA.cpp
ProjectNonLandmarkPoints.cpp
SketchMapBase.cpp
SketchMapConjGrad.cpp
SketchMap.cpp
SketchMapPointwise.cpp
SketchMapRead.cpp
SketchMapSmacof.cpp
SMACOF.cpp
SmacoffMDS.cpp
)
ADDMODULENEEDS(${module_name}
core tools reference gridtools analysis
)
