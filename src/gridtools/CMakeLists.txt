#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "gridtools")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
ActionWithGrid.cpp
ActionWithInputGrid.cpp
ActionWithIntegral.cpp
AverageOnGrid.cpp
ContourFindingBase.cpp
ConvertToFES.cpp
DumpCube.cpp
DumpGrid.cpp
FindContour.cpp
FindContourSurface.cpp
FindSphericalContour.cpp
FourierTransform.cpp
GridPrintingBase.cpp
GridSearch.cpp
GridToXYZ.cpp
GridVessel.cpp
HistogramOnGrid.cpp
IntegrateGrid.cpp
InterpolateGrid.cpp
)
ADDMODULENEEDS(${module_name}
	core tools vesselbase
)
ADDMODULEDEPENDENCIES(${module_name}
	core tools vesselbase
)
