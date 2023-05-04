#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "multicolvar")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
ActionVolume.cpp
AlphaBeta.cpp
Angles.cpp
AtomValuePack.cpp
Bridge.cpp
BridgedMultiColvarFunction.cpp
CatomPack.cpp
CenterOfMultiColvar.cpp
CoordinationNumbers.cpp
Density.cpp
DihedralCorrelation.cpp
DistanceFromContour.cpp
Distances.cpp
DumpMultiColvar.cpp
FilterBetween.cpp
FilterLessThan.cpp
FilterMoreThan.cpp
InPlaneDistances.cpp
LocalAverage.cpp
MultiColvarBase.cpp
MultiColvarCombine.cpp
MultiColvarDensity.cpp
MultiColvarFilter.cpp
MultiColvarProduct.cpp
NumberOfLinks.cpp
Torsions.cpp
VolumeAround.cpp
VolumeBetweenContours.cpp
VolumeCavity.cpp
VolumeGradientBase.cpp
VolumeInCylinder.cpp
VolumeInSphere.cpp
VolumeTetrapore.cpp
XAngle.cpp
XDistances.cpp
XYDistances.cpp
XYTorsion.cpp
)
ADDMODULENEEDS(${module_name}
core tools vesselbase gridtools
)
