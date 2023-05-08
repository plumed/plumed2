#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "function")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" ON)
ADDMODULETOKERNEL(${module_name}
Combine.cpp
Custom.cpp
Ensemble.cpp
FuncPathGeneral.cpp
FuncPathMSD.cpp
FuncSumHills.cpp
Function.cpp
LocalEnsemble.cpp
Piecewise.cpp
Sort.cpp
Stats.cpp
Target.cpp
)
ADDMODULENEEDS(${module_name}
	core reference tools lepton
)
ADDMODULEDEPENDENCIES(${module_name}
	reference #tools core
	#tools #lepton
)
