#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "cltools")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module cltools" ON)
ADDMODULETOKERNEL(${module_name}
CLTool.cpp
Completion.cpp
Driver.cpp
DriverDouble.cpp
DriverFloat.cpp
GenExample.cpp
GenJson.cpp
GenTemplate.cpp
Info.cpp
kT.cpp
Manual.cpp
PdbRenumber.cpp
pesmd.cpp
SimpleMD.cpp
SumHills.cpp
)
ADDMODULEDEPENDENCIES(${module_name}
core config tools molfile xdrfile
)
