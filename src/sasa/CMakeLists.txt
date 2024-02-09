#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "sasa")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
Sasa.cpp
sasa_HASEL.cpp
sasa_LCPO.cpp
)
ADDMODULENEEDS(${module_name}
	core tools config  
)
ADDMODULEDEPENDENCIES(${module_name}
	core #tools #config  
)
