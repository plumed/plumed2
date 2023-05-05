# set(module_name wrapper)
# set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")
add_library(PlumedStaticWrapper OBJECT PlumedStatic.cpp)
add_library(PlumedWrapper OBJECT Plumed.c)

#add here the lib so and the lib static for Plumed.so and PlumedKernel.so
#find a way to add all the activated libraries here or add them to tha main CMakeList.txt?
#Then cmake should short all the mess for you, hopefully