set(module_name wrapper)
#this is not a module, this is set up to make easier the dependency interface
set(module_${module_name} ON CACHE INTERNAL "always active core module ${module_name}")
add_library(PlumedStaticWrapper OBJECT PlumedStatic.cpp)
add_library(PlumedDynamicWrapper OBJECT Plumed.c)
foreach(lib PlumedStaticWrapper PlumedDynamicWrapper)
    target_link_libraries(${lib} PUBLIC config)
endforeach(lib PlumedStaticWrapper PlumedDynamicWrapper)


#add here the lib so and the lib static for Plumed.so and PlumedKernel.so
#find a way to add all the activated libraries here or add them to tha main CMakeList.txt?
#Then cmake should short all the mess for you, hopefully