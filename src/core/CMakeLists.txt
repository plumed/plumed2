#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "core")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module core")
ADDMODULETOKERNEL(${module_name}
ActionAnyorder.cpp
ActionAtomistic.cpp
Action.cpp
ActionPilot.cpp
ActionRegister.cpp
ActionSet.cpp
ActionSetup.cpp
ActionShortcut.cpp
ActionWithArguments.cpp
ActionWithValue.cpp
ActionWithVirtualAtom.cpp
Atoms.cpp
CLTool.cpp
CLToolMain.cpp
CLToolRegister.cpp
Colvar.cpp
DataFetchingObject.cpp
ExchangePatterns.cpp
FlexibleBin.cpp
GenericMolInfo.cpp
GREX.cpp
MDAtoms.cpp
PlumedMain.cpp
PlumedMainInitializer.cpp
TargetDist.cpp
Value.cpp
WithCmd.cpp
)
ADDMODULEDEPENDENCIES(${module_name}
config tools lepton
)
