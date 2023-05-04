#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "core")
#Note that the macros here require this directory added as a subdir of plumed/src
set(module_${module_name} ON CACHE INTERNAL "always active module ${module_name}")
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
#file generated during the build process
PlumedMainMap.inc PlumedMainEnum.inc
GREXMap.inc GREXEnum.inc
CLToolMainMap.inc CLToolMainEnum.inc
)
ADDMODULEDEPENDENCIES(${module_name}
config tools lepton
)

#execute_process(
#        COMMAND bash ${PLUMED_MAKETOOLS}/makecmd_cmake map "${CMAKE_CURRENT_SOURCE_DIR}/PlumedMain.cpp"  PlumedMainMap.inc
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)

MACRO(CREATEINC cppFile type outputFile)
    add_custom_command(OUTPUT ${outputFile}
        DEPENDS ${cppFile}
        COMMAND bash ${PLUMED_MAKETOOLS}/makecmd_cmake map "${CMAKE_CURRENT_SOURCE_DIR}/${cppFile}"  ${outputFile}
        COMMENT "Generating ${outputFile}"
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDMACRO()

CREATEINC(PlumedMain.cpp map PlumedMainMap.inc)
CREATEINC(PlumedMain.cpp enum PlumedMainEnum.inc)
CREATEINC(GREX.cpp map GREXMap.inc)
CREATEINC(GREX.cpp enum GREXEnum.inc)
CREATEINC(CLToolMain.cpp map CLToolMainMap.inc)
CREATEINC(CLToolMain.cpp enum CLToolMainEnum.inc)



