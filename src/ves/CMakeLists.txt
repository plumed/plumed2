#automatically generated CMakeLists.txt, may not work!
#the variable module_name is set up as a sugar to reduce "copy-paste" errors
set (module_name "ves")
#Note that the macros here require this directory added as a subdir of plumed/src
option(module_${module_name} "activate module ${module_name}" OFF)
ADDMODULETOKERNEL(${module_name}
BasisFunctions.cpp
BF_Chebyshev.cpp
BF_Combined.cpp
BF_Cosine.cpp
BF_CubicBsplines.cpp
BF_Custom.cpp
BF_Fourier.cpp
BF_Gaussians.cpp
BF_Legendre.cpp
BF_Powers.cpp
BF_Sine.cpp
BF_Wavelets.cpp
CoeffsBase.cpp
CoeffsMatrix.cpp
CoeffsVector.cpp
FermiSwitchingFunction.cpp
GridIntegrationWeights.cpp
GridLinearInterpolation.cpp
GridProjWeights.cpp
LinearBasisSetExpansion.cpp
MD_LinearExpansionPES.cpp
Opt_Adam.cpp
Opt_BachAveragedSGD.cpp
Opt_Dummy.cpp
Optimizer.cpp
Opt_RobbinsMonroSGD.cpp
OutputBasisFunctions.cpp
OutputFesBias.cpp
OutputTargetDistribution.cpp
TargetDistModifer.cpp
TargetDistribution.cpp
TD_Chi.cpp
TD_ChiSquared.cpp
TD_Custom.cpp
TD_Exponential.cpp
TD_ExponentiallyModifiedGaussian.cpp
TD_Gaussian.cpp
TD_GeneralizedExtremeValue.cpp
TD_GeneralizedNormal.cpp
TD_Grid.cpp
TD_LinearCombination.cpp
TD_Multicanonical.cpp
TD_MultithermalMultibaric.cpp
TD_ProductCombination.cpp
TD_ProductDistribution.cpp
TD_Uniform.cpp
TD_VonMises.cpp
TD_WellTempered.cpp
VesBias.cpp
VesDeltaF.cpp
VesLinearExpansion.cpp
VesTools.cpp
WaveletCoeffs.cpp
WaveletGrid.cpp
)
ADDMODULEDEPENDENCIES(${module_name}
bias cltools colvar config core tools lepton
)
