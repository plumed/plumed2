\page VES Variationally Enhanced Sampling (VES code) 

<!-- 
description: Module that implements enhanced sampling methods based on Variationally Enhanced Sampling
authors: Omar Valsson
reference: \cite Valsson-PRL-2014
-->

The VES code is a module for PLUMED that implements enhanced sampling methods
based on _Variationally Enhanced Sampling_ (VES) \cite Valsson-PRL-2014.
The VES code is developed by [Omar Valsson](http://www.valsson.info), 
see the [homepage of the VES code](http://www.ves-code.org) for further information.

The VES code is an optional module that needs to be enabled when configuring the
compilation of PLUMED by using the '\-\-enable-modules=ves' 
(or '\-\-enable-modules=all') flag when running the 'configure' script. 

In the \ref ves_tutorials "tutorials" you can learn how to use the methods 
implemented in the VES code.

The various components of the VES code module are listed and described in the following sections 

- \subpage ves_biases
- \subpage ves_basisf
- \subpage ves_targetdist
- \subpage ves_optimizer
- \subpage ves_utils
- \subpage ves_cltools
- \subpage ves_tutorials


\page ves_biases Biases

The following list contains the biases available in the VES code.

@VES_BIAS@

\page ves_basisf Basis functions

The following list contains the one-dimensional basis functions available in the VES code.

@VES_BASISF@


\page ves_targetdist Target Distributions

The following list contains the target distributions available in the VES code.

@VES_TARGETDIST@


\page ves_optimizer Optimizers

The following list contains the optimizers available in the VES code.

@VES_OPTIMIZER@


\page ves_utils Utilities

The following list contains various utilities available in the VES code. 

@VES_UTILS@



\page ves_cltools Command Line Tools

The following list contains the command line tools available in the VES code.

@VES_TOOLS@



\page ves_tutorials Tutorials

The following tutorials are available for the VES code. 

\subpage ves_tutorial_lugano_2017

@VES_TUTORIALS@




\page ves_tutorial_lugano_2017 MARVEL-VES School February 2017

\image html ves-lugano2017-logo.png  width=800px

Tutorials from the [MARVEL School on Variationally Enhanced Sampling]
(https://sites.google.com/site/vesschool2017/home) that was held in
Lugano, February 14-17, 2017.

\par Suggested readings

Metadynamics:

[Enhancing Important Fluctuations: Rare Events and Metadynamics from a Conceptual Viewpoint](https://doi.org/10.1146/annurev-physchem-040215-112229), Annual Reviews in Physical Chemistry 2016



Variationally Enhanced Sampling:

[Variational Approach to Enhanced Sampling and Free Energy Calculations](https://doi.org/10.1103/PhysRevLett.113.090601), Physical Review Letters 2014

[Variationally Optimized Free-Energy Flooding for Rate Calculation](https://doi.org/10.1103/PhysRevLett.115.070601), Physical Review Letters 2015



\par Tuesday February 14

\ref marvel-1 "Tutorial 1": Introduction to PLUMED and analyzing molecular simulations

\par Wednesday February 15

\ref ves-lugano2017-metad "Tutorial 2": Biasing with metadynamics

\ref ves-lugano2017-ves1 "Tutorial 3": Biasing with variationally enhanced sampling

\par Thursday February 16

\ref ves-lugano2017-ves2 "Tutorial 4": Further on variationally enhanced sampling

Tutorial 5: Advanced collective variables
- \ref marvel-2 "Path CVs"
- \ref belfast-10 "Multicolvar"
- \ref belfast-3 "Dimensionality reduction"

\par Friday February 17

\ref ves-lugano2017-kinetics "Tutorial 6": Obtaining kinetics from molecular simulations
