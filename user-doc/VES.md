\page VES VES code

<!-- 
description: Module that implements enhanced sampling methods based on Variationally Enhanced Sampling
authors: Omar Valsson
reference: \cite Valsson-PRL-2014
-->

The VES code is a module for PLUMED that implements enhanced sampling methods
based on _Variationally Enhanced Sampling_ (VES) \cite Valsson-PRL-2014.
The VES code is developed by [Omar Valsson](http://www.valsson.info), 
see the [homepage of the VES code](http://www.ves-code.org) for further information.

The theory of VES is briefly explained \subpage ves_theory "here".

The VES code is an optional module that needes to be enabled when configuring the
compilation of PLUMED by using the '\-\-enable-modules=ves' 
(or '\-\-enable-modules=all') flag when running the 'configure' script. 

The various components of the VES code module are listed and described in the following sections 

- \subpage ves_tutorials
- \subpage ves_biases
- \subpage ves_basisf
- \subpage ves_targetdist
- \subpage ves_optimizer
- \subpage ves_utils
- \subpage ves_cltools


\page ves_theory Theory of VES


\par Variational Principle

In Variationally Enhanced Sampling \cite Valsson-PRL-2014 an external biasing potential \f$V(\mathbf{s})\f$ that acts in the space spanned by some set of collective variables (CVs) \f$\mathbf{s}=(s_1,s_2,\ldots,s_d)\f$ is constructed by minimizing the following convex functional
\f[
\Omega[V] = \frac{1}{\beta} \log
\frac{\int d\mathbf{s} \, e^{-\beta [F(\mathbf{s}+V(\mathbf{s}))]}}
{\int d\mathbf{s} \, e^{-\beta F(\mathbf{s})}}
+ \int d\mathbf{s} \, p(\mathbf{s}) V(\mathbf{s}),
\f]
where \f$F(\mathbf{s})\f$ is the free energy surface (FES) associated to the CVs at temperature \f$T\f$
(and \f$ \beta^{-1}=k_{\mathrm{B}}T \f$) and \f$p(\mathbf{s})\f$ is a so-called target distribution that is assumed to be normalized. It can be shown that the minimum of \f$\Omega[V]\f$ is given by
\f[
V(\mathbf{s}) = - F(\mathbf{s}) - \frac{1}{\beta} \log p(\mathbf{s}).
\f]
Under the influence of this minimum bias potential the biased equilibrium CV distribution \f$P_{V}(\mathbf{s})\f$ is equal to the target distribution \f$p(\mathbf{s})\f$,
\f[
P_{V}(\mathbf{s}) =
\frac{e^{-\beta\left[F(\mathbf{s})+V(\mathbf{s})\right]}}
{\int d\mathbf{s}\, e^{-\beta\left[F(\mathbf{s})+V(\mathbf{s})\right]}}
= p(\mathbf{s}).
\f]
The role of the target distribution \f$p(\mathbf{s})\f$ is therefore to determine
the sampling of the CVs that is achieved when minimizing \f$\Omega[V]\f$.

\par Minimization

The minimization is performed by assuming some given functional form for the bias potential \f$V(\mathbf{s};\boldsymbol{\alpha})\f$ that depends on some set of variational parameters \f$\boldsymbol{\alpha}=(\alpha_{1},\alpha_{2},\ldots,\alpha_{n})\f$. The convex function \f$\Omega(\boldsymbol{\alpha})=\Omega[V(\boldsymbol{\alpha})]\f$ is then minimized through
a gradient-based optimization technique.
The elements of the gradient \f$\nabla\Omega(\boldsymbol{\alpha})\f$ are defined as
\f[
\frac{\partial \Omega(\boldsymbol{\alpha})}{\partial \alpha_{i}} =
-<\frac{\partial V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{i}}>_{V\boldsymbol{\alpha})}
+<\frac{\partial V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{i}}>_{p},
\f]
where the first term is an average obtained in the biased simulation under the influence of the bias \f$V(\mathbf{s};\boldsymbol{\alpha})\f$ and the second term is an average over the target distribution \f$p(\mathbf{s})\f$.
Similarly the elements of the Hessian \f$H(\boldsymbol{\alpha})\f$ are defined as
\f[
\frac{\partial^2 \Omega(\boldsymbol{\alpha})}{\partial \alpha_{i} \partial \alpha_{j}} =
-<\frac{\partial^2 V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{i} \partial \alpha_{j}}>_{V(\boldsymbol{\alpha})}
+<\frac{\partial^2 V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{i} \partial \alpha_{j}}>_{p}
+ \, \beta \, \mathrm{Cov}\left[
\frac{\partial V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{i}},
\frac{\partial V(\mathbf{s};\boldsymbol{\alpha})}{\partial \alpha_{j}}
\right]_{V\boldsymbol{\alpha})}
\f]
where the covariance
\f$\mathrm{Cov} \left[a,b\right]_{V(\boldsymbol{\alpha})}=<a\,b>_{V(\boldsymbol{\alpha})} -
<a>_{V(\boldsymbol{\alpha})} <b>_{V(\boldsymbol{\alpha})} \f$
is obtained in a biased simulation under the influence of the bias \f$V(\mathbf{s};\boldsymbol{\alpha})\f$.

The gradient (and the Hessian) are inherently noisy as they need to be sampled in a biased simulation.
Therefore it is generally better to employ stochastic optimization
methods to perform the minimization of \f$\Omega(\boldsymbol{\alpha})\f$.

\par Linear Expansion

Most general, and normally most convenient, is to consider a linear expansion in some set of basis functions,
\f[
V(\mathbf{s};\boldsymbol{\alpha})=\sum_{\mathbf{k}} f_{\mathbf{k}}(\mathbf{s}),
\f]
where the basis functions \f$f_{\mathbf{k}}(\mathbf{s})\f$ are for example Fourier series (i.e. plane waves), or a product of some orthogonal polynomials like Legendre or Chebyshev. For such a linear expansion the gradients simplifies to being averages of the basis functions,
\f[
\frac{\partial \Omega(\boldsymbol{\alpha})}{\partial \alpha_{\mathbf{k}}} =
-<f_{\mathbf{k}}(\mathbf{s})>_{V\boldsymbol{\alpha})}
+<f_{\mathbf{k}}(\mathbf{s})>_{p},
\f]
while the Hessian is just the covariance of the basis functions,
\f[
\frac{\partial^2 \Omega(\boldsymbol{\alpha})}{\partial \alpha_{\mathbf{k}} \partial \alpha_{\mathbf{l}}} =
 \beta \, \mathrm{Cov}\left[
f_{\mathbf{k}}(\mathbf{s}),
f_{\mathbf{l}}(\mathbf{s})
\right]_{V(\boldsymbol{\alpha})}.
\f]





\page ves_biases Bias

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

[Enhancing Important Fluctuations: Rare Events and Metadynamics from a Conceptual Viewpoint](https://doi.org/10.1146/annurev-physchem-040215-112229), Annu. Rev. Phys. Chem. 2016



Variationally Enhanced Sampling:

[Variational Approach to Enhanced Sampling and Free Energy Calculations](https://doi.org/10.1103/PhysRevLett.113.090601), Phys. Rev. Lett. 2014

[Variationally Optimized Free-Energy Flooding for Rate Calculation](https://doi.org/10.1103/PhysRevLett.115.070601), Phys. Rev. Lett. 2015



\par Tuesday February 14

\ref lugano-1 "Tutorial 1": Introduction to PLUMED and analyzing molecular simulations

\par Wednesday February 15

\ref ves-lugano2017-metad "Tutorial 2": Biasing with metadynamics

\ref ves-lugano2017-ves1 "Tutorial 3": Biasing with variationally enhanced sampling

\par Thursday February 16

\ref ves-lugano2017-ves2 "Tutorial 4": Further on variationally enhanced sampling

Tutorial 5: Advanced collective variables
- \ref lugano-2 "Path CVs"
- \ref belfast-10 "Multicolvar"
- \ref belfast-3 "Dimensionality reduction"

\par Friday February 17

\ref ves-lugano2017-kinetics "Tutorial 6": Obtaining kinetics from molecular simulations
