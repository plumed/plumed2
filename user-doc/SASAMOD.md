\page SASAMOD SASA collective variable

<!-- 
description: Solvent Accessible Surface Area collective variable (SASA)
authors: Andrea Arsiccio
reference: \cite Hasel1988 \cite Weiser1999 \cite Arsiccio-C-SASA-2022 \cite Arsiccio-T-SASA-2021 \cite Arsiccio-P-SASA-2021
-->

\section Overview

This SASA module contains methods for the calculation of the solvent accessible surface area (SASA) of proteins using either the fast algorithm by Hasel et al. \cite Hasel1988 or the LCPO algorithm \cite Weiser1999. This module can be used to include the SASA as a collective variable in metadynamics simulations, and also for implicit solvent simulations as described in \cite Arsiccio-C-SASA-2022, \cite Arsiccio-T-SASA-2021 and \cite Arsiccio-P-SASA-2021.

\section Installation 
This module is not installed by default. Add '\-\-enable-modules=sasa' to your './configure' command when building PLUMED to enable these features.

\section Usage
Currently, all features of the SASA module are included in two SASA functions: \ref SASA_HASEL \ref SASA_LCPO

\section Module Contents
- \subpage SASAMODColvar

\page SASAMODColvar CVs Documentation

The following list contains descriptions of biases developed for the PLUMED-SASA module. They can be used in combination with other biases outside of the SASA module.

@SASAMOD_COLVAR@
