Membrane fusion process, when two separate membranes merge, is crucial in life. The fusion of lipid bilayers follows a series of discrete steps with two relevant intermediates: hemifusion structures and fusion pores. The hemifusion structures mix lipids from the involved membranes without cargo exchange, 
while the fusion pore require an aqueous channel to connect the contents.

To study the hemifusion stage computationally, Hub and Awasthi developed a CV that initially nucleated hydrophilic pores in lipid bilayers (see second paper cited below) and later extended it to induce the hemifusion stalks (see fourth paper cited below). Di Bartolo and Masone implemented that CV in PLUMED (see first paper cited below).

Then, to nucleate and expand the fusion pore, based on Hub's work  in single lipid bilayers (see third paper cited below), Di Bartolo and Masone implemented two others CVs to nucleate and expand fusion pores.

## Installation 

This module is not installed by default. Add `--enable-modules=membranefusion` to your './configure' command when building PLUMED to enable these features.

## Usage

This module contains three CVs to:

- Induce a hemifusion stalk: [MEMFUSIONP](MEMFUSIONP.md)
- Nucleate a fusion pore: [FUSIONPORENUCLEATIONP](FUSIONPORENUCLEATIONP.md)
- Expand a fusion pore: [FUSIONPOREEXPANSIONP](FUSIONPOREEXPANSIONP.md)

