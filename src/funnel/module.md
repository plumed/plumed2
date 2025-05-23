FM is a combination of Metadynamics bias potential with a funnel-shape restraint potential applied to the target structure of a binding interaction. 
The latter is composed of a cone restraint, which covers the ligand binding site, and a cylindrical one that heads towards the solvent. 
When inside the funnel volume, the ligand does not feel any restraint potential, proceeding as regular Metadynamics.
Upon reaching the boundaries of the funnel, a repulsive bias is applied forcing the ligand to remain in the allowed funnel space. 
The result is an acceleration in the sampling of the binding/unbinding process, leading to a swift convergence of the calculation and a well-defined binding free-energy surface.

## Installation 

This module is not installed by default. Add `--enable-modules=funnel` to your './configure' command when building PLUMED to enable these features.

## Usage

This module is a direct evolution of the original FM that is discussed in the first paper cited below since it incorporates an alignment function that removes the necessity to block the target macromolecule in the simulation box.

The user can follow a comprehensive protocol that is provided in the second paper cited below, which will help in all stages of the simulation, including pre- and post-processing.
An example of an input file can be found on <a href="https://www.plumed-nest.org/eggs/19/039/">FUNNEL-NEST's webpage</a>

## Module Contents

The funnel module is composed of a collective variable that calculates the position of a ligand with respect to a line and a potential that creates a funnel-shape restraint centered on the line ([FUNNEL_PS](FUNNEL_PS.md) and [FUNNEL](FUNNEL.md), respectively).
