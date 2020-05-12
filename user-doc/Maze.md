\page MAZE MAZE

<!-- 
description: Module that implements enhanced sampling methods for ligand 
unbinding from protein tunnels.
authors: Jakub Rydzewski
reference: \cite RydzewskiMaze 
-->

maze is a module for PLUMED2, which implements enhanced sampling methods for 
ligand unbinding from protein tunnels. The maze module is developed and 
maintained by [Jakub Rydzewski](http://www.fizyka.umk.pl/~jr) at the Institute 
of Physics, Nicolaus Copernicus University, Torun, Poland. See this 
[link](https://www.fizyka.umk.pl/~jr/maze.html) for additional information.

The maze module is an optional module for PLUMED2 that needs to be enabled when 
configuring the compilation of PLUMED2. You can either pass a flag
'\-\-enable-modules=maze' or a '\-\-enable-modules=all' when running the 
configure script. 

See the following sections for further information:

- \subpage maze_loss
- \subpage maze_optimizer
- \subpage maze_bias

\page maze_loss Loss

The following list contains the loss functions available in the maze module.

@MAZE_LOSS@

\page maze_optimizer Optimizers

The following list contains the optimizers available in the maze module.

@MAZE_OPTIMIZER@

\page maze_bias Biases

The following list contains the biases available in the maze module.

@MAZE_BIAS@
