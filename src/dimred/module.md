This module contains implementations of various dimensionality reduction algorithms.
This [tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html) provides an introduction 
to the ways in which these tools have been used in the papers outlined below. 

In general, all these actions take a collection of vectors in input. You can use values from any action that outputs a 
vector as input. Most commonly, however, the input to these actions is the output from a [COLLECT_FRAMES](COLLECT_FRAMES.md) 
shortcut or a subset of the points that were collected with a [COLLECT_FRAMES](COLLECT_FRAMES.md) shortcut that have been selected 
with one of the tools from the [landmarks](module_landmarks.md) module.
