This module contains implementations of various dimensionality reduction algorithms.
This [tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html) provides an introduction 
to the ways in which these tools have been used in the papers outlined below. 

In general, all these actions take a collection of vectors in input. You can use values from any action that outputs a 
vector as input. Most commonly, however, the input to these actions is the output from a [COLLECT_FRAMES](COLLECT_FRAMES.md) 
shortcut or a subset of the points that were collected with a [COLLECT_FRAMES](COLLECT_FRAMES.md) shortcut that have been selected 
with one of the tools from the [landmarks](module_landmarks.md) module.

Many dimensionality reduction algorithms then work in a manner similar to the way we use when we make maps. You start with distances
between London, Belfast, Paris and Dublin and then you try to arrange points on a piece of paper so that the (suitably transformed)
distances between the points in your map representing each of those cities are related to the true distances between the cities.
Stating this more mathematically dimensionality reduction algorithms endeavor to find an <a href="http://en.wikipedia.org/wiki/Isometry">isometry</a>
between points distributed in a high-dimensional space and a set of points distributed in a low-dimensional plane.
In other words, if we have $M$ $D$-dimensional points, $\mathbf{X}$,
and we can calculate dissimilarities between pairs them, $D_{ij}$, we can, with an MDS calculation, try to create $M$ projections,
$\mathbf{x}$, of the high dimensionality points in a $d$-dimensional linear space by trying to arrange the projections so that the
Euclidean distances between pairs of them, $d_{ij}$, resemble the dissimilarities between the high dimensional points.  In short we minimize:

$$
\chi^2 = \sum_{i \ne j} w_i w_j \left( F(D_{ij}) - f(d_{ij}) \right)^2
$$

where $F(D_{ij})$ is some transformation of the distance between point $X^{i}$ and point $X^{j}$ and $f(d_{ij})$ is some transformation
of the distance between the projection of $X^{i}$, \f$x^i\f$, and the projection of $X^{j}$, $x^j$.  $w_i$ and $w_j$ are the weights
of configurations $X^i$ and $^j$ respectively. 
