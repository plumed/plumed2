The actions implemented in this module were inspired by concept of the Willard-Chandler surface that is described in the first
two papers referenced below. The way that this idea and these implementations has been used for a range of applications is then discussed
in the other paper referenced below.

Willard and Chander's method for finding the position of the water surface starts by introducing an instantaneous density field:

$$
\rho(x,y,z) = \sum_{i=1}^N K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right)
$$

In this expression the sum runs over a collection of atoms.  Each of these atoms is located at a position $(x_i,y_i,z_i)$ and the $K$ is a smooth, 3-dimensional Kernel 
function (e.g. a Gaussian) that integrates to one.  Each atom thus contributes a function that is peaked at the atom's position and which decays over some characteristic 
length scale $\lambda$ to the final field $\rho(x,y,z)$. The value of $\rho(x,y,z)$ is thus large 
if this function is evaluated at a point that near to the specified atoms and small if it is evaluated at some position that is far from them.  Consequently,
we can find interfaces that separate the region where the atoms of interest are located and their surroundings by finding the 2-dimensional manifold containing points for which:    

$$
\varphi(x,y,z) - \varphi_0 = 0    
$$

In this expression $\varphi_0$ is a parameter and $\varphi(x,y,z)$ is calculated using the first equation above. This procedures thus converts a particle based representation of the 
atoms and their surroundings into a coarse grained representation that describes the extents region containing atoms and the extent of the region that does not contain atoms.
In other words, if we have performed a simulation of coexisting liquid and gas phases we can use the expressions above to locate the interface between the two phases.

Similarly, if we search for a contour using the second expression above in a phase field that is defined as follows:

$$
\rho'(x,y,z) = \frac{ \sum_{i=1}^N s_i K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }{ \sum_{i=1}^N K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }
$$

where $s_i$ is the value of one of the symmetry functions described in the documentation for the [symfunc](module_symfunc.md) that can be used to distinguish solid-like atoms from liquid-like atoms evaluated 
based on the environment around atom $i$, we can find solid-liquid interfaces. This works becuase $\rho'(x,y,z)$ is large in regions of the box where there are many solid-like atoms and small in regions of the 
simulation box that are filled with liquid-like atoms.
