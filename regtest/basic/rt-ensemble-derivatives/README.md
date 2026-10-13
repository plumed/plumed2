# ENSEMBLE derivatives

Finite-difference both the sampled observable and its reweighting energy on
each replica independently. Cover uniform and unequal weights, raw and central
moments of order 2 and 3, powers 1 and 2, and both mean and moment outputs.

The default is two replicas and two MPI ranks per replica. The executable
accepts the replica count as its first argument; the MPI size must be divisible
by that count. Three replicas test odd layouts and nonzero third central
moments. The force oracle is the derivative of the single reported ensemble
value, not the sum of replicated values.
