# First-step inactive domain exchange

Two MPI ranks own one atom each. The first MD step activates only a box-only
VOLUME output; an atomic RESTRAINT becomes active on the next step through
STRIDE=2. The atoms move on every step.

Check the harmonic bias and both local forces against analytic values, with
the existing 1e-12 absolute tolerance, for synchronous and asynchronous data
sharing. STRIDE=2 intentionally doubles the force at active steps.

The first-step domain exchange initializes constant atom properties even if
the domain action is inactive. It must be received before firststep is
cleared. Otherwise the asynchronous path consumes stale coordinates on later
active steps and leaves unmatched messages. The synchronous path is a control.
No reference, CV definition, variance floor, or force clipping is changed.
