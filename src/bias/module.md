This module contains various methods for adding bias potentials.  
These bias potentials are the core of the enhanced sampling algorithms
that allow you to force rare events to occur on the short timescales that 
are accessible on the molecular dynamics timescale.

In addition, the module also contains tools such as [REWEIGHT_BIAS](REWEIGHT_BIAS.md),
[REWEIGHT_METAD](REWEIGHT_METAD.md) and [REWEIGHT_TEMP_PRESS](REWEIGHT_TEMP_PRESS.md).
These actions allow you to extract "unbiased" free energy surfaces from biased simulations
and to thus extract free energy differences from your enhanced sampling calculations. 
