METAD with ADAPTIVE=DIFF and replica exchange (GREX interface).

Four replicas, partners 0<->1 and 2<->3. Each replica follows one of two "configurations"
(a distance oscillating around 1 or around 3). At step 100 all pairs attempt an exchange that is
accepted (the configurations are swapped), at step 200 only the pair 0<->1 attempts an exchange
that is rejected.

- after the accepted exchange the hills keep a width of the order of the fluctuations
  (the running average of the CV is restarted, the variance is kept)
- the rejected exchange does not change anything: HILLS.0 and HILLS.1 must be identical to
  HILLS.2 and HILLS.3 (checked in output)
