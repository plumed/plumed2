Exercise constant-only ENSEMBLE during input reading, the trigger for #1606.
Check means and raw/central second moments, with and without powers and
reweighting. One replica reproduces the minimal crash; two replicas use
constants 1 and 3 and unequal weights 1:3. Four ranks also exercise broadcasts
within each replica. The expected values are analytic and checked through BIASVALUE.
