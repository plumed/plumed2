# GPU parallelism in PLUMED

For certain actions in PLUMED you can use the USEGPU flag. This flag turns on an experimental GPU parallized version of the 
command. GPU parallelism in PLUMED has been implemented using [openACC](https://www.openacc.org) and is currently experimental. We are actively working
on these features at the moment. __There is thus no guarantee that the GPU accelerated versions of actions are any faster than 
the CPU versions.__ If you have experimented with these features on your own calculations we would love to hear from you (even 
if your experience was negative.)
