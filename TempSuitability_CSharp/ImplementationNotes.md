In the original implementation, we iterated once (at each cell location) through the 2-hour slices comprising the overall modelled period. 
At each step a new cohort of mosquitoes was modelled as emerging and then an inner loop was run over the life of that cohort 
in 2-hour slices to add the infectiousness contribution of that cohort to an overall sum-of-infectiousness for every 
time slice day in which it was alive. 

In this re-implementation, we iterate a single time (at each cell location) through the 2-hour slices comprising the overall runtime, and maintain a population 
of mosquitoes using a FIFO queue structure. At each step the oldest cohort is removed, a new one added, and 
the contribution of all living cohorts is summarised. 

The number of calculations is the same, but much less data has to be held in memory and access to the data is more linear, 
allowing the code to be better optimised. For example we do not need to maintain (or access) arrays of daylight hours or survival probability, but 
calculate these values only when they are required and then discard them. We haven't yet tested whether 
these improvements outweigh the overhead of running in C# rather than C.

Additionally, parallelisation is now cell-based and multithreaded, rather than multiprocessing on blocks of cells. 

