This is a re-implementation of the model described in: 
Air temperature suitability for Plasmodium falciparum malaria transmission in Africa 2000-2012: a high-resolution spatiotemporal prediction (Weiss et al, 2014, DOI: 10.1186/1475-2875-13-171)

The model described in the above-referenced paper was written in IDL with the core numerical calculation code written in C. The IDL code 
dealt with reading and writing the data image files, and created temporary data files for the C code to use.
This implementation is written in C# throughout and all I/O will be handled within the program. 

At present the code is still in development.

