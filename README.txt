README

----------------------------------------------------------------------------------------
	MCMC inversions of trace and/or major elements from mantle-derived melts 
	for the thermochemical state of the source (mantle).
         
 	Reference: Oliveira, B., Afonso, J.C., Klocing, M. (2021), Melting 
            	conditions and mantle source composition from probabilistic 
            	joint inversion of major and rare earth element concentrations,
            	Geochim. Cosmochim. Acta, doi:
            	
	Contact: oliveira.bravo.b@gmail.com
----------------------------------------------------------------------------------------

Main files: 	(1) mcmc_wrapper.m  	===>  
		(2) amrum.m  		===>  
		(3) melting_master.m  	===>  
		(4) perplex_lookup_isentrope.m  ===>
		(5) run_nearest_diseqm.m ===>
		(6) diseqm_trace_parallel_fractional


(1) mcmc_wrapper.m
Main mcmc wrapper

(2) amrum.m
Inversion, calls the melting model

(3) melting_master.m
Gets Tp (mantle potential temperature), ztop (final depth of melting) and source composition and 
calls the forward melting problem for both trace and major element composition.

(4) perplex_lookup_isentrope.m
It identifies the Tp-Al2O3-Na2O values that fall closest to trial set of Tp-Al2O3-Na2O parameters 
i.e the closest 6 values within the used Tp-Al2O3-Na2O discretization space (a prism).

(5) run_nearest_diseqm.m
With each of the six vertexes, it either
	retrieves precomputed PTC paths for given chemical composition 
	(option isen = 1 in model_parameters.txt)
	or,
	 computes isentropic PTC paths on the fly
	(option isen = 2 in model_parameters.txt)
And calls the solver for traces (diseqm_trace_fractional.m) for each set of Tp-Al2O3-Na2O values.
Once it retrieves the composition of pooled magmas from each isentrope, it computes an averaged
composition using linear interpolation.

(6) diseqm_trace_fractional.m
Computes both major and trace elment composition of pooled magma for the given trace parameters 
(trace_elements_parameters.txt) and PTC isentropic path. It uses,
	Experimental or infinite diffusion coefficients (<diff> in model_parameters.txt).
	PTC dependant or constant partition coefficients (<part> in model_parameters.txt).


Default parameters reproduce the results in Section 5, where both major and trace element 
compositional data have been used to infer mantle potential temperature (Tp), final depth of 
melting (ztop), and source composition (Al2O3, Na2O and REE's). 
	

----------------------
DISCLAIMER
----------------------

isen = 2 compiles perplex in MATLAB using mex files. The current version is only valid for Linux.
Contact the authors with any enquiry.


----------------------
OTHER ADVANCED OPTIONS
----------------------

To include more oxides in the inversion modify/expand:
	<Sig_sol> in mcmc_wrapper.m
	<params.par0> in mcmc_wrapper.m
	<params.bounds> in mcmc_wrapper.m
And select isen = 2 in model_parameters.txt
	
To modify the misfit function:
	misfit.m
	
To increas/decrease the resolution of the Tp-Al2O3-Na2O discretization space where interpolations
are performed, modify:
	<*_stencil> in mcmc_wrapper.m

Contant the authors with any question or comment.

-------------------------------------------------------------------------------------------------
Copyright (C) 2021  B. Oliveira, J.C. Afonso & M. Klocking

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 Reference: Oliveira, B., Afonso, J.C., Klocing, M. (2021), Melting 
            conditions and mantle source composition from probabilistic 
            joint inversion of major and rare earth element concentrations,
            Geochim. Cosmochim. Acta, doi:



  
