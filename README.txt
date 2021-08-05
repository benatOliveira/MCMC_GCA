README

----------------------------------------------------------------------------------------
         MCMC Two-Phase system: Trace and Major element composition inversion
         Created by Beñat Oliveira Bravo
                 oliveira.bravo.b@gmail.com
		 benat.oliveira-bravo@mq.edu.au
----------------------------------------------------------------------------------------

Main files: 	(1) mcmc_wrapper.m  	===>  
		(2) dramrum.m  		===>  
		(3) melting_master.m  	===>  
		(4) perplex_lookup_isentrope.m  ===>
		(5) diseqm_trace_fractional.m


(1) mcmc_wrapper.m
Main mcmc wrapper

(2) dramrum.m
Inversion, calls the melting model

(3) melting_master.m
Gets Tp (mantle potential temperature) and ztop (final depth of melting) and computes PTC paths
Three options (<method> in model_parameters.txt):
	MELTS (not tested)
	PERPLEX (working, perplex_lookup_isentrope)
	ZOU (working, before using it contact Benat)

(4) perplex_lookup_isentrope.m
Different options (<isen> in model_parameters.txt):
if isen == 1
	Computes PTC path along an isentrope for the given model parameters (model_parameters.txt)
	It uses Asimow 1997 approximations
elseif isen == 2
	Computes PTC path along an isentrope for the given model parameters (model_parameters.txt)
	It uses Connolly's law, and iterates the temperature to ensure isentropic path.
elseif isen == 3
	Uses precomputed PTC paths for given chemical composition
end
Once PTC is obtained, it calls the solver for traces (diseqm_trace_fractional.m)

(5) diseqm_trace_fractional.m
Computes trace composition for the given trace parameters (trace_elements_parameters.txt).
Experimental and infinite diffusion coefficients (<diff> in model_parameters.txt):
It also computes major composition


  