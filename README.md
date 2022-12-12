Uploaded: Dec. 8, 2022

Uploaded by: Tim Goudge (tgoudge@jsg.utexas.edu)

This repository contains the model scripts used for Goudge et al. "Assessing Controls on the Incomplete Draining of Martian Open-Basin Lakes". This modeling approach is a modification of the approach initially outlined by Fassett and Goudge (2021, doi: 10.1029/2021JE006979). The original modeling approach scripts are available here: https://github.com/cfassett/ANUGA-erosion.

All the parameters are set via the command line (see example call below), with the exception of the model domain grain size, which is set in model_params.py ("grainD").

The code call from the command line (so long as ANUGA is also installed correctly) is:

	python TAG_implementation_of_FG21_lakeerodeR.py 20000 500 2.8 150 0.01


The parameters set on the command line (in order), with values given in the example above, are:

	Lake Radius [m] = 20000

	Lake Depth [m] = 500
	
	Topographic Exponent (for basin topography) [unitless] = 2.8
	
	Rim Height [m] = 150
	
	Regional Slope [unitless] = 0.01
