# RANS_stableOF50


## Introduction to the problem and solution
In past simulations of free-surface waves (both breaking and non-breaking) using RANS models, there have been a marked a collective tendency to over-estimate the turbulence levels. In cases involving breaking waves, this has even been most pronounced prior to breaking with turbulence levels pre-breaking of same order of magnitude as in the surf-zone. The underlying cause was identified by Mayer and Madsen (2000) who showed that turbulence production exists in potential flow waves. They further showed that if omega falls below a certain threshold, the model becomes unstable, leading to an exponential growth in both the eddy viscosity and turbulent kinetic energy.


This code is still being developed and is not yet finished. Stablized turbulence models for OpenFOAM-5.0

Standard two-equation RANS models are unstable when applied to free-surface waves as documented in Larsen and Fuhrman 2018. In Larsen and Fuhrman 2018 it was likewise shown how standard two-equation turbulence models can be formally stabilized. This repository contains formally stabilized versions of standard two-equation turbulence models for foam-extend 3.1.

## Installation
Source OpenFOAM-5.0:

	source /home/$USER/OpenFOAM/OpenFOAM-5.0/etc/bashrc

In a linux terminal download the package with git by typing:

	git clone https://github.com/BjarkeEltardLarsen/RANS_stableOF50.git
	
Create folder for turbulence model (if the folders already exist skip this part)

	mkdir $WM_PROJECT_USER_DIR/src $WM_PROJECT_USER_DIR/src/turbulence $WM_PROJECT_USER_DIR/src/turbulence/incompressible

Move the folder to the user source code

	mv RANS_stableFE31 $WM_PROJECT_USER_DIR/src/turbulence/incompressible/
	
Go to the directory and compile the turbulence models

	cd $WM_PROJECT_USER_DIR/src/turbulence/incompressible/RANS_stableFE31
	
	wmake libso
	
Move the tutorials to the desired folder e.g FOAM_RUN

	mv Tutorials $FOAM_RUN
	
## Usage
Include the libary of the stabilized turbulence models in the system/controlDict folder

	libs
	(
    	"libMyStableRASModels.so"
	);

Change the constant/RASproperties
Each of the four models can be chosen by uncommenting the desired model.
If lambda2=0 the models default to their standard OpenFOAM implementation (but with the buoyancy production term added). 

	RASModel        kOmegaStab;
	//RASModel        kOmegaSSTStab;
	//RASModel        kEpsilonStab;
	//RASModel        RNGkEpsilonStab;


	turbulence      on;

	printCoeffs     on;

	kOmegaStabCoeffs
	{
  	lambda2 0.05;
	}

	kOmegaSSTStabCoeffs
	{
  	lambda2 0.05;
	}

	kEpsilonStabCoeffs
	{
  	lambda2 0.05;
	}

	RNGkEpsilonStabCoeffs
	{
  	lambda2 0.05;
	}

Change the system/fvSchemes from (if they are specified individually)

	ddt(k)
	ddt(epsilon)
	ddt(omega)
	...
	div(phi,omega)
	div(phi,epsilon)
	div(phi,k)
	
to 
	
	ddt(rho,k)
	ddt(rho,epsilon)
	ddt(rho,omega)
	...
	div(rho*phi,omega)
	div(rho*phi,epsilon)
	div(rho*phi,k)

## Tutorials

The tutorials consist of eight different version of the compiled turbulence models. In each of these tutorials a wave is inialized in a domain which is exactly one wave length long with cyclic boundaries and a slip condition at the bed. This in other words an idealized case where the flow can considered very close to potential flow. 

The cases ending with "Stab" have lambda2=0.05 and are stabilized version of the new turbulence models. In these cases the initial turbulence levels will decay in time. 

Cases not ending with Stab corresponds to the standard models (with buoyancy production included) and here growth rather than decay in the turbulence level can be seen.

To run a tutorial go to the folder of the tutorials and choose one of the eight cases. 
Then run the RunScript by typing (the tutorials will take a few ours depending on speed of the machine)
	
	RunScript
	
	
