/*
Program:     ContaminationFlow
Description: Monte Carlo simulator for satellite contanimation studies
Authors:     Rudolf Schönmann / Hoai My Van
Copyright:   TU Munich
Forked from: Molflow (CERN) (https://cern.ch/molflow)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

/*
 *  This file contains calculations for the contamination
 */

#include "SimulationLinux.h"
#include "GLApp/MathTools.h"
#include <math.h>

// Global handles
extern Simulation* sHandle; //Declared at molflowSub.cpp

// calculation of used values
double calcNmono(SubprocessFacet *iFacet)
{
	return iFacet->sh.area/Sqr(76E-12);
}

double calcdNsurf(){
	return sHandle->wp.gasMass/12.011;
}

double calcKrealvirt(int m){ //TODO not sure yet
	double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
	return timeCorrection;
}

double calcRealCovering(SubprocessFacet *iFacet){ //TODO not sure yet
	double covering=0.0;

	size_t nbMoments = sHandle->moments.size();
	for (size_t m = 0; m <= nbMoments; m++) {
		covering += (double)iFacet->tmpCounter[m].hit.covering*calcKrealvirt(m);

	}
	return covering;
}


// calculations for simulation
double calcCoveringUpdate(SubprocessFacet *iFacet)
{
	//TODO: adapt units, this may not yet be the correct result
	double N_mono= calcNmono(iFacet);
	double dN_surf=calcdNsurf();

	return dN_surf/N_mono;

}

void calcStickingnew(SubprocessFacet *iFacet) {
	double s1 = 0.1;
	double s2 = 0.2;
	double E_ad = pow(10, -21);
	//double E_de = 1.5*pow(10, -21);
	double kb = 1.38 * pow(10, -23);

	double temperature;
	double covering=calcRealCovering(iFacet);

	if(covering>0.0){
		temperature=iFacet->sh.temperature;
		if (covering < 1) {
			iFacet->sh.sticking = (s1*(1 - covering) + s2 * covering)*(1 - exp(-E_ad / (kb*temperature)));
		}
		else
		{
			iFacet->sh.sticking  = s2 * (1 - exp(-E_ad / (kb*temperature)));
		}

	}
	else{
		iFacet->sh.sticking =0.0;
	}

}

double calcDesorption(SubprocessFacet *iFacet){
	double tau=1;
	double d=1;
	double E_de= 2.5E-21;
	double kb = 1.38E-23;

	double temperature;
	double covering=calcRealCovering(iFacet);

	double desorption=0.0;

	if(covering>0.0){
		temperature=iFacet->sh.temperature;
		double N_mono= calcNmono(iFacet);
		double dN_surf=calcdNsurf();

		desorption= 1/tau * pow(covering,d) *exp(-E_de/(kb*temperature)) * N_mono/dN_surf;

	}

	return desorption;
}


