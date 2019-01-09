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
#include <assert.h>
// Global handles
extern Simulation* sHandle; //Declared at molflowSub.cpp

// calculation of used values
double calcNmono(SubprocessFacet *iFacet)
{
	return (iFacet->sh.area*1E-4)/Sqr(76E-12);
}

double calcdNsurf(){
	return sHandle->wp.gasMass/12.011;
}

std::tuple<double, double> calctotalDesorption(){ //adapted from totaloutgassingworker in worker.cpp
	double desrate, totaldes=0.0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					double facetdes=calcDesorption(&f);
					desrate+=facetdes/ (1.38E-23*f.sh.temperature);
					totaldes+=sHandle->wp.latestMoment * facetdes / (1.38E-23*f.sh.temperature);;
			}
	}
	return {std::make_tuple(desrate, totaldes)};
}
// TODO is this correct?
double calcKrealvirt(SubprocessFacet *iFacet, int m){ //TODO not sure yet
	double desrate, totaldes=0.0;
	std::tie( desrate,  totaldes)=calctotalDesorption();
	double timeCorrection = m == 0 ? (sHandle->wp.finalOutgassingRate+desrate) : (sHandle->wp.totalDesorbedMolecules + totaldes) / sHandle->wp.timeWindowSize;
	//std::cout <<sHandle->wp.finalOutgassingRate <<"\t" <<(sHandle->wp.totalDesorbedMolecules + calcDesorption(iFacet)) / sHandle->wp.timeWindowSize <<std::endl;
	assert(timeCorrection>0.0);
	return timeCorrection;

}


double calcRealCovering(SubprocessFacet *iFacet){ //TODO not sure yet
	double covering=0.0;

	size_t nbMoments = sHandle->moments.size();
	for (size_t m = 0; m <= nbMoments; m++) {
		covering += ((double)iFacet->tmpCounter[m].hit.covering)*calcKrealvirt(iFacet,m);
	}
	return covering;
}

double calcCovering(SubprocessFacet *iFacet){ //TODO not sure yet
	double covering=0.0;

	size_t nbMoments = sHandle->moments.size();
	for (size_t m = 0; m <= nbMoments; m++) {
		covering += (double)iFacet->tmpCounter[m].hit.covering;

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
	double covering=calcRealCovering(iFacet); // double covering=calcCovering(iFacet);

	if(covering>0.0){
		temperature=iFacet->sh.temperature;
		if (covering < 1) {
			iFacet->sh.sticking = (s1*(1.0 - covering) + s2 * covering)*(1.0 - exp(-E_ad / (kb*temperature)));
		}
		else
		{
			iFacet->sh.sticking  = s2 * (1.0 - exp(-E_ad / (kb*temperature)));
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
	double covering=calcCovering(iFacet); // double covering=calcRealCovering(iFacet);

	double desorption=0.0;

	if(covering>0.0){
		temperature=iFacet->sh.temperature;
		double N_mono= calcNmono(iFacet);
		double dN_surf=calcdNsurf();

		desorption= 1.0/tau * pow(covering,d) *exp(-E_de/(kb*temperature)) * N_mono/dN_surf;

	}

	return desorption;
}


