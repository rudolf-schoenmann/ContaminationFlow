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
extern CoveringHistory* covhistory;

int getFacetIndex(SubprocessFacet *iFacet){
	int idx = 0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				//std::cout <<iFacet << '\t' <<&f <<std::endl;
					if(iFacet==&f){
						return idx;
					}
					idx+=1;
			}
	}
	return -1;
}

// calculation of used values
double calcNmono(SubprocessFacet *iFacet){//Calculates the Number of (carbon equivalent) particles of one monolayer
	return (iFacet->sh.area*1E-4)/(pow(76E-12, 2));
}

double calcdNsurf(){//Calculates the (carbon equivalent) mass
	return sHandle->wp.gasMass/12.011;
}

std::tuple<double, double> calctotalDesorption(){ //adapted from totaloutgassingworker in worker.cpp
	//Since this old Molflow code 'calctotalDesorption' calculates the Number of particles which have left all facets due to OUTGASSING;
	//Therefore, in order to get he Number of particles which have left all facets due to OUTGASSING and due to the desorption of the adsorbate
	//you have to add the desrate (= Number of particles which have left all facets due to desorption of the adsorbate);
	double desrate, totaldes=0.0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					double facetdes=calcDesorptionRate(&f);
					desrate+=facetdes/ (1.38E-23*f.sh.temperature);
					totaldes+=sHandle->wp.latestMoment * facetdes / (1.38E-23*f.sh.temperature);;
			}
	}
	return {std::make_tuple(desrate, totaldes)};
}
/*
// TODO is this correct? => Definitively not! Use GetMolecules per TP
double calcKrealvirt(SubprocessFacet *iFacet, int moment){ //TODO not sure yet; c.f. Worker::GetMoleculesPerTP
	double desrate, totaldes=0.0;
	std::tie( desrate,  totaldes)=calctotalDesorption();
	double timeCorrection = moment == 0 ? (sHandle->wp.finalOutgassingRate+desrate) : (sHandle->wp.totalDesorbedMolecules + totaldes) / sHandle->wp.timeWindowSize;
	//Hier entspricht jetzt ein einziger Hit der (Ausgas + Desorptions)-Rate!
	//Rudi: Lass uns den Code von GetMoleculesPerTP nutzen. Dann müssen wir nur noch schauen, wie wir das mit dem globalHitCache machen
	assert(timeCorrection>0.0);
	return timeCorrection;

}
*/
double GetMoleculesPerTP(size_t moment) // alternative for calcKrealvirt
//Returns how many physical molecules one test particle represents
{
	if (sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed == 0) return 0; //avoid division by 0
	double desrate, totaldes=0.0;
	std::tie( desrate,  totaldes)=calctotalDesorption();
	if (moment == 0) {
		//Constant flow
		//Each test particle represents a certain real molecule influx per second
		return (sHandle->wp.finalOutgassingRate+desrate) / sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
	}
	else {
		//Time-dependent mode
		//Each test particle represents a certain absolute number of real molecules
		return ((sHandle->wp.totalDesorbedMolecules+totaldes) / sHandle->wp.timeWindowSize) / sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
	}
}

double calcRealCovering(SubprocessFacet *iFacet){ //TODO not sure yet

	double covering= ((double)iFacet->tmpCounter[0].hit.covering)*GetMoleculesPerTP(0); // only one moment used; one moment means stationary simulation, first moment is moment 0
	return covering;
}

double calcCovering(SubprocessFacet *iFacet){ //TODO not sure yet
	double covering = (double)iFacet->tmpCounter[0].hit.covering; // only one moment used
	return covering;
}


// calculations for simulation
double calcCoveringUpdate(SubprocessFacet *iFacet)
{
	double N_mono= calcNmono(iFacet);
	double dN_surf=calcdNsurf();

	return dN_surf/N_mono;

}

void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer) {//Calculates sticking coefficient dependent on covering.
	double s1 = 0.1;
	double s2 = 0.2;
	double E_ad = 1E-21;
	//double E_de = 1.5E-21;
	double kb = 1.38E-23;
	double covering;
	double temperature;
	//int facetidx = getFacetIndex(iFacet);
	//std::cout <<facetidx <<std::endl<<std::endl;
	BYTE *buffer;
	buffer = hitbuffer->buff;
	/*
	GlobalHitBuffer *gHits;
	double globalCovering;
	//Should wie introduce globalCovering? It does some kind of average covering over all facets.
	//That is not really helpful, but could be a measure of how many particles are leaving the geometry through a hole/pump.
	//globalCovering = sHandle->tmpGlobalResult.globalHits.hit.covering;
	*/
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);
	covering = facetHitBuffer->hit.covering;
	/*
	 * Ich würde das covering lieber aus dem Hitbuffer neu berechnen lassen.
	 * covhistory können wir ja behalten. Vielleicht nur im Hauptprozess, da hat es am meisten Sinn...
	 * Wie schaut covhistory aus, wenn wir eine Texture haben?
	 *
	assert(facetidx < covhistory->pointintime_list.back().second.size());
	double covering = covhistory->pointintime_list.back().second[facetidx];
	*/
	//std::cout <<facetidx <<'\t' <<covering<<std::endl;
	//double covering=calcCovering(iFacet); // double covering=calcRealCovering(iFacet);

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

double calcDesorption(SubprocessFacet *iFacet){//This returns ((d'covering')/dt)de. So to speak desorption rate in units of [1/s]
	double tau=1E-13;
	double d=1;
	double E_de= 1.5E-21;
	double kb = 1.38E-23;

	double temperature;

	int facetidx = getFacetIndex(iFacet);
	//std::cout <<facetidx <<std::endl<<std::endl;
	assert(facetidx < covhistory->pointintime_list.back().second.size());
	double covering = covhistory->pointintime_list.back().second[facetidx];
	//std::cout <<facetidx <<'\t' <<covering<<std::endl;
	//double covering=calcCovering(iFacet); // double covering=calcRealCovering(iFacet);

	double desorption=0.0;

	if(covering>0.0){
		temperature=iFacet->sh.temperature;
		desorption= 1.0/tau * pow(covering,d) *exp(-E_de/(kb*temperature));
	}

	return desorption;
}

double calcDesorptionRate(SubprocessFacet *iFacet) {//This returns ((d'covering')/dt)de * Nmono * kb*T. So to speak desorption rate in units of [Pa m³/s]
	double desorption = calcDesorption(iFacet);
	double desorptionRate = desorption * calcNmono(iFacet) * 1.38E-23* iFacet->sh.temperature;
	return desorptionRate;
}
