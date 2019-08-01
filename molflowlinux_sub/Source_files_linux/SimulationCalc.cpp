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
extern Simulation* sHandle; //Declared at molflowlinux_main.cpp
extern ProblemDef* p;

//get values from buffer/handle
int getFacetIndex(SubprocessFacet *iFacet){ // finds index of facet. index used for CoveringHistory class
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

llong getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns covering from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.covering;
}

llong getHits(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number hits from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbHitEquiv; //nbMCHit or nbHitEquiv?
}

llong getnbDesorbed(Databuff *hitbuffer_sum){
	BYTE *buffer;
	buffer = hitbuffer_sum->buff;
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)buffer;

	return gHits->globalHits.hit.nbDesorbed;
}

//-----------------------------------------------------------

// calculation of useful intermediate values
double calcNmono(SubprocessFacet *iFacet){//Calculates the Number of (carbon equivalent) particles of one monolayer
	return (iFacet->sh.area*1E-4)/(pow(carbondiameter, 2));
}

double calcdNsurf(){//Calculates the (carbon equivalent relative) mass factor
	return sHandle->wp.gasMass/12.011;
}

long double calcCoverage(SubprocessFacet *iFacet, Databuff *hitbuffer){ // calculates coverage depending on covering (number particles on facet)
	llong covering;
	long double coverage;

	covering = getCovering( iFacet, hitbuffer);

	coverage = (long double)covering /(long double)(calcNmono(iFacet)/calcdNsurf());
	return coverage;
}

std::tuple<double, double> calctotalDesorption(){// desorptionrate as well as total Number of desorbed particles
	double desrate =0.0, totaldes=0.0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					double facetdes = f.sh.desorption;
					//std::cout<< "f.sh.desorption = " << f.sh.desorption << std::endl;
					desrate+=facetdes/ (1.38E-23*f.sh.temperature);
					//std::cout<< "desrate = " << desrate << std::endl;
					totaldes+=sHandle->wp.latestMoment * facetdes / (1.38E-23*f.sh.temperature);;
			}
	}
	return {std::make_tuple(desrate, totaldes)};
}


//Brauchen wir nicht, weil wir CaclTotalOutgassing haben.
/*
std::tuple<double, double> calctotalOutgassing(){
	double outgassrate, totalout=0.0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					double facetoutgassrate = f.sh.outgassing;
					outgassrate+=facetoutgassrate/ (1.38E-23*f.sh.temperature);
					totalout+=sHandle->wp.latestMoment * facetoutgassrate / (1.38E-23*f.sh.temperature);;
			}
	}
	return {std::make_tuple(outgassrate, totalout)};
}
*/


//Brauchen wir wahrscheinlich nicht.
/*
std::tuple<double, double> calctotalParticles_out(Databuff *hitbuffer){
	std::tuple<double,double> particles_outDes = calctotalDesorption(hitbuffer);
	double desrate = std::get<0>(particles_outDes);
	double totaldes = std::get<1>(particles_outDes);
	std::tuple<double,double> particles_outOut = calctotalOutgassing();
	double outgassrate = std::get<0>(particles_outOut);
	double totalout = std::get<1>(particles_outOut);
	double sum_rate = desrate + outgassrate;
	double sum_total = totaldes + totalout;

	return{std::make_tuple(sum_rate, sum_total)};
}
*/

/*
// Use GetMolecules per TP
double calcKrealvirt(SubprocessFacet *iFacet, int moment){
	double desrate, totaldes=0.0;
	std::tie( desrate,  totaldes)=calctotalDesorption();
	double timeCorrection = moment == 0 ? (sHandle->wp.finalOutgassingRate+desrate) : (sHandle->wp.totalDesorbedMolecules + totaldes) / sHandle->wp.timeWindowSize;
	//Hier entspricht jetzt ein einziger Hit der (Ausgas + Desorptions)-Rate!
	//Rudi: Lass uns den Code von GetMoleculesPerTP nutzen. Dann müssen wir nur noch schauen, wie wir das mit dem globalHitCache machen
	assert(timeCorrection>0.0);
	return timeCorrection;

}
*/

//-----------------------------------------------------------
// calculation of used values

double GetMoleculesPerTP(Databuff *hitbuffer_sum, llong nbDesorbed_old) // Calculation of Krealvirt
//Returns how many physical molecules one test particle represents
{
	llong nbDesorbed = getnbDesorbed(hitbuffer_sum)-nbDesorbed_old;
	if (nbDesorbed == 0) return 0; //avoid division by 0

	double desrate, totaldes =0.0;
	std::tie(desrate,  totaldes)=calctotalDesorption();
	CalcTotalOutgassingWorker();
	//Constant flow
	//Each test particle represents a certain real molecule influx per second
	/*
	std::cout << "wp.finalOutgassingRate [Pa*m^3/s] = " << sHandle->wp.finalOutgassingRate << std::endl;
	std::cout << "desrate [1/s] = " << desrate << std::endl;
	std::cout << "wp.finalOutgassingRate + desrate [1/s] = " << sHandle->wp.finalOutgassingRate + desrate << std::endl;
	std::cout << "gHits->globalHits.hit.nbDesorbed [1] = " << nbDesorbed<< std::endl;
	std::cout << "(wp.finalOutgassingRate + desrate)/gHits->globalHits.hit.nbDesorbed [1/s] = "<< (sHandle->wp.finalOutgassingRate + desrate)/nbDesorbed << std::endl;
	*/
	return (sHandle->wp.finalOutgassingRate+desrate) / nbDesorbed;
	}
/* //Wahrscheinlich brauchen wir das nicht.
double calcRealCovering(SubprocessFacet *iFacet){

	double covering= ((double)iFacet->tmpCounter[0].hit.covering)*GetMoleculesPerTP(0); // only one moment used; one moment means stationary simulation, first moment is moment 0
	return covering;
}

*/

// Brauchen wir eigentlich auch nicht.
/*
double calcCoveringUpdate(SubprocessFacet *iFacet)
{
	double N_mono= calcNmono(iFacet);
	double dN_surf=calcdNsurf();
	//return dN_surf/N_mono;
	return 1; //Das ist natürlich falsch und nur zu Testzwecken eingebaut.

}*/

void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer) {//Calculates sticking coefficient dependent on covering.

	llong covering=getCovering(iFacet,hitbuffer);

	if (covering>=100){
		iFacet->sh.sticking = 1;}
	else{
		iFacet->sh.sticking = 0.0;
	}
}

long double calcDesorption(SubprocessFacet *iFacet, Databuff *hitbuffer){//This returns ((d'coverage')/dt)de. So to speak desorption rate in units of [1/s]
	long double coverage;
	double temperature;
	long double desorption=0.0;

	coverage = calcCoverage(iFacet,hitbuffer);
	llong covering=getCovering(iFacet,hitbuffer);

	temperature=iFacet->sh.temperature;

	if(coverage==0||covering<100){
		return 0.0;
	}
	long double tau=(long double)(h/(kb*temperature));

	desorption= (long double)(1.0/tau) * powl(coverage,(long double)p->d) *expl(-(long double)p->E_de/(long double)(kb*temperature));
	//if (Desorption Energy/Temperature) >~ 1.02E-20J/K, desorption will be zero
	return desorption;
}

long double calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	long double desorption = calcDesorption(iFacet, hitbuffer);
	long double desorptionRate = desorption * (long double)(calcNmono(iFacet) /calcdNsurf()) * (long double)(kb* iFacet->sh.temperature);
	return desorptionRate;
}

double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f){
	double scaleY = 1.0 / (f->sh.area  /*/(double)PROFILE_SIZE*/*1E-4); //0.01: Pa->mbar
	//TODO is this correct?
	return scaleY *GetMoleculesPerTP(hitbuffer_sum,0) * f->tmpCounter[0].hit.sum_1_per_ort_velocity;
}

double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f){
	double scaleY = 1.0 / (f->sh.area  /*/(double)PROFILE_SIZE*/*1E-4)* sHandle->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
	return f->tmpCounter[0].hit.sum_1_per_ort_velocity*scaleY * GetMoleculesPerTP(hitbuffer_sum,0);
}
