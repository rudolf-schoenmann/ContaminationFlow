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
extern SimulationHistory* simHistory;

//-----------------------------------------------------------

//----get values from buffer/handle
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

boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet){ // returns covering from simHistory
	return simHistory->coveringList.getCurrent(iFacet);
}

double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of hits from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbHitEquiv; //TODO nbMCHit or nbHitEquiv?
}

llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of desorbed testparticles from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbDesorbed;
}
llong getnbAdsorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of adsorbed testparticles from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbAbsEquiv;
}

llong getnbDesorbed(Databuff *hitbuffer_sum){
	BYTE *buffer;
	buffer = hitbuffer_sum->buff;
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)buffer;

	return gHits->globalHits.hit.nbDesorbed;
}

//-----------------------------------------------------------

//----calculation of useful intermediate values
double calcNmono(SubprocessFacet *iFacet){//Calculates the Number of (carbon equivalent) particles of one monolayer
	return (iFacet->sh.area*1E-4)/(pow(p->particleDia, 2));
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

boost::multiprecision::float128 calcCoverage(SubprocessFacet *iFacet){ // calculates coverage depending on covering (number particles on facet)
	boost::multiprecision::uint128_t covering = getCovering(iFacet);

	boost::multiprecision::float128 coverage = boost::multiprecision::float128(covering) /boost::multiprecision::float128(calcNmono(iFacet)/calcdNsurf());
	return coverage;
}

boost::multiprecision::float128 calctotalDesorption(){// calculates the desorbed particles of all facets
	boost::multiprecision::float128 desrate(0.0);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				if(f.sh.temperature==0) {continue;}

				boost::multiprecision::float128 facetdes = f.sh.desorption;
				desrate+=facetdes;
			}
	}
	return desrate;
}

double calcStep(long double variable, double start, double end, double inflection_point, double Wtr){

	if(start==end)
		return start;
	else
		return (double)tanh(((variable - inflection_point)/Wtr) *2*tuneE) * (end - start)/2 +(start+end)/2; //tanh(adjust width) * adjust height + adjust bias
		//formula remains the same. But written like this, seems to be easier to be understood!
		// factor of 2, since tanh(x*2*tuneE) is then 0,99 for x = 1/2 as well as -0,99 for x = -1/2; here the normalized transition width is 1/2 - (-1/2) = 1;
	/* wrong:
	else if(start>end)
		return (double)tanh((inflection_point-variable) * (2*tuneE)/Wtr) * (start - end)/2 +(start+end)/2; //tanh(adjust width) * adjust height + adjust bias
	else
		return (-1.0)*(double)tanh((inflection_point-variable) * (2*tuneE)/Wtr) * (start - end)/2 +(start+end)/2; //-tanh(adjust width) * adjust height + adjust bias
		//in the 'else' case an additional minus sign (it's the absolute value of '(start-end)/2') has been forgotten!
	*/
}

double calcEnergy(SubprocessFacet *iFacet){ //TODO verify
	long double coverage=calcCoverage(iFacet).convert_to<long double>();
	return calcStep(coverage, p->E_de, p->H_vap,1, p->W_tr);
}


//-----------------------------------------------------------
// calculation of used values


boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum, llong smallCoveringFactor){ // Calculation of the new Krealvirt
//Returns how many physical molecules one test particle represents

	llong nbDesorbed = getnbDesorbed(hitbuffer_sum);
	if (nbDesorbed == 0) return 0; //avoid division by 0

	boost::multiprecision::float128 desrate(0.0);
	desrate=calctotalDesorption();

	CalcTotalOutgassingWorker();


	return (boost::multiprecision::float128(sHandle->wp.finalOutgassingRate) +desrate) / (boost::multiprecision::float128(nbDesorbed)/smallCoveringFactor);
	}


void calcStickingnew(SubprocessFacet *iFacet) {//Calculates sticking coefficient dependent on covering.

	llong covering=getCovering(iFacet).convert_to<llong>();

	if (covering>=0){
		iFacet->sh.sticking = p->sticking;}
	else{
		iFacet->sh.sticking = 0.0;
	}

}

boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet){//This returns Delta'covering' in units of [1].
	boost::multiprecision::float128 coverage;
	double temperature;
	boost::multiprecision::float128 desorption(0.0);
	double time_step = simHistory->stepSize;
	coverage = calcCoverage(iFacet);
	temperature=iFacet->sh.temperature;
	boost::multiprecision::float128 tau_0=static_cast<boost::multiprecision::float128>(h/(kb*temperature));
	boost::multiprecision::float128 energy_de=static_cast<boost::multiprecision::float128>(p->E_de);
	boost::multiprecision::float128 enthalpy_vap=static_cast<boost::multiprecision::float128>(p->H_vap);
	boost::multiprecision::float128 tau_subst = tau_0 * boost::multiprecision::exp(energy_de/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the substrate

	if(coverage==0 || temperature==0){
		return 0.0;
	}
	else if (coverage <= 1){
		desorption = coverage *(1 - boost::multiprecision::exp(-time_step/tau_subst));//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
	}
	else{//coverage > 1
		boost::multiprecision::float128 tau_ads = tau_0 * boost::multiprecision::exp(enthalpy_vap/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the adsorbate
		if ((coverage - 1) >= (time_step/tau_ads)){//There are more layers (excluding the first monolayer), than desorbing while the iteration time.
			desorption = time_step/tau_ads;//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
		}
		else{//(coverage - 1) < (time_step/tau_ads): There are less layers (excluding the first monolayer), than desorbing while the iteration time.
			double time_step_ads;//time while particles desorbe form multilayer until one single monolayer is reached
			time_step_ads = (double)(tau_ads*(coverage - 1));
			double time_step_subst;//time while particles desorbe form single monolayer
			time_step_subst = time_step - time_step_ads;
			desorption = coverage -1 + (1 - boost::multiprecision::exp(-time_step_subst/tau_subst));//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
		}
	}
	desorption = desorption *(calcNmono(iFacet)/calcdNsurf());//This returns Delta'covering' in units of [1]. 1 means one particle.
	return desorption;
}
/* Don't needed in the new Krealvirt approach...
boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	boost::multiprecision::float128 desorption = calcDesorption(iFacet);
	boost::multiprecision::float128 desorptionRate = desorption * boost::multiprecision::float128(kb* iFacet->sh.temperature * calcNmono(iFacet) /calcdNsurf());
	return desorptionRate;
}
*/

double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f, llong smallCoveringFactor){
	double scaleY = 1.0 / (f->sh.area * 1E-4); //1E4 is conversion from m2 to cm2
	return scaleY *GetMoleculesPerTP(hitbuffer_sum, smallCoveringFactor).convert_to<double>() * f->tmpCounter[0].hit.sum_1_per_ort_velocity;
}

double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f, llong smallCoveringFactor){//calculates Pressure of facet. Output value's unit is mbar.
	double scaleY = 1.0 / (f->sh.area  * 1E-4)* sHandle->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
	return f->tmpCounter[0].hit.sum_v_ort*scaleY* GetMoleculesPerTP(hitbuffer_sum, smallCoveringFactor).convert_to<double>();
}


//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
double calcEnergy(SubprocessFacet *iFacet, Databuff *hitbuffer){ //TODO verify

	long double coverage=calcCoverage(iFacet,hitbuffer);
	return calcStep(coverage, p->E_de, p->H_vap,1, p->W_tr);
}
void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer) {//Calculates sticking coefficient dependent on covering.

	llong covering=getCovering(iFacet,hitbuffer);

	if (covering>=0){
		iFacet->sh.sticking = p->sticking;}
	else{
		iFacet->sh.sticking = 0.0;
	}
}

boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet, Databuff *hitbuffer){//This returns ((d'coverage')/dt)de. So to speak desorption rate in units of [1/s]
	boost::multiprecision::float128 coverage;
	llong covering=getCovering(iFacet,hitbuffer);
	double temperature;
	boost::multiprecision::float128 desorption(0.0);

	coverage = boost::multiprecision::float128 (calcCoverage(iFacet,hitbuffer));
	temperature=iFacet->sh.temperature;

	if(coverage==0 || temperature==0){
		return 0.0;
	}
	boost::multiprecision::float128 factor(1.0);

	boost::multiprecision::float128 d = boost::multiprecision::float128(calcStep((long double)(coverage), 1, 0, 1, p->W_tr));
	boost::multiprecision::float128 tau_1=static_cast<boost::multiprecision::float128>(1.0)/static_cast<boost::multiprecision::float128>(h/(kb*temperature));

	boost::multiprecision::float128 energy_de=static_cast<boost::multiprecision::float128>(calcEnergy(iFacet,hitbuffer));

	desorption = tau_1 * boost::multiprecision::pow(coverage,static_cast<boost::multiprecision::float128>(d)) *boost::multiprecision::exp(-energy_de/static_cast<boost::multiprecision::float128>(kb*temperature));
	//if (Desorption Energy/Temperature) >~ 1.02E-20J/K, desorption will be zero

	return factor * desorption;
}

boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	boost::multiprecision::float128 desorption = calcDesorption(iFacet, hitbuffer);
	boost::multiprecision::float128 desorptionRate = desorption * boost::multiprecision::float128(kb* iFacet->sh.temperature * calcNmono(iFacet) /calcdNsurf());
	return desorptionRate;
}*/
