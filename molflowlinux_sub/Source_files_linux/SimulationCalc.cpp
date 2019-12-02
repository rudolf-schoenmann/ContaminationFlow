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

double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number hits from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbHitEquiv; //TODO nbMCHit or nbHitEquiv?
}

llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number hits from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.nbDesorbed;
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

//std::tuple<boost::multiprecision::float128, boost::multiprecision::float128> calctotalDesorption(){// desorptionrate as well as total Number of desorbed particles
boost::multiprecision::float128 calctotalDesorption(){// desorptionrate as well as total Number of desorbed particles
	boost::multiprecision::float128 desrate(0.0);
	//boost::multiprecision::float128 totaldes(0.0);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				if(f.sh.temperature==0) {continue;}

				boost::multiprecision::float128 facetdes = f.sh.desorption;
				//std::cout<< "f.sh.desorption = " << f.sh.desorption << std::endl;
				desrate+=facetdes/ boost::multiprecision::float128(1.38E-23*f.sh.temperature);
				//std::cout<< "desrate = " << desrate << std::endl;
				//totaldes+= facetdes * boost::multiprecision::float128(sHandle->wp.latestMoment/(1.38E-23*f.sh.temperature));
			}
	}
	//return {std::make_tuple(desrate, totaldes)};
	return desrate;
}

double calcEnergy(SubprocessFacet *iFacet, Databuff *hitbuffer){ //TODO verify
	long double coverage=calcCoverage(iFacet,hitbuffer);

	return (double)tanh((1-coverage) * (2*tuneE)/p->W_tr) * (p->E_de - p->H_vap)/2 +(p->E_de+p->H_vap)/2; //tanh(adjust width) * adjust height + adjust bias

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

//-----------------------------------------------------------
// calculation of used values

boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum, llong nbDesorbed_old) // Calculation of Krealvirt
//Returns how many physical molecules one test particle represents per time
{
	llong nbDesorbed = getnbDesorbed(hitbuffer_sum)-nbDesorbed_old;
	if (nbDesorbed == 0) return 0; //avoid division by 0

	boost::multiprecision::float128 desrate(0.0);
	//boost::multiprecision::float128 totaldes(0.0);
	desrate=calctotalDesorption();

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
	return (boost::multiprecision::float128(sHandle->wp.finalOutgassingRate) +desrate) / boost::multiprecision::float128(nbDesorbed);
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
	if(covering<llong(p->coveringLimit)){
		factor= boost::multiprecision::float128(0.8);
		//return boost::multiprecision::float128(0.0);
		//return boost::multiprecision::pow(boost::multiprecision::float128(10.),boost::multiprecision::float128(-100.));
	}

	boost::multiprecision::float128 tau_1=static_cast<boost::multiprecision::float128>(1.0/(h/(kb*temperature)));

	boost::multiprecision::float128 energy_de=static_cast<boost::multiprecision::float128>(calcEnergy(iFacet,hitbuffer));

	desorption = tau_1 * boost::multiprecision::pow(coverage,static_cast<boost::multiprecision::float128>(p->d)) *boost::multiprecision::exp(-energy_de/static_cast<boost::multiprecision::float128>(kb*temperature));
	//if (Desorption Energy/Temperature) >~ 1.02E-20J/K, desorption will be zero

	/*if(desorption.convert_to<long double>()==0.0 || desorption.convert_to<double>()==0.0){
		std::cout <<"Facet "<<getFacetIndex(iFacet)<<": Desorption double = 0.0. Desorption long double = "<<desorption.convert_to<long double>()<<". Desorption float128 = "<<desorption <<std::endl;
		p->outFile <<"Facet "<<getFacetIndex(iFacet)<<": Desorption double = 0.0. Desorption long double = "<<desorption.convert_to<long double>()<<". Desorption float128 = "<<desorption <<std::endl;
	}*/

	return factor * desorption;
}

boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	boost::multiprecision::float128 desorption = calcDesorption(iFacet, hitbuffer);
	boost::multiprecision::float128 desorptionRate = desorption * boost::multiprecision::float128(kb* iFacet->sh.temperature * calcNmono(iFacet) /calcdNsurf());
	return desorptionRate;
}

double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f){
	double scaleY = 1.0 / (f->sh.area * 1E-4); //0.01: Pa->mbar
	//TODO is this correct?
	return scaleY *GetMoleculesPerTP(hitbuffer_sum,0).convert_to<double>() * f->tmpCounter[0].hit.sum_1_per_ort_velocity;
}

double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f){
	double scaleY = 1.0 / (f->sh.area  * 1E-4)* sHandle->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
	return f->tmpCounter[0].hit.sum_1_per_ort_velocity*scaleY * GetMoleculesPerTP(hitbuffer_sum,0).convert_to<double>();
}
