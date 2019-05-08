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
//extern CoveringHistory* covhistory;


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


// calculation of used values
double calcNmono(SubprocessFacet *iFacet){//Calculates the Number of (carbon equivalent) particles of one monolayer
	return (iFacet->sh.area*1E-4)/(pow(carbondiameter, 2));
}

double calcdNsurf(){//Calculates the (carbon equivalent relative) mass factor
	return sHandle->wp.gasMass/12.011;
}


llong getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns covering from hitbuffer
	BYTE *buffer;
	buffer = hitbuffer->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return facetHitBuffer->hit.covering;
}
/*
void testprintcovering(Databuff *hitbuffer){
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
					std::cout <<'\t' <<getCovering(&f, hitbuffer);
				}
		}
		std::cout <<std::endl;
}*/

double calcCoverage(SubprocessFacet *iFacet, Databuff *hitbuffer){ // calculates coverage depending on covering (number particles on facet)
	llong covering;
	double coverage;

	covering = getCovering( iFacet, hitbuffer);

	coverage = covering /(calcNmono(iFacet)/calcdNsurf());
	return coverage;
}

std::tuple<double, double> calctotalDesorption(){// desorptionrate as well as total Number of desorbed particles
	double desrate =0.0, totaldes=0.0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					double facetdes = f.sh.desorption;
					std::cout<< "f.sh.desorption = " << f.sh.desorption << std::endl;
					desrate+=facetdes/ (1.38E-23*f.sh.temperature);
					std::cout<< "desrate = " << desrate << std::endl;
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

llong getnbDesorbed(Databuff *hitbuffer_sum){
	BYTE *buffer;
	buffer = hitbuffer_sum->buff;
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)buffer;

	return gHits->globalHits.hit.nbDesorbed;
}

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
	std::cout << "wp.finalOutgassingRate [Pa*m^3/s] = " << sHandle->wp.finalOutgassingRate << std::endl;
	std::cout << "desrate [1/s] = " << desrate << std::endl;
	std::cout << "wp.finalOutgassingRate + desrate [1/s] = " << sHandle->wp.finalOutgassingRate + desrate << std::endl;
	std::cout << "gHits->globalHits.hit.nbDesorbed [1] = " << nbDesorbed<< std::endl;
	std::cout << "(wp.finalOutgassingRate + desrate)/gHits->globalHits.hit.nbDesorbed [1/s] = "<< (sHandle->wp.finalOutgassingRate + desrate)/nbDesorbed << std::endl;
	return (sHandle->wp.finalOutgassingRate+desrate) / nbDesorbed;
	}
/* //Wahrscheinlich brauchen wir das nicht.
double calcRealCovering(SubprocessFacet *iFacet){ //TODO not sure yet

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
	double coverage;
	double temperature;

	temperature=iFacet->sh.temperature;
	coverage = calcCoverage(iFacet,hitbuffer);
	if (coverage < 1) {
		iFacet->sh.sticking = (p->s1*(1.0 - coverage) + p->s2 * coverage)*(1.0 - exp(-p->E_ad / (kb*temperature)));
	}
	else
	{
		iFacet->sh.sticking  = p->s2 * (1.0 - exp(-p->E_ad / (kb*temperature)));
	}



}

double calcDesorption(SubprocessFacet *iFacet, Databuff *hitbuffer){//This returns ((d'coverage')/dt)de. So to speak desorption rate in units of [1/s]
	double coverage;
	double temperature;
	double desorption=0.0;

	coverage = calcCoverage(iFacet,hitbuffer);
	temperature=iFacet->sh.temperature;
	desorption= 1.0/tau * pow(coverage,p->d) *exp(-p->E_de/(kb*temperature)); //what if coverage ==0??

	return desorption;
}

double calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	double desorption = calcDesorption(iFacet, hitbuffer);
	double desorptionRate = desorption * (calcNmono(iFacet) / calcdNsurf()) * 1.38E-23* iFacet->sh.temperature;
	return desorptionRate;
}

void UpdateCovering(Databuff *hitbuffer_phys, Databuff *hitbuffer_sum, double time_step){//Updates Covering after one Iteration using Krealvirt, resets other counters
	//If one wants to read out pressure and particle density, this must be done before calling UpdateCovering.
	//Calculates with the summed up counters of hitbuffer_sum how many test particles are equivalent to one physical particle.
	//Then the physical values are stored in the hitbuffer.
	double Krealvirt = GetMoleculesPerTP(hitbuffer_sum,0);
	std::cout <<"Krealvirt = " << Krealvirt << std::endl;
	llong covering_phys;
	llong covering_sum;
	double covering_check;
	BYTE *buffer_phys;
	buffer_phys = hitbuffer_phys->buff;
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;
	double test_time_step = pow(10,-14);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer_phys = (FacetHitBuffer *)(buffer_phys + f.sh.hitOffset);
				covering_phys = facetHitBuffer_phys->hit.covering;
				FacetHitBuffer *facetHitBuffer_sum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_sum = facetHitBuffer_sum->hit.covering;
				std::cout<<std::endl << "Facet " << &f << std::endl;
				std::cout << "covering_sum = " << covering_sum << std::endl;
				std::cout<< "covering_phys_before = " << covering_phys << std::endl;
				if (covering_sum > covering_phys){
					llong covering_delta = static_cast < llong > ((covering_sum - covering_phys)*Krealvirt*test_time_step); //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
					covering_phys += covering_delta;
					std::cout << "covering rises"<< std::endl;
				}
				else{
					covering_check = covering_phys + (covering_phys - covering_sum)*Krealvirt*(-1)*test_time_step; //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
					std::cout <<"covering_check = " << covering_check << std::endl;
					if(!(covering_check<0)){
						llong covering_delta = static_cast < llong > ((covering_phys - covering_sum)*Krealvirt*test_time_step); //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
						covering_phys -= covering_delta;
						std::cout << "covering decreases but remains positive" << std::endl;
					}
					else {
						std::cout<<"Upps! Covering darf nicht negativ sein. Iteration wird nicht upgedated."<<std::endl;
						std::cout<<"test: "<<(double)(covering_phys/(covering_phys - covering_sum)*Krealvirt) <<std::endl;
						//std::cout << covering_check << std::endl;
						//nichts updaten
						//iteration neu starten mit weniger nbSteps; Wie viel weniger? 1/10 der vorigen Anzahl?
					}
				}
				std::cout<< "covering_phys_after = " << covering_phys << std::endl;
				facetHitBuffer_phys->hit.covering = covering_phys;
				//Reset of Hitbuffer_phys for the next Iteration Step
				facetHitBuffer_phys->hit.nbAbsEquiv = 0;
				facetHitBuffer_phys->hit.nbDesorbed = 0;
				facetHitBuffer_phys->hit.nbMCHit = 0;
				facetHitBuffer_phys->hit.nbHitEquiv = 0;
				facetHitBuffer_phys->hit.sum_1_per_ort_velocity = 0;
				facetHitBuffer_phys->hit.sum_v_ort = 0;
				facetHitBuffer_phys->hit.sum_1_per_velocity = 0;
				//Reset of Hitbuffer_sum for the next Iteration Step
				facetHitBuffer_sum->hit.nbAbsEquiv = 0;
				facetHitBuffer_sum->hit.nbDesorbed = 0;
				facetHitBuffer_sum->hit.nbMCHit = 0;
				facetHitBuffer_sum->hit.nbHitEquiv = 0;
				facetHitBuffer_sum->hit.sum_1_per_ort_velocity = 0;
				facetHitBuffer_sum->hit.sum_v_ort = 0;
				facetHitBuffer_sum->hit.sum_1_per_velocity = 0;
		}
	}
	if(covering_check){ //TODO always true?
	//Reset GlobalHitBuffer
	GlobalHitBuffer *gHits_phys;
	gHits_phys = (GlobalHitBuffer *)buffer_phys;
	gHits_phys->globalHits.hit.nbMCHit = 0;
	gHits_phys->globalHits.hit.nbHitEquiv = 0;
	gHits_phys->globalHits.hit.nbAbsEquiv = 0;
	gHits_phys->globalHits.hit.nbDesorbed = 0;
	GlobalHitBuffer *gHits_sum;
	gHits_sum = (GlobalHitBuffer *)buffer_sum;
	gHits_sum->globalHits.hit.nbMCHit = 0;
	gHits_sum->globalHits.hit.nbHitEquiv = 0;
	gHits_sum->globalHits.hit.nbAbsEquiv = 0;
	gHits_sum->globalHits.hit.nbDesorbed = 0;
	}
}

//---------------------test--------------------------------------
double preTestTimeStep(SimulationHistory *history, Databuff *hitbuffer_sum, double Krealvirt){
	double test_time_step = pow(10,-14);
	llong covering_phys;
	llong covering_sum;
	double covering_check;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					//FacetHitBuffer *facetHitBuffer_sum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
					covering_phys = history->coveringList.getCurrent(&f);
					covering_sum = getCovering(&f, hitbuffer_sum);


					if (covering_sum < covering_phys){
						covering_check = covering_phys + (covering_phys - covering_sum)*Krealvirt*(-1)*test_time_step; //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
						if(covering_check<0){
							test_time_step=(double)(covering_phys/((covering_phys - covering_sum)*Krealvirt));
							std::cout<<"Change of test_time_step: "<<test_time_step <<std::endl;
							//std::cout << covering_check << std::endl;
							//nichts updaten
							//iteration neu starten mit weniger nbSteps; Wie viel weniger? 1/10 der vorigen Anzahl?
						}
					}
			}
		}

	return test_time_step;

}

void UpdateCovering(SimulationHistory *history, Databuff *hitbuffer_sum, double time_step,llong *nbDesorbed_old){//Updates Covering after one Iteration using Krealvirt, resets other counters
	//If one wants to read out pressure and particle density, this must be done before calling UpdateCovering.
	//Calculates with the summed up counters of hitbuffer_sum how many test particles are equivalent to one physical particle.
	//Then the physical values are stored in the hitbuffer.
	double Krealvirt = GetMoleculesPerTP(hitbuffer_sum, *nbDesorbed_old);
	std::cout <<"nbDesorbed before and after:\t" << *nbDesorbed_old <<'\t';
	*nbDesorbed_old = getnbDesorbed(hitbuffer_sum);
	std::cout << *nbDesorbed_old <<std::endl;
	std::cout <<"Krealvirt = " << Krealvirt << std::endl;

	llong covering_phys;
	llong covering_sum;
	double covering_check;
	//BYTE *buffer_sum;
	//buffer_sum = hitbuffer_sum->buff;
	double test_time_step = preTestTimeStep(history, hitbuffer_sum,  Krealvirt);

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
				//FacetHitBuffer *facetHitBuffer_sum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_phys = history->coveringList.getCurrent(&f);
				covering_sum = getCovering(&f, hitbuffer_sum);

				std::cout<<std::endl << "Facet " << &f << std::endl;
				std::cout << "covering_sum = " << covering_sum << std::endl;
				std::cout<< "covering_phys_before = " << covering_phys << std::endl;

				if (covering_sum > covering_phys){
					llong covering_delta = static_cast < llong > ((covering_sum - covering_phys)*Krealvirt*test_time_step); //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
					covering_phys += covering_delta;
					std::cout << "covering rises by " <<covering_delta << std::endl;
				}
				else{
					covering_check = covering_phys + (covering_phys - covering_sum)*Krealvirt*(-1)*test_time_step; //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
					std::cout <<"covering_check = " << covering_check << std::endl;
					if(!(covering_check<0)){
						llong covering_delta = static_cast < llong > ((covering_phys - covering_sum)*Krealvirt*test_time_step); //Fehlt noch mal Delta_t (timestep)! [...] (Sekunden) als Test!
						covering_phys -= covering_delta;
						std::cout << "covering decreases but remains positive" << std::endl;
					}
					else {
						std::cout<<"Upps! Covering darf nicht negativ sein. Iteration wird nicht upgedated."<<std::endl;
						std::cout<<(covering_phys - covering_sum) <<std::endl;
						//std::cout << covering_check << std::endl;
						//nichts updaten
						//iteration neu starten mit weniger nbSteps; Wie viel weniger? 1/10 der vorigen Anzahl?
					}
				}
				std::cout<< "covering_phys_after = " << covering_phys << std::endl;
				std::cout<< "coveringThreshhold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << std::endl;
				history->coveringList.setCurrent(&f, covering_phys);
		}
	}
}

void UpdateCoveringphys(SimulationHistory *history, Databuff *hitbuffer_sum, Databuff *hitbuffer){
	llong covering_phys;
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;

	BYTE *buffer;
	buffer = hitbuffer->buff;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				FacetHitBuffer *facetHitSum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_phys = history->coveringList.getCurrent(&f);
				facetHitBuffer->hit.covering=covering_phys;
				facetHitSum->hit.covering=covering_phys;
			}
	}
}
