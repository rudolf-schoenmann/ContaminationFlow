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
#include <cmath>
#include <cassert>
#include <Random.h>

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

FacetHitBuffer* getFacetHitBuffer(SubprocessFacet *iFacet, Databuff *hitbuffer){
	return (FacetHitBuffer *)(hitbuffer->buff + iFacet->sh.hitOffset);
}

llong getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns covering from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->hit.covering;
}

boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet){ // returns covering from simHistory
	return simHistory->coveringList.getCurrent(iFacet);
}

double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of hits from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->hit.nbHitEquiv; //TODO nbMCHit or nbHitEquiv?
}

llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of desorbed testparticles from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->hit.nbDesorbed;
}
llong getnbAdsorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of adsorbed testparticles from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->hit.nbAbsEquiv;
}

llong getnbDesorbed(Databuff *hitbuffer_sum){
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)hitbuffer_sum->buff;

	return gHits->globalHits.hit.nbDesorbed;
}
/*
std::tuple<double, double, double> getVelocities(SubprocessFacet *iFacet, Databuff *hitbuffer_sum){
	BYTE *buffer;
	buffer = hitbuffer_sum->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return {std::make_tuple(facetHitBuffer->hit.sum_1_per_ort_velocity,facetHitBuffer->hit.sum_v_ort, facetHitBuffer->hit.sum_1_per_velocity)};
}*/

//-----------------------------------------------------------

//----calculation of useful intermediate values
double calcNmono(SubprocessFacet *iFacet){//Calculates the Number of (carbon equivalent) particles of one monolayer
	// area facet / area particle
	return (iFacet->sh.area*1E-4)/(pow(p->particleDia, 2));
}

boost::multiprecision::float128 calcCoverage(SubprocessFacet *iFacet){ // calculates coverage depending on covering (number particles on facet)
	boost::multiprecision::uint128_t covering = getCovering(iFacet);

	return boost::multiprecision::float128(covering) /boost::multiprecision::float128(calcNmono(iFacet));
}

boost::multiprecision::float128 calctotalDesorption(){// calculates the desorbed particles of all facets
	boost::multiprecision::float128 desrate(0.0);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(f.sh.temperature==0) {continue;}
			desrate+=f.sh.desorption;
		}
	}
	return desrate;
}

double calcOutgassingFactor(SubprocessFacet *iFacet){
	return simHistory->stepSize_outgassing/(kb * iFacet->sh.temperature);
}

//-----------------------------------------------------------
// calculation of used values


boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum){ // Calculation of the new Krealvirt
//Returns how many physical molecules one test particle represents

	llong nbDesorbed = getnbDesorbed(hitbuffer_sum);
	if (nbDesorbed == 0) return 0; //avoid division by 0

	boost::multiprecision::float128 des=calctotalDesorption();
	CalcTotalOutgassingWorker();

	return (boost::multiprecision::float128(sHandle->wp.totalOutgassingParticles) +des) / (boost::multiprecision::float128(nbDesorbed)/boost::multiprecision::float128(simHistory->smallCoveringFactor));
	}


void calcSticking(SubprocessFacet *iFacet) {//Calculates sticking coefficient dependent on covering.

	llong covering=getCovering(iFacet).convert_to<llong>();
	// sticking constant, zero if no covering
	if (covering>0){
		iFacet->sh.sticking = p->sticking;}
	else{
		iFacet->sh.sticking = 0.0;
	}

}

boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet){//This returns Delta'covering' in units of [1].
	boost::multiprecision::float128 desorption(0.0);

	boost::multiprecision::float128 coverage= calcCoverage(iFacet);
	double temperature = iFacet->sh.temperature;
	boost::multiprecision::float128 time_step = boost::multiprecision::float128(simHistory->stepSize);

	boost::multiprecision::float128 tau_0=static_cast<boost::multiprecision::float128>(h/(kb*temperature));
	boost::multiprecision::float128 energy_de=static_cast<boost::multiprecision::float128>(p->E_de);
	boost::multiprecision::float128 enthalpy_vap=static_cast<boost::multiprecision::float128>(p->H_vap);
	boost::multiprecision::float128 tau_subst = tau_0 * boost::multiprecision::exp(energy_de/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the substrate

	if(coverage==boost::multiprecision::float128(0) || temperature==0){
		return 0.0;
	}
	else if (coverage <= boost::multiprecision::float128(1)){
		desorption = coverage *(boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step/tau_subst));//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
	}
	else{//coverage > 1
		boost::multiprecision::float128 tau_ads = tau_0 * boost::multiprecision::exp(enthalpy_vap/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the adsorbate
		if(tau_ads == tau_subst){//special case, where desorption rate is constant until a monolayer is reached.
			if ((coverage - boost::multiprecision::float128(1)) >= (time_step/tau_ads)){//There are more layers (excluding the first monolayer), than desorbing while the iteration time.
						desorption = time_step/tau_ads;//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
			}
			else{//(coverage - 1) < (time_step/tau_ads): There are less layers (excluding the first monolayer), than desorbing while the iteration time.
						//time while particles desorbe form multilayer until one single monolayer is reached
						boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(1));
						//time while particles desorbe form single monolayer
						boost::multiprecision::float128 time_step_subst = time_step - time_step_ads;
						desorption = coverage - boost::multiprecision::float128(1) + (boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step_subst/tau_subst));//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
						//Anyway: Here tau_ads == tau_subst! So, it is distinguished to provide a better understanding, what is happening physically
			}
		}
		else{//desorption rate is constant until coverage equals two.
			boost::multiprecision::float128 a = (-1/tau_ads)+(1/tau_subst);
			boost::multiprecision::float128 b = (1/tau_ads)-(2/tau_subst);
			boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(2));
			boost::multiprecision::float128 time_step_mixed = time_step - time_step_ads;
			if((2+a/b)*boost::multiprecision::exp(b*time_step_mixed)-(a/b)>1){//case, where still more than one monolayer is left
				desorption = coverage - boost::multiprecision::float128(2) + boost::multiprecision::float128(2) - (2+a/b)*boost::multiprecision::exp(b*time_step_mixed)+(a/b);
				}
			else{//case, where coverage will be smaller than one.
				time_step_mixed =(1/b)* boost::multiprecision::log((1+a/b)/(2+a/b));
				boost::multiprecision::float128 time_step_subst = time_step - time_step_ads - time_step_mixed;
				desorption = coverage - boost::multiprecision::float128(1) + (boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step_subst/tau_subst));
			}
		}
	}
	return desorption * boost::multiprecision::float128(calcNmono(iFacet));//This returns Delta'covering' in units of [1]. 1 means one particle.
}

double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f){

	if(f->tmpCounter[0].hit.sum_1_per_ort_velocity==std::numeric_limits<double>::infinity()){
		return 0.0;
		}
	else{
		//Correction for double-density effect (measuring density on desorbing/absorbing facets):

		//Normally a facet only sees half of the particles (those moving towards it). So it multiplies the "seen" density by two.
		//However, in case of desorption or sticking, the real density is not twice the "seen" density, but a bit less, therefore this reduction factor
		//If only desorption, or only absorption, the correction factor is 0.5, if no des/abs, it's 1.0, and in between, see below
		double densityCorrection=1.0;
		if (f->tmpCounter[0].hit.nbMCHit > 0 || f->tmpCounter[0].hit.nbDesorbed > 0) {
			if (f->tmpCounter[0].hit.nbAbsEquiv > 0.0 || f->tmpCounter[0].hit.nbDesorbed > 0) {//otherwise save calculation time
				densityCorrection-= (f->tmpCounter[0].hit.nbAbsEquiv + (double)f->tmpCounter[0].hit.nbDesorbed) / (f->tmpCounter[0].hit.nbHitEquiv + (double)f->tmpCounter[0].hit.nbDesorbed) / 2.0;
			}
		}


		double scaleY = 1.0 / (f->getArea() * 1E-4); //1E4 is conversion from m2 to cm2
		double scaleTime=1.0/(p->counterWindowPercent * simHistory->stepSize);
		return scaleTime * scaleY * densityCorrection * GetMoleculesPerTP(hitbuffer_sum).convert_to<double>() * f->tmpCounter[0].hit.sum_1_per_ort_velocity;
	}
}

double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f){//calculates Pressure of facet. Output value's unit is mbar.
	double scaleY = 1.0 / (f->getArea()  * 1E-4)* sHandle->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
	double scaleTime=1.0/(p->counterWindowPercent * simHistory->stepSize);
	return scaleTime * scaleY * GetMoleculesPerTP(hitbuffer_sum).convert_to<double>() * f->tmpCounter[0].hit.sum_v_ort ;
}

double calcStartTime(SubprocessFacet *iFacet, bool desorbed_b, bool printWarning){

	if(desorbed_b){// if desorption
		if(p->desWindowPercent==0.0) return 0.0;

		boost::multiprecision::float128 t_start(0.0);
		boost::multiprecision::float128 rand_t=boost::multiprecision::float128(rnd());
		//if(rand_t>boost::multiprecision::float128(0.99999)){rand_t=boost::multiprecision::float128(1.0);}

		boost::multiprecision::float128 time_step = boost::multiprecision::float128(p->desWindowPercent*simHistory->stepSize);

		boost::multiprecision::float128 coverage = calcCoverage(iFacet);
		double temperature=iFacet->sh.temperature;

		boost::multiprecision::float128 tau_0=static_cast<boost::multiprecision::float128>(h/(kb*temperature));
		boost::multiprecision::float128 energy_de=static_cast<boost::multiprecision::float128>(p->E_de);
		boost::multiprecision::float128 enthalpy_vap=static_cast<boost::multiprecision::float128>(p->H_vap);
		boost::multiprecision::float128 tau_subst = tau_0 * boost::multiprecision::exp(energy_de/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the substrate

		if (coverage <= boost::multiprecision::float128(1)){
			t_start= - tau_subst * boost::multiprecision::log(boost::multiprecision::float128(1)-rand_t*(boost::multiprecision::float128(1)-boost::multiprecision::exp(-time_step/tau_subst))) ;
		}
		else{//coverage > 1
				boost::multiprecision::float128 tau_ads = tau_0 * boost::multiprecision::exp(enthalpy_vap/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the adsorbate
				if(tau_ads == tau_subst){//special case, where desorption rate is constant until a monolayer is reached.
					if ((coverage - boost::multiprecision::float128(1)) >= (time_step/tau_ads)){//There are more layers (excluding the first monolayer), than desorbing while the iteration time.
								t_start = rand_t * time_step;
					}
					else{//(coverage - 1) < (time_step/tau_ads): There are less layers (excluding the first monolayer), than desorbing while the iteration time.
						if(printWarning){	std::ostringstream tmpstream (std::ostringstream::app);
											tmpstream << "!!! Warning: Facet "<<getFacetIndex(iFacet) <<" is predicted to reach monolayer this iteration. Coverage = " <<coverage  <<" !!!" << std::endl;
											printStream(tmpstream.str());
										}
						boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(1));
						boost::multiprecision::float128 time_step_subst = time_step - time_step_ads;
						if(rand_t<(coverage-boost::multiprecision::float128(1))/(coverage-boost::multiprecision::exp(-time_step_subst/tau_subst))){
										t_start = rand_t * time_step_ads;
										}
						else{
										t_start=time_step_ads - tau_subst * boost::multiprecision::log(boost::multiprecision::float128(1)-rand_t*(boost::multiprecision::float128(1)-boost::multiprecision::exp(-time_step_subst/tau_subst))) ;
										}
					//Anyway: Here tau_ads == tau_subst! So, it is distinguished to provide a better understanding, what is happening physically
					}
				}
				else{//desorption rate is constant until coverage equals two.
					boost::multiprecision::float128 a = (-1/tau_ads)+(1/tau_subst);
					boost::multiprecision::float128 b = (1/tau_ads)-(2/tau_subst);
					boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(2));
					if((coverage - boost::multiprecision::float128(2)>=time_step/tau_ads)){//There are more layers (excluding the first two layers), than desorbing while the iteration time.
						t_start = rand_t * time_step;
						}
					else{//(coverage - 2) < (time_step/tau_ads): There are less layers (excluding the first two layers), than desorbing while the iteration time.
						boost::multiprecision::float128 time_step_mixed = time_step - time_step_ads;
						if((2+a/b)*boost::multiprecision::exp(b*time_step_mixed)-a/b>1){//case, where still more than one monolayer is left
						//desorption = coverage - boost::multiprecision::float128(2) + boost::multiprecision::float128(2) - (2+a/b)*boost::multiprecision::exp(b*time_step_mixed);
							if(rand_t<(coverage-boost::multiprecision::float128(2))/(coverage-(2+a/b)*boost::multiprecision::exp(b*time_step_mixed)+a/b)){
								t_start = rand_t * time_step_ads;
								}
							else{
								t_start=time_step_ads + (1/b)* boost::multiprecision::log((1/(2+(a/b)))*(boost::multiprecision::float128(2)-rand_t*(boost::multiprecision::float128(2)-(2+(a/b))*boost::multiprecision::exp(b*time_step_mixed)+a/b)+a/b));
								}
							}
						else{//case, where coverage will be smaller than one.
							if(printWarning){	std::ostringstream tmpstream (std::ostringstream::app);
											tmpstream << "!!! Warning: Facet "<<getFacetIndex(iFacet) <<" is predicted to reach monolayer this iteration. Coverage = " <<coverage  <<" !!!" << std::endl;
											printStream(tmpstream.str());
											}
							time_step_mixed =(1/b)* boost::multiprecision::log((1+a/b)/(2+a/b));
							boost::multiprecision::float128 time_step_subst = time_step - time_step_ads - time_step_mixed;
							if(rand_t<(coverage-boost::multiprecision::float128(1))/(coverage-boost::multiprecision::exp(-time_step_subst/tau_subst))){
								if(rand_t<(coverage-boost::multiprecision::float128(2))/(coverage-boost::multiprecision::float128(1))){
									t_start = rand_t * time_step_ads;
								}
								else{
									t_start=time_step_ads + (1/b)* boost::multiprecision::log((1/(2+(a/b)))*(boost::multiprecision::float128(2)-rand_t*(boost::multiprecision::float128(2)-(2+(a/b))*boost::multiprecision::exp(b*time_step_mixed)+a/b)+a/b));
									}
								}
							else{
							t_start=time_step_ads + time_step_mixed - tau_subst * boost::multiprecision::log(boost::multiprecision::float128(1)-rand_t*(boost::multiprecision::float128(1)-boost::multiprecision::exp(-time_step_subst/tau_subst)));
								}
							}
						}
					}
			}
		/*else{//coverage > 1
			boost::multiprecision::float128 tau_ads = tau_0 * boost::multiprecision::exp(enthalpy_vap/static_cast<boost::multiprecision::float128>(kb*temperature));//tau for particles desorbing on the adsorbate
			if ((coverage - boost::multiprecision::float128(1)) >= (time_step/tau_ads)){//There are more layers (excluding the first monolayer), than desorbing while the iteration time.
				t_start = rand_t * time_step;
			}
			else{//(coverage - 1) < (time_step/tau_ads): There are less layers (excluding the first monolayer), than desorbing while the iteration time.
				if(printWarning){
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream << "!!! Warning: Facet "<<getFacetIndex(iFacet) <<" is predicted to reach monolayer this iteration. Coverage = " <<coverage  <<" !!!" << std::endl;
					printStream(tmpstream.str());
				}

				boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(1));
				boost::multiprecision::float128 time_step_subst = time_step - time_step_ads;

				if(rand_t<(coverage-boost::multiprecision::float128(1))/(coverage-boost::multiprecision::exp(-time_step_subst/tau_subst))){
					t_start = rand_t * time_step_ads;
				}
				else{
					t_start=time_step_ads - tau_subst * boost::multiprecision::log(boost::multiprecision::float128(1)-rand_t*(boost::multiprecision::float128(1)-boost::multiprecision::exp(-time_step_subst/tau_subst))) ;
				}
			}*/
		//if(t_start.convert_to<double>()>24*3600*7)
		//	std::cout << t_start.convert_to<double>()<<"\t"<<rand_t<<std::endl;
		return t_start.convert_to<double>();
	}
	else{//if outgassing
		return rnd() * simHistory->stepSize_outgassing;
	}

}

//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
long double calcCoverage(SubprocessFacet *iFacet, Databuff *hitbuffer){ // calculates coverage depending on covering (number particles on facet)
	llong covering;
	long double coverage;

	covering = getCovering( iFacet, hitbuffer);

	coverage = (long double)covering /(long double)(calcNmono(iFacet)/calcdNsurf());
	return coverage;
}

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


//----------deprecated functions because of new K_real/virt approach

/*
double calcStep(long double variable, double start, double end, double inflection_point, double Wtr){

	if(start==end)
		return start;
	else
		return (double)tanh(((variable - inflection_point)/Wtr) *2*tuneE) * (end - start)/2 +(start+end)/2; //tanh(adjust width) * adjust height + adjust bias
}

double calcEnergy(SubprocessFacet *iFacet){ //TODO verify
	long double coverage=calcCoverage(iFacet).convert_to<long double>();
	return calcStep(coverage, p->E_de, p->H_vap,1, p->W_tr);
}

boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet) {//This returns ((d'coverage')/dt)de * (Nmono/dNSurf) * kb*T. So to speak desorption rate in units of [Pa m³/s]
	boost::multiprecision::float128 desorption = calcDesorption(iFacet);
	boost::multiprecision::float128 desorptionRate = desorption * boost::multiprecision::float128(kb* iFacet->sh.temperature * calcNmono(iFacet) /calcdNsurf());
	return desorptionRate;
}
*/
