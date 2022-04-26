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

boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns covering from hitbuffer
	boost::multiprecision::uint128_t cov = getFacetHitBuffer(iFacet,hitbuffer)->covering;
	cov = (cov.backend().size() == 0) ? boost::multiprecision::uint128_t(0) : cov; //Not necessary but used as precaution
	return cov;
}

boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet){ // returns covering from simHistory
	return simHistory->coveringList.getCurrent(iFacet);
}

boost::multiprecision::uint128_t getPredictedCovering(SubprocessFacet *iFacet){ // returns facet's covering from predictList.
	return simHistory->coveringList.getPredict(iFacet);
}

double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of hits from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->nbHitEquiv; //TODO nbMCHit or nbHitEquiv?
}

llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of desorbed testparticles from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->nbDesorbed;
}
llong getnbOutgassed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of outgassed testparticles from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->nbOutgassed;
}
llong getnbAdsorbed(SubprocessFacet *iFacet, Databuff *hitbuffer){ // returns number of adsorbed testparticles from hitbuffer
	return getFacetHitBuffer(iFacet,hitbuffer)->nbAbsEquiv;
}

llong getnbDesorbed(Databuff *hitbuffer_sum){
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)hitbuffer_sum->buff;

	return gHits->globalHits.nbDesorbed;
}

llong getnbOutgassed(Databuff *hitbuffer_sum){
	GlobalHitBuffer *gHits;
	gHits = (GlobalHitBuffer *)hitbuffer_sum->buff;

	return gHits->globalHits.nbOutgassed;
}
/*
std::tuple<double, double, double> getVelocities(SubprocessFacet *iFacet, Databuff *hitbuffer_sum){
	BYTE *buffer;
	buffer = hitbuffer_sum->buff;
	FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + iFacet->sh.hitOffset);

	return {std::make_tuple(facetHitBuffer->sum_1_per_ort_velocity,facetHitBuffer->sum_v_ort, facetHitBuffer->sum_1_per_velocity)};
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

boost::multiprecision::float128 calcPredictedCoverage(SubprocessFacet *iFacet){ //Calculates coverage from predictList
	return boost::multiprecision::float128(getPredictedCovering(iFacet))/boost::multiprecision::float128(calcNmono(iFacet));
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
	llong nbOutgassed = getnbOutgassed(hitbuffer_sum);
	if (nbDesorbed == 0 && nbOutgassed == 0) return 0; //avoid division by 0

	boost::multiprecision::float128 des=calctotalDesorption();
	CalcTotalOutgassingWorker();
	/*std::ostringstream tmpstream (std::ostringstream::app);
	tmpstream << "total desorption =  "<< des <<std::endl;
	tmpstream << "total outgassing =  "<< sHandle->wp.totalOutgassingParticles <<std::endl;
	tmpstream << "total outgassing/total desorption =  "<< sHandle->wp.totalOutgassingParticles/des <<std::endl;
	printStream(tmpstream.str());*/
	return (boost::multiprecision::float128(sHandle->wp.totalOutgassingParticles) +des) / (boost::multiprecision::float128(nbDesorbed+nbOutgassed)/boost::multiprecision::float128(simHistory->smallCoveringFactor));
	}


void calcSticking(SubprocessFacet *iFacet) {//Calculates sticking coefficient dependent on covering.

	boost::multiprecision::uint128_t covering=getCovering(iFacet);
	// sticking constant, zero if no covering
	if (covering>0){
		iFacet->sh.sticking = p->sticking;}
	else{
		iFacet->sh.sticking = 0.0;
	}

}

boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet){//This returns Delta'covering' in units of [1].
	//std::ostringstream tmpstream (std::ostringstream::app);

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
						//time while particles desorb form multilayer until one single monolayer is reached
						boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(1));
						//time while particles desorb form single monolayer
						boost::multiprecision::float128 time_step_subst = time_step - time_step_ads;
						desorption = coverage - boost::multiprecision::float128(1) + (boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step_subst/tau_subst));//This returns Delta'coverage' in units of [1]. 1 means one monolayer.
						//Anyway: Here tau_ads == tau_subst! So, it is distinguished to provide a better understanding, what is happening physically
			}
		}
		else{//desorption rate is constant until coverage equals two.
			if(coverage >= 2){
				if((coverage - boost::multiprecision::float128(2)) >= (time_step/tau_ads)){//There are more layers (excluding the first two layers), than desorbing while the iteration time.
					desorption = time_step/tau_ads;//This returns Delta'coverage' in units of [1]. 1 means one monolayer.

				}
				else{//(coverage - 2) < (time_step/tau_ads): There are less layers (excluding the first two layers), than desorbing while the iteration time.
					boost::multiprecision::float128 b = (-1/tau_ads)+(1/tau_subst);
					boost::multiprecision::float128 a = (1/tau_ads)-(2/tau_subst);
					boost::multiprecision::float128 time_step_ads = tau_ads*(coverage - boost::multiprecision::float128(2));
					boost::multiprecision::float128 time_step_mixed = time_step - time_step_ads;

					if((2+a/b)*boost::multiprecision::exp(b*time_step_mixed)-(a/b)>1){//case, where still more than one monolayer is left
						desorption = coverage - boost::multiprecision::float128(2) + boost::multiprecision::float128(2) - (2+a/b)*boost::multiprecision::exp(b*time_step_mixed)+(a/b);
						/*tmpstream << "coverage of Facet "<< getFacetIndex(iFacet) << " = "<< coverage <<std::endl;
						tmpstream << "time_step_ads of Facet "<< getFacetIndex(iFacet) << " = "<< time_step_ads <<std::endl;
						tmpstream << "time_step of Facet "<< getFacetIndex(iFacet) << " = "<< time_step <<std::endl;
						tmpstream << "time_step_mixed of Facet "<< getFacetIndex(iFacet) << " = "<< time_step_mixed <<std::endl;
						tmpstream << "Desorption of facet "<< getFacetIndex(iFacet) <<" will be so low, that more than one monolayer will be left. " <<std::endl;
						tmpstream << "tau_subst"<<  " = "<< tau_subst <<std::endl;
						tmpstream << "tau_ads"<<  " = "<< tau_ads <<std::endl;
						tmpstream << "a"<<  " = "<< a <<std::endl;
						tmpstream << "b"<<  " = "<< b <<std::endl;
						tmpstream << "a/b"<<  " = "<< a/b <<std::endl;
						tmpstream << "(2+a/b)*boost::multiprecision::exp(b*time_step_mixed) - (a/b) of Facet "<< getFacetIndex(iFacet) << " = "<< (2+a/b)*boost::multiprecision::exp(b*time_step_mixed)-(a/b) <<std::endl;
						*/}
					else{//case, where coverage will be smaller than one.
						time_step_mixed =(1/b)* boost::multiprecision::log((1+a/b)/(2+a/b));
						boost::multiprecision::float128 time_step_subst = time_step - time_step_ads - time_step_mixed;
						desorption = coverage - boost::multiprecision::float128(1) + (boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step_subst/tau_subst));
						//tmpstream << "Desorption of facet "<< getFacetIndex(iFacet) <<" will be so high, that less than one monolayer will be left. " <<std::endl;
						}
				}
			}
			else{// 1 <= coverage <= 2
				boost::multiprecision::float128 b = (-1/tau_ads)+(1/tau_subst);
				boost::multiprecision::float128 a = (1/tau_ads)-(2/tau_subst);
				if((2+a/b)*boost::multiprecision::exp(b*time_step)-(a/b)>1){//case, where still more than one monolayer is left
					desorption = coverage - boost::multiprecision::float128(2) + boost::multiprecision::float128(2) - (2+a/b)*boost::multiprecision::exp(b*time_step)+(a/b);
					//tmpstream << "Desorption of facet "<< getFacetIndex(iFacet) <<" will be so low, that more than one monolayer will be left. " <<std::endl;
					//tmpstream << "(2+a/b)*boost::multiprecision::exp(b*time_step) - (a/b) of Facet "<< getFacetIndex(iFacet) << " = "<< (2+a/b)*boost::multiprecision::exp(b*time_step)-(a/b) <<std::endl;
					}
				else{//case, where coverage will be smaller than one.
					boost::multiprecision::float128 time_step_mixed =(1/b)* boost::multiprecision::log((1+a/b)/(2+a/b));
					boost::multiprecision::float128 time_step_subst = time_step - time_step_mixed;
					desorption = coverage - boost::multiprecision::float128(1) + (boost::multiprecision::float128(1) - boost::multiprecision::exp(-time_step_subst/tau_subst));
					//tmpstream << "Desorption of facet "<< getFacetIndex(iFacet) <<" will be so high, that less than one monolayer will be left. " <<std::endl;
					}
			}
		}
	}
	//tmpstream << "Desorption [layers] of facet "<< getFacetIndex(iFacet) <<" = " << desorption <<std::endl;
	//tmpstream << "Desorption [particles] of facet "<< getFacetIndex(iFacet) <<" = " << desorption * boost::multiprecision::float128(calcNmono(iFacet)) <<std::endl;
	//printStream(tmpstream.str());
	if(desorption <0){
		desorption = 0;
	}
	return desorption * boost::multiprecision::float128(calcNmono(iFacet));//This returns Delta'covering' in units of [1]. 1 means one particle.
}

double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f){

	if(f->tmpCounter[0].sum_1_per_ort_velocity==std::numeric_limits<double>::infinity()){
		return 0.0;
		}
	else{
		//Correction for double-density effect (measuring density on desorbing/absorbing facets):

		//Normally a facet only sees half of the particles (those moving towards it). So it multiplies the "seen" density by two.
		//However, in case of desorption or sticking, the real density is not twice the "seen" density, but a bit less, therefore this reduction factor
		//If only desorption, or only absorption, the correction factor is 0.5, if no des/abs, it's 1.0, and in between, see below
		double densityCorrection=1.0;
		if (f->tmpCounter[0].nbMCHit > 0 || f->tmpCounter[0].nbDesorbed > 0 || f->tmpCounter[0].nbOutgassed > 0) {
			if (f->tmpCounter[0].nbAbsEquiv > 0.0 || f->tmpCounter[0].nbDesorbed > 0 || f->tmpCounter[0].nbOutgassed > 0) {//otherwise save calculation time
				densityCorrection-= (f->tmpCounter[0].nbAbsEquiv + (double)f->tmpCounter[0].nbDesorbed  + (double)f->tmpCounter[0].nbOutgassed) / (f->tmpCounter[0].nbHitEquiv + (double)f->tmpCounter[0].nbDesorbed+ (double)f->tmpCounter[0].nbOutgassed) / 2.0;
			}
		}


		double scaleY = 1.0 / (f->getArea() * 1E-4); //1E4 is conversion from m2 to cm2
		double scaleTime=1.0/(p->counterWindowPercent * simHistory->stepSize);
		return scaleTime * scaleY * densityCorrection * GetMoleculesPerTP(hitbuffer_sum).convert_to<double>() * f->tmpCounter[0].sum_1_per_ort_velocity;
	}
}

double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f){//calculates Pressure of facet. Output value's unit is mbar.
	double scaleY = 1.0 / (f->getArea()  * 1E-4)* sHandle->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar;  //1E4 is conversion from m2 to cm2, 0.01: Pa->mbar
	double scaleTime=1.0/(p->counterWindowPercent * simHistory->stepSize);
	std::cout << "f->tmpCounter[0].sum_v_ort of Facet " << f->tmpCounter[0].sum_v_ort <<std::endl;
	return scaleTime * scaleY * GetMoleculesPerTP(hitbuffer_sum).convert_to<double>() * f->tmpCounter[0].sum_v_ort ;
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
					boost::multiprecision::float128 b = (-1/tau_ads)+(1/tau_subst);
					boost::multiprecision::float128 a = (1/tau_ads)-(2/tau_subst);
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
		//regular outgassing
		if (iFacet->sh.outgassing_paramId == -1) { //constant outgassing
			return rnd() * simHistory->stepSize_outgassing;
			}
		else {
			//time-dependent outgassing
			double time_step = simHistory->stepSize; //length of the current iteration
			double t_start =simHistory->lastTime; //start of iteration
			double t_stop =t_start + time_step;	//end of iteration
			double end_of_outgassing = sHandle->IDs[iFacet->sh.IDid].back().first; //last point of the defined and loaded outgassing table
			double start_of_outgassing = sHandle->IDs[iFacet->sh.IDid].front().first; //first point of the defined and loaded outgassing table
			double outgassing_start = 0;//start of outgassing within the iteration step
			double outgassing_end = 0;//end of outgassing within the iteration step
			double facet_outgassing = 0;//over time integrated outgassing of facet during the iteration step
			double outgassing_before_step = 0;
			if(t_start >= end_of_outgassing){//case, when t_start is after the last point in time, where an outgassing is defined
				facet_outgassing =0;
			}
			else if(t_start <= end_of_outgassing && t_stop >= end_of_outgassing){//case, when t_start is before and
				//t_stopp is after the last point in time, where an outgassing is defined
					if (t_start <= start_of_outgassing){
						outgassing_start = start_of_outgassing;
						outgassing_end = end_of_outgassing;
					}
					else{
						outgassing_start = t_start;
						outgassing_end = end_of_outgassing;
						}
					facet_outgassing = InterpolateY(outgassing_end, sHandle->IDs[iFacet->sh.IDid], false, true) - InterpolateY(outgassing_start, sHandle->IDs[iFacet->sh.IDid], false, true);
					}
			else{
				//same as =>else if(t_start <= end_of_outgassing && t_stop <= end_of_outgassing){//case, when t_start is before and
				//t_stopp is before the last point in time, where an outgassing is defined
				if (t_start <= start_of_outgassing){
					if(t_stop <= start_of_outgassing){
						facet_outgassing =0;
					}
					else{//t_stop >= start_of_outgassing
						outgassing_start = start_of_outgassing;
						outgassing_end = t_stop;
						facet_outgassing = InterpolateY(outgassing_end, sHandle->IDs[iFacet->sh.IDid], false, true) - InterpolateY(outgassing_start, sHandle->IDs[iFacet->sh.IDid], false, true);
					}
				}
				else{//t_start >= start_of_outgassing
					outgassing_start = t_start;
					outgassing_end = t_stop;
					facet_outgassing = InterpolateY(outgassing_end, sHandle->IDs[iFacet->sh.IDid], false, true) - InterpolateY(outgassing_start, sHandle->IDs[iFacet->sh.IDid], false, true);
				}
			}
			outgassing_before_step = InterpolateY(outgassing_start, sHandle->IDs[iFacet->sh.IDid], false, true);
			return InterpolateX(outgassing_before_step + rnd() * facet_outgassing, sHandle->IDs[iFacet->sh.IDid], false, true);
		}
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
