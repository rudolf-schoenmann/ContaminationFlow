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
 * This file contains the functions that adds up Databuff structs
 */

#include "SimulationLinux.h"
#include "GLApp/MathTools.h"
#include <math.h>

extern Simulation *sHandle;
extern ProblemDef* p;
extern SimulationHistory* simHistory;

// Step size for intervals
double getStepSize(){
	double T_min = p->t_min;//set minimal time resolution to 1E-4 seconds.

	//Dynamical calculation of min_time is not straight forward, since 'manageTimeStep()' can change it.
	//Dynamical calculation can be done later, if it is regarded as useful.
	double t_start = T_min*exp((double)simHistory->currentStep*(log(p->maxTimeS/(T_min))/(double)p->iterationNumber));
	double t_stop = T_min*exp((double)(simHistory->currentStep+1)*(log(p->maxTimeS/(T_min))/(double)p->iterationNumber));
	double Delta = t_stop - t_start;
	double Delta_final = p->t_max - t_start;
	if(t_stop > p->t_max){
		return Delta_final;
	}
	else{
		return Delta;
	}
	/*if(simHistory->currentStep==0){
		double T_min = estimateAverageFlightTime();
		return T_min*exp((double)simHistory->currentStep*(log(p->maxTimeS/T_min)/(double)p->iterationNumber));
		}
	else{
		return exp((double)simHistory->currentStep*(log(p->maxTimeS/simHistory->coveringList.pointintime_list[1].first)/(double)p->iterationNumber));
		}*/

}

double manageStepSize(){
	double step_size;
	bool steptoolong; //If step is too long, we will decrease the step size.
	bool needforCheck = true;

	while(needforCheck){
		step_size = getStepSize();
		steptoolong=false;
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				if(f.sh.desorption==0||f.sh.temperature==0)continue;

				boost::multiprecision::uint128_t covering_phys = simHistory->coveringList.getLast(&f);

				if ((boost::multiprecision::uint128_t)((f.sh.desorption/boost::multiprecision::float128(kb* f.sh.temperature))*boost::multiprecision::float128(step_size) +boost::multiprecision::float128(0.5))>covering_phys){
					steptoolong=true;
				}
			}
		}
		if(steptoolong){
			if(simHistory->currentStep == 0){
				needforCheck = false;
			}
			else{
			simHistory->currentStep-=1;
			std::cout<<"Decrease simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
			p->outFile<<"Decrease simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
			}
		}
		else {
			needforCheck = false;
		}
	}
	return step_size;
}
//-----------------------------------------------------------
//Update Covering
//Simhistory version
void UpdateCovering(Databuff *hitbuffer_sum, llong smallCoveringFactor){//Updates Covering after one Iteration using Krealvirt,
	//the current time step (of the iteration step) and the smallCoveringFactor,
	//Calculates with the summed up counters of hitbuffer_sum how many test particles are equivalent to one physical particle.
	//simTime in ms

	boost::multiprecision::float128 Krealvirt = GetMoleculesPerTP(hitbuffer_sum);
	//llong nbDesorbed = getnbDesorbed(hitbuffer_sum)-simHistory->nbDesorbed_old;
	//std::cout <<"nbDesorbed before and after:\t" << history->nbDesorbed_old <<'\t';
	//simHistory->nbDesorbed_old = getnbDesorbed(hitbuffer_sum); //Not needed anymore.
	//std::cout << history->nbDesorbed_old <<std::endl;

	boost::multiprecision::uint128_t covering_phys;
	boost::multiprecision::uint128_t covering_sum;//covering as it is summed up of all subprocesses. In case, it is multiplied by smallCoveringFactor
	boost::multiprecision::float128 covering_sum_netto ;//used to devide covering_sum by the smallCoveringFactor; float;
	//boost::multiprecision::float128 covering_check;

	double time_step = getStepSize();
	if(Krealvirt==boost::multiprecision::float128(0.0)){ //if no Krealvirt(no desorption), do not increase currentStep, time_step=0
		time_step=0;
	}
	else{
		double error_event=0.0;
		double error_covering=0.0;
		double area=0.0;

		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				if(simHistory->errorList_event.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0 || f.sh.isVipFacet)//ignore facet if no hits (=inf error)
					continue;
				if(simHistory->errorList_covering.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0 || f.sh.isVipFacet)//ignore facet if no hits (=inf error)
					continue;

				error_event+=simHistory->errorList_event.getCurrent(&f)*f.sh.area;
				error_covering+=simHistory->errorList_covering.getCurrent(&f)*f.sh.area;
				area+=f.sh.area;
			}
		}
		// Print total error and error per facet of this iteration
		std::ostringstream tmpstream (std::ostringstream::app);
		tmpstream <<"Total Error (event) averaged over facets "<<error_event/area <<std::endl;
		simHistory->errorList_event.printCurrent(tmpstream);
		tmpstream << std::endl<<"Total Error (covering) averaged over facets "<<error_covering/area <<std::endl;
		simHistory->errorList_covering.printCurrent(tmpstream);

		if(!p->vipFacets.empty()){
			tmpstream <<"Vip Facets:"<<std::endl;
			for(unsigned int i = 0; i < p->vipFacets.size(); i++){
				tmpstream <<"\t"<<p->vipFacets[i].first <<"\t" << simHistory->errorList_event.getCurrent(p->vipFacets[i].first)<<std::endl;
			}
			tmpstream <<std::endl;
		}
		std::cout <<tmpstream.str();
		p->outFile<<tmpstream.str();


		//if targetError not reached: do not update currentstep
		if(checkErrorSub(p->targetError, error_covering/area, 1.0))
			{simHistory->currentStep += 1;}

	}

	std::cout <<"Krealvirt = " << Krealvirt << std::endl;
	std::cout << "Covering difference will be multiplied by Krealvirt*(time step): " << Krealvirt*boost::multiprecision::float128(time_step) << std::endl;

	p->outFile <<"Krealvirt = " << Krealvirt << std::endl;
	p->outFile << "Covering difference will be multiplied by Krealvirt*(time step): " << Krealvirt*boost::multiprecision::float128(time_step) << std::endl;
	//std::cout <<"testing timestep: " <<time_step <<'\t' <<estimateAverageFlightTime() <<std::endl;

	//double rounding=1/simHistory->numFacet;
	boost::multiprecision::float128 rounding(0.5);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {

				covering_phys = simHistory->coveringList.getLast(&f);
				covering_sum = boost::multiprecision::uint128_t(getCovering(&f, hitbuffer_sum));
				//std::cout << "covering_sum_brutto " << covering_sum << std::endl;
				//p->outFile << "covering_sum_brutto " << covering_sum << std::endl;
				covering_sum_netto = boost::multiprecision::float128( (static_cast < boost::multiprecision::float128 >(covering_sum))/smallCoveringFactor);
				//std::cout << "covering_sum_netto " << covering_sum_netto << std::endl;
				//p->outFile << "covering_sum_netto " << covering_sum_netto << std::endl;
				covering_sum = boost::multiprecision::uint128_t (static_cast <boost::multiprecision::uint128_t>(covering_sum_netto));

				std::cout<<std::endl << "Facet " << getFacetIndex(&f)<< std::endl;
				std::cout << "covering_sum = " << covering_sum  << " = "<< boost::multiprecision::float128(covering_sum) << std::endl;
				std::cout<< "covering_phys_before = " << covering_phys << " = "<< boost::multiprecision::float128(covering_phys) << std::endl;

				p->outFile<<std::endl << "Facet " << getFacetIndex(&f) << std::endl;
				p->outFile<< "covering_sum = " << covering_sum  << " = "<< boost::multiprecision::float128(covering_sum) << std::endl;
				p->outFile<< "covering_phys_before = " << covering_phys << " = "<< boost::multiprecision::float128(covering_phys) << std::endl;

				if (covering_sum > covering_phys){
					//+0.5 for rounding
					boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding + boost::multiprecision::float128(covering_sum - covering_phys)*Krealvirt*boost::multiprecision::float128(time_step));
					covering_phys += covering_delta;
					std::cout << "covering rises by " <<covering_delta << " = "<<boost::multiprecision::float128(covering_delta) << std::endl;
					p->outFile << "covering rises by " <<covering_delta << " = "<<boost::multiprecision::float128(covering_delta) << std::endl;
				}
				else{
					//+0.5 for rounding
					boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding + boost::multiprecision::float128(covering_phys - covering_sum)*Krealvirt*boost::multiprecision::float128(time_step));
					if(covering_phys+1==covering_delta){
						std::cout<<"!!!-----Correct covering_delta by 1: "<<covering_phys<<" - "<<covering_delta <<" because of rounding-----!!!"<<std::endl;
						p->outFile<<"!!!-----Correct covering_delta by 1: "<<covering_phys<<" - "<<covering_delta <<" because of rounding-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}
					else if(covering_phys<covering_delta && covering_phys+10>=covering_delta){
						std::cout<<"!!!-----Correct covering_delta: "<<covering_phys<<" - "<<covering_delta <<" because of bad statistic-----!!!"<<std::endl;
						p->outFile<<"!!!-----Correct covering_delta: "<<covering_phys<<" - "<<covering_delta <<" because of bad statistic-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}
					else if(covering_phys<covering_delta){
						std::cout<<"!!!-----Covering gets negative: "<<covering_phys<<" - "<<covering_delta <<". Correct Covering=0-----!!!"<<std::endl;
						p->outFile<<"!!!-----Covering gets negative: "<<covering_phys<<" - "<<covering_delta <<". Correct Covering=0-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}

					covering_phys -= covering_delta;
					std::cout << "covering decreases by "<<covering_delta << " = " << boost::multiprecision::float128(covering_delta) << std::endl;
					p->outFile << "covering decreases by "<<covering_delta << " = " << boost::multiprecision::float128(covering_delta) << std::endl;

				}
				std::cout<< "covering_phys_after = " << covering_phys << " = " << boost::multiprecision::float128(covering_phys) << std::endl;
				std::cout<< "coveringThreshhold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << " = " << boost::multiprecision::float128(sHandle->coveringThreshold[getFacetIndex(&f)]) << std::endl;
				p->outFile<< "covering_phys_after = " << covering_phys << " = " << boost::multiprecision::float128(covering_phys) << std::endl;
				p->outFile<< "coveringThreshhold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << " = " << boost::multiprecision::float128(sHandle->coveringThreshold[getFacetIndex(&f)]) << std::endl;

				simHistory->coveringList.setCurrent(&f, covering_phys);
		}
	}

	simHistory->coveringList.appendCurrent(simHistory->lastTime+time_step);
	simHistory->lastTime+=time_step;
	simHistory->hitList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->errorList_event.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->errorList_covering.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->desorbedList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->stepSize=time_step;//For what do I need this? Maybe could be uncommented...
}

void UpdateErrorMain(Databuff *hitbuffer_sum){

	double num_hit_it=0;
	double num_des_ad_it =0;

	//count all hits and all desorption events for all facets in num_hit_it
	//count all adsorption and all desorption events for all facets in num_des_ad_it
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			/*
			num_hit_it+=f.sh.opacity * (getHits(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) ); // I think, we can replace the 'f.sh.opacity' by 1,0 or typecasting
			num_des_ad_it += f.sh.opacity * (getnbAdsorbed(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );// In case of a opacity being not 1 hits, adsorbs and desorbs happen
			//randomly and therefore counters will be raised or not. So there should be no need for multiplying with this factor 'f.sh.opacity' here.
			 */
			num_hit_it+= (double) (getHits(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );
			num_des_ad_it += (double)(getnbAdsorbed(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );
		}
	}

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			//double num_hit_f=f.sh.opacity * ( getHits(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );
			double num_hit_f= (double)( getHits(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );
			//double num_des_ad_f=f.sh.opacity * ( getHits(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );
			double num_des_ad_f= (double)( getnbAdsorbed(&f,hitbuffer_sum) + getnbDesorbed(&f, hitbuffer_sum) );

			if(num_hit_f/num_hit_it<p->hitRatioLimit){// threshold. If reached, small number of hits neglected
				num_hit_it-=num_hit_f;
				num_hit_f=0;
			}
			if(num_des_ad_f/num_des_ad_it<p->hitRatioLimit){// threshold. If reached, small number of hits neglected
				num_des_ad_it-=num_des_ad_f;
				num_des_ad_f=0;
			}

			if(f.sh.opacity==0){
				simHistory->errorList_event.setCurrent(&f, 0.0);
				simHistory->errorList_covering.setCurrent(&f, 0.0);
			}
			else{
				double error_event=pow((1/num_hit_f)*(1-num_hit_f/num_hit_it),0.5);
				simHistory->errorList_event.setCurrent(&f, error_event);
				double error_covering=pow((1/num_des_ad_f)*(1-num_des_ad_f/num_des_ad_it),0.5);
				simHistory->errorList_covering.setCurrent(&f, error_covering);
			}

			simHistory->hitList.setLast(&f,getHits(&f,hitbuffer_sum));
			simHistory->desorbedList.setLast(&f,getnbDesorbed(&f,hitbuffer_sum));
		}
	}

	simHistory->errorList_event.appendCurrent(simHistory->lastTime);
	simHistory->errorList_covering.appendCurrent(simHistory->lastTime);
	//simHistory->hitList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateCovering before UpdateErrorMain
	//simHistory->desorbedList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateCovering before UpdateErrorMain

}

std::tuple<std::vector<double>,std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration(){//calculates statistical uncertainties of error_event and error_covering at the end of
	// the simulation for writigin these in the output file. While simulating only the error values of the subprocesses are used to decide, if the targeted error level has been reached.
	std::vector<double> errorPerIt_event;
	errorPerIt_event =std::vector<double> ();
	std::vector<double> errorPerIt_covering;
	errorPerIt_covering =std::vector<double> ();

	std::vector<boost::multiprecision::uint128_t> covPerIt;
	covPerIt =std::vector<boost::multiprecision::uint128_t> ();

	for(unsigned int it=0; it<simHistory->errorList_event.pointintime_list.size();it++){
		// Total error/covering for each iteration
		double error=0.0;
		double area=0.0;
		boost::multiprecision::uint128_t covering=0;

		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				int idx=getFacetIndex(&f);
				covering+=simHistory->coveringList.pointintime_list[it].second[idx];

				double err=simHistory->errorList_event.pointintime_list[it].second[idx];
				if(err== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
					continue;

				error+=err*f.sh.area;
				area+=f.sh.area;
			}
		}
		errorPerIt_event.push_back(error/area);
		covPerIt.push_back(covering);
	}

	for(unsigned int it=0; it<simHistory->errorList_covering.pointintime_list.size();it++){
			// Total error/covering for each iteration
			double error=0.0;
			double area=0.0;
			boost::multiprecision::uint128_t covering=0;

			for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
				for (SubprocessFacet& f : sHandle->structures[s].facets) {
					int idx=getFacetIndex(&f);
					covering+=simHistory->coveringList.pointintime_list[it].second[idx];

					double err=simHistory->errorList_covering.pointintime_list[it].second[idx];
					if(err== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
						continue;

					error+=err*f.sh.area;
					area+=f.sh.area;
				}
			}
			errorPerIt_covering.push_back(error/area);
		}

	return std::make_tuple(errorPerIt_event,errorPerIt_covering,covPerIt);
}

// Copy covering to buffer
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer){
	boost::multiprecision::uint128_t covering_phys;
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;

	BYTE *buffer;
	buffer = hitbuffer->buff;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				FacetHitBuffer *facetHitSum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_phys = simHistory->coveringList.getLast(&f);
				facetHitBuffer->hit.covering=covering_phys.convert_to<llong>();
				facetHitSum->hit.covering=covering_phys.convert_to<llong>();
			}
	}

	simHistory->flightTime=0.0;
	simHistory->nParticles=0;
}

//-----------------------------------------------------------
// Update Buffer of main process using buffer of sub process
//Simhistory version

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank, llong smallCoveringFactor) {
	BYTE *buffer, *subbuff;
	GlobalHitBuffer *gHits, *subHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif


	size_t nbMoments = (size_t)sHandle->moments.size();

	buffer = mainbuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	//added subbuffer that contains simulation results from a subprocess, to be added to mainbuffer
	subbuff=subbuffer->buff;
	subHits=(GlobalHitBuffer *)subbuff;

/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl <<std::endl;

	std::cout <<subHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<subHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<subHits->globalHits.hit.covering <<std::endl;
	std::cout <<subHits->distTraveled_total  <<std::endl;
	std::cout <<subHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;*/

	// Global hits and leaks: adding local hits to shared memory
	sHandle->tmpGlobalResult.globalHits.hit.nbMCHit= gHits->globalHits.hit.nbMCHit += subHits->globalHits.hit.nbMCHit;
	sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv=gHits->globalHits.hit.nbHitEquiv += subHits->globalHits.hit.nbHitEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv=gHits->globalHits.hit.nbAbsEquiv += subHits->globalHits.hit.nbAbsEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed=gHits->globalHits.hit.nbDesorbed += subHits->globalHits.hit.nbDesorbed;
	sHandle->tmpGlobalResult.distTraveled_total=gHits->distTraveled_total += subHits->distTraveled_total;
	sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly=gHits->distTraveledTotal_fullHitsOnly += subHits->distTraveledTotal_fullHitsOnly;

	/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;
	std::cout <<gHits->hitCacheSize <<std::endl;*/

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	//sHandle->wp.sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < subHits->leakCacheSize; leakIndex++)
		sHandle->tmpGlobalResult.leakCache[leakIndex]=gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = subHits->leakCache[leakIndex];
	sHandle->tmpGlobalResult.nbLeakTotal=gHits->nbLeakTotal += subHits->nbLeakTotal;
	sHandle->tmpGlobalResult.lastLeakIndex=gHits->lastLeakIndex = (gHits->lastLeakIndex + subHits->leakCacheSize) % LEAKCACHESIZE;
	sHandle->tmpGlobalResult.leakCacheSize=gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + subHits->leakCacheSize);


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	if (rank == 0) {
			for (size_t hitIndex = 0; hitIndex < subHits->hitCacheSize; hitIndex++)
				sHandle->tmpGlobalResult.hitCache[hitIndex]=gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = subHits->hitCache[hitIndex];

			if (subHits->hitCacheSize > 0) {
				sHandle->tmpGlobalResult.lastHitIndex = gHits->lastHitIndex = (gHits->lastHitIndex + subHits->hitCacheSize) % HITCACHESIZE;
				sHandle->tmpGlobalResult.hitCache[gHits->lastHitIndex].type = gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
				sHandle->tmpGlobalResult.hitCacheSize = gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + subHits->hitCacheSize);
			}
		}
	/*
	std::cout <<gHits->hitCacheSize <<std::endl;
	std::cout <<gHits->leakCacheSize <<std::endl;
	std::cout <<gHits->nbLeakTotal <<std::endl<<std::endl;*/


	//Global histograms
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {
		BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
		BYTE *subhist = subbuff + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
		if (sHandle->wp.globalHistogramParams.recordBounce) {
			double* nbHitsHistogram = (double*)histCurrentMoment;
			double* nbHitsSub=(double*)subhist;
			for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
				sHandle->tmpGlobalHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
			}
		}

		if (sHandle->wp.globalHistogramParams.recordDistance) {
			double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			double* distanceSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
				sHandle->tmpGlobalHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
			}
		}
		if (sHandle->wp.globalHistogramParams.recordTime) {
			double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
			double* timeSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
			for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
				sHandle->tmpGlobalHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
			}
		}

	}



	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	//int num=0;
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {

		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//if (f.hitted) {
				//std::cout <<"Facet " <<num <<std::endl; // Da wird immer "0" angezeigt. Was ist der Sinn?
				for (unsigned int m = 0; m < (1 + nbMoments); m++) { // Add hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitSub = (FacetHitBuffer *)(subbuff + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					llong covering_phys= smallCoveringFactor * history->coveringList.getLast(&f).convert_to<llong>();
					llong covering_sum = facetHitSub->hit.covering;
/*
					std::cout <<sizeof(GlobalHitBuffer) <<std::endl;
					std::cout <<f.sh.hitOffset  <<std::endl;
					std::cout <<sizeof(FacetHitBuffer) <<std::endl;
					std::cout <<facetHitSub->hit.covering <<std::endl;*/

					/*if(f.globalId == 1){
					std::cout <<"facet number\t\t\t"<< f.globalId <<std::endl;
					std::cout <<"buffer before\t\t\t" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/

					f.tmpCounter[m].hit.nbAbsEquiv = facetHitBuffer->hit.nbAbsEquiv += facetHitSub->hit.nbAbsEquiv;
					f.tmpCounter[m].hit.nbDesorbed = facetHitBuffer->hit.nbDesorbed += facetHitSub->hit.nbDesorbed;
					f.tmpCounter[m].hit.nbMCHit = facetHitBuffer->hit.nbMCHit += facetHitSub->hit.nbMCHit;
					f.tmpCounter[m].hit.nbHitEquiv = facetHitBuffer->hit.nbHitEquiv += facetHitSub->hit.nbHitEquiv;;
					f.tmpCounter[m].hit.sum_1_per_ort_velocity = facetHitBuffer->hit.sum_1_per_ort_velocity += facetHitSub->hit.sum_1_per_ort_velocity;
					f.tmpCounter[m].hit.sum_v_ort = facetHitBuffer->hit.sum_v_ort += facetHitSub->hit.sum_v_ort;
					f.tmpCounter[m].hit.sum_1_per_velocity = facetHitBuffer->hit.sum_1_per_velocity += facetHitSub->hit.sum_1_per_velocity;
					//facetHitBuffer->hit.covering += facetHitSub->hit.covering; //We do that in another way.

					if (facetHitSub->hit.covering > covering_phys){
					facetHitBuffer->hit.covering += (covering_sum - covering_phys);
					}
					else{
						if(facetHitBuffer->hit.covering > (covering_phys - covering_sum)){
							facetHitBuffer->hit.covering -= (covering_phys- covering_sum);
							}
						else{
							facetHitBuffer->hit.covering = 0;//Counter cannot be negative! Maybe we could interrupt the iteration here?
						}
					}
					f.tmpCounter[m].hit.covering = facetHitBuffer->hit.covering;

					/*
					if(f.globalId == 1){
					std::cout <<"buffer afterwards" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/
				}

				if (f.sh.isProfile) { //(MY) Add profiles
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						ProfileSlice *shProfileSub = (ProfileSlice *)(subbuff + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							f.profile[m][j] = shProfile[j] += shProfileSub[j];
						}
					}
				}

				if (f.sh.isTextured) {// Add texture
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						TextureCell *shTextureSub = (TextureCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)
						//std::cout <<"timecorrection: " <<timeCorrection <<std::endl<<std::endl;

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Add temporary hit counts
								f.texture[m][add] = shTexture[add] += shTextureSub[add];

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f.textureCellIncrements[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f.largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}
							}
						}
					}
				}

				if (f.sh.countDirection) {
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						DirectionCell *shDirSub = (DirectionCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;
								f.direction[m][add].dir.x = shDir[add].dir.x += shDirSub[add].dir.x;
								f.direction[m][add].dir.y = shDir[add].dir.y += shDirSub[add].dir.y;
								f.direction[m][add].dir.z = shDir[add].dir.z += shDirSub[add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								f.direction[m][add].count = shDir[add].count += shDirSub[add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					size_t *shAngleMapSub = (size_t *)(subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							f.angleMap.pdf[add] = shAngleMap[add] += shAngleMapSub[add];
						}
					}
				}

				//Facet histograms

				for (unsigned int m = 0; m < (1 + nbMoments); m++) {
					BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
					BYTE *histSub = subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
					if (f.sh.facetHistogramParams.recordBounce) {
						double* nbHitsHistogram = (double*)histCurrentMoment;
						double* nbHitsSub = (double*)histSub;
						for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
							f.tmpHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
						}
					}
					if (f.sh.facetHistogramParams.recordDistance) {
						double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
						double* distanceSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize());
						for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
							f.tmpHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
						}
					}
					if (f.sh.facetHistogramParams.recordTime) {
						double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
						double* timeSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
						for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
							f.tmpHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
						}
					}
				}

			//} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
		sHandle->tmpGlobalResult.texture_limits[v]=gHits->texture_limits[v];
	}

	//ResetTmpCounters(); //This will reset counters -> but needed for Tmin

	//extern char* GetSimuStatus();
	//SetState(NULL, GetSimuStatus(), false, true); // (Rudi) Don't need that.

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}
