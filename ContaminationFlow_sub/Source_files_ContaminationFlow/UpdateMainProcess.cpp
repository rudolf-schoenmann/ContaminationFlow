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

#include "SimulationContaminationFlow.h"
#include "GLApp/MathTools.h"
#include <cmath>

extern Simulation *sHandle;
extern ProblemDef* p;
extern SimulationHistory* simHistory;

// Step size for intervals
double getStepSize(){
	double t_min = p->t_min;
	double t_max = p->maxTimeS;
	double last_step_size=0;
	double test_step_size=0;
	double iterationNumber = (double)p->iterationNumber - 1;// steps: 0, 1, 2, ..., 'p->iterationNumber - 1' => these are 'p->iterationNumber' steps
		if(simHistory->currentStep == 0){
			return t_min;
		}
		else{
			last_step_size=simHistory->stepSize;//This is the old stepsize, which will be updated here.
			//exponentially growing points in time (t_i);time step = t_(i+1)- t_i
			double t_start = t_min*exp((double)(simHistory->currentStep-1)*(log(t_max/t_min)/iterationNumber));
			double t_stop = t_min*exp((double)simHistory->currentStep*(log(t_max/t_min)/iterationNumber));
			test_step_size = t_stop - t_start;
			for(int i = 40; i > 0; i--){
				if(test_step_size<last_step_size){
					//ith order polynomial growing points in time (t_i);time step = t_(i+1)- t_i
					t_start = t_min + (t_max - t_min)*pow(((simHistory->currentStep-1)/iterationNumber),i*0.25);
					t_stop = t_min + (t_max - t_min)*pow((simHistory->currentStep/iterationNumber),i*0.25);
					test_step_size = t_stop - t_start;
				}
				if(i == 1 && test_step_size<last_step_size){
					return t_min;
				}
			}
			return t_stop-t_start<p->t_max?t_stop-t_start:p->t_max;
		}
	}

//-----------------------------------------------------------
//Update Covering
//Simhistory version

void UpdateCovering(Databuff *hitbuffer_sum){//Updates Covering after an Iteration using Krealvirt and the smallCoveringFactor.
	//Calculates with the summed up counters of hitbuffer_sum, how many test particles are equivalent to one physical particle.

	boost::multiprecision::float128 Krealvirt = GetMoleculesPerTP(hitbuffer_sum);
	boost::multiprecision::uint128_t covering_phys;
	boost::multiprecision::uint128_t covering_sum;//covering as it is summed up of all subprocesses. In case, it is multiplied by smallCoveringFactor
	std::ostringstream tmpstream (std::ostringstream::app);

	if (!(p->usePCMethod==1&&simHistory->pcStep == 0)) {
		// Calculate total error of this iteration
		double total_error_event=0.0;
		double total_error_covering=0.0;
		std::tie(total_error_event,total_error_covering)=CalcErrorAll();

		std::string monitoredFacets=(!p->doFocusGroupOnly||p->focusGroup.second.size()==simHistory->numFacet)?"all":"selected";

		// Print total error and error per facet of this iteration
		tmpstream <<"Target Error (only "+p->errorMode+" monitored) "<<p->targetError <<std::endl <<std::endl;
		tmpstream <<"Total Error (event) averaged over "+monitoredFacets+" facets "<<total_error_event <<std::endl;
		simHistory->errorList_event.printCurrent(tmpstream);
		tmpstream << std::endl<<"Total Error (covering) averaged over "+monitoredFacets+" facets "<<total_error_covering<<std::endl;
		simHistory->errorList_covering.printCurrent(tmpstream);

		HistoryList<double> *listptr;
		listptr = getErrorList(p->errorMode);
		if(!p->vipFacets.empty()){
			tmpstream <<"Vip Facets:"<<std::endl;
			for(unsigned int i = 0; i < p->vipFacets.size(); i++){
				tmpstream <<"\t"<<p->vipFacets[i].first <<"\t" << listptr->getCurrent(p->vipFacets[i].first)<<std::endl;
			}
			tmpstream <<std::endl;
		}

		//if targetError not reached: do not update currentstep
		if(checkError(p->targetError, (p->errorMode=="covering")?total_error_covering:total_error_event, 1.0, p->errorMode)){simHistory->currentStep += 1;}

	}
	tmpstream <<"Krealvirt = " << Krealvirt << std::endl;
	tmpstream << "Covering difference will be multiplied by Krealvirt: " << Krealvirt << std::endl;
	printStream(tmpstream.str());

	boost::multiprecision::float128 rounding(0.5);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			std::ostringstream tmpstream (std::ostringstream::app);

			if (!(p->usePCMethod==2&&simHistory->pcStep > 0)) {
				covering_phys = simHistory->coveringList.getLast(&f);
			}
			else{//p->usePCMethod==2&&simHistory->pcStep > 0
				covering_phys = simHistory->coveringList.getForelast(&f);
			}
			covering_sum = getCovering(&f, hitbuffer_sum);

			tmpstream <<std::endl << "Facet " << getFacetIndex(&f)<< std::endl;
			tmpstream << "covering_sum  = " <<covering_sum << " / "<< simHistory->smallCoveringFactor<< " = " << (static_cast < boost::multiprecision::float128 >(covering_sum))/(static_cast < boost::multiprecision::float128 >(simHistory->smallCoveringFactor)) << std::endl;
			tmpstream << "covering_phys_before = " << covering_phys << " = "<< boost::multiprecision::float128(covering_phys) << std::endl;

			covering_phys = static_cast < boost::multiprecision::uint128_t >(covering_phys*simHistory->smallCoveringFactor);//scale up covering_phys; covering_sum (in all buffers) is already scaled up before the iteration!
			if(!std::isinf(simHistory->errorList_covering.getCurrent(&f))){
				if(simHistory->errorList_covering.getCurrent(&f) < p->noupdateError){
							if (covering_sum > covering_phys){
								boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding +(boost::multiprecision::float128(covering_sum - covering_phys)/simHistory->smallCoveringFactor)*Krealvirt);
								covering_phys = static_cast < boost::multiprecision::uint128_t >(covering_phys/simHistory->smallCoveringFactor)+ covering_delta;//scale down covering_phys here; covering_delta already scaled down!
								tmpstream << "covering rises by " <<covering_delta << " = "<<boost::multiprecision::float128(covering_delta) << std::endl;
							}
							else{
								boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding + (boost::multiprecision::float128(covering_phys - covering_sum)/simHistory->smallCoveringFactor)*Krealvirt);
								covering_phys = static_cast < boost::multiprecision::uint128_t >(covering_phys/simHistory->smallCoveringFactor);//scale down covering_phys here; covering_delta already scaled down!
								// Check, if covering got negative

								if(covering_phys+1==covering_delta){
									tmpstream<<"!!!-----Correct covering_delta by 1: "<<covering_phys<<" - "<<covering_delta <<" because of rounding-----!!!"<<std::endl;
									covering_delta=covering_phys;
								}
								else if(covering_phys<covering_delta && covering_phys+10>=covering_delta){
									tmpstream <<"!!!-----Correct covering_delta: "<<covering_phys<<" - "<<covering_delta <<" because of bad statistics-----!!!"<<std::endl;
									covering_delta=covering_phys;
								}
								else if(covering_phys<covering_delta){
									tmpstream <<"!!!-----Covering gets negative: "<<covering_phys<<" - "<<covering_delta <<". Correct: Covering=0-----!!!"<<std::endl;
									covering_delta=covering_phys;
								}
								covering_phys -= covering_delta;
								tmpstream << "covering decreases by "<<covering_delta << " = " << boost::multiprecision::float128(covering_delta) << std::endl;
							}
				}
				else{tmpstream << "covering is not updated! error_covering = "<<simHistory->errorList_covering.getCurrent(&f) << std::endl;
				}
			}
			else{tmpstream << "covering is not updated! error_covering = inf" << std::endl;

			}
			tmpstream << "covering_phys_after = " << covering_phys << " = " << boost::multiprecision::float128(covering_phys) << std::endl;
			tmpstream << "coveringThreshold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << " = " << boost::multiprecision::float128(sHandle->coveringThreshold[getFacetIndex(&f)]) << std::endl;

			if (p->usePCMethod== 0||(p->usePCMethod==1&&simHistory->pcStep == 1)){
				simHistory->coveringList.setCurrent(&f, covering_phys);
			}
			else{
				simHistory->coveringList.setPredict(&f, covering_phys);
			} 

			printStream(tmpstream.str());
		}
	}
	// Save covering to simHistory
	simHistory->lastTime+= simHistory->stepSize;//update the value of lastTime
	if (p->usePCMethod != 2){//normal mode or normal Predictor-Corrector mode
		if (simHistory->pcStep == (p->usePCMethod?1:0)) {//update simHistory with current values and time value
			simHistory->coveringList.appendCurrent(simHistory->lastTime);
		}
	}
	else{//TimeStepControl Mode
		if(simHistory->pcStep > 0){// in this case we overwrite (erase, then append) the last point in time,
			//which we set at pcStep zero.
			simHistory->coveringList.erase(simHistory->coveringList.getlastindex());//erase last covering values as well as last time value
			simHistory->coveringList.appendPredict(simHistory->lastTime);
		}
		else{//simHistory->pcStep == 0
			simHistory->coveringList.appendPredict(simHistory->lastTime);
		}

	}
}

// Copy covering to hitbuffers
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer, bool step_size_change){
	boost::multiprecision::uint128_t covering_phys;
	if (!step_size_change){
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
					covering_phys = simHistory->coveringList.getLast(&f);
					getFacetHitBuffer(&f,hitbuffer)->covering=covering_phys;
					getFacetHitBuffer(&f,hitbuffer_sum)->covering=covering_phys;
				}
			}
	}
	else{
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
					covering_phys = simHistory->coveringList.getForelast(&f);
					getFacetHitBuffer(&f,hitbuffer)->covering=covering_phys;
					getFacetHitBuffer(&f,hitbuffer_sum)->covering=covering_phys;
				}
			}
	}
	simHistory->nParticles=0;
}

void UpdateErrorMain(){
	if(simHistory->pcStep > 0){// in this case we overwrite (erase, then append) the last point in time
		simHistory->errorList_event.erase(simHistory->errorList_event.getlastindex());
		simHistory->errorList_covering.erase(simHistory->errorList_covering.getlastindex());
		}
	simHistory->errorList_event.appendCurrent(simHistory->lastTime);
	simHistory->errorList_covering.appendCurrent(simHistory->lastTime);
}

void UpdateParticleDensityAndPressure(Databuff *hitbuffer_sum){
	//std::ostringstream tmpstream (std::ostringstream::app);
	//tmpstream<<std::endl;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			// Calculate particle density and pressure of facet
			simHistory->particleDensityList.setCurrent(&f, calcParticleDensity(hitbuffer_sum , &f));
			simHistory->pressureList.setCurrent(&f, calcPressure(hitbuffer_sum , &f));
		}
	}
	//tmpstream<<std::endl;
	//printStream(tmpstream.str());

	// Update history lists for particle density and pressure
	if(simHistory->pcStep > 0){// in this case we overwrite (erase, then append) the last point in time
		simHistory->particleDensityList.erase(simHistory->particleDensityList.getlastindex());
		simHistory->pressureList.erase(simHistory->pressureList.getlastindex());
	}
	simHistory->particleDensityList.appendCurrent(simHistory->lastTime);
	simHistory->pressureList.appendCurrent(simHistory->lastTime);
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration(){
	// Calculates statistical uncertainties of error_event and error_covering at the end of the simulation
	// While simulating only the error values of the subprocesses are used to decide, if the targeted error level has been reached.
	std::vector<double> errorPerIt_event;
	errorPerIt_event =std::vector<double> ();
	std::vector<double> errorPerIt_covering;
	errorPerIt_covering =std::vector<double> ();

	std::vector<boost::multiprecision::uint128_t> covPerIt;
	covPerIt =std::vector<boost::multiprecision::uint128_t> ();

	for(unsigned int it=0; it<simHistory->errorList_event.historyList.first.size();it++){
		// Total error/covering for each iteration
		double total_error_event=0.0;
		double total_error_covering=0.0;
		std::tie(total_error_event,total_error_covering)=CalcErrorAll(it);

		boost::multiprecision::uint128_t covering=boost::multiprecision::uint128_t(0);
		for (unsigned int idx=0; idx<simHistory->numFacet; idx++) {
			covering+=simHistory->coveringList.historyList.second[idx][it];
		}

		errorPerIt_event.push_back(total_error_event);
		errorPerIt_covering.push_back(total_error_covering);
		covPerIt.push_back(covering);
	}

	return std::make_tuple(errorPerIt_event,errorPerIt_covering,covPerIt);
}

//-----------------------------------------------------------
// Update Buffer of main process using buffer of sub process
//Simhistory version

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank) {
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

	// Global hits and leaks: adding local hits to shared memory
	sHandle->tmpGlobalResult.globalHits.nbMCHit= gHits->globalHits.nbMCHit += subHits->globalHits.nbMCHit;
	sHandle->tmpGlobalResult.globalHits.nbHitEquiv=gHits->globalHits.nbHitEquiv += subHits->globalHits.nbHitEquiv;
	sHandle->tmpGlobalResult.globalHits.nbAbsEquiv=gHits->globalHits.nbAbsEquiv += subHits->globalHits.nbAbsEquiv;
	sHandle->tmpGlobalResult.globalHits.nbDesorbed=gHits->globalHits.nbDesorbed += subHits->globalHits.nbDesorbed;
	sHandle->tmpGlobalResult.globalHits.nbOutgassed=gHits->globalHits.nbOutgassed += subHits->globalHits.nbOutgassed;
	sHandle->tmpGlobalResult.distTraveled_total=gHits->distTraveled_total += subHits->distTraveled_total;
	sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly=gHits->distTraveledTotal_fullHitsOnly += subHits->distTraveledTotal_fullHitsOnly;

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
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {

		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//if (f.hitted) {
				//std::cout <<"Facet " <<num <<std::endl; // Da wird immer "0" angezeigt. Was ist der Sinn?
				for (unsigned int m = 0; m < (1 + nbMoments); m++) { // Add hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitSub = (FacetHitBuffer *)(subbuff + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					boost::multiprecision::uint128_t covering_phys= simHistory->smallCoveringFactor * history->coveringList.getLast(&f);
					boost::multiprecision::uint128_t covering_sum = facetHitSub->covering;
					f.tmpCounter[m].nbAbsEquiv = facetHitBuffer->nbAbsEquiv += facetHitSub->nbAbsEquiv;
					f.tmpCounter[m].nbDesorbed = facetHitBuffer->nbDesorbed += facetHitSub->nbDesorbed;
					f.tmpCounter[m].nbOutgassed = facetHitBuffer->nbOutgassed += facetHitSub->nbOutgassed;
					f.tmpCounter[m].nbMCHit = facetHitBuffer->nbMCHit += facetHitSub->nbMCHit;
					f.tmpCounter[m].nbHitEquiv = facetHitBuffer->nbHitEquiv += facetHitSub->nbHitEquiv;;
					f.tmpCounter[m].sum_1_per_ort_velocity = facetHitBuffer->sum_1_per_ort_velocity += facetHitSub->sum_1_per_ort_velocity;
					f.tmpCounter[m].sum_v_ort = facetHitBuffer->sum_v_ort += facetHitSub->sum_v_ort;
					f.tmpCounter[m].sum_1_per_velocity = facetHitBuffer->sum_1_per_velocity += facetHitSub->sum_1_per_velocity;
					//facetHitBuffer->covering += facetHitSub->covering; //We do that in another way.

					if (facetHitSub->covering > covering_phys){
						facetHitBuffer->covering += (covering_sum - covering_phys);
					}
					else{
						if(facetHitBuffer->covering > (covering_phys - covering_sum)){
							facetHitBuffer->covering -= (covering_phys- covering_sum);
							}
						else{
							facetHitBuffer->covering = 0;//Counter cannot be negative! Maybe we could interrupt the iteration here?
						}
					}
					f.tmpCounter[m].covering = facetHitBuffer->covering;
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
						//double dCoef = gHits->globalHits.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

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

	//ResetTmpCounters(); //This will reset counters -> but needed for t_min

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}


//----------deprecated functions because of new K_real/virt approach
/*
double manageStepSize(){
	double step_size;
	bool steptoolong; //If step is too long, we will decrease the step size.
	bool needforCheck = true;

	while(needforCheck){
		step_size = simHistory->stepSize;
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
				std::ostringstream tmpstream (std::ostringstream::app);
				tmpstream<<"Decrease simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
				printStream(tmpstream.str());
			}
		}
		else {
			needforCheck = false;
		}
	}
	return step_size;
}
*/

//----------Not used anymore
/*
void printVelocities(Databuff *hitbuffer){
	std::ostringstream tmpstream (std::ostringstream::app);
	tmpstream<<std::endl<<"Velocity counter after iteration"<<std::endl;
	tmpstream<<std::setw(7+12)<<std::right <<"1_ort_v"<<std::setw(12)<<std::right <<"ort_v"<<std::setw(12)<<std::right <<"1_v"<<std::endl;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			double sum_1_per_ort_velocity, sum_v_ort, sum_1_per_velocity=0.0;
			std::tie(sum_1_per_ort_velocity, sum_v_ort, sum_1_per_velocity)=getVelocities(&f, hitbuffer);
			tmpstream<<"Facet "<<getFacetIndex(&f) <<std::setw(12)<<std::right<<sum_1_per_ort_velocity <<std::setw(12)<<std::right<<sum_v_ort <<std::setw(12)<<std::right<<sum_1_per_velocity <<std::endl;
		}
	}
	tmpstream<<std::endl;
	printStream(tmpstream.str());
}
*/
