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
 * This file contains functions that saves simulation results in Databuff structs
 */


#include "SimulationLinux.h"
extern Simulation *sHandle;
extern ProblemDef* p;
extern SimulationHistory* simHistory;

//-----------------------------------------------------------
//----Update values of subprocess
//sticking
void UpdateSticking(){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			calcStickingnew(&f);
		}
	}
}

//desorption
bool UpdateDesorptionRate (){
	boost::multiprecision::float128 totaldes(0.0);
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.desorption = calcDesorption(&f); // TODO long double -> double here. Change f.sh.desorption to long double? bit need to adapt in Windows and new buffers
			//f.sh.desorption is now a number of desorbed particles per time step!
				if(f.sh.temperature==0) {continue;}
				totaldes+= f.sh.desorption * boost::multiprecision::float128(sHandle->wp.latestMoment/ (1.38E-23*f.sh.temperature));
		}
	}

	if((boost::multiprecision::float128(sHandle->wp.totalDesorbedMolecules)+totaldes)<boost::multiprecision::pow(boost::multiprecision::float128(10.),boost::multiprecision::float128(-50))){return false;}
	else{return true;}
}

void UpdateSojourn(){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.enableSojournTime = true;
			//f.sh.sojournFreq = (kb*f.sh.temperature)/h; //Don't need that at the moment. If temperature won't be constant later, reactivation is necessary.
			//f.sh.sojournE = calcEnergy(&f); //Don't need that. We decide in each case of residence at random, if the particle will adsorb at substrate or adsorbate.
		}
	}
}

void UpdateErrorSub(){

	double num_hit_it=0;
	double num_des_ad_it=0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) { //save current num total hits in currentList, add difference current-old to num_hit_it
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			num_hit_it+=f.sh.opacity * (f.tmpCounter[0].hit.nbHitEquiv + (double)f.tmpCounter[0].hit.nbDesorbed);
			num_des_ad_it+=f.sh.opacity * (f.tmpCounter[0].hit.nbAbsEquiv + (double)f.tmpCounter[0].hit.nbDesorbed);
		}
	}

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			double num_hit_f=f.sh.opacity * ( f.tmpCounter[0].hit.nbHitEquiv + (double)f.tmpCounter[0].hit.nbDesorbed);
			double num_des_ad_f=f.sh.opacity * ( f.tmpCounter[0].hit.nbAbsEquiv + (double)f.tmpCounter[0].hit.nbDesorbed);

			//neglect hits if very small compared to total hits
			if(num_hit_f/num_hit_it<(p->hitRatioLimit)/pow(simHistory->numSubProcess,0.5)){// To be consistent with the ignored facets for calculating the error after summation
				//over all subprocesses, here the hitRationLimit must be reduced with the correction factor due to multiple subprocesses.
				num_hit_it-=num_hit_f;
				num_hit_f=0;
			}
			if(num_des_ad_f/num_des_ad_it<(p->hitRatioLimit)/pow(simHistory->numSubProcess,0.5)){// To be consistent with the ignored facets for calculating the error after summation
				//over all subprocesses, here the hitRationLimit must be reduced with the correction factor due to multiple subprocesses.
				num_des_ad_it-=num_des_ad_f;
				num_des_ad_f=0;
			}



			if(f.sh.opacity==0){
				simHistory->errorList_event.setCurrent(&f, 0.0);
				simHistory->errorList_covering.setCurrent(&f, 0.0);
			}
			else{
				double error_event=pow((1/num_hit_f)*(1-num_hit_f/num_hit_it),0.5);
				double error_covering=pow((1/num_des_ad_f)*(1-num_des_ad_f/num_des_ad_it),0.5);
				simHistory->errorList_event.setCurrent(&f, error_event);
				simHistory->errorList_covering.setCurrent(&f, error_covering);
			}
			simHistory->hitList.setCurrent(&f,f.tmpCounter[0].hit.nbHitEquiv);
			simHistory->desorbedList.setCurrent(&f,f.tmpCounter[0].hit.nbDesorbed);
		}
	}

}

double UpdateError(std::string mode){//calculates the averaged total error weighted with the facets area to decide, if the desired uncertainty level is reached
	UpdateErrorSub();

	double error=0.0;
	double area=0.0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if(mode=="covering"){
				if(simHistory->errorList_covering.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
					continue;
				error+=simHistory->errorList_covering.getCurrent(&f)*f.sh.area;
			}
			else if(mode =="event"){
				if(simHistory->errorList_event.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
					continue;
				error+=simHistory->errorList_event.getCurrent(&f)*f.sh.area;
			}
			else{
				std::cout<<"------------! Error mode '"<<mode <<"' not implemented !------------\n";
				return 0.0;
			}
			area+=f.sh.area;
		}
	}

	return error/area;
}

std::tuple<double,double> UpdateErrorAll(){//calculates the averaged total error weighted with the facets area to decide, if the desired uncertainty level is reached
	double error_event=0.0; double error_covering=0.0;
	double area=0.0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if(simHistory->errorList_event.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
				continue;
			if(simHistory->errorList_covering.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
				continue;
			error_event+=simHistory->errorList_event.getCurrent(&f)*f.sh.area;
			error_covering+=simHistory->errorList_covering.getCurrent(&f)*f.sh.area;

			area+=f.sh.area;
		}
	}

	return std::make_tuple(error_event/area, error_covering/area);
}

bool checkErrorSub(double targetError, double currentError, double factor, std::string mode){
	bool vipCheck = currentError<=targetError;
	if(!p->vipFacets.empty()){
		for(unsigned int i = 0; i < p->vipFacets.size(); i++){
			if(mode=="covering"){
				vipCheck = vipCheck && (simHistory->errorList_covering.getCurrent(p->vipFacets[i].first)== std::numeric_limits<double>::infinity() || simHistory->errorList_covering.getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
			}
			else if(mode=="event"){
				vipCheck = vipCheck && (simHistory->errorList_event.getCurrent(p->vipFacets[i].first)== std::numeric_limits<double>::infinity() || simHistory->errorList_event.getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
			}
			else{
				std::cout<<"------------! Error mode '"<<mode <<"' not implemented !------------\n";
				return true;
			}
		}
	}
	return vipCheck;
}

//-----------------------------------------------------------
// copy sHandle values to buffer of sub process
void UpdateMCSubHits(Databuff *databuffer, int rank) {
	BYTE *buffer;
	GlobalHitBuffer *gHits;
	//TEXTURE_MIN_MAX texture_limits_old[3];
	int j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif
	/* (Rudi) Don't need that.
	SetState(NULL, "Waiting for 'hits' dataport access...", false, true);
	sHandle->lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
	SetState(NULL, "Updating MC hits...", false, true);
	if (!sHandle->lastHitUpdateOK) return; //Timeout, will try again later
	*/

	buffer = databuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	size_t nbMoments=(size_t)sHandle->moments.size();

	// Global hits and leaks: save local hits to shared memory
	gHits->globalHits.hit.nbMCHit = sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	gHits->globalHits.hit.nbHitEquiv = sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;
	gHits->globalHits.hit.nbAbsEquiv = sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv;
	gHits->globalHits.hit.nbDesorbed = sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
	//gHits->globalHits.hit.covering=0;
	gHits->distTraveled_total = sHandle->tmpGlobalResult.distTraveled_total;
	gHits->distTraveledTotal_fullHitsOnly = sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly;

	//Memorize current limits, then do a min/max search //(My) not needed for subprocesses


	// Leak saved
	for (size_t leakIndex = 0; leakIndex < sHandle->tmpGlobalResult.leakCacheSize; leakIndex++)//TODO which one correct?
		gHits->leakCache[(leakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
		//gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
	gHits->nbLeakTotal = sHandle->tmpGlobalResult.nbLeakTotal;
	gHits->lastLeakIndex = sHandle->tmpGlobalResult.leakCacheSize;
	gHits->leakCacheSize = sHandle->tmpGlobalResult.leakCacheSize;


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	//if (rank == 1) {// (MY) save hitcache
		for (size_t hitIndex = 0; hitIndex < sHandle->tmpGlobalResult.hitCacheSize; hitIndex++)//TODO which one correct?
			gHits->hitCache[(hitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];
			//gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];

		if (sHandle->tmpGlobalResult.hitCacheSize > 0) {
			gHits->lastHitIndex = sHandle->tmpGlobalResult.hitCacheSize;
			gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			gHits->hitCacheSize = sHandle->tmpGlobalResult.hitCacheSize;
		}
	//}


	//Global histograms saved

		for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
			BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();

			double* nbHitsHistogram = (double*)histCurrentMoment;
			double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());

			if (sHandle->wp.globalHistogramParams.recordBounce) {
				for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
					nbHitsHistogram[i] = sHandle->tmpGlobalHistograms[m].nbHitsHistogram[i];
				}
			}

			if (sHandle->wp.globalHistogramParams.recordDistance) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] = sHandle->tmpGlobalHistograms[m].distanceHistogram[i];
				}
			}

			if (sHandle->wp.globalHistogramParams.recordTime) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] = sHandle->tmpGlobalHistograms[m].timeHistogram[i];
				}
			}
		}

	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if (f.hitted) {

				for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) save hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));

					facetHitBuffer->hit.nbAbsEquiv = f.tmpCounter[m].hit.nbAbsEquiv;
					facetHitBuffer->hit.nbDesorbed = f.tmpCounter[m].hit.nbDesorbed;
					facetHitBuffer->hit.nbMCHit = f.tmpCounter[m].hit.nbMCHit;
					facetHitBuffer->hit.nbHitEquiv = f.tmpCounter[m].hit.nbHitEquiv;
					facetHitBuffer->hit.sum_1_per_ort_velocity = f.tmpCounter[m].hit.sum_1_per_ort_velocity;
					facetHitBuffer->hit.sum_v_ort = f.tmpCounter[m].hit.sum_v_ort;
					facetHitBuffer->hit.sum_1_per_velocity = f.tmpCounter[m].hit.sum_1_per_velocity;
					facetHitBuffer->hit.covering= f.tmpCounter[m].hit.covering;
				}

				if (f.sh.isProfile) {//(MY) save profile
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							shProfile[j] = f.profile[m][j];
						}
					}
				}

				if (f.sh.isTextured) {//(MY) save texture
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						//double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Temporary hit counts
								shTexture[add] = f.texture[m][add]; //(My) removed +
								// Autoscaling will be done in SimulationMCmain.cpp

							}
						}
					}
				}

				if (f.sh.countDirection) {//(MY) save this
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;
								shDir[add].dir.x = f.direction[m][add].dir.x;
								shDir[add].dir.y = f.direction[m][add].dir.y;
								shDir[add].dir.z = f.direction[m][add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								shDir[add].count = f.direction[m][add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {//(MY) save this
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							shAngleMap[add] = f.angleMap.pdf[add];
						}
					}
				}

				//Facet histograms
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) Save histogrms
						BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();

						if (f.sh.facetHistogramParams.recordBounce) {
							double* nbHitsHistogram = (double*)histCurrentMoment;
							for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
								nbHitsHistogram[i] = f.tmpHistograms[m].nbHitsHistogram[i];
							}
						}


						if (f.sh.facetHistogramParams.recordDistance) {
							double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
								distanceHistogram[i] = f.tmpHistograms[m].distanceHistogram[i];
							}
						}

						if (f.sh.facetHistogramParams.recordTime) {
							double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
								timeHistogram[i] = f.tmpHistograms[m].timeHistogram[i];
							}
						}

					}

			}
		} // End nbFacet
	} // End nbSuper

	//if there were no textures: //(My) not needed for subprocesses


	//ResetTmpCounters();

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

	return;
}

//-----------------------------------------------------------
//reset of counters/buffers
/*
void initcounterstozero(Databuff *databuffer){
	BYTE *buffer;
	GlobalHitBuffer *gHits;

	buffer=NULL;
	buffer = databuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	size_t nbMoments=(size_t)sHandle->moments.size();
	//size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);

	gHits->globalHits.hit.nbMCHit = 0;
	gHits->globalHits.hit.nbHitEquiv = 0.0;
	gHits->globalHits.hit.nbAbsEquiv = 0.0;
	gHits->globalHits.hit.nbDesorbed = 0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			for (unsigned int m = 0; m < (1 + nbMoments); m++) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
				facetHitBuffer->hit.nbAbsEquiv = 0.0;
				facetHitBuffer->hit.nbDesorbed = 0;
				facetHitBuffer->hit.nbMCHit = 0;
				facetHitBuffer->hit.nbHitEquiv = 0.0;
				facetHitBuffer->hit.sum_1_per_ort_velocity = 0.0;
				facetHitBuffer->hit.sum_v_ort = 0.0;
				facetHitBuffer->hit.sum_1_per_velocity = 0.0;

			}
		}
	}

}*/

void initbufftozero(Databuff *databuffer){

	BYTE *buffer;
	GlobalHitBuffer *gHits;

	buffer=NULL;
	buffer = databuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	// Global hits and leaks: save local hits to shared memory
	gHits->globalHits.hit.nbMCHit = 0;
	gHits->globalHits.hit.nbHitEquiv = 0.0;
	gHits->globalHits.hit.nbAbsEquiv = 0.0;
	gHits->globalHits.hit.nbDesorbed = 0;
	//gHits->globalHits.hit.covering=0;
	gHits->distTraveled_total = 0.0;
	gHits->distTraveledTotal_fullHitsOnly = 0.0;

	size_t nbMoments=(size_t)sHandle->moments.size();
	//Memorize current limits, then do a min/max search //(My) not needed for subprocesses


	// Leak saved
	/*
	for (size_t leakIndex = 0; leakIndex < sHandle->tmpGlobalResult.leakCacheSize; leakIndex++)//TODO which one correct?
		gHits->leakCache[(leakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
		//gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
	gHits->nbLeakTotal = sHandle->tmpGlobalResult.nbLeakTotal;
	gHits->lastLeakIndex = sHandle->tmpGlobalResult.leakCacheSize;
	gHits->leakCacheSize = sHandle->tmpGlobalResult.leakCacheSize;*/


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	//if (rank == 1) {// (MY) save hitcache
	/*
		for (size_t hitIndex = 0; hitIndex < sHandle->tmpGlobalResult.hitCacheSize; hitIndex++)//TODO which one correct?
			gHits->hitCache[(hitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];
			//gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];

		if (sHandle->tmpGlobalResult.hitCacheSize > 0) {
			gHits->lastHitIndex = sHandle->tmpGlobalResult.hitCacheSize;
			gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			gHits->hitCacheSize = sHandle->tmpGlobalResult.hitCacheSize;
		}*/
	//}

	//Global histograms saved

		for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
			BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();

			double* nbHitsHistogram = (double*)histCurrentMoment;
			double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());

			if (sHandle->wp.globalHistogramParams.recordBounce) {
				for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
					nbHitsHistogram[i] = 0.0;
				}
			}

			if (sHandle->wp.globalHistogramParams.recordDistance) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] = 0.0;
				}
			}

			if (sHandle->wp.globalHistogramParams.recordTime) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] = 0.0;
				}
			}
		}


	//here: init values with zero
	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	int j, x, y, s;

	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {

				//if (f.hitted) {

					for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
						FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
						facetHitBuffer->hit.nbAbsEquiv = 0.0;
						facetHitBuffer->hit.nbDesorbed = 0;
						facetHitBuffer->hit.nbMCHit = 0;
						facetHitBuffer->hit.nbHitEquiv = 0.0;
						facetHitBuffer->hit.sum_1_per_ort_velocity = 0.0;
						facetHitBuffer->hit.sum_v_ort = 0.0;
						facetHitBuffer->hit.sum_1_per_velocity = 0.0;
					}

					if (f.sh.isProfile) {//(MY) removed +
						for (unsigned int m = 0; m < (1 + nbMoments); m++) {
							ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
							for (j = 0; j < (int)PROFILE_SIZE; j++) {
								shProfile[j].countEquiv=0.0; shProfile[j].sum_1_per_ort_velocity=0.0; shProfile[j].sum_v_ort=0.0;
							}
						}
					}

					if (f.sh.isTextured) {//(MY)
						for (unsigned int m = 0; m < (1 + nbMoments); m++) {
							TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));

							for (y = 0; y < (int)f.sh.texHeight; y++) {
								for (x = 0; x < (int)f.sh.texWidth; x++) {
									size_t add = x + y * f.sh.texWidth;

									//Add temporary hit counts
									shTexture[add].countEquiv=0.0; shTexture[add].sum_1_per_ort_velocity=0.0; shTexture[add].sum_v_ort_per_area=0.0;

								}
							}
						}
					}

					if (f.sh.countDirection) {//(MY) removed +
						for (unsigned int m = 0; m < (1 + nbMoments); m++) {
							DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
							for (y = 0; y < (int)f.sh.texHeight; y++) {
								for (x = 0; x < (int)f.sh.texWidth; x++) {
									size_t add = x + y * f.sh.texWidth;
									shDir[add].dir.x = 0.0;
									shDir[add].dir.y = 0.0;
									shDir[add].dir.z = 0.0;
									//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
									shDir[add].count = 0;
								}
							}
						}
					}

					if (f.sh.anglemapParams.record) {//(MY) removed +
						size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
						for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
							for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
								size_t add = x + y * f.sh.anglemapParams.phiWidth;
								shAngleMap[add] = 0;
							}
						}
					}

					//Facet histograms

						for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
							BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();

							if (f.sh.facetHistogramParams.recordBounce) {
								double* nbHitsHistogram = (double*)histCurrentMoment;
								for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
									nbHitsHistogram[i] = 0.0;
								}
							}


							if (f.sh.facetHistogramParams.recordDistance) {
								double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
								for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
									distanceHistogram[i] = 0.0;
								}
							}

							if (f.sh.facetHistogramParams.recordTime) {
								double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
								for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
									timeHistogram[i] = 0.0;
								}
							}

						}

				//} // End if(hitted)

			} // End nbFacet
		} // End nbSuper

	return;
}

//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
void UpdateSticking(Databuff *hitbuffer){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			calcStickingnew(&f, hitbuffer);
		}
	}
}

bool UpdateDesorptionRate (Databuff *hitbuffer){
	boost::multiprecision::float128 totaldes(0.0);
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.desorption = calcDesorptionRate(&f, hitbuffer); // TODO long double -> double here. Change f.sh.desorption to long double? bit need to adapt in Windows and new buffers
				if(f.sh.temperature==0) {continue;}
				totaldes+= f.sh.desorption * boost::multiprecision::float128(sHandle->wp.latestMoment/ (1.38E-23*f.sh.temperature));
		}
	}

	if((boost::multiprecision::float128(sHandle->wp.totalDesorbedMolecules)+totaldes)<boost::multiprecision::pow(boost::multiprecision::float128(10.),boost::multiprecision::float128(-50))){return false;}
	else{return true;}
}

void UpdateSojourn(Databuff *hitbuffer){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.enableSojournTime = true;
			f.sh.sojournFreq = (kb*f.sh.temperature)/h;
			f.sh.sojournE = calcEnergy(&f, hitbuffer);
		}
	}
}*/
