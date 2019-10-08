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
void UpdateSticking(Databuff *hitbuffer){
	//std::cout <<"    Facet information:" <<std::endl;
	int i = 0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			calcStickingnew(&f, hitbuffer);
			/*if(i==1){
			std::cout <<"    (Just for Facet " << i << " to have less output in console.)" << std::endl;
			std::cout <<"    Facet " <<i <<":" <<std::endl;
			std::cout <<"\t Facet area [cm²]\t\t\t"<< f.sh.area << std::endl;
			std::cout <<"\t Facet temperature [K] \t\t\t" << f.sh.temperature << std::endl;
			std::cout <<"\t Sticking [1]\t\t\t\t" <<f.sh.sticking <<std::endl;
			std::cout <<"\t Outgassing [Pa m³/s]\t\t\t" <<f.sh.outgassing <<std::endl;
			std::cout <<"\t Desorption rate [1/s]\t\t\t" <<calcDesorption(&f, hitbuffer) <<std::endl;
			std::cout <<"\t Desorption rate [Pa m³/s]\t\t" <<calcDesorptionRate(&f, hitbuffer) << std::endl;
			std::cout <<"\t Outgassing + Desorption rate [Pa m³/s]\t" <<f.sh.outgassing + calcDesorptionRate(&f, hitbuffer)<<std::endl;
			std::cout <<"\t Covering [Number of particles]\t\t" <<calcCovering(&f) <<std::endl;
			std::cout <<"\t Coverage [1 = Monolayer]\t\t" << calcCovering(&f) / (calcNmono(&f)/calcdNsurf())<< std::endl;
			//std::cout <<"\t real covering\t\t\t\t" <<calcRealCovering(&f) <<std::endl; // Das macht hier gar keinen Sinn. Erst wenn alle Subprozesse ihre Counter
			// zum Prozess 0 geschickt haben und alle Counter addiert wurden, wird Krealvirt berechnet. Dann kann man ein real Covering ausgeben.
			}*/
			i+=1;
		}
	}
	//std::cout <<std::endl;
}
//desorption
void UpdateDesorptionRate (Databuff *hitbuffer){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.desorption = (double)calcDesorptionRate(&f, hitbuffer); // TODO long double -> double here. Change f.sh.desorption to long double? bit need to adapt in Windows and new buffers
		}
	}

}

void UpdateSojourn(Databuff *hitbuffer){
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			f.sh.enableSojournTime = true;
			f.sh.sojournFreq = (kb*f.sh.temperature)/h;
			f.sh.sojournE = calcEnergy(&f, hitbuffer);
		}
	}
}

void UpdateErrorSub(){

	double num_hit_it=0;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) { //save current num total hits in currentList, add difference current-old to num_hit_it
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			num_hit_it+=f.sh.opacity * (f.tmpCounter[0].hit.nbHitEquiv + f.tmpCounter[0].hit.nbDesorbed);
		}
	}

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			double num_hit_f=f.sh.opacity * ( f.tmpCounter[0].hit.nbHitEquiv + f.tmpCounter[0].hit.nbDesorbed);

			//neglect hits if very small compared to total hits
			if(num_hit_f/num_hit_it<p->hitRatioLimit){
				num_hit_it-=num_hit_f; //TODO also adapt facet counters??
				num_hit_f=0;
			}

			if(f.sh.opacity==0){simHistory->errorList.setCurrentList(&f, 0.0);}
			else{
				double error=pow((1/num_hit_f)*(1-num_hit_f/num_hit_it),0.5);

				simHistory->errorList.setCurrentList(&f, error);
			}
			simHistory->hitList.setCurrentList(&f,f.tmpCounter[0].hit.nbHitEquiv);
			simHistory->desorbedList.setCurrentList(&f,f.tmpCounter[0].hit.nbDesorbed);
		}
	}

}

double UpdateError(){
	UpdateErrorSub();

	double error=0.0;
	double area=0.0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if(simHistory->errorList.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
				continue;

			error+=simHistory->errorList.getCurrent(&f)*f.sh.area;
			area+=f.sh.area;
		}
	}

	return error/area;
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

			} // End if(hitted)
			/*
			else // if not hitted, initialize to zero
			{
				for (unsigned int m = 0; m < (1 + nbMoments); m++) {//hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					facetHitBuffer->hit.nbAbsEquiv = 0.0;
					facetHitBuffer->hit.nbDesorbed = 0;
					facetHitBuffer->hit.nbMCHit = 0;
					facetHitBuffer->hit.nbHitEquiv = 0.0;
					facetHitBuffer->hit.sum_1_per_ort_velocity = 0.0;
					facetHitBuffer->hit.sum_v_ort = 0.0;
					facetHitBuffer->hit.sum_1_per_velocity = 0.0;
					//facetHitBuffer->hit.covering= 0.0; //If a facet is not hit, its covering stays the same instead ob being zero!
				}

				if (f.sh.isProfile) {//profile
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							shProfile[j].countEquiv=0.0; shProfile[j].sum_1_per_ort_velocity=0.0; shProfile[j].sum_v_ort=0.0;
						}
					}
				}

				if (f.sh.isTextured) {//texture
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//temporary hit counts
								shTexture[add].countEquiv=0.0; shTexture[add].sum_1_per_ort_velocity=0.0; shTexture[add].sum_v_ort_per_area=0.0;

							}
						}
					}
				}

				if (f.sh.countDirection) {//
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

				if (f.sh.anglemapParams.record) {//
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							shAngleMap[add] = 0;
						}
					}
				}

				//Facet histograms

					for (unsigned int m = 0; m < (1 + nbMoments); m++) {//facet histograms
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
			}//end else
			*/
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
//Wahrscheinlich müssten hier auch noch alle Profiles und Textures resetet werden!
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
