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

extern Simulation *sHandle; //delcared in molflowSub.cpp

void UpdateMainHits(Databuff *databuffer,Databuff *subbuffer, int rank) {//Schreibt (=Summiert) die Daten aller Hitbuffer aus den Subprozessen in den Hitbuffer_original des Hauptprozesses.
	switch (sHandle->wp.sMode) {
	case MC_MODE:
	{
		 // std::cout <<"shandle size " <<(double)sHandle->moments.size() << std::endl;
		UpdateMCmainHits(databuffer, subbuffer, rank, (size_t)sHandle->moments.size());
		//if (dpLog) UpdateLog(dpLog, timeout);
	}
		break;
	case AC_MODE:

		//UpdateACHits(dpHit, prIdx, timeout);
		break;
	}

}

void UpdateMCmainHits(Databuff *mainbuffer, Databuff *subbuffer,int rank, size_t nbMoments) {
	BYTE *buffer, *subbuff;
	GlobalHitBuffer *gHits, *subHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif


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
	gHits->globalHits.hit.nbMCHit += subHits->globalHits.hit.nbMCHit;
	gHits->globalHits.hit.nbHitEquiv += subHits->globalHits.hit.nbHitEquiv;
	gHits->globalHits.hit.nbAbsEquiv += subHits->globalHits.hit.nbAbsEquiv;
	gHits->globalHits.hit.nbDesorbed += subHits->globalHits.hit.nbDesorbed;
	gHits->distTraveled_total += subHits->distTraveled_total;
	gHits->distTraveledTotal_fullHitsOnly += subHits->distTraveledTotal_fullHitsOnly;
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
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = subHits->leakCache[leakIndex];
	gHits->nbLeakTotal += subHits->nbLeakTotal;
	gHits->lastLeakIndex = (gHits->lastLeakIndex + subHits->leakCacheSize) % LEAKCACHESIZE;
	gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + subHits->leakCacheSize);


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	if (rank == 0) {
		for (size_t hitIndex = 0; hitIndex < subHits->hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = subHits->hitCache[hitIndex];

		if (subHits->hitCacheSize > 0) {
			gHits->lastHitIndex = (gHits->lastHitIndex + subHits->hitCacheSize) % HITCACHESIZE;
			gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + subHits->hitCacheSize);
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
					nbHitsHistogram[i] += nbHitsSub[i];
				}
			}

			if (sHandle->wp.globalHistogramParams.recordDistance) {
				double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
				double* distanceSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] += distanceSub[i];
				}
			}
			if (sHandle->wp.globalHistogramParams.recordTime) {
				double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
				double* timeSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] += timeSub[i];
				}
			}

		}


	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	int num=0;
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {

		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//if (f.hitted) {
				//std::cout <<"Facet " <<num <<std::endl; // Da wird immer "0" angezeigt. Was ist der Sinn?
				for (unsigned int m = 0; m < (1 + nbMoments); m++) { // Add hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitSub = (FacetHitBuffer *)(subbuff + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
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

					facetHitBuffer->hit.nbAbsEquiv += facetHitSub->hit.nbAbsEquiv;
					facetHitBuffer->hit.nbDesorbed += facetHitSub->hit.nbDesorbed;
					facetHitBuffer->hit.nbMCHit += facetHitSub->hit.nbMCHit;
					facetHitBuffer->hit.nbHitEquiv += facetHitSub->hit.nbHitEquiv;;
					facetHitBuffer->hit.sum_1_per_ort_velocity += facetHitSub->hit.sum_1_per_ort_velocity;
					facetHitBuffer->hit.sum_v_ort += facetHitSub->hit.sum_v_ort;
					facetHitBuffer->hit.sum_1_per_velocity += facetHitSub->hit.sum_1_per_velocity;
					facetHitBuffer->hit.covering += facetHitSub->hit.covering;

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
							shProfile[j] += shProfileSub[j];
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
								shTexture[add] += shTextureSub[add];

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
								shDir[add].dir.x += shDirSub[add].dir.x;
								shDir[add].dir.y += shDirSub[add].dir.y;
								shDir[add].dir.z += shDirSub[add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								shDir[add].count += shDirSub[add].count;
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
							shAngleMap[add] += shAngleMapSub[add];
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
								nbHitsHistogram[i] += nbHitsSub[i];
							}
						}
						if (f.sh.facetHistogramParams.recordDistance) {
							double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
							double* distanceSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
								distanceHistogram[i] += distanceSub[i];
							}
						}
						if (f.sh.facetHistogramParams.recordTime) {
							double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							double* timeSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
								timeHistogram[i] += timeSub[i];
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
	}


	ResetTmpCounters();
	//extern char* GetSimuStatus();
	//SetState(NULL, GetSimuStatus(), false, true); // (Rudi) Don't need that.

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}
