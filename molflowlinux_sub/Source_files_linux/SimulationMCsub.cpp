#include "SimulationLinux.h"


extern Simulation *sHandle; //delcared in molflowSub.cpp



void UpdateSubHits(Databuff *databuffer, int rank) {
	switch (sHandle->wp.sMode) {
	case MC_MODE:
	{
		UpdateSubMCHits(databuffer, rank, sHandle->moments.size());
		//if (dpLog) UpdateLog(dpLog, timeout);
	}
		break;
	case AC_MODE:

		//UpdateACHits(dpHit, prIdx, timeout);
		break;
	}

}

void UpdateSubMCHits(Databuff *databuffer, int rank, size_t nbMoments) {
	BYTE *buffer;
	GlobalHitBuffer *gHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
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

	// Global hits and leaks: adding local hits to shared memory (My) removed +, etc
	gHits->globalHits.hit.nbMCHit = sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	gHits->globalHits.hit.nbHitEquiv = sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;
	gHits->globalHits.hit.nbAbsEquiv = sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv;
	gHits->globalHits.hit.nbDesorbed = sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
	gHits->distTraveled_total = sHandle->tmpGlobalResult.distTraveled_total;
	gHits->distTraveledTotal_fullHitsOnly = sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly;

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	//sHandle->wp.sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak (MY) removed +, etc
	for (size_t leakIndex = 0; leakIndex < sHandle->tmpGlobalResult.leakCacheSize; leakIndex++)
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
	gHits->nbLeakTotal = sHandle->tmpGlobalResult.nbLeakTotal;
	gHits->lastLeakIndex = sHandle->tmpGlobalResult.leakCacheSize;
	gHits->leakCacheSize = sHandle->tmpGlobalResult.leakCacheSize;

	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	if (rank == 1) {// (MY) removed +, etc//(MY) TODO Init for else? or not?
		for (size_t hitIndex = 0; hitIndex < sHandle->tmpGlobalResult.hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];

		if (sHandle->tmpGlobalResult.hitCacheSize > 0) {
			gHits->lastHitIndex = sHandle->tmpGlobalResult.hitCacheSize;
			gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			gHits->hitCacheSize = sHandle->tmpGlobalResult.hitCacheSize;
		}
	}

	//Global histograms (MY) change that?

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
			else{
				for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
					nbHitsHistogram[i] = 0.0;
				}
			}

			if (sHandle->wp.globalHistogramParams.recordDistance) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] = sHandle->tmpGlobalHistograms[m].distanceHistogram[i];
				}
			}
			else{
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] = 0.0;
				}
			}

			if (sHandle->wp.globalHistogramParams.recordTime) {//(MY) removed +
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] = sHandle->tmpGlobalHistograms[m].timeHistogram[i];
				}
			}
			else{
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] = 0.0;
				}
			}
		}


	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			initbufftozero(nbMoments, buffer,&f);

			if (f.hitted) {

				for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					facetHitBuffer->hit.nbAbsEquiv = f.tmpCounter[m].hit.nbAbsEquiv;
					facetHitBuffer->hit.nbDesorbed = f.tmpCounter[m].hit.nbDesorbed;
					facetHitBuffer->hit.nbMCHit = f.tmpCounter[m].hit.nbMCHit;
					facetHitBuffer->hit.nbHitEquiv = f.tmpCounter[m].hit.nbHitEquiv;
					facetHitBuffer->hit.sum_1_per_ort_velocity = f.tmpCounter[m].hit.sum_1_per_ort_velocity;
					facetHitBuffer->hit.sum_v_ort = f.tmpCounter[m].hit.sum_v_ort;
					facetHitBuffer->hit.sum_1_per_velocity = f.tmpCounter[m].hit.sum_1_per_velocity;
				}

				if (f.sh.isProfile) {//(MY) removed +
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							shProfile[j] = f.profile[m][j];
						}
					}
				}

				if (f.sh.isTextured) {//(MY)
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						//double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Add temporary hit counts
								shTexture[add] = f.texture[m][add]; //(My) removed +
								/* Will be done in SimulationMCmain.cpp
								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								//TODO (MY) adapt this part of code;

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f.textureCellIncrements[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale (MY) how to solve f.largeEnough?
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f.largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only) (MY)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}*/
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
								shDir[add].dir.x = f.direction[m][add].dir.x;
								shDir[add].dir.y = f.direction[m][add].dir.y;
								shDir[add].dir.z = f.direction[m][add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								shDir[add].count = f.direction[m][add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {//(MY) removed +
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							shAngleMap[add] = f.angleMap.pdf[add];
						}
					}
				}

				//Facet histograms

					for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +
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

		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
	}

	//ReleaseDataport(dpHit); // (Rudi) Don't need that.

	ResetTmpCounters();
	//extern char* GetSimuStatus();
	//SetState(NULL, GetSimuStatus(), false, true); // (Rudi) Don't need that.

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}

void initbufftozero(size_t nbMoments, BYTE *buffer,SubprocessFacet *f){
	//here: init values with zero
	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	int j, x, y;

	for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +, initializing for else?
						FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
						facetHitBuffer->hit.nbAbsEquiv = 0.0;
						facetHitBuffer->hit.nbDesorbed = 0;
						facetHitBuffer->hit.nbMCHit = 0;
						facetHitBuffer->hit.nbHitEquiv = 0.0;
						facetHitBuffer->hit.sum_1_per_ort_velocity = 0.0;
						facetHitBuffer->hit.sum_v_ort = 0.0;
						facetHitBuffer->hit.sum_1_per_velocity = 0.0;
					}
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {
							ProfileSlice *shProfile = (ProfileSlice *)(buffer + f->sh.hitOffset + facetHitsSize + m * f->profileSize);
							for (j = 0; j < (int)PROFILE_SIZE; j++) {
								shProfile[j].countEquiv=0.0; shProfile[j].sum_1_per_ort_velocity=0.0; shProfile[j].sum_v_ort=0.0;
							}
						}
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {
		TextureCell *shTexture = (TextureCell *)(buffer + (f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + m * f->textureSize));
		for (y = 0; y < (int)f->sh.texHeight; y++) {
			for (x = 0; x < (int)f->sh.texWidth; x++) {
				size_t add = x + y * f->sh.texWidth;
				shTexture[add].countEquiv=0.0; shTexture[add].sum_1_per_ort_velocity=0.0; shTexture[add].sum_v_ort_per_area=0.0;
			}
		}
	}
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {
		DirectionCell *shDir = (DirectionCell *)(buffer + (f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*m));
		for (y = 0; y < (int)f->sh.texHeight; y++) {
			for (x = 0; x < (int)f->sh.texWidth; x++) {
				size_t add = x + y * f->sh.texWidth;
				shDir[add].dir.x = 0.0;
				shDir[add].dir.y = 0.0;
				shDir[add].dir.z = 0.0;
				//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
				shDir[add].count = 0;
			}
		}
	}
	size_t *shAngleMap = (size_t *)(buffer + f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*(1 + nbMoments));
	for (y = 0; y < (int)(f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes); y++) {
		for (x = 0; x < (int)f->sh.anglemapParams.phiWidth; x++) {
			size_t add = x + y * f->sh.anglemapParams.phiWidth;
			shAngleMap[add] = 0;
		}
	}
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {//(MY) removed +, initialized for else
		BYTE *histCurrentMoment = buffer + f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*(1 + nbMoments) + f->sh.anglemapParams.GetRecordedDataSize() + m * f->sh.facetHistogramParams.GetDataSize();

		double* nbHitsHistogram = (double*)histCurrentMoment;
		double* distanceHistogram = (double*)(histCurrentMoment + f->sh.facetHistogramParams.GetBouncesDataSize());
		double* timeHistogram = (double*)(histCurrentMoment + f->sh.facetHistogramParams.GetBouncesDataSize() + f->sh.facetHistogramParams.GetDistanceDataSize());


			for (size_t i = 0; i < f->sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
				nbHitsHistogram[i] = 0.0;
			}

			for (size_t i = 0; i < (f->sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
				distanceHistogram[i] = 0.0;
			}

			for (size_t i = 0; i < (f->sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
				timeHistogram[i] = 0.0;
			}

	}

	return;
}