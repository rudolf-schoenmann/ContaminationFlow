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
 * This file contains calculations for the iterative algorithm
 */

#include "SimulationLinux.h"
#include "levmar.h"
extern SimulationHistory* simHistory;
extern Simulation *sHandle;

//Estimation of Tmin
double estimateTmin(Databuff *hitbuffer){ //not ready yet => finish //TODO
BYTE *buffer;
buffer = hitbuffer->buff;
//Ich muss der Funktion noch einen Hitbuffer übergeben. Ich brauche ja 'covering'.
	double tmin=0;
	double sum_v_avg = 0;
	double normalization_factor_v = 0;
	//llong covering;
	//avarage of <v> (<v> is the average velocity on a single facet depending on temperature)
    // over all facets weighted with the rate of outgoing particles (outgassing + desorption) per facet
	double tmin_particles_out = 0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				llong covering = facetHitBuffer->hit.covering;
				double v_avg_therm = 3*pow((3.14159265359/8),0.5)* pow((8.314472* f.sh.temperature/(sHandle->wp.gasMass*0.001)),0.5); //0.001 to convert MolarMass form g to kg
				sum_v_avg += v_avg_therm * (f.sh.outgassing + f.sh.desorption.convert_to<double>());
				normalization_factor_v += f.sh.outgassing + f.sh.desorption.convert_to<double>();
				if ((f.sh.outgassing + f.sh.desorption.convert_to<double>()) > 0){ //avoid division by 0
					double ttemp= (double)(covering/((f.sh.outgassing + f.sh.desorption.convert_to<double>())/ (1.38E-23*f.sh.temperature)));
					if (!tmin_particles_out){
						tmin_particles_out = (ttemp);
					}
					if (tmin_particles_out > ttemp){
						tmin_particles_out = ttemp;}
				}				
			 }
	}
	//if (normalization_factor_v == 0) std::cout << "divide by zero!" << std::endl;
	double v_avg = sum_v_avg/normalization_factor_v;
	std::cout << "v_avg = sum_v_avg/normalization_factor_v = " << sum_v_avg << "/" << normalization_factor_v << std::endl;
	double av_path_length = sHandle->tmpGlobalResult.distTraveled_total*0.01/sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;// The 0.01 is because the distance should be saved in cm. ToDo: Check up, if that's right.
	std::cout <<  "av_path_length = sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbMCHit = " << sHandle->tmpGlobalResult.distTraveled_total << "/" << sHandle->tmpGlobalResult.globalHits.hit.nbMCHit << std::endl;
	std::cout << "tmin = av_path_length /v_avg = " << av_path_length << " / " << v_avg << std::endl;
	tmin = av_path_length /v_avg;

	//Time step information for developing

	std::cout << "Time step information for developing" << std::endl;
	std::cout << "estimateTmin: tmin = " <<tmin<< "ms"<< std::endl;
	std::cout << "estimateTmin: tmin_particles_out = " <<tmin_particles_out<< "ms"<< std::endl;
	/*We don't need that anymore because we solved the problem of too much desorbed particles via a threshold value.
	if (tmin < tmin_particles_out){
		std::cout << "estimateTmin: tmin = " <<tmin<< "ms is chosen."<< std::endl;
	 	return tmin;
		}
	else{
		std::cout << "estimateTmin: tmin_particles_out = " <<tmin_particles_out<< "ms is chosen."<< std::endl;
		return tmin_particles_out;
		}
	*/
	return tmin;
}


double estimateAverageFlightTime(){
	return simHistory->flightTime/simHistory->nParticles;
}

//-----------------------------------------------------------
//Functions for covering threshold

void allocateCovering(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong temp=size+1;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{facetHitBuffer->hit.covering=facetHitBuffer->hit.covering/temp + facetHitBuffer->hit.covering%temp;}
				else
					{facetHitBuffer->hit.covering=facetHitBuffer->hit.covering/temp;}
			}
	}
}


void initCoveringThresh(){
	sHandle->coveringThreshold=std::vector<llong> ();
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for(uint i=0; i<sHandle->structures[j].facets.size(); i++){
					sHandle->coveringThreshold.push_back(0);
				}
		}
}


void setCoveringThreshold(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{cov_sim=facetHitBuffer->hit.covering/num_sim + facetHitBuffer->hit.covering%num_sim;}
				else
					{cov_sim=facetHitBuffer->hit.covering/num_sim;}

				sHandle->coveringThreshold[getFacetIndex(&f)]=facetHitBuffer->hit.covering-cov_sim;
			}
	}
}

void setCoveringThreshold(int size, int rank){
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				if(rank==1)
					{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim + simHistory->coveringList.getCurrent(&f).convert_to<llong>()%num_sim;}
				else
					{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim;}

				sHandle->coveringThreshold[getFacetIndex(&f)]=simHistory->coveringList.getCurrent(&f).convert_to<llong>()-cov_sim;
			}
	}
}

