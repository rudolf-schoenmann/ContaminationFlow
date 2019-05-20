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

double estimateTmin_RudiTest(Databuff *hitbuffer){ //not ready yet => finish //TODO
BYTE *buffer;
buffer = hitbuffer->buff;
//Ich muss der Funktion noch einen Hitbuffer übergeben. Ich brauche ja 'covering'.
	double tmin=0;
	double sum_v_avg = 0;
	double normalization_factor_v = 0;
	llong covering;
	//avarage of <v> (<v> is the average velocity on a single facet depending on temperature)
    // over all facets weighted with the rate of outgoing particles (outgassing + desorption) per facet
	double tmin_particles_out = 0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				llong covering = facetHitBuffer->hit.covering;
				double v_avg_therm = 3*pow((3.14159265359/8),0.5)* pow((8.314472* f.sh.temperature/(sHandle->wp.gasMass*0.001)),0.5); //0.001 to convert MolarMass form g to kg
				sum_v_avg += v_avg_therm * (f.sh.outgassing + f.sh.desorption);
				normalization_factor_v += f.sh.outgassing + f.sh.desorption;
				if ((f.sh.outgassing + f.sh.desorption) > 0){ //avoid division by 0
					double ttemp= (double)(covering/((f.sh.outgassing + f.sh.desorption)/ (1.38E-23*f.sh.temperature)));
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
	//std::cout << "v_avg = sum_v_avg/normalization_factor_v = " << sum_v_avg << "/" << normalization_factor_v << std::endl;
	double av_path_length = sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;//Muss geändert werden, da diese Werte Null sind im Prozess Null.
	//std::cout <<  "av_path_length = Handle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbMCHit = " << sHandle->tmpGlobalResult.distTraveled_total << "/" << sHandle->tmpGlobalResult.globalHits.hit.nbMCHit << std::endl;
	//std::cout << "tmin = av_path_length /v_avg = " << av_path_length << " / " << v_avg << std::endl;
	tmin = av_path_length /v_avg;

	//Time step information for developing

	std::cout << "Time step information for developing" << std::endl;
	std::cout << "estimateTmin_RudiTest: tmin = " <<tmin<< "ms"<< std::endl;
	std::cout << "estimateTmin_RudiTest: tmin_particles_out = " <<tmin_particles_out<< "ms"<< std::endl;
	std::cout << "estimateTmin: tmin would be " << estimateTmin() << "ms" << std::endl;


	if (tmin < tmin_particles_out){
		std::cout << "estimateTmin_RudiTest: tmin = " <<tmin<< "ms is chosen."<< std::endl;
	 	return tmin;
		}
	else{
		std::cout << "estimateTmin_RudiTest: tmin_particles_out = " <<tmin_particles_out<< "ms is chosen."<< std::endl;
		return tmin_particles_out;
		}
	return tmin;
}


#include <fstream>
#include <sstream>

extern Simulation *sHandle;

double estimateTminFlightTime(){
	return simHistory->flightTime/simHistory->nParticles;
}


double estimateTmin(){
	double sum_v_ort=0;
	double sum_1_v_ort=0;
	double facetcounter=0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
				//sum_1_v_orth
				//sum_1_v_ort+=f.tmpCounter[m].hit.sum_1_per_ort_velocity*f.sh.area;
				sum_v_ort+=f.tmpCounter[0].hit.sum_v_ort;
				sum_1_v_ort+=f.tmpCounter[0].hit.sum_1_per_velocity;
				facetcounter++;
		}
	}
	double dist_total=0.01*(double)sHandle->tmpGlobalResult.distTraveled_total; //factor 0.01 to convert distance from cm to m.
	//Muss geändert werden, da diese Werte Null sind im Prozess Null.
	double hits= (double)sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;//Muss geändert werden, da diese Werte Null sind im Prozess Null.
	double hits2= pow((double)sHandle->tmpGlobalResult.globalHits.hit.nbMCHit,2);//Muss geändert werden, da diese Werte Null sind im Prozess Null.
	double particlenumber = (double)sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;//Muss geändert werden, da diese Werte Null sind im Prozess Null.
	/*
	std::cout << "_______________________________________________________________________________________________________"<< std::endl;
	std::cout <<"Current version of estimate Tmin"<< std::endl;
	std::cout <<"dist total [m]\t\t" << dist_total << std::endl;
	std::cout <<"hits \t\t\t" << hits << std::endl;
	std::cout <<"hits^2 \t\t\t" << hits2 << std::endl;
	std::cout <<"sum_v_ort [m/s]\t\t" << sum_v_ort<< std::endl;
	std::cout <<"<v_ort>[m/s]\t\t"<< sum_v_ort/hits << std::endl;
	std::cout <<"1/<v_ort>[s/m]\t\t"<< hits/sum_v_ort<< std::endl;
	std::cout <<"sum_1_v_ort [s/m]\t" << sum_1_v_ort << std::endl;
	std::cout <<"<1/v_ort>[s/m]\t\t"<<sum_1_v_ort/hits << std::endl;
	std::cout <<"Alternative 1 for Tmin\t (dist_total*1000)/sum_v_ort [ms]\t" <<dist_total*1000/sum_v_ort <<std::endl;
	std::cout <<"Alternative 3 for Tmin\t(dist_total/hits^2)*1000*sum_1_v_ort [ms]\t" <<(dist_total/hits2)*1000*sum_1_v_ort <<std::endl;
	std::cout << "Currently used: Alternative 3" << std::endl;
	std::cout << "_______________________________________________________________________________________________________"<< std::endl;
	*/
	return (dist_total/hits2)*sum_1_v_ort*1000;
}

SimulationHistory::SimulationHistory(){
	numFacet=0;
	nParticles=-1;
	flightTime=0.0;
}

SimulationHistory::SimulationHistory(Databuff *hitbuffer){
	std::vector<llong> currentstep;
	currentstep =std::vector<llong> ();
	numFacet=0;
	nParticles=-1;
	flightTime=0.0;


	llong covering;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				covering = getCovering(&f, hitbuffer);
				currentstep.push_back(covering);
				f.tmpCounter[0].hit.covering=covering;
				numFacet+=1;
			}
	}
	coveringList.appendList(currentstep, 0);
	coveringList.initCurrent(numFacet);
}


void SimulationHistory::appendList(Databuff *hitbuffer, double time){

	std::vector<llong> currentstep;
	currentstep =std::vector<llong> ();

	if(time==-1.0) //TODO improve: here steps instead of time
		time=coveringList.pointintime_list.back().first+1.0;

	llong covering;

	//int i=0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering=getCovering(&f, hitbuffer);
			currentstep.push_back(covering);
			//i+=1;
		}
	}
	coveringList.appendList(currentstep, -1);

}

void SimulationHistory::print(){
	coveringList.print();
}


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
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
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

				sHandle-> coveringThreshold[getFacetIndex(&f)]=facetHitBuffer->hit.covering-cov_sim;
			}
	}
}

void TestMinCovering(Databuff *hitbuffer){
	BYTE *buffer;
	buffer = hitbuffer->buff;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				facetHitBuffer->hit.covering=100;
				//f.tmpCounter[0].hit.covering=facetHitBuffer->hit.covering;
			}
	}
}

