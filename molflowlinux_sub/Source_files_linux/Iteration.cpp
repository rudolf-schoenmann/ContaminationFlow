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
#include <math.h>
//#include "levmar.h"
extern SimulationHistory* simHistory;
extern Simulation *sHandle;
extern ProblemDef *p;


//Update Error lists
std::tuple<double,double,double,double> getErrorVariables(SubprocessFacet* f, Databuff *hitbuffer_sum){
	if(hitbuffer_sum==NULL){
		return std::make_tuple(f->tmpCounter[0].nbHitEquiv, (double)f->tmpCounter[0].nbDesorbed, (double)f->tmpCounter[0].nbOutgassed, f->tmpCounter[0].nbAbsEquiv);
	}
	else{
		return std::make_tuple(getHits(f,hitbuffer_sum),(double)getnbDesorbed(f, hitbuffer_sum),(double)getnbOutgassed(f, hitbuffer_sum),getnbAdsorbed(f,hitbuffer_sum));
	}
}

void UpdateErrorList(Databuff *hitbuffer_sum){ // hitbuffer_sum==NULL: subprocess error, else: main error
	double num_hit_it=0;
	double num_des_ad_it=0;
	double nbhits=0.0; double nbdes=0.0; double nbout=0.0;double nbads=0.0;

	double factor=hitbuffer_sum==NULL?pow(simHistory->numSubProcess,0.5):1.0; // To be consistent with the ignored facets for calculating the error after summation over all subprocesses, here the hitRationLimit must be reduced with the correction factor due to multiple subprocesses.

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) { //save current num total hits in currentList, add difference current-old to num_hit_it
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			std::tie(nbhits,nbdes,nbout,nbads)=getErrorVariables(&f, hitbuffer_sum);
			num_hit_it+=(nbhits + nbdes + nbout);
			num_des_ad_it+=(nbads + nbdes);
		}
	}

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			std::tie(nbhits,nbdes,nbout,nbads)=getErrorVariables(&f, hitbuffer_sum);
			double num_hit_f=(nbhits + nbdes + nbout);
			double num_des_ad_f=(nbads + nbdes);

			//neglect events/covering change if very small compared to total hits
			if(num_hit_f/num_hit_it<(p->hitRatioLimit)/factor){
				num_hit_it-=num_hit_f;
				num_hit_f=0;
			}
			if(num_des_ad_f/num_des_ad_it<(p->hitRatioLimit)/factor){
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
			if(hitbuffer_sum==NULL){ // Sub process: set currentList values
				simHistory->hitList.setCurrent(&f,nbhits);
				simHistory->desorbedList.setCurrent(&f,nbdes);
				simHistory->outgassedList.setCurrent(&f,nbout);
				simHistory->adsorbedList.setCurrent(&f,nbads);
			}
			else{ // Main process: set historyList values
				simHistory->hitList.setLast(&f,nbhits);
				simHistory->desorbedList.setLast(&f,nbdes);
				simHistory->outgassedList.setLast(&f,nbout);
				simHistory->adsorbedList.setLast(&f,nbads);
			}
		}
	}
}

// Calculate total error
std::tuple<double,double> CalcErrorAll(int it){//calculates the averaged total error weighted with the facets area to decide, if the desired uncertainty level is reached
	double error_event=0.0; double error_covering=0.0;
	double area_event=0.0; double area_covering=0.0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			int idx=getFacetIndex(&f);
			if(f.sh.opacity==0||(p->doFocusGroupOnly && std::find(std::begin(p->focusGroup.second),std::end(p->focusGroup.second),idx)==std::end(p->focusGroup.second))) {continue;} //ignore facet if no opacity

			double error_event_facet=(it==-1)?simHistory->errorList_event.getCurrent(idx):simHistory->errorList_event.historyList.second[idx][it];
			if(!std::isinf(error_event_facet)){//ignore facet if no hits (=inf error)
				error_event+=error_event_facet*pow(f.sh.area,3/2);
				area_event+=f.sh.area;
			}

			double error_covering_facet=(it==-1)?simHistory->errorList_covering.getCurrent(idx):simHistory->errorList_covering.historyList.second[idx][it];
			if(!std::isinf(error_covering_facet)){//ignore facet if no hits (=inf error)
				error_covering+=error_covering_facet*pow(f.sh.area,3/2);
				area_covering+=f.sh.area;
			}
		}
	}
	return std::make_tuple(error_event/(pow(area_event,3/2)), error_covering/(pow(area_covering,3/2)));
}

// Check if error reached targets
bool checkError(double targetError, double currentError, double factor, std::string mode){
	bool vipCheck = currentError<=targetError;
	if(!p->vipFacets.empty()){
		HistoryList<double> *listptr;
		listptr = getErrorList(mode);
		if(listptr==NULL){
			return true;
		}

		for(unsigned int i = 0; i < p->vipFacets.size(); i++){
			//vipCheck = vipCheck && (listptr->getCurrent(p->vipFacets[i].first)== std::numeric_limits<double>::infinity() || listptr->getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
			vipCheck = vipCheck && (listptr->getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
		}
	}
	return vipCheck;
}

HistoryList<double>* getErrorList(std::string mode){
	if(mode=="covering"){
		return &simHistory->errorList_covering;
	}
	else if(mode=="event"){
		return &simHistory->errorList_event;
	}
	else{
		std::cout<<"------------! Error mode '"<<mode <<"' not implemented !------------\n";
		return NULL;
	}
}

//-----------------------------------------------------------
//Functions for covering threshold
void initCoveringThresh(){
	// Initializes size of cveringThreshold
	sHandle->coveringThreshold=std::vector<llong> ();
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for(uint i=0; i<sHandle->structures[j].facets.size(); i++){
			sHandle->coveringThreshold.push_back(0);
		}
	}
}

void setCoveringThreshold(int size, int rank){
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet per subprocess
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(rank==1)
				{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim + simHistory->coveringList.getCurrent(&f).convert_to<llong>()%num_sim;}
			else
				{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim;}

			// Covering threshold: current covering of face - particles that can be desorbed
			sHandle->coveringThreshold[getFacetIndex(&f)]=simHistory->coveringList.getCurrent(&f).convert_to<llong>()-cov_sim;
		}
	}
}

//-----------------------------------------------------------

//Estimation of Tmin
double estimateTmin(Databuff *hitbuffer){ //not ready yet => finish //TODO
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
				boost::multiprecision::uint128_t covering = getFacetHitBuffer(&f,hitbuffer)->covering;
				double v_avg_therm = 3*pow((3.14159265359/8),0.5)* pow((8.314472* f.sh.temperature/(sHandle->wp.gasMass*0.001)),0.5); //0.001 to convert MolarMass form g to kg
				sum_v_avg += v_avg_therm * (f.sh.outgassing + f.sh.desorption.convert_to<double>());
				normalization_factor_v += f.sh.outgassing + f.sh.desorption.convert_to<double>();
				if ((f.sh.outgassing + f.sh.desorption.convert_to<double>()) > 0){ //avoid division by 0
					double ttemp= (double)((double) covering/((f.sh.outgassing + f.sh.desorption.convert_to<double>())/ (kb*f.sh.temperature)));
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
	double av_path_length = sHandle->tmpGlobalResult.distTraveled_total*0.01/sHandle->tmpGlobalResult.globalHits.nbMCHit;// The 0.01 is because the distance should be saved in cm. ToDo: Check up, if that's right.
	std::cout <<  "av_path_length = sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.nbMCHit = " << sHandle->tmpGlobalResult.distTraveled_total << "/" << sHandle->tmpGlobalResult.globalHits.nbMCHit << std::endl;
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


//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
void setCoveringThreshold(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{cov_sim=facetHitBuffer->covering/num_sim + facetHitBuffer->covering%num_sim;}
				else
					{cov_sim=facetHitBuffer->covering/num_sim;}

				sHandle->coveringThreshold[getFacetIndex(&f)]=facetHitBuffer->covering-cov_sim;
			}
	}
}
*/

//----------Not used anymore
/*
void allocateCovering(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong temp=size+1;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{facetHitBuffer->covering=facetHitBuffer->covering/temp + facetHitBuffer->covering%temp;}
				else
					{facetHitBuffer->covering=facetHitBuffer->covering/temp;}
			}
	}
}
*/
