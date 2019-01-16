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
#include <fstream>
#include <sstream>

extern Simulation *sHandle;
/*
double estimateTmin(){//TODO something is wrong here
	double tmin=1;
	int facetcounter=0;
	double sum_1_v_ort=0.0;
	double sum_abs=0.0;
	//double sum_hits=0.0;

	//mittelwert 1/v_ort
	size_t nbMoments = sHandle->moments.size();
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			for (size_t m = 0; m <= nbMoments; m++) {
				//sum_1_v_orth
				//sum_1_v_ort+=f.tmpCounter[m].hit.sum_1_per_ort_velocity*f.sh.area; //(s/m)*m^2
				sum_1_v_ort+=f.tmpCounter[m].hit.sum_v_ort*f.sh.area*1E-4; //(s/m)*m^2
				facetcounter++;

				//covering
				double coveringtemp=f.tmpCounter[m].hit.covering;
				sum_abs+=calcNmono(&f)*coveringtemp; //mass
				//sum_hits+=f.tmpCounter[m].hit.nbHitEquiv;//or nbMChit
			}
		}
	}
	sum_1_v_ort=1/sum_1_v_ort;
	sum_1_v_ort /=facetcounter;//(s*m)
	//double Nmean=sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv/sum_abs;//(1/mass) //nbHitEquiv or nbMCHit
	double Nmean=1/sum_abs;//(1/mass) //nbHitEquiv or nbMCHit
	std::cout <<"test " <<sum_1_v_ort <<std::endl;
	std::cout <<"test " <<Nmean <<std::endl;

	//TODO durchschnittstrecke
	//double dmean=sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;//m
	double dmean=sHandle->tmpGlobalResult.distTraveled_total;//m
	//double dmean= sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly;
	std::cout <<"test " <<dmean <<std::endl;

	tmin=dmean*sum_1_v_ort*Nmean;//s
	return tmin*1000;//ms
}*/

double estimateTmin(){
	double sum_v_ort=0;
	double sum_1_v_ort=0;
	double facetcounter=0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
				//sum_1_v_orth
				//sum_1_v_ort+=f.tmpCounter[m].hit.sum_1_per_ort_velocity*f.sh.area; //(s/m)*m^2
				sum_v_ort+=f.tmpCounter[0].hit.sum_v_ort; //(s/m)*m^2
				sum_1_v_ort+=f.tmpCounter[0].hit.sum_1_per_velocity;
				facetcounter++;

		}
	}

	double dist_total=(double)sHandle->tmpGlobalResult.distTraveled_total;
	double hits2= pow((double)sHandle->tmpGlobalResult.globalHits.hit.nbMCHit,2);
	std::cout <<"dist_total, hits^2, sum_v_ort, sum_1_v_ort: \t\t\t" <<dist_total <<'\t'<<hits2<<'\t'<<sum_v_ort<<'\t' <<sum_1_v_ort <<std::endl;
	std::cout <<"Alternative 1 for Tmin\t (dist_total*1000)/sum_v_ort [ms]\t" <<dist_total*1000/sum_v_ort <<std::endl;
	std::cout <<"Alternative 1 for Tmin\t(dist_total/hits^2)*1000/sum_v_ort [ms]\t" <<(dist_total/hits2)*1000/sum_v_ort <<std::endl;
	return (dist_total/hits2)*sum_1_v_ort*1000;
}

CoveringHistory::CoveringHistory(){
	pointintime_list=std::vector< std::pair<double,std::vector<double>> >();
}

CoveringHistory::CoveringHistory(Databuff *hitbuffer){
	pointintime_list=std::vector< std::pair<double,std::vector<double>> >();
	std::vector<double> currentstep;
	currentstep =std::vector<double> ();

	BYTE *buffer;
	buffer = hitbuffer->buff;

	double covering;
	std::cout <<"Reading covering values from buffer\t";
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				covering = facetHitBuffer->hit.covering;
				std::cout <<covering <<"\t";
				currentstep.push_back(covering);
				f.tmpCounter[0].hit.covering=covering;
			}
	}
	std::cout <<std::endl;
	pointintime_list.push_back(std::make_pair(0.0,currentstep));
}

void CoveringHistory::appendList(double time){

	std::vector<double> currentstep;
	currentstep =std::vector<double> ();

	double covering;

	int i=0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//covering=calcRealCovering(&f);
			covering=calcCovering(&f);
			currentstep.push_back(covering);
			i+=1;
		}
	}
	pointintime_list.push_back(std::make_pair(time,currentstep));

}

void::CoveringHistory::print(){

	std::cout <<"time";
	for(int i=0;i<pointintime_list.size();i++)
	{
		if(i==0){
			for(int j=0; j<pointintime_list[i].second.size();j++)
					{
						std::cout <<"\t Covering for Facet " <<j;
					}
		}
		std::cout<<std::endl;
		std::cout <<pointintime_list[i].first;

		for(int j=0; j<pointintime_list[i].second.size();j++)
		{
			std::cout <<"\t\t" <<pointintime_list[i].second[j];
			if(pointintime_list[i].second[j] == 0.0)
				std::cout <<"\t";
		}

	}
	std::cout<<std::endl<<std::endl;
}

void::CoveringHistory::write(std::string filename){
	//std::string write = "/home/van/history"+std::to_string(num)+".txt";
	std::ofstream outfile(filename,std::ofstream::out|std::ios::trunc);

	for(int i=0;i<pointintime_list.size();i++)
	{
		outfile <<pointintime_list[i].first;

		for(int j=0; j<pointintime_list[i].second.size();j++)
		{
			outfile <<'\t' <<pointintime_list[i].second[j];
		}

		outfile <<'\n';

	}
	outfile.close();
}

void::CoveringHistory::read(std::string filename){
	pointintime_list.clear();
	//pointintime_list_read.clear();
	//std::string read = "/home/van/history"+std::to_string(num)+".txt";
	std::string line;

	std::ifstream input(filename,std::ifstream::in);
	std::cout <<"Reading in covering history from " <<filename <<std::endl;
	while(std::getline(input,line)){
		std::vector<double> currentstep;
		currentstep =std::vector<double> ();

		double covering;
		double time;
		std::istringstream is( line );

		is >> time;
		while(!is.eof()){
			is >> covering;
			currentstep.push_back(covering);

		}
		pointintime_list.push_back(std::make_pair(time,currentstep));
	}
	input.close();

	int i=0;
	double num_mol=0.0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
				num_mol=pointintime_list.back().second[i]/calcCoveringUpdate(&f);
				f.tmpCounter[0].hit.covering = pointintime_list.back().second[i];
				calcStickingnew(&f); // calculate new sticking for new covering value

				std::cout <<"Facet "<<i <<"\t covering: " <<f.tmpCounter[0].hit.covering <<"\t Corresponding number of particles: " <<num_mol <<std::endl;
				i+=1;
		}
	}
	std::cout <<std::endl;

}
