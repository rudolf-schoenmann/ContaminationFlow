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

extern Simulation *sHandle;

double estimateTmin_RudiTest(Databuff *hitbuffer){ //not ready yet => finish
BYTE *buffer;
buffer = hitbuffer->buff;
//Ich muss der Funktion noch einen Hitbuffer übergeben. Ich brauche ja 'covering'.
	double tmin=0;
	double sum_v_avg = 0;
	double normalization_factor_v = 0;
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
					if (!tmin_particles_out){
						tmin_particles_out = (covering/(f.sh.outgassing + f.sh.desorption)/ (1.38E-23*f.sh.temperature));
					}
					if (tmin_particles_out > (covering/(f.sh.outgassing +f.sh.desorption)/ (1.38E-23*f.sh.temperature))){
						tmin_particles_out = (covering/(f.sh.outgassing +f.sh.desorption)/ (1.38E-23*f.sh.temperature));}
				}				
			 }
	}
	double v_avg = sum_v_avg/normalization_factor_v;
	double av_path_length = sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	tmin = av_path_length /v_avg;

	if (tmin < tmin_particles_out)
	 	return tmin;
	else
		return tmin_particles_out;
	
	
	std::cout << "estimateTmin_RudiTest: tmin = " <<tmin<< "ms"<< std::endl;
	//std::cout << "estimateTmin_RudiTest: tmin_particles_out = " <<tmin_particles_out<< "ms"<< std::endl;
	std::cout << "_______________________________________________________________________________________________________"<< std::endl<<std::endl;
	
	return tmin;
}


#include <fstream>
#include <sstream>

extern Simulation *sHandle;


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
	double hits= (double)sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	double hits2= pow((double)sHandle->tmpGlobalResult.globalHits.hit.nbMCHit,2);
	double particlenumber = (double)sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
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

CoveringHistory::CoveringHistory(){
	pointintime_list=std::vector< std::pair<double,std::vector<double>> >();
}

CoveringHistory::CoveringHistory(Databuff *hitbuffer){
	pointintime_list=std::vector< std::pair<double,std::vector<double>> >();
	std::vector<double> currentstep;
	currentstep =std::vector<double> ();

	BYTE *buffer;
	buffer = hitbuffer->buff;

	llong covering;
	//std::cout <<"Reading covering values from buffer\t";
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				covering = facetHitBuffer->hit.covering;
				//std::cout <<covering <<"\t";
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

	std::cout <<"time\t";
	for(int i=0;i<pointintime_list.size();i++)
	{
		if(i==0){
			for(int j=0; j<pointintime_list[i].second.size();j++)
					{
					std::cout <<"\tCovering for Facet " <<j;
					}
		}
		std::cout<<std::endl;
		std::cout <<pointintime_list[i].first;

		for(int j=0; j<pointintime_list[i].second.size();j++)
		{
			std::cout <<"\t\t" <<pointintime_list[i].second[j]<<"\t";
			/*if(pointintime_list[i].second[j] == 0.0)
				std::cout <<"\t";*/
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

void::CoveringHistory::read(std::string filename, Databuff *hitbuffer){//Rudi: Not ready yet.
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
				num_mol=pointintime_list.back().second[i]; //Rudi: Maybe wrong, since we changed covering and introduced coverage.
				f.tmpCounter[0].hit.covering = pointintime_list.back().second[i];
				//calcStickingnew(&f, hitbuffer); // calculate new sticking for new covering value
				// 1) Update the hitbuffer with the last covering value
				// 2) calcStickingnew(&f, hitbuffer);
				std::cout <<"Facet "<<i <<"\t covering: " <<f.tmpCounter[0].hit.covering <<"\t Corresponding number of particles: " <<num_mol <<std::endl;
				i+=1;
		}
	}
	std::cout <<std::endl;

}
