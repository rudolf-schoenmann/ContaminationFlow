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

double estimateTmin_RudiTest(){
	double tmin=1;
	// TODO Code eine Sekunde lang simulieren lassen.
	//const double distribution_factor = 3*(3.14159265359/8)^0.5*(8,314472/MolareGasMasse in (kg/mol))^0.5;
	double sum_v_avg = 0;
	double normalization_factor_v = 0;
	double sum_sc = 0;
	double normalization_factor_sc = sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	//avarage of <v> (<v> is the average velocity on a single facet depending on temperature)
	// over all facets weighted with the rate of outgoing particles (outgassing + desorption) per facet
	size_t nbMoments = sHandle->moments.size();
	/*
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			for (size_t m = 0; m <= nbMoments; m++) {

				//avarage of <v>
				sum_v_avg += distribution_factor*(f.sh.temperature)^0.5 * (f.sh.outgassing + Desorptionsrate);
				normalization_factor_v += (f.sh.outgassing + Desorptionsrate)

				//average sticking coefficient
				sum_sc += f.sh.sticking*f.tmpCounter.hit.nbMCHit;
			}
		}
	}*/
	//double avg_v_avg = sum_v_avg/normalization_factor_v;
	//double avg_sc = sum_sc/normalization_factor_sc;
	//double avg_path_length = sHandle->tmpGlobalResult.distTraveled_total / sHandle->tmpGlobalResult.globalHits.hit.nbMCHit/ avg_sc;
	//tmin = vg_path_length /avg_v_avg;
	return tmin;

}

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
	//double dmean=c/sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;//m
	double dmean=sHandle->tmpGlobalResult.distTraveled_total;//m
	//double dmean= sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly;
	std::cout <<"test " <<dmean <<std::endl;

	tmin=dmean*sum_1_v_ort*Nmean;//s
	return tmin*1000;//ms
}

TimeTest::TimeTest(){
	pointintime_list=std::vector< std::pair<double,std::vector<double>> >();
}

void TimeTest::appendList(double time){

	std::vector<double> currentstep;
	currentstep =std::vector<double> ();

	double covering;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering=calcRealCovering(&f);
			currentstep.push_back(covering);
		}
	}
	pointintime_list.push_back(std::make_pair(time,currentstep));

}

void::TimeTest::print(){
	for(int i=0;i<pointintime_list.size();i++)
	{
		std::cout <<pointintime_list[i].first;

		for(int j=0; j<pointintime_list[i].second.size();j++)
		{
			std::cout <<'\t' <<pointintime_list[i].second[j];
		}

		std::cout <<std::endl;

	}
}
