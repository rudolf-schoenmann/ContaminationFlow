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
 * This file contains functions of the reduced worker class. Functions needed for iterative algorithm.
 */

#include "worker.h"
#include "Simulation.h"
#include "SimulationContaminationFlow.h"
#include <vector>
#include "MathTools.h"
#include "Random.h"
extern Simulation *sHandle; //delcared in molflowSub.cpp
extern SimulationHistory* simHistory;

void CalcTotalOutgassingWorker() {
	// Compute the outgassing of all source facet
	sHandle->wp.totalDesorbedMolecules = sHandle->wp.finalOutgassingRate_Pa_m3_sec = sHandle->wp.finalOutgassingRate = sHandle->wp.totalOutgassingParticles= 0.0;
	//double time_step = simHistory->stepSize;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			//if(f.sh.temperature==0) {continue;}

			if (f.sh.desorbType != DES_NONE) { //there is a kind of desorption
				if (f.sh.useOutgassingFile) { //outgassing file
					//this is the Molflow+ original code block. Modifications for ContaminationFlow necessary???
					for (unsigned int l = 0; l < (f.sh.outgassingMapWidth*f.sh.outgassingMapHeight); l++) {
						sHandle->wp.totalDesorbedMolecules += sHandle->wp.latestMoment * f.outgassingMap[l] / (kb*f.sh.temperature);
						sHandle->wp.finalOutgassingRate += f.outgassingMap[l] / (kb*f.sh.temperature);
						sHandle->wp.finalOutgassingRate_Pa_m3_sec += f.outgassingMap[l];

						sHandle->wp.totalOutgassingParticles += f.outgassingMap[l] * calcOutgassingFactor(&f); // Number of particles outgassing in time stepSize_outgassing

						//Modifications like in the regular outgassing case necessary!?!
					}
				}
				else { //regular outgassing
					if (f.sh.outgassing_paramId == -1) { //constant outgassing
						//This following three lines are still the old code.
						sHandle->wp.totalDesorbedMolecules += sHandle->wp.latestMoment * f.sh.outgassing / (kb*f.sh.temperature);
						sHandle->wp.finalOutgassingRate += f.sh.outgassing / (kb*f.sh.temperature);  //Outgassing molecules/sec
						sHandle->wp.finalOutgassingRate_Pa_m3_sec += f.sh.outgassing;
						sHandle->wp.totalOutgassingParticles +=f.sh.outgassing * calcOutgassingFactor(&f); // Number of particles outgassing in time stepSize_outgassing
						}
					else {
						//time-dependent outgassing
						double time_step = simHistory->stepSize; //length of the current iteration
						double t_start =simHistory->lastTime; //start of iteration
						double t_stop =t_start + time_step;	//end of iteration
						double end_of_outgassing = sHandle->IDs[f.sh.IDid].back().first; //last point of the defined and loaded outgassing table
						double start_of_outgassing = sHandle->IDs[f.sh.IDid].front().first; //first point of the defined and loaded outgassing table
						double outgassing_start = 0;//start of outgassing within the iteration step
						double outgassing_end = 0;//end of outgassing within the iteration step
						double facet_outgassing = 0;//over time integrated outgassing of facet during the iteration step
						if(t_start >= end_of_outgassing){//case, when t_start is after the last point in time, where an outgassing is defined
							continue;
						}
						else if(t_start <= end_of_outgassing && t_stop >= end_of_outgassing){//case, when t_start is before and
							//t_stopp is after the last point in time, where an outgassing is defined
							if (t_start <= start_of_outgassing){
								outgassing_start = start_of_outgassing;
								outgassing_end = end_of_outgassing;
							}
							else{
								outgassing_start = t_start;
								outgassing_end = end_of_outgassing;
							}
						}
						else{
						//same as =>else if(t_start <= end_of_outgassing && t_stop <= end_of_outgassing){//case, when t_start is before and
							//t_stopp is before the last point in time, where an outgassing is defined
							if (t_start <= start_of_outgassing){
								if(t_stop <= start_of_outgassing){
									continue;
								}
								else{//t_stop >= start_of_outgassing
									outgassing_start = start_of_outgassing;
									outgassing_end = t_stop;
								}
							}
							else{//t_start >= start_of_outgassing
								outgassing_start = t_start;
								outgassing_end = t_stop;
							}
						}
						facet_outgassing = InterpolateY(outgassing_end, sHandle->IDs[f.sh.IDid], false, true) - InterpolateY(outgassing_start, sHandle->IDs[f.sh.IDid], false, true);
						sHandle->wp.totalOutgassingParticles +=facet_outgassing / (kb*f.sh.temperature);
						sHandle->wp.totalDesorbedMolecules +=facet_outgassing / (kb*f.sh.temperature);
					}
				}
			}
		}
	}
}

//----------Not used anymore
/*
std::vector<std::pair<double, double>> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size){
	std::vector<std::pair<double, double>> cdf; cdf.reserve(size);
	double Kb = 1.38E-23;
	double R = 8.3144621;
	double a = sqrt(Kb*gasTempKelvins / (gasMassGramsPerMol*1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

	//Generate cumulative distribution function
	double mostProbableSpeed = sqrt(2 * R*gasTempKelvins / (gasMassGramsPerMol / 1000.0));
	double binSize = 4.0*mostProbableSpeed / (double)size; //distribution generated between 0 and 4*V_prob

	for (size_t i = 0; i < size; i++) {
		double x = (double)i*binSize;
		double x_square_per_2_a_square = pow(x, 2) / (2 * pow(a, 2));
		cdf.push_back(std::make_pair(x, 1 - exp(-x_square_per_2_a_square)*(x_square_per_2_a_square + 1)));

	}


	return cdf;
}

int GetCDFId(double temperature) {

	int i;
	for (i = 0; i<(int)sHandle->temperatures.size() && (abs(temperature - (double)(sHandle->temperatures[i]))>1E-5); i++); //check if we already had this temperature
	if (i >= (int)sHandle->temperatures.size()) i = -1; //not found
	return i;
}

int GenerateNewCDF(double temperature){
	size_t i = sHandle->temperatures.size();
	sHandle->temperatures.push_back(temperature);
	sHandle->CDFs.push_back(Generate_CDF(temperature, sHandle->wp.gasMass, CDF_SIZE));
	return (int)i;
}

std::vector<std::pair<double, double>> Generate_ID(int paramId){
	std::vector<std::pair<double, double>> ID;
	//First, let's check at which index is the latest moment
	size_t indexBeforeLastMoment;
	for (indexBeforeLastMoment = 0; indexBeforeLastMoment < sHandle->parameters[paramId].GetSize() &&
		(sHandle->parameters[paramId].GetX(indexBeforeLastMoment) < sHandle->wp.latestMoment); indexBeforeLastMoment++);
		if (indexBeforeLastMoment >= sHandle->parameters[paramId].GetSize()) indexBeforeLastMoment = sHandle->parameters[paramId].GetSize() - 1; //not found, set as last moment

	//Construct integral from 0 to latest moment
	//Zero
	ID.push_back(std::make_pair(0.0, 0.0));

	//First moment
	ID.push_back(std::make_pair(sHandle->parameters[paramId].GetX(0),
			sHandle->parameters[paramId].GetX(0)*sHandle->parameters[paramId].GetY(0)*0.100)); //for the first moment (0.1: mbar*l/s -> Pa*m3/s)

	//Intermediate moments
	for (size_t pos = 1; pos <= indexBeforeLastMoment; pos++) {
		if (IsEqual(sHandle->parameters[paramId].GetY(pos) , sHandle->parameters[paramId].GetY(pos-1))) //two equal values follow, simple integration by multiplying
			ID.push_back(std::make_pair(sHandle->parameters[paramId].GetX(pos),
			ID.back().second +
			(sHandle->parameters[paramId].GetX(pos) - sHandle->parameters[paramId].GetX(pos-1))*sHandle->parameters[paramId].GetY(pos)*0.100));
		else { //difficult case, we'll integrate by dividing to 20 equal sections
			for (double delta = 0.05; delta < 1.0001; delta += 0.05) {
				double delta_t = sHandle->parameters[paramId].GetX(pos) - sHandle->parameters[paramId].GetX(pos-1);
				double time = sHandle->parameters[paramId].GetX(pos-1) + delta*delta_t;
				double avg_value = (sHandle->parameters[paramId].InterpolateY(time - 0.05*delta_t,false) + sHandle->parameters[paramId].InterpolateY(time,false))*0.100 / 2.0;
				ID.push_back(std::make_pair(time,
					ID.back().second +
					0.05*delta_t*avg_value));
			}
		}
	}

	//wp.latestMoment
	double valueAtlatestMoment = sHandle->parameters[paramId].InterpolateY(sHandle->wp.latestMoment,false);
	if (IsEqual(valueAtlatestMoment , sHandle->parameters[paramId].GetY(indexBeforeLastMoment))) //two equal values follow, simple integration by multiplying
		ID.push_back(std::make_pair(sHandle->wp.latestMoment,
		ID.back().second +
		(sHandle->wp.latestMoment - sHandle->parameters[paramId].GetX(indexBeforeLastMoment))*sHandle->parameters[paramId].GetY(indexBeforeLastMoment)*0.100));
	else { //difficult case, we'll integrate by dividing two 5equal sections
		for (double delta = 0.0; delta < 1.0001; delta += 0.05) {
			double delta_t = sHandle->wp.latestMoment - sHandle->parameters[paramId].GetX(indexBeforeLastMoment);
			double time = sHandle->parameters[paramId].GetX(indexBeforeLastMoment) + delta*delta_t;
			double avg_value = (sHandle->parameters[paramId].GetY(indexBeforeLastMoment)*0.100 + sHandle->parameters[paramId].InterpolateY(time, false)*0.100) / 2.0;
			ID.push_back(std::make_pair(time,
				ID.back().second +
				0.05*delta_t*avg_value));
		}
	}

	return ID;

}

int GenerateNewID(int paramId){
	size_t i = sHandle->desorptionParameterIDs.size();
	sHandle->desorptionParameterIDs.push_back(paramId);
	sHandle->IDs.push_back(Generate_ID(paramId));
	return (int)i;
}

int GetParamId(const std::string name) {
	int foundId = -1;
	for (int i = 0; foundId == -1 && i < (int)sHandle->parameters.size(); i++)
		if (name.compare(sHandle->parameters[i].name) == 0) foundId = i;
	return foundId;
}
*/
