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
 * This file contains the serialization of the Simulation
 */

#include "SimulationLinux.h"
#include <array>
//extern Simulation *sHandle;
extern CoveringHistory* covhistory;


bool simulateSub(Databuff *hitbuffer, int rank, int simutime){

	covhistory = new CoveringHistory (hitbuffer);
	double timestep=1000; // desired length per iteration for simulation, here hardcoded to 1 second
	double realtimestep; // actual time elapsed for iteration step

	// Set end of simulation flag
	bool eos=false;

	// Read covering list, saves list in covhistory
	// std::string name1 = "/home/van/simcovering.txt";
	//covhistory->read(name1); //TODO parameterübergabe in kommandezeile, nicht hardcoden!

	UpdateSticking(hitbuffer);
	UpdateDesorptionRate(hitbuffer);

	// Start Simulation = create first particle
	StartSimulation();


	// Run Simulation for simutime steps. One step ~ 1 seconds
	for(double i=0; i<(double)(simutime) && !eos;i+=realtimestep){
		if(i!=0||covhistory->pointintime_list.empty()){
			covhistory->appendList(hitbuffer,i); //append list with current time and covering //TODO offset if covering list does not end with timestep 0
			}

		if(i+timestep>=(double)(simutime)){ //last timestep
			std::tie(eos,realtimestep) = SimulationRun((double)simutime-i, hitbuffer); // Some additional simulation, as iteration step  does not run for exactly timestep ms
			covhistory->appendList(hitbuffer,i+realtimestep); // append list with last entry
			break;
			}
		std::tie(eos, realtimestep) = SimulationRun(timestep, hitbuffer);      // Run during timestep ms, performs MC steps
	}

	// Save simulation results in hitbuffer
	UpdateMCSubHits(hitbuffer, rank);

	// Update quantaties for contamination

	//UpdateSticking();
	estimateTmin();
	estimateTmin_RudiTest(hitbuffer);


	//Save history to new file
	std::string name0 = "/home/van/history"+std::to_string(rank)+".txt";
	covhistory->print();
	covhistory->write(name0);

	ResetTmpCounters(); //resets counter in sHandle
	return !eos;
}

double convertunit(double simutime, std::string unit){
	for(int i=0; i<6;i++){
		// day to seconds
		if(unit==day[i]) return (simutime*3600.0*24.0*1000.0);
	}
	for(int i=0; i<8;i++){
		// hour to seconds
			if(unit==hour[i]) return (simutime*3600.0*1000.0);
		}
	for(int i=0; i<8;i++){
		// minute to seconds
			if(unit==min[i]) return (simutime*60.0*1000.0);
		}
	//default: seconds
	return simutime*1000.0;

}

