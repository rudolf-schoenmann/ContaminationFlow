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

bool simulateSub(Databuff *hitbuffer, int rank, int simutime){

	TimeTest test;
	double timestep=1000;
	double realtimestep;

	// Set end of simulation flag
	bool eos=false;

	// Start Simulation = create first particle
	StartSimulation();

	// Run Simulation for simutime steps. One step ~ 1 seconds
	for(double i=0; i<(double)(simutime) && !eos;i+=realtimestep){
		test.appendList(i);
		if(i+timestep>=(double)(simutime)){
			std::tie(eos,realtimestep) = SimulationRun((double)simutime-i); // Some additional simulation, as Simulation does not run for exactly timestep ms
			test.appendList(i+realtimestep);
			break;
		}
		std::tie(eos, realtimestep) = SimulationRun(timestep);      // Run during timestep ms, performs MC steps


		//calc new timestep?

	}

	// Save simulation results in hitbuffer
	UpdateSubHits(hitbuffer, rank);

	// Update quantaties for contamination
	UpdadeSticking();

	//std::cout <<"test " <<sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv  <<std::endl;
	std::cout <<"test of estimatetTmin() [unit = milliseconds]: " <<estimateTmin() <<std::endl;



	std::string name0 = "/home/van/history"+std::to_string(rank)+".txt";
	test.print();
	test.write(name0);

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

