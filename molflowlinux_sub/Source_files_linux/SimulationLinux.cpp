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

	// Set end of simulation flag
	bool eos=false;

	// Start Simulation = create first particle
	StartSimulation();

	// Run Simulation for simutime steps. One step ~ 1 seconds
	for(int i=0; i<simutime && !eos;i++){
		eos = SimulationRun();      // Run during 1 sec, performs MC steps
	}

	// Save simulation results in hitbuffer
	UpdateSubHits(hitbuffer, rank);
	UpdadeSticking();

	//std::cout <<"test " <<sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv  <<std::endl;
	std::cout <<"test " <<estimateTmin() <<std::endl;

	ResetTmpCounters(); //resets counter in sHandle
	return !eos;
}

double convertunit(double simutime, std::string unit){
	for(int i=0; i<6;i++){
		// day to seconds
		if(unit==day[i]) return (simutime*3600.0*24.0);
	}
	for(int i=0; i<8;i++){
		// hour to seconds
			if(unit==hour[i]) return (simutime*3600.0);
		}
	for(int i=0; i<8;i++){
		// minute to seconds
			if(unit==min[i]) return (simutime*60.0);
		}
	//default: seconds
	return simutime;

}

