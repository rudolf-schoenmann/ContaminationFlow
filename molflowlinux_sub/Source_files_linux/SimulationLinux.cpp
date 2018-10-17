/*
 * Buffer.cpp
 *
 *  Created on: 17.09.2018
 *      Author: Van
 */

#include "SimulationLinux.h"
#include <array>

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

