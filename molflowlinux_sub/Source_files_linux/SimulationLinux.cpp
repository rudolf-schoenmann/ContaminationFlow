/*
 * Buffer.cpp
 *
 *  Created on: 17.09.2018
 *      Author: Van
 */

#include "SimulationLinux.h"
#include <array>

extern Simulation *sHandle; //delcared in molflowSub.cpp

bool simulateSub(Databuff *hitbuffer, int rank, double simutime, std::string unit){

	unsigned int newsimutime = (int)(convertunit(simutime, unit)+0.5);
	std::cout << "Simulation time " << simutime << unit <<" converted to " <<newsimutime <<"s" <<std::endl;

	bool eos=false;

	InitSimulation();//creates simulation handle
	SetReady(); //checks simulation status?
	StartSimulation();
	for(unsigned int i=0; i<newsimutime && !eos;i++){
		eos = SimulationRun();      // Run during 1 sec, performs MC steps
	}
	UpdateSubHits(hitbuffer, rank);

	return true;
}

double convertunit(double simutime, std::string unit){
	for(int i=0; i<6;i++){
		if(unit==day[i]) return (simutime*3600.0*24.0);
	}
	for(int i=0; i<8;i++){
			if(unit==hour[i]) return (simutime*3600.0);
		}
	for(int i=0; i<8;i++){
			if(unit==min[i]) return (simutime*60.0);
		}
	return simutime;

}

