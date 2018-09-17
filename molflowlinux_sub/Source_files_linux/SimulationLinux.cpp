/*
 * Buffer.cpp
 *
 *  Created on: 17.09.2018
 *      Author: Van
 */

#include "SimulationLinux.h"

bool simulateSub(double simutime, Databuff *hitbuffer, int rank){

	bool eos=false;

	InitSimulation();
	SetReady();
	StartSimulation();
	eos = SimulationRun(simutime);      // Run during simutime sec
	UpdateHits(hitbuffer, rank);

	return !eos;
}
