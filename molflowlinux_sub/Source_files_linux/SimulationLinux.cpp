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
#include <fstream>
#include <sstream>

extern SimulationHistory* simHistory;


std::tuple<bool, std::vector<int> > simulateSub(Databuff *hitbuffer, int rank, int simutime){

	simHistory = new SimulationHistory (hitbuffer);
	double timestep=1000; // desired length per iteration for simulation, here hardcoded to 1 second
	double realtimestep; // actual time elapsed for iteration step

	// Set end of simulation flag
	bool eos=false;
	std::vector<int> facetNum;
	facetNum =std::vector<int> ();

	// Read covering list, saves list in covhistory
	// std::string name1 = "/home/van/simcovering.txt";
	//covhistory->read(name1); //TODO parameterübergabe in kommandezeile, nicht hardcoden!

	UpdateSticking(hitbuffer);
	UpdateDesorptionRate(hitbuffer);

	// Start Simulation = create first particle
	if(!StartSimulation())
		return {std::make_tuple(true,facetNum)};

	// Run Simulation for simutime steps. One step ~ 1 seconds
	for(double i=0; i<(double)(simutime) && !eos;i+=realtimestep){
		//if(i!=0||covhistory->pointintime_list.empty()){
		if(simHistory->coveringList.empty()){
			simHistory->appendList(hitbuffer,i); //append list with current time and covering //TODO offset if covering list does not end with timestep 0
			}

		if(i+timestep>=(double)(simutime)){ //last timestep
			std::tie(eos,realtimestep) = SimulationRun((double)simutime-i, hitbuffer, rank); // Some additional simulation, as iteration step  does not run for exactly timestep ms
			//covhistory->appendList(hitbuffer,i+realtimestep); // append list with last entry
			break;
			}
		std::tie(eos, realtimestep) = SimulationRun(timestep, hitbuffer, rank);      // Run during timestep ms, performs MC steps
	}

	// Save simulation results in hitbuffer
	UpdateMCSubHits(hitbuffer, rank);

	// Update quantaties for contamination

	//UpdateSticking();
	estimateTmin();
	estimateTmin_RudiTest(hitbuffer);


	//Save history to new file
	//std::string name0 = "/home/van/history"+std::to_string(rank)+".txt";
	if(rank==1)
		simHistory->print();
	//covhistory->write(name0);



	int num;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
					num=getFacetIndex(&f);
					if(f.tmpCounter[0].hit.covering==sHandle->coveringThreshold[num])
						{facetNum.push_back(num);}
				}
		}

	ResetTmpCounters(); //resets counter in sHandle
	delete[] simHistory;simHistory=NULL;
	return {std::make_tuple(eos,facetNum)};
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


ProblemDef::ProblemDef(int argc, char *argv[]){
	loadbufferPath= argc > 1 ? argv[1] :"";
	hitbufferPath=argc > 2 ? argv[2]:"";
	resultbufferPath=argc > 3 ? argv[3]:"";

	iterationNumber = 43200;
	s1=0.99;
	s2=0.2;
	E_de=1E-21;
	E_ad=1.5E-21;
	d=1;

	simulationTime = argc > 4? std::atof(argv[4]): 10.0;
	unit = argc > 5? argv[5]:"s";

	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	std::cout << "Simulation time " << simulationTime << unit << " converted to " << simulationTimeMS << "ms" << std::endl;

}

ProblemDef::ProblemDef(){
	loadbufferPath= "/home/van/Buffer/loadbuffer_alle_RT";
	hitbufferPath="/home/van/Buffer/hitbuffer_allee-6";
	resultbufferPath="/home/van/resultbuffer";

	iterationNumber = 43200;
	s1=0.99;
	s2=0.2;
	E_de=1E-21;
	E_ad=1.5E-21;
	d=1;

	simulationTime = 10.0;
	unit = "s";

	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);

}

void ProblemDef::readArg(int argc, char *argv[], int rank){
	loadbufferPath= argc > 1 ? argv[1] :loadbufferPath;
	hitbufferPath=argc > 2 ? argv[2]:hitbufferPath;
	resultbufferPath=argc > 3 ? argv[3]:resultbufferPath;

	simulationTime = argc > 4? std::atof(argv[4]): simulationTime;
	unit = argc > 5? argv[5]:unit;

	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	if (rank==0) {std::cout << "Simulation time " << simulationTime << unit << " converted to " << simulationTimeMS << "ms" << std::endl;}

}

void ProblemDef::readInputfile(std::string filename, int rank){
	std::string line;
	std::ifstream input(filename,std::ifstream::in);
	if (rank==0) {std::cout <<"Test to read input arguments from " <<filename <<std::endl;}

	while(std::getline(input,line)){
		if(rank==0) {std::cout <<line <<std::endl;}
		std::string stringIn;
		double doubleIn;
		int intIn;
		std::istringstream is( line );

		is >> stringIn;

		if(stringIn == "loadbufferPath") {if(rank==0) {std::cout <<line <<std::endl;}is >> stringIn; loadbufferPath=stringIn;}
		else if(stringIn == "hitbufferPath") {is >> stringIn; hitbufferPath=stringIn;}
		else if(stringIn=="resultbufferPath"){is >> stringIn; resultbufferPath=stringIn;}
		else if(stringIn == "simulationTime") {is >>doubleIn; simulationTime = doubleIn;}
		else if(stringIn == "unit"){is >> stringIn; unit=stringIn;}

		else if(stringIn =="iterationNumber"){is >> intIn; iterationNumber=intIn;}
		else if(stringIn =="s2"){is >> doubleIn; s2=doubleIn;}
		else if(stringIn =="s1"){is >> doubleIn; s1=doubleIn;}
		else if(stringIn =="d"){is >> doubleIn; d=doubleIn;}
		else if(stringIn =="E_de"){is >> doubleIn; E_de=doubleIn;}
		else if(stringIn =="E_ad"){is >> doubleIn; E_ad=doubleIn;}

		else{std::cout <<stringIn <<" not a valid argument." <<std::endl;}



	}
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	if (rank==0) {std::cout << "New Simulation time " << simulationTimeMS << "ms" << std::endl;}
}

void ProblemDef::writeInputfile(std::string filename, int rank){
	std::ofstream outfile(filename,std::ofstream::out|std::ios::trunc);
	if (rank==0) {std::cout <<"Test to write input arguments to " <<filename <<std::endl;}

	outfile <<"loadbufferPath" <<'\t' <<loadbufferPath <<std::endl;
	outfile <<"hitbufferPath" <<'\t' <<hitbufferPath <<std::endl;
	outfile <<"resultbufferPath" <<'\t' <<resultbufferPath <<std::endl;
	outfile <<"simulationTime" <<'\t' <<simulationTime <<std::endl;
	outfile <<"unit" <<'\t' <<unit <<std::endl;

	outfile <<"iterationNumber" <<'\t' <<iterationNumber<<std::endl;
	outfile <<"s1" <<'\t' <<s1<<std::endl;
	outfile <<"s2" <<'\t' <<s2<<std::endl;
	outfile <<"d" <<'\t' <<d<<std::endl;
	outfile <<"E_de" <<'\t' <<E_de<<std::endl;
	outfile <<"E_ad" <<'\t' <<E_ad<<std::endl;

}



