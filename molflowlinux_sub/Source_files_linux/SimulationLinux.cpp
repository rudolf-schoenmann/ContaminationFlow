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
#include <unistd.h>
#include <libgen.h>
#include <time.h>
#include <sys/stat.h>

extern SimulationHistory* simHistory;
extern Simulation *sHandle;
extern ProblemDef* p;

// Simulation on subprocess
std::tuple<bool, std::vector<int> > simulateSub(Databuff *hitbuffer, int rank, int simutime){

	//Replaced consrtuctor with update function
	simHistory->updateHistory(hitbuffer);

	//timesteps
	double timestep=500; // desired length per iteration for simulation, here hardcoded to 1 second
	double realtimestep; // actual time elapsed for iteration step


	// Set end of simulation flag
	bool eos=false;

	// Facets that have reached the covering threshold
	std::vector<int> facetNum;
	facetNum =std::vector<int> ();

	//Update values of subprocess using hitbuffer
	UpdateSticking(hitbuffer);
	UpdateDesorptionRate(hitbuffer);

	// Start Simulation = create first particle
	if(!StartSimulation())
		return {std::make_tuple(true,facetNum)};

	// Run Simulation for timestep milliseconds
	for(double i=0; i<(double)(simutime) && !eos;i+=realtimestep){
		if(i>=(double(simutime)*0.99)){break;}
		if(simHistory->coveringList.empty()){
			simHistory->appendList(hitbuffer,i); //append list with initial covering
			}

		if(i+timestep>=(double)(simutime)){ //last timestep
			std::tie(eos,realtimestep) = SimulationRun((double)simutime-i, hitbuffer); // Some additional simulation, as iteration step  does not run for exactly timestep ms
			}
		else{
			std::tie(eos, realtimestep) = SimulationRun(timestep, hitbuffer);      // Run during timestep ms, performs MC steps
		}
		std::cout <<"  Elapsed calculation time for step (substep of one iteration step) for process " <<rank <<": "  <<realtimestep <<"ms" <<std::endl;
		p->outFile <<"  Elapsed calculation time for step (substep of one iteration step) for process " <<rank <<": "  <<realtimestep <<"ms" <<std::endl;
	}

	// Save simulation results in hitbuffer
	UpdateMCSubHits(hitbuffer, rank);


	// Update quantaties for contamination
	//UpdateSticking();
	//estimateTmin();
	//estimateTmin_RudiTest(hitbuffer);


	//Print simHistory example
	if(rank==1)
		simHistory->print();

	//std::cout <<simHistory->flightTime <<'\t' <<simHistory->nParticles <<std::endl;


	// Find Facets that have reached the covering threshold
	int num;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
					num=getFacetIndex(&f);
					if(f.tmpCounter[0].hit.covering==sHandle->coveringThreshold[num])
						{facetNum.push_back(num);}
				}
		}

	ResetTmpCounters(); //resets counter in sHandle
	return {std::make_tuple(eos,facetNum)};
}
//-----------------------------------------------------------
//helpful functions
double convertunit(double simutime, std::string unit){
	for(int i=0; i<8;i++){
		// year to seconds to MS
			if(unit==year[i]) return (simutime*365.25*24.0*3600.0*1000.0);
		}
	for(int i=0; i<8;i++){
		// month to seconds to MS
			if(unit==month[i]) return (simutime*30.4375*24.0*3600.0*1000.0);
		}
	for(int i=0; i<6;i++){
		// day to seconds to MS
		if(unit==day[i]) return (simutime*3600.0*24.0*1000.0);
	}
	for(int i=0; i<8;i++){
		// hour to seconds to MS
			if(unit==hour[i]) return (simutime*3600.0*1000.0);
		}
	for(int i=0; i<8;i++){
		// minute to seconds to MS
			if(unit==min[i]) return (simutime*60.0*1000.0);
		}
	//default: seconds to MS
	return simutime*1000.0;

}
//-----------------------------------------------------------
void printConsole(std::string str,std::ofstream outFile){
	std::cout <<str;
	outFile <<str;
}


//ProblemDef class

std::string get_path( )
{
        char arg1[20];
        char exepath[256] = {0};

        sprintf( arg1, "/proc/%d/exe", getpid() );
        readlink( arg1, exepath, sizeof(exepath) );
        return std::string( exepath );
}

ProblemDef::ProblemDef(int argc, char *argv[]){
	loadbufferPath= argc > 1 ? argv[1] :"";
	hitbufferPath=argc > 2 ? argv[2]:"";
	resultbufferPath=argc > 3 ? argv[3]:"";

	iterationNumber = 43200;
	//s1=0.99;
	//s2=0.2;
	E_de=1E-21;
	//E_ad=1E-21;
	d=1;
	H_vap=0.8E-19;
	W_tr=1.0;

	simulationTime = argc > 4? std::atof(argv[4]): 10.0;
	unit = argc > 5? argv[5]:"s";

	maxTime=10.0;
	maxUnit="y";

	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);

	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

}

ProblemDef::ProblemDef(){
	std::string path=get_path();
	std::cout <<path <<std::endl;
	char *test=&path[0u];
	std::string test2(dirname(dirname(test)));
	resultpath=test2+"/results/"+std::to_string(time(0));
	mkdir(resultpath.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	loadbufferPath= "/home/van/Buffer/loadbuffer_alle_RT";
	hitbufferPath="/home/van/Buffer/hitbuffer_allee-6";
	resultbufferPath=resultpath+"/resultbuffer";

	outFile.open(resultpath+"/console.txt", std::fstream::app);

	iterationNumber = 43200;
	//s1=1;
	//s2=0.2;
	E_de=1E-21;
	//E_ad=1E-21;
	d=1;
	H_vap=0.8E-19;
	W_tr=1.0;

	maxTime=10.0;
	maxUnit="y";
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

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

	writeInputfile(resultpath+"/InputFile.txt",rank);

}

void ProblemDef::readInputfile(std::string filename, int rank){
	std::string line;
	std::ifstream input(filename,std::ifstream::in);
	//if (rank==0) {std::cout <<"Test to read input arguments from " <<filename <<std::endl;}

	while(std::getline(input,line)){
		//if(rank==0) {std::cout <<line <<std::endl;}
		std::string stringIn;
		double doubleIn;
		int intIn;
		std::istringstream is( line );

		is >> stringIn;

		if(stringIn == "loadbufferPath") {if(rank==0) {std::cout <<line <<std::endl;}is >> stringIn; loadbufferPath=stringIn;}
		else if(stringIn == "hitbufferPath") {is >> stringIn; hitbufferPath=stringIn;}
		//else if(stringIn=="resultbufferPath"){is >> stringIn; resultbufferPath=stringIn;}
		else if(stringIn == "simulationTime") {is >>doubleIn; simulationTime = doubleIn;}
		else if(stringIn == "unit"){is >> stringIn; unit=stringIn;}

		else if(stringIn =="iterationNumber"){is >> intIn; iterationNumber=intIn;}
		else if(stringIn == "maxTime") {is >>doubleIn; maxTime = doubleIn;}
		else if(stringIn == "maxUnit"){is >> stringIn; maxUnit=stringIn;}

		//else if(stringIn =="s2"){is >> doubleIn; s2=doubleIn;}
		//else if(stringIn =="s1"){is >> doubleIn; s1=doubleIn;}
		else if(stringIn =="d"){is >> doubleIn; d=doubleIn;}
		else if(stringIn =="E_de"){is >> doubleIn; E_de=doubleIn;}
		//else if(stringIn =="E_ad"){is >> doubleIn; E_ad=doubleIn;}
		else if(stringIn =="H_vap"){is >> doubleIn; H_vap=doubleIn;}
		else if(stringIn =="W_tr"){is >> doubleIn; W_tr=doubleIn;}

		else{std::cout <<stringIn <<" not a valid argument." <<std::endl;}



	}
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

	writeInputfile(resultpath+"/InputFile.txt",rank);
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
	outfile <<"maxTime" <<'\t' <<maxTime <<std::endl;
	outfile <<"maxUnit" <<'\t' <<maxUnit <<std::endl;

	//outfile <<"s1" <<'\t' <<s1<<std::endl;
	//outfile <<"s2" <<'\t' <<s2<<std::endl;
	outfile <<"d" <<'\t' <<d<<std::endl;
	outfile <<"E_de" <<'\t' <<E_de<<std::endl;
	//outfile <<"E_ad" <<'\t' <<E_ad<<std::endl;

	outfile <<"H_vap" <<'\t' <<H_vap <<std::endl;
	outfile <<"W_tr" <<'\t' <<W_tr <<std::endl;

}

void ProblemDef::printInputfile(std::ostream& out){ //std::cout or p->outFile

	out  <<std::endl<<"Print input arguments"<<std::endl;

	out  <<"resultPath" <<'\t' <<resultpath <<std::endl;
	out  <<"loadbufferPath" <<'\t' <<loadbufferPath <<std::endl;
	out  <<"hitbufferPath" <<'\t' <<hitbufferPath <<std::endl;
	out  <<"resultbufferPath" <<'\t' <<resultbufferPath <<std::endl<<std::endl;

	out  <<"simulationTime" <<'\t' <<simulationTime <<std::endl;
	out  <<"unit" <<'\t' <<unit <<std::endl;

	out  <<"iterationNumber" <<'\t' <<iterationNumber<<std::endl;
	out  <<"maxTime" <<'\t' <<maxTime <<std::endl;
	out  <<"maxUnit" <<'\t' <<maxUnit <<std::endl<<std::endl;

	//out  <<"s1" <<'\t' <<s1<<std::endl;
	//out  <<"s2" <<'\t' <<s2<<std::endl;
	out  <<"d" <<'\t' <<d<<std::endl;
	out  <<"E_de" <<'\t' <<E_de<<std::endl;
	//out  <<"E_ad" <<'\t' <<E_ad<<std::endl<<std::endl;
	out <<"H_vap" <<'\t' <<H_vap <<std::endl;
	out <<"W_tr" <<'\t' <<W_tr <<std::endl;

	out  << "Simulation time " << simulationTime << unit << " converted to " << simulationTimeMS << "ms" << std::endl;
	out  << "Maximum simulation time " << maxTime << maxUnit << " converted to " << maxTimeS << "s" << std::endl<<std::endl;

}

//-----------------------------------------------------------
//SimulationHistory

SimulationHistory::SimulationHistory(){
	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	nbDesorbed_old=0;
	lastTime=0.0;
	currentStep=0;
}

SimulationHistory::SimulationHistory(Databuff *hitbuffer){
	std::vector<llong> currentstep;
	currentstep =std::vector<llong> ();

	std::vector<llong> currentHits;
	currentHits =std::vector<llong> ();

	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	nbDesorbed_old= getnbDesorbed(hitbuffer);

	llong numHit;
	llong covering;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				covering = getCovering(&f, hitbuffer);
				numHit=getHits(&f, hitbuffer);
				currentstep.push_back(covering);
				currentHits.push_back(numHit);
				f.tmpCounter[0].hit.covering=covering;
				numFacet+=1;
			}
	}

	coveringList.appendList(currentstep, 0);
	coveringList.initCurrent(numFacet);
	hitList.appendList(currentHits, 0);
	hitList.initCurrent(numFacet);
	errorList.initCurrent(numFacet);
	errorList.appendCurrent(0.0);
	currentStep=0;
}

void SimulationHistory::updateHistory(Databuff *hitbuffer){
	std::vector<llong> currentstep;
	currentstep =std::vector<llong> ();
	std::vector<llong> currentHits;
	currentHits =std::vector<llong> ();
	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	nbDesorbed_old= getnbDesorbed(hitbuffer);

	llong numHit;
	llong covering;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				covering = getCovering(&f, hitbuffer);
				numHit=getHits(&f, hitbuffer);
				currentstep.push_back(covering);
				currentHits.push_back(numHit);
				f.tmpCounter[0].hit.covering=covering;
				numFacet+=1;
			}
	}
	coveringList.reset();
	coveringList.appendList(currentstep, 0);
	coveringList.initCurrent(numFacet);
	hitList.reset();
	hitList.appendList(currentHits, 0);
	hitList.initCurrent(numFacet);


}


void SimulationHistory::appendList(Databuff *hitbuffer, double time){

	std::vector<llong> currentstep;
	currentstep =std::vector<llong> ();

	if(time==-1.0) //One step
		time=coveringList.pointintime_list.back().first+1.0;

	llong covering;

	//int i=0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering=getCovering(&f, hitbuffer);
			currentstep.push_back(covering);
			//i+=1;
		}
	}
	coveringList.appendList(currentstep, time);

}

void SimulationHistory::print(bool write){
	coveringList.print(std::cout, "Accumulative covering");
	hitList.print(std::cout, "Accumulative number hits");
	errorList.print(std::cout, "Error per iteration");
	if(write){
		coveringList.print(p->outFile, "Accumulative covering");
		hitList.print(p->outFile, "Accumulative number hits");
		errorList.print(p->outFile, "Error per iteration");
	}
}

void SimulationHistory::write(std::string path){
	coveringList.write(path+"/covering.txt");
}


