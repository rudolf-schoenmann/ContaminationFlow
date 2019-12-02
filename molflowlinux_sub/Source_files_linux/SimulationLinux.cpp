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
#include <limits>

extern SimulationHistory* simHistory;
extern Simulation *sHandle;
extern ProblemDef* p;

// Simulation on subprocess
std::tuple<bool, std::vector<int> > simulateSub(Databuff *hitbuffer, int rank, int simutime){

	//Calculate target values for error and desorbed particles
	int targetParticles=p->targetParticles/simHistory->numSubProcess;
	double targetError=p->targetError*pow(simHistory->numSubProcess,0.5);

	//Replaced constructor with update function
	bool smallCovering; llong smallCoveringFactor;
		std::tie(smallCovering,smallCoveringFactor) = simHistory->updateHistory(hitbuffer);
	if(rank==1){
		std::cout <<std::endl <<"Currentstep: " << simHistory->currentStep <<". Step size: " <<simHistory->stepSize <<std::endl;
		std::cout <<"Target Particles: " << targetParticles <<". Target Error: " <<targetError <<std::endl <<std::endl;

		p->outFile <<std::endl <<"Currentstep: " << simHistory->currentStep <<". Step size: " <<simHistory->stepSize <<std::endl;
		p->outFile <<"Target Particles: " << targetParticles <<". Target Error: " <<targetError <<std::endl <<std::endl;
	}

	//Values for simulation
	double timestep=1000; // desired length per iteration for simulation, here hardcoded to 1 second
	double realtimestep; // actual time elapsed for iteration step
	double i;			// Time elapsed between checking of targets
	double totalTime=0.0;	//Total simulated time
	bool eos=false;		// Set end of simulation flag
	double totalError=1.;//Total error

	// Facets that have reached the covering threshold
	std::vector<int> facetNum;
	facetNum =std::vector<int> ();

	//----Simulation

	// Start Simulation = create first particle
	if(!sHandle->currentParticle.lastHitFacet){
			if(!StartSimulation())
				return {std::make_tuple(true,facetNum)};
		}


	// Run Simulation for timestep milliseconds
	for(int j=0; j<p->maxSimPerIt && !(j>0 && simHistory->nParticles>targetParticles && totalError<targetError)&& !eos; j++){

		for(i=0; i<(double)(simutime) && !eos;i+=realtimestep){

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

		}
		totalTime+=i;

		totalError=UpdateError();

		if(j%(int)(30000/simutime)==0 || (simHistory->nParticles>targetParticles && totalError<targetError)|| eos){
			// Print information every 30s or if target reached
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"  "<<rank<<": Step "<<std::setw(4)<<std::right <<j <<"    &    Total time " <<std::setw(10)<<std::right <<totalTime <<"ms    &    Desorbed particles "<<std::setw(10)<<std::right<<simHistory->nParticles <<"    &    Total error "  <<std::setw(10)<<std::left<<totalError<<std::endl;

			if(totalError>targetError){
				simHistory->hitList.printCurrent(tmpstream, std::to_string(rank)+": hitlist");
				simHistory->desorbedList.printCurrent(tmpstream, std::to_string(rank)+": desorbedlist");
				simHistory->errorList.printCurrent(tmpstream, std::to_string(rank)+": errorlist");
				tmpstream <<std::endl;
			}

			std::cout <<tmpstream.str();
			p->outFile <<tmpstream.str();
		}
	}

	// Find Facets that have reached the covering threshold
		int num;
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				if(f.sh.desorption==0) continue;

				num=getFacetIndex(&f);
				if(smallCovering){
					//std::cout <<"\t  "<<rank << ": Facet " <<num <<": Covering " <<f.tmpCounter[0].hit.covering;
					f.tmpCounter[0].hit.covering/=smallCoveringFactor;
					sHandle->coveringThreshold[num]/=smallCoveringFactor;
					//std::cout <<"  to  " <<f.tmpCounter[0].hit.covering <<std::endl;
				}

				if(f.tmpCounter[0].hit.covering<=sHandle->coveringThreshold[num])
					{facetNum.push_back(num);}
			}
		}

	// Save simulation results in hitbuffer
	UpdateMCSubHits(hitbuffer, rank);

	if(!facetNum.empty())
			eos=true;

	ResetTmpCounters(); //resets counter in sHandle
	return {std::make_tuple(eos,facetNum)};
}

std::tuple<bool, std::vector<int> > simulateSub2(Databuff *hitbuffer,int rank, int simutime){

	//Calculate target values for error and desorbed particles
	int targetParticles=p->targetParticles/simHistory->numSubProcess;
	double targetError=p->targetError*pow(simHistory->numSubProcess,0.5);

	//Replaced constructor with update function
	bool smallCovering; llong smallCoveringFactor;
	std::tie(smallCovering,smallCoveringFactor) = simHistory->updateHistory();
	if(rank==1){
		std::cout <<std::endl <<"Currentstep: " << simHistory->currentStep <<". Step size: " <<simHistory->stepSize <<std::endl;
		std::cout <<"Target Particles: " << targetParticles <<". Target Error: " <<targetError <<std::endl <<std::endl;

		p->outFile <<std::endl <<"Currentstep: " << simHistory->currentStep <<". Step size: " <<simHistory->stepSize <<std::endl;
		p->outFile <<"Target Particles: " << targetParticles <<". Target Error: " <<targetError <<std::endl <<std::endl;
	}

	//Values for simulation
	double timestep=1000; // desired length per iteration for simulation, here hardcoded to 1 second
	double realtimestep; // actual time elapsed for iteration step
	double i;			// Time elapsed between checking of targets
	double totalTime=0.0;	//Total simulated time
	bool eos=false;		// Set end of simulation flag
	double totalError=1.;//Total error

	// Facets that have reached the covering threshold
	std::vector<int> facetNum;
	facetNum =std::vector<int> ();

	//----Simulation

	// Start Simulation = create first particle
	if(!sHandle->currentParticle.lastHitFacet){
		if(!StartSimulation())
			return {std::make_tuple(true,facetNum)};
	}


	// Run Simulation for timestep milliseconds
	for(int j=0; j<p->maxSimPerIt && !(j>0 && simHistory->nParticles>targetParticles && totalError<targetError)&& !eos; j++){

		for(i=0; i<(double)(simutime) && !eos;i+=realtimestep){

			if(i>=(double(simutime)*0.99)){break;}
			if(simHistory->coveringList.empty()){
				simHistory->appendList(i); //append list with initial covering
				}

			if(i+timestep>=(double)(simutime)){ //last timestep
				std::tie(eos,realtimestep) = SimulationRun((double)simutime-i); // Some additional simulation, as iteration step  does not run for exactly timestep ms
				}
			else{
				std::tie(eos, realtimestep) = SimulationRun(timestep);      // Run during timestep ms, performs MC steps
			}

		}
		totalTime+=i;

		totalError=UpdateError();

		if(j%(int)(30000/simutime)==0 || (simHistory->nParticles>targetParticles && totalError<targetError)|| eos || j >= p->maxSimPerIt-1){
			// Print information every 30s or if target reached
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"  "<<rank<<": Step "<<std::setw(4)<<std::right <<j <<"    &    Total time " <<std::setw(10)<<std::right <<totalTime <<"ms    &    Desorbed particles "<<std::setw(10)<<std::right<<simHistory->nParticles <<"    &    Total error "  <<std::setw(10)<<std::left<<totalError<<std::endl;

			if(totalError>targetError){
				simHistory->hitList.printCurrent(tmpstream, std::to_string(rank)+": hitlist");
				simHistory->desorbedList.printCurrent(tmpstream, std::to_string(rank)+": desorbedlist");
				simHistory->coveringList.printCurrent(tmpstream, std::to_string(rank)+": coveringlist");
				simHistory->errorList.printCurrent(tmpstream, std::to_string(rank)+": errorlist");
				tmpstream <<std::endl;
			}

			std::cout <<tmpstream.str();
			p->outFile <<tmpstream.str();
		}
	}

	// Find Facets that have reached the covering threshold
	int num;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(f.sh.desorption==0) continue;

			num=getFacetIndex(&f);
			if(smallCovering){
				//std::cout <<"\t  "<<rank << ": Facet " <<num <<": Covering " <<f.tmpCounter[0].hit.covering;
				f.tmpCounter[0].hit.covering/=smallCoveringFactor;
				sHandle->coveringThreshold[num]/=smallCoveringFactor;
				//std::cout <<"  to  " <<f.tmpCounter[0].hit.covering <<std::endl;
			}

			if(f.tmpCounter[0].hit.covering<=sHandle->coveringThreshold[num])
				{facetNum.push_back(num);}
		}
	}

	// Save simulation results in hitbuffer
	UpdateMCSubHits(hitbuffer, rank);

	if(!facetNum.empty())
		eos=true;

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

//----get path of executable
std::string get_path( )
{
        char arg1[20];
        char exepath[256] = {0};

        sprintf( arg1, "/proc/%d/exe", getpid() );
        readlink( arg1, exepath, sizeof(exepath) );
        return std::string( exepath );
}

//----ProblemDef class
ProblemDef::ProblemDef(){
	loadbufferPath= "/home/van/Buffer/loadbuffer_alle_RT";
	hitbufferPath="/home/van/Buffer/hitbuffer_allee-6";

	iterationNumber = 43200;
	//s1=1;
	//s2=0.2;
	E_de=1E-21;
	//E_ad=1E-21;
	d=1;
	H_vap=0.8E-19;
	W_tr=1.0;
	sticking=0.0;
	coveringLimit=0;

	maxTime=10.0;
	maxUnit="y";
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

	simulationTime = 10.0;
	unit = "s";
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);

	saveResults=true;

	targetParticles=1000;
	targetError=0.001;
	hitRatioLimit=1E-5;
	Tmin=1E-4;


	coveringMinThresh=10000;
	//coveringMinFactor=100.0;
	//coveringMaxFactor=1000.0;

	maxStepSize=std::numeric_limits<double>::infinity();
	maxSimPerIt=std::numeric_limits<int>::infinity();

}

void ProblemDef::createOutput(int save){
	if(save!=0){
		std::string path=get_path();
		//std::cout <<path <<std::endl;
		char *test=&path[0u];
		std::string test2(dirname(dirname(test)));
		resultpath=test2+"/results/"+std::to_string(time(0));
		mkdir(resultpath.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

		resultbufferPath=resultpath+"/resultbuffer";
		outFile.open(resultpath+"/console.txt", std::fstream::app);
		saveResults=true;
	}
	else{
		outFile.setstate(std::ios_base::badbit);
		saveResults=false;
	}
}

void ProblemDef::readArg(int argc, char *argv[], int rank){
	int tmpsave=argc > 3 ? (int)atof(argv[3]):1;
	createOutput(tmpsave);

	loadbufferPath= argc > 1 ? argv[1] :loadbufferPath;
	hitbufferPath=argc > 2 ? argv[2]:hitbufferPath;


	simulationTime = argc > 4? std::atof(argv[4]): simulationTime;
	unit = argc > 5? argv[5]:unit;

	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);

	if(saveResults)
		writeInputfile(resultpath+"/InputFile.txt",rank);

}

void ProblemDef::readInputfile(std::string filename, int rank, int save){
	createOutput(save);

	std::string line;
	std::ifstream input(filename,std::ifstream::in);
	//if (rank==0) {std::cout <<"Test to read input arguments from " <<filename <<std::endl;}

	while(std::getline(input,line)){
		//if(rank==0) {std::cout <<line <<std::endl;}
		std::string stringIn;
		double doubleIn;
		int intIn;
		llong llongIn;
		std::istringstream is( line );

		is >> stringIn;

		if(stringIn == "loadbufferPath") {if(rank==0) {std::cout <<line <<std::endl;}is >> stringIn; loadbufferPath=stringIn;}
		else if(stringIn == "hitbufferPath") {is >> stringIn; hitbufferPath=stringIn;}
		//else if(stringIn=="resultbufferPath"){is >> stringIn; resultbufferPath=stringIn;} //deprecated
		else if(stringIn == "simulationTime") {is >>doubleIn; simulationTime = doubleIn;}
		else if(stringIn == "unit"){is >> stringIn; unit=stringIn;}

		else if(stringIn =="iterationNumber"){is >> intIn; iterationNumber=intIn;}
		else if(stringIn == "maxTime") {is >>doubleIn; maxTime = doubleIn;}
		else if(stringIn == "maxUnit"){is >> stringIn; maxUnit=stringIn;}

		//else if(stringIn =="s2"){is >> doubleIn; s2=doubleIn;} //deprecated
		//else if(stringIn =="s1"){is >> doubleIn; s1=doubleIn;} //deprecated
		else if(stringIn =="d"){is >> doubleIn; d=doubleIn;}
		else if(stringIn =="E_de"){is >> doubleIn; E_de=doubleIn;}
		//else if(stringIn =="E_ad"){is >> doubleIn; E_ad=doubleIn;} //deprecated
		else if(stringIn =="H_vap"){is >> doubleIn; H_vap=doubleIn;}
		else if(stringIn =="W_tr"){is >> doubleIn; W_tr=doubleIn;}
		else if(stringIn =="sticking"){is >> doubleIn; sticking=doubleIn;}

		else if(stringIn =="targetParticles"){is >> intIn; targetParticles=intIn;}
		else if(stringIn == "targetError") {is >>doubleIn; targetError = doubleIn;}
		else if(stringIn == "hitRatioLimit") {is >>doubleIn; hitRatioLimit = doubleIn;}
		else if(stringIn == "Tmin") {is >>doubleIn; Tmin = doubleIn;}

		else if(stringIn =="coveringLimit"){is >> llongIn; coveringLimit=llongIn;}
		else if(stringIn == "coveringMinThresh") {is >>llongIn; coveringMinThresh = llongIn;}
		//else if(stringIn == "coveringMinFactor") {is >>doubleIn; coveringMinFactor = doubleIn;}
		//else if(stringIn == "coveringMaxFactor") {is >>doubleIn; coveringMaxFactor = doubleIn;}

		else if(stringIn =="maxStepSize"){is >>doubleIn; maxStepSize=doubleIn;}
		else if(stringIn =="maxSimPerIt"){is >> intIn; maxSimPerIt=intIn;}

		else{std::cout <<stringIn <<" not a valid argument." <<std::endl;}



	}
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

	if(saveResults)
		writeInputfile(resultpath+"/InputFile.txt",rank);
}

void ProblemDef::writeInputfile(std::string filename, int rank){

	if(rank==0){
	std::ofstream outfile(filename,std::ofstream::out|std::ios::trunc);
	std::cout <<"Write input arguments to " <<filename <<std::endl;

	outfile <<"loadbufferPath" <<'\t' <<loadbufferPath <<std::endl;
	outfile <<"hitbufferPath" <<'\t' <<hitbufferPath <<std::endl;
	//outfile <<"resultbufferPath" <<'\t' <<resultbufferPath <<std::endl; //deprecated
	outfile <<"simulationTime" <<'\t' <<simulationTime <<std::endl;
	outfile <<"unit" <<'\t' <<unit <<std::endl;

	outfile <<"iterationNumber" <<'\t' <<iterationNumber<<std::endl;
	outfile <<"maxTime" <<'\t' <<maxTime <<std::endl;
	outfile <<"maxUnit" <<'\t' <<maxUnit <<std::endl;

	outfile <<"sticking" <<'\t' <<sticking<<std::endl;
	//outfile <<"s1" <<'\t' <<s1<<std::endl; //deprecated
	//outfile <<"s2" <<'\t' <<s2<<std::endl; //deprecated
	outfile <<"d" <<'\t' <<d<<std::endl;
	outfile <<"E_de" <<'\t' <<E_de<<std::endl;
	//outfile <<"E_ad" <<'\t' <<E_ad<<std::endl; //deprecated

	outfile <<"H_vap" <<'\t' <<H_vap <<std::endl;
	outfile <<"W_tr" <<'\t' <<W_tr <<std::endl;

	outfile <<"targetError" <<'\t' <<targetError <<std::endl;
	outfile <<"targetParticles" <<'\t' <<targetParticles <<std::endl;
	outfile <<"hitRatioLimit" <<'\t' <<hitRatioLimit <<std::endl;
	outfile <<"Tmin" <<Tmin <<'\t' <<std::endl;

	outfile <<"coveringLimit" <<"\t" <<coveringLimit<<std::endl;
	outfile <<"coveringMinThresh" <<"\t" <<coveringMinThresh<<std::endl;

	outfile <<"maxStepSize" <<"\t" <<maxStepSize<<std::endl;
	outfile <<"maxSimPerIt" <<"\t" <<maxSimPerIt<<std::endl;
	//outfile <<"coveringMinFactor" <<"\t" <<coveringMinFactor;
	//outfile <<"coveringMaxFactor" <<"\t" <<coveringMaxFactor;

	}

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

	//out  <<"s1" <<'\t' <<s1<<std::endl; //deprecated
	//out  <<"s2" <<'\t' <<s2<<std::endl; //deprecated
	out <<"sticking" <<'\t' <<sticking <<std::endl;
	out  <<"d" <<'\t' <<d<<std::endl;
	out  <<"E_de" <<'\t' <<E_de<<std::endl;
	//out  <<"E_ad" <<'\t' <<E_ad<<std::endl<<std::endl; //deprecated
	out <<"H_vap" <<'\t' <<H_vap <<std::endl;
	out <<"W_tr" <<'\t' <<W_tr <<std::endl;

	out <<"targetError" <<'\t' <<targetError <<std::endl;
	out <<"targetParticles" <<'\t' <<targetParticles <<std::endl;
	out <<"hitRatioLimit" <<'\t' <<hitRatioLimit <<std::endl;
	out <<"Tmin" <<'\t' <<Tmin <<std::endl;

	out <<"coveringLimit" <<"\t" <<coveringLimit <<std::endl;
	out <<"coveringMinThresh" <<"\t" <<coveringMinThresh <<std::endl;
	//out <<"coveringMinFactor" <<"\t" <<coveringMinFactor;
	//out <<"coveringMaxFactor" <<"\t" <<coveringMaxFactor;
	out <<"maxStepSize" <<"\t" <<maxStepSize<<std::endl;
	out <<"maxSimPerIt" <<"\t" <<maxSimPerIt<<std::endl;

	out  << "Simulation time " << simulationTime << unit << " converted to " << simulationTimeMS << "ms" << std::endl;
	out  << "Maximum simulation time " << maxTime << maxUnit << " converted to " << maxTimeS << "s" << std::endl<<std::endl;

}

//-----------------------------------------------------------
//----SimulationHistory class

SimulationHistory::SimulationHistory(int world_size){
	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	nbDesorbed_old=0;
	lastTime=0.0;
	currentStep=0;
	stepSize=0.0;
	numSubProcess=world_size-1;
	startNewParticle=false;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			numFacet+=1;
		}
	}
	coveringList.initCurrent(numFacet);
	hitList.initCurrent(numFacet);
	desorbedList.initCurrent(numFacet);
	errorList.initCurrent(numFacet);


}

SimulationHistory::SimulationHistory(Databuff *hitbuffer, int world_size){
	/*
	std::vector<boost::multiprecision::uint128_t> currentCov;
	currentCov =std::vector<boost::multiprecision::uint128_t> ();

	std::vector<llong> currentDes;
	currentDes =std::vector<llong> ();

	std::vector<double> currentHits;
	currentHits =std::vector<double> ();
	*/

	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	nbDesorbed_old= getnbDesorbed(hitbuffer);

	double numHit;
	llong numDes;
	boost::multiprecision::uint128_t covering;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				covering = boost::multiprecision::uint128_t(getCovering(&f, hitbuffer));
				numHit=getHits(&f, hitbuffer);
				numDes=getnbDesorbed(&f, hitbuffer);

				coveringList.currentList.push_back(covering);
				hitList.currentList.push_back(numHit);
				desorbedList.currentList.push_back(numDes);

				f.tmpCounter[0].hit.covering=(llong)covering;
				numFacet+=1;
			}
	}

	coveringList.appendCurrent(0);
	//coveringList.initCurrent(numFacet);
	hitList.appendCurrent(0);
	//hitList.initCurrent(numFacet);
	desorbedList.appendCurrent(0);
	//desorbedList.initCurrent(numFacet);
	errorList.initCurrent(numFacet);
	errorList.appendCurrent(0.0);

	currentStep=0;
	stepSize=0.0;

	numSubProcess=world_size-1;
	startNewParticle=false;

}

std::tuple<bool, llong > SimulationHistory::updateHistory(Databuff *hitbuffer){
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	startNewParticle=false;

	nbDesorbed_old= getnbDesorbed(hitbuffer);

	double numHit;
	llong numDes;
	llong smallCoveringFactor=0;
	boost::multiprecision::uint128_t covering;

	long double avg=0;
	long double avgFacet=0.0;

	bool smallCovering=false;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				covering = boost::multiprecision::uint128_t(getCovering(&f, hitbuffer));
				numHit=getHits(&f, hitbuffer);
				numDes=getnbDesorbed(&f, hitbuffer);

				coveringList.setCurrentList(&f, covering);
				hitList.setCurrentList(&f, numHit);
				desorbedList.setCurrentList(&f, numDes);

				f.tmpCounter[0].hit.covering=covering.convert_to<llong>();

				if(llong(covering)<p->coveringMinThresh && f.sh.desorption>0.0){
				avg+=covering.convert_to<long double>();
				avgFacet+=1.0;
				smallCovering=true;
				//std::cout <<"Facet "<<numFacet-1 <<": "<<llong(covering) <<" smaller than " <<p->coveringMinThresh <<std::endl;
				}
			}
	}
	coveringList.reset();
	coveringList.appendCurrent(0);


	hitList.reset();
	hitList.appendCurrent(0);

	desorbedList.reset();
	desorbedList.appendCurrent(0);

	errorList.reset();
	errorList.initCurrent(numFacet);

	stepSize=manageStepSize(false);

	if(avgFacet>0.0){
		if(smallCovering && llong(avg/avgFacet)<=p->coveringMinThresh){
			//long double minCov(p->coveringMinFactor); long double threshCov(p->coveringMinThresh); long double maxCov(p->coveringMaxFactor);
			//smallCoveringFactor=p->coveringMaxFactor - llong(((maxCov-minCov)/threshCov)*(avg/avgFacet));
			//smallCoveringFactor=llong(calcStep((long double)(avg/avgFacet), 2*p->coveringMaxFactor-p->coveringMinFactor, p->coveringMinFactor, 0, 0.1*double(p->coveringMinThresh)));

			smallCoveringFactor=llong(1.0+1.1*double(p->coveringMinThresh)/(double(avg/avgFacet)));

			std::cout <<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for avg "<< avg/avgFacet<<std::endl;
			p->outFile<<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for avg "<< avg/avgFacet<<std::endl;

			for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
				for (SubprocessFacet& f : sHandle->structures[s].facets) {
					if(f.sh.desorption==0) continue;
					f.tmpCounter[0].hit.covering*=smallCoveringFactor;
					sHandle->coveringThreshold[getFacetIndex(&f)]*=smallCoveringFactor;
				}
			}
		}
	}

	return {std::make_tuple(smallCovering,smallCoveringFactor)};
}

std::tuple<bool, llong > SimulationHistory::updateHistory(){

	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	startNewParticle=false;

	//nbDesorbed_old= getnbDesorbed(hitbuffer);

	llong smallCoveringFactor=0;
	boost::multiprecision::uint128_t covering;

	long double avg=0;
	long double avgFacet=0.0;

	bool smallCovering=false;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering = getCovering(&f);

			f.tmpCounter[0].hit.covering=covering.convert_to<llong>();

			if(llong(covering)<p->coveringMinThresh && f.sh.desorption>0.0){
				avg+=covering.convert_to<long double>();
				avgFacet+=1.0;
				smallCovering=true;
				//std::cout <<"Facet "<<numFacet-1 <<": "<<llong(covering) <<" smaller than " <<p->coveringMinThresh <<std::endl;
			}
		}
	}
	coveringList.pointintime_list.clear();
	coveringList.appendCurrent(0);

	hitList.reset();
	hitList.initCurrent(numFacet);
	hitList.appendCurrent(0);

	desorbedList.reset();
	desorbedList.initCurrent(numFacet);
	desorbedList.appendCurrent(0);

	errorList.reset();
	errorList.initCurrent(numFacet);

	stepSize=manageStepSize(false);

	if(avgFacet>0.0){
			if(smallCovering && llong(avg/avgFacet)<=p->coveringMinThresh){
				//long double minCov(p->coveringMinFactor); long double threshCov(p->coveringMinThresh); long double maxCov(p->coveringMaxFactor);
				//smallCoveringFactor=p->coveringMaxFactor - llong(((maxCov-minCov)/threshCov)*(avg/avgFacet));
				//smallCoveringFactor=llong(calcStep((long double)(avg/avgFacet), 2*p->coveringMaxFactor-p->coveringMinFactor, p->coveringMinFactor, 0, 0.1*double(p->coveringMinThresh)));

				smallCoveringFactor=llong(1.0+1.1*double(p->coveringMinThresh)/(double(avg/avgFacet)));

				std::cout <<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for avg "<< avg/avgFacet<<std::endl;
				p->outFile<<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for avg "<< avg/avgFacet<<std::endl;

				for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
					for (SubprocessFacet& f : sHandle->structures[s].facets) {
						if(f.sh.desorption==0) continue;
						f.tmpCounter[0].hit.covering*=smallCoveringFactor;
						sHandle->coveringThreshold[getFacetIndex(&f)]*=smallCoveringFactor;
					}
				}
			}
		}

	return {std::make_tuple(smallCovering,smallCoveringFactor)};
}


void SimulationHistory::appendList(Databuff *hitbuffer, double time){

	std::vector<boost::multiprecision::uint128_t> currentCov;
	currentCov =std::vector<boost::multiprecision::uint128_t> ();

	if(time==-1.0) //One step
		time=coveringList.pointintime_list.back().first+1.0;

	boost::multiprecision::uint128_t covering;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering=boost::multiprecision::uint128_t(getCovering(&f, hitbuffer));
			currentCov.push_back(covering);
		}
	}
	coveringList.appendList(currentCov, time);
}

void SimulationHistory::appendList(double time){

	if(time==-1.0) //One step
		time=coveringList.pointintime_list.back().first+1.0;

	coveringList.appendCurrent(time);

}

void SimulationHistory::print(bool write){
	std::vector<double> errorPerIt;
	std::vector<boost::multiprecision::uint128_t> covPerIt;
	std::tie(errorPerIt,covPerIt) = CalcPerIteration();

	coveringList.print(std::cout,covPerIt, "Accumulative covering");
	hitList.print(std::cout, "Accumulative number hits");
	desorbedList.print(std::cout, "Accumulative number desorbed");
	errorList.print(std::cout,errorPerIt, "Error per iteration");

	if(write){
		coveringList.print(p->outFile,covPerIt, "Accumulative covering");
		hitList.print(p->outFile, "Accumulative number hits");
		desorbedList.print(p->outFile, "Accumulative number desorbed");
		errorList.print(p->outFile,errorPerIt, "Error per iteration");
	}
}

void SimulationHistory::write(std::string path){
	coveringList.write(path+"/covering.txt");
}


