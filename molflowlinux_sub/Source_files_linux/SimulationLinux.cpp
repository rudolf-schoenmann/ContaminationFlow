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
 * This file contains the serialization of the simulation in the sub processes
 */

#include "SimulationLinux.h"
#include <array>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <libgen.h>
#include <time.h>
#include <sys/stat.h>
#include <mpi.h>

extern SimulationHistory* simHistory;
extern Simulation *sHandle;
extern ProblemDef* p;

// Simulation on subprocess
std::tuple<bool, std::vector<int> > simulateSub2(Databuff *hitbuffer,int rank, int simutime){

	//Calculate target values for error and number of desorbed particles
	int targetParticles=p->targetParticles/simHistory->numSubProcess;
	double targetError=p->targetError*pow(simHistory->numSubProcess,0.5);

	if(rank==1){
		std::ostringstream tmpstream (std::ostringstream::app);
		tmpstream <<std::endl <<"Currentstep: " << simHistory->currentStep <<". Step size: " <<simHistory->stepSize <<std::endl;
		tmpstream<<simHistory->lastTime <<" + " <<simHistory->stepSize << " = " <<simHistory->lastTime+simHistory->stepSize << " < " <<p->outgassingTimeWindow <<" ? => stepSize_outgassing = " <<simHistory->stepSize_outgassing <<std::endl;
		tmpstream <<"Target Particles: " << targetParticles <<". Target Error: " <<targetError <<std::endl <<std::endl;
		printStream(tmpstream.str());

	}

	//Values for simulation
	double timestep=1000; // Desired length per iteration step for simulation, here hardcoded to 1 second
	double realtimestep; // Actual time elapsed for iteration step
	double i;			// Time elapsed between checking of targets (i.e., if each iteration step)
	double totalTime=0.0;	//Total simulation time
	bool eos=false;		// End of simulation flag
	double totalError=1.;//Total error

	int j_old=0;
	bool j_print=false;

	// Facets that have reached the covering threshold
	std::vector<int> facetNum;
	facetNum =std::vector<int> ();

	//----Simulation

	// Start Simulation = create first particle
	if(!StartSimulation())
		return {std::make_tuple(true,facetNum)};



	// Run Simulation for timestep milliseconds
	// Termination condition: maximum number of iteration reached or target error & number of desorbed particles reached or maximum desoption/covering threshold reached
	for(int j=0; j<p->maxSimPerIt && !(j>0 && simHistory->nParticles>targetParticles && checkErrorSub(targetError, totalError, pow(simHistory->numSubProcess,0.5)))&& !eos; j++){
		for(i=0; i<(double)(simutime) && !eos;i+=realtimestep){
			// Until simutime is reached, do simulation
			if(i>=(double(simutime)*0.99)){break;}
			if(simHistory->coveringList.empty()){
				simHistory->appendList(i); //append list with initial covering
				}

			if(i+timestep>=(double)(simutime)){ //last timestep
				std::tie(eos,realtimestep) = SimulationRun((double)simutime-i); // Some additional simulation to reach desired simutime, as iteration step  does not run for exactly timestep ms
				}
			else{
				std::tie(eos, realtimestep) = SimulationRun(timestep);      // Run for timestep ms, performs MC steps
			}

		}
		totalTime+=i;
		j_old=j;
		if(i/simutime>=2.0*simutime/1000.0){
			//j+=int(i/simutime -1); // correct number of steps for significantly longer simulation MCSteps: for small sojourn time, sHandle->stepPerSec can be smaller than 1 => 1 MCStep can take significantly longer than 1 second
			j=int((totalTime-simutime)/1000.0);
		}

		j_print=false;
		for(; j_old<=j;j_old++){
			if(j_old%(int)(30000/simutime)==0){
				j_print=true;
				break;
			}
		}

		// Calculate error for this iteration step
		totalError=UpdateError("covering");

		if(j_print || (simHistory->nParticles>targetParticles && checkErrorSub(targetError, totalError, pow(simHistory->numSubProcess,0.5)))|| eos || j >= p->maxSimPerIt-1){
			// Print current history lists every 30s or if target reached
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<" Subprocess "<<rank<<": Step "<<std::setw(4)<<std::right <<j <<"    &    Total time " <<std::setw(10)<<std::right <<totalTime <<"ms    &    Adsorbed particles "<<std::setw(10)<<std::right<<simHistory->nParticles <<"    &    Total error "  <<std::setw(10)<<std::left<<totalError<<std::endl;
			tmpstream << std::endl;

			tmpstream << "Facet" << std::setw(23)<<std::right<< " ";
			int num;
			for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
				for (SubprocessFacet& f : sHandle->structures[j].facets) {
						num=getFacetIndex(&f);
						tmpstream <<"\t"<< std::setw(12)<<std::right << num;
				}
			}
			tmpstream << std::endl;
			simHistory->hitList.printCurrent(tmpstream, std::to_string(rank)+": hitlist");
			simHistory->desorbedList.printCurrent(tmpstream, std::to_string(rank)+": desorbedlist");
			simHistory->coveringList.printCurrent(tmpstream, std::to_string(rank)+": coveringlist");
			//simHistory->errorList_event.printCurrent(tmpstream, std::to_string(rank)+": errorlist_event");
			simHistory->errorList_covering.printCurrent(tmpstream, std::to_string(rank)+": errorlist_covering");
			tmpstream <<std::endl;

			printStream(tmpstream.str());
		}
	}

	// Find Facets that have reached the covering threshold
	int num;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(f.sh.desorption==0) continue;
			num=getFacetIndex(&f);

			/* Just output stuff
			llong cv = f.tmpCounter[0].hit.covering;
			std::cout << "At the end of the calculation rank " << rank << " f.tmpCounter[0].hit.covering = " << cv << std::endl;
			p->outFile << "At the end of the calculation rank " << rank << " f.tmpCounter[0].hit.covering = " << cv << std::endl;
			*/

			if(f.tmpCounter[0].hit.covering<=sHandle->coveringThreshold[num]) {facetNum.push_back(num);}
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

void printStream(std::string string){
	std::cout <<string;
	p->outFile <<string;
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

//-----------------------------------------------------------
//----ProblemDef class
ProblemDef::ProblemDef(){
	loadbufferPath= "/home/van/Buffer/loadbuffer_alle_RT";
	hitbufferPath="/home/van/Buffer/hitbuffer_allee-6";

	iterationNumber = 43200;
	particleDia=carbondiameter;
	E_de=1E-21;
	H_vap=0.8E-19;
	W_tr=1.0;
	sticking=0.0;

	maxTime=10.0;
	maxUnit="y";
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

	simulationTime = 10.0;
	unit = "s";
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);

	saveResults=true;

	targetParticles=1000;
	targetError=0.001;
	hitRatioLimit=0;
	t_min=1E-4;

	t_max=std::numeric_limits<double>::max();
	maxSimPerIt=std::numeric_limits<int>::max();
	histSize=std::numeric_limits<int>::max();

	outgassingTimeWindow=0.0;
	counterWindowPercent=0.1;
	desWindowPercent=1.0;

	rollingWindowSize=10;

	coveringMinThresh=1000000;

	vipFacets = std::vector< std::pair<int,double> >();

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

	while(std::getline(input,line)){
		std::string stringIn;
		double doubleIn;
		int intIn;
		llong llongIn;
		std::istringstream is( line );

		is >> stringIn;

		if(stringIn == "loadbufferPath") {is >> stringIn; loadbufferPath=stringIn;}
		else if(stringIn == "hitbufferPath") {is >> stringIn; hitbufferPath=stringIn;}
		else if(stringIn == "simulationTime") {is >>doubleIn; simulationTime = doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn == "unit"){is >> stringIn; unit=stringIn;}

		else if(stringIn =="iterationNumber"){is >> intIn; iterationNumber=intIn>0?intIn:0;}
		else if(stringIn == "maxTime") {is >>doubleIn; maxTime = doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn == "maxUnit"){is >> stringIn; maxUnit=stringIn;}

		else if(stringIn =="particleDia"){is >> doubleIn; particleDia=doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="E_de"){is >> doubleIn; E_de=doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="H_vap"){is >> doubleIn; H_vap=doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="W_tr"){is >> doubleIn; W_tr=doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="sticking"){is >> doubleIn; sticking=doubleIn>0.0?doubleIn:0.0;}

		else if(stringIn =="targetParticles"){is >> intIn; targetParticles=intIn>0?intIn:0;}
		else if(stringIn == "targetError") {is >>doubleIn; targetError = doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn == "hitRatioLimit") {is >>doubleIn; hitRatioLimit = doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn == "t_min") {is >>doubleIn; t_min = doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="maxStepSize" || stringIn =="t_max"){is >>doubleIn; t_max=doubleIn>0.0?doubleIn:0.0;}

		else if(stringIn == "coveringMinThresh") {is >>llongIn; coveringMinThresh = llongIn;}

		else if(stringIn =="maxSimPerIt"){is >> intIn; maxSimPerIt=intIn>1?intIn:1;}
		else if(stringIn =="histSize"){is >> intIn; histSize=intIn>1?intIn:1;}
		else if(stringIn =="counterWindowPercent"){is >>doubleIn; doubleIn=doubleIn<1.0?doubleIn:1.0; counterWindowPercent=doubleIn>0.0?doubleIn:0.0;}
		else if(stringIn =="desWindowPercent"){is >>doubleIn; doubleIn=doubleIn<1.0?doubleIn:1.0; desWindowPercent=doubleIn>0.0?doubleIn:0.0; }
		else if(stringIn == "outgassingTimeWindow"){is >>doubleIn; outgassingTimeWindow=doubleIn>0.0?doubleIn:0.0; }

		else if(stringIn =="rollingWindowSize"){is >> intIn; rollingWindowSize=intIn>1?intIn:1;}

		else if(stringIn=="vipFacets"){
			int vipf = 0; double vipe=0.0;
			unsigned int vipfacet=0;
			while(is>>doubleIn){
				if(vipfacet%2==0){
					if(rint(doubleIn)==doubleIn){
						vipf=rint(doubleIn);
					}
					else{
						break;
					}
				}
				else{
					vipe = doubleIn;
					vipFacets.push_back(std::make_pair(vipf,vipe));
				}


				vipfacet+=1;
			}
			if(rank==0) {
				std::cout <<vipFacets.size()<<" vip facet(s)" <<std::endl;
				for(vipfacet=0; vipfacet<vipFacets.size(); vipfacet++){
					std::cout <<"\t"<<vipFacets[vipfacet].first <<"\t" <<vipFacets[vipfacet].second <<std::endl;
				}
			}
		}

		else{if(rank==0) std::cout <<stringIn <<" not a valid argument." <<std::endl;}

	}
	simulationTimeMS = (int) (convertunit(simulationTime, unit) + 0.5);
	maxTimeS=convertunit(maxTime, maxUnit)/1000.0;

	if(saveResults)
		writeInputfile(resultpath+"/InputFile.txt",rank);
}

void ProblemDef::writeInputfile(std::string filename, int rank){
	if(rank==0){
		std::ofstream outfile(filename,std::ofstream::out|std::ios::trunc);
		printInputfile(outfile, false);
	}
}

void ProblemDef::printInputfile(std::ostream& out, bool printConversion){ //std::cout or p->outFile
	if(printConversion) out  <<std::endl<<"Print input arguments"<<std::endl;

	if(printConversion) out  <<"resultPath" <<'\t' <<resultpath <<std::endl;
	out  <<"loadbufferPath" <<'\t' <<loadbufferPath <<std::endl;
	out  <<"hitbufferPath" <<'\t' <<hitbufferPath <<std::endl;
	if(printConversion) out  <<"resultbufferPath" <<'\t' <<resultbufferPath <<std::endl;
	if(printConversion) out <<std::endl;

	out  <<"simulationTime" <<'\t' <<simulationTime <<std::endl;
	out  <<"unit" <<'\t' <<unit <<std::endl;

	out  <<"iterationNumber" <<'\t' <<iterationNumber<<std::endl;
	out  <<"maxTime" <<'\t' <<maxTime <<std::endl;
	out  <<"maxUnit" <<'\t' <<maxUnit <<std::endl;
	if(printConversion) out <<std::endl;

	out <<"sticking" <<'\t' <<sticking <<std::endl;
	out <<"particleDia" <<'\t' <<particleDia<<std::endl;
	out  <<"E_de" <<'\t' <<E_de<<std::endl;
	out <<"H_vap" <<'\t' <<H_vap <<std::endl;
	out <<"W_tr" <<'\t' <<W_tr <<std::endl;

	out <<"targetError" <<'\t' <<targetError <<std::endl;
	out <<"targetParticles" <<'\t' <<targetParticles <<std::endl;
	out <<"hitRatioLimit" <<'\t' <<hitRatioLimit <<std::endl;
	out <<"t_min" <<'\t' <<t_min <<std::endl;

	out <<"t_max" <<"\t" <<t_max<<std::endl;
	out <<"maxSimPerIt" <<"\t" <<maxSimPerIt<<std::endl;

	out <<"coveringMinThresh" <<"\t" <<coveringMinThresh <<std::endl;

	out <<"histSize" <<"\t" <<histSize<<std::endl;
	out <<"counterWindowPercent" <<"\t" <<counterWindowPercent<<std::endl;
	out <<"desWindowPercent" <<"\t" <<desWindowPercent<<std::endl;
	out <<"outgassingTimeWindow" <<"\t" <<outgassingTimeWindow <<std::endl;

	out <<"rollingWindowSize" <<"\t" <<rollingWindowSize <<std::endl;


	if(!vipFacets.empty()){
		out <<"vipFacets";
		for(unsigned int i=0; i<vipFacets.size();i++){
			out <<"\t" <<vipFacets[i].first <<"\t" <<vipFacets[i].second;
		}
		out <<std::endl;
	}
	if(printConversion){
		out  << "Simulation time " << simulationTime << unit << " converted to " << simulationTimeMS << "ms" << std::endl;
		out  << "Maximum simulation time " << maxTime << maxUnit << " converted to " << maxTimeS << "s" << std::endl<<std::endl;
	}
}

//-----------------------------------------------------------
void checkSmallCovering(int rank, Databuff *hitbuffer_sum){
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;

	llong smallCoveringFactor=1;
	boost::multiprecision::uint128_t covering;
	boost::multiprecision::uint128_t mincov = boost::multiprecision::uint128_t(p->coveringMinThresh);

	bool smallCovering=false;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering = getCovering(&f);

			if(llong(covering)<p->coveringMinThresh && f.sh.desorption>0.0 && covering >0){
				smallCovering=true;
				if(covering<mincov){
					mincov=covering;
				}
				//std::cout <<"Facet "<<numFacet-1 <<": "<<llong(covering) <<" smaller than " <<p->coveringMinThresh <<std::endl;
			}
		}
	}

	if(smallCovering){
		smallCoveringFactor=llong(1.0+1.1*double(p->coveringMinThresh)/(double(mincov)));
		/*
		std::cout <<"Small covering found for rank " << rank << ": multiply covering and threshold by " <<smallCoveringFactor <<" for mincov "<< mincov <<std::endl;
		p->outFile<<"Small covering found for rank " << rank << ": multiply covering and threshold by " <<smallCoveringFactor <<" for mincov "<< mincov <<std::endl;
		*/
		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				FacetHitBuffer *facetHitSum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				facetHitSum->hit.covering *=smallCoveringFactor;
				f.tmpCounter[0].hit.covering *=smallCoveringFactor;
				sHandle->coveringThreshold[getFacetIndex(&f)] *= smallCoveringFactor;
			}
		}
	}
	simHistory->smallCoveringFactor=smallCoveringFactor;
}

//----SimulationHistory class

SimulationHistory::SimulationHistory(int world_size){
	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;
	currentStep=0;
	stepSize=0.0;
	stepSize_outgassing=0.0;
	numSubProcess=world_size-1;
	smallCoveringFactor=1;


	//normalFacets = std::vector<unsigned int>();

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if(!p->vipFacets.empty()){
				for(unsigned int i =0; i< p->vipFacets.size(); i++ ){
					if(int(numFacet)==p->vipFacets[i].first){
						f.sh.isVipFacet=true;
						break;
					}
				}
			}
			//if(!f.sh.isVipFacet)
			//	normalFacets.push_back(numFacet);
			numFacet+=1;
		}
	}
	coveringList.initList(numFacet);
	coveringList.initCurrent(numFacet);

	hitList.initList(numFacet);
	hitList.initCurrent(numFacet);

	desorbedList.initList(numFacet);
	desorbedList.initCurrent(numFacet);

	errorList_event.initList(numFacet);
	errorList_event.initCurrent(numFacet);

	errorList_covering.initList(numFacet);
	errorList_covering.initCurrent(numFacet);

	particleDensityList.initList(numFacet);
	particleDensityList.initCurrent(numFacet);

	pressureList.initList(numFacet);
	pressureList.initCurrent(numFacet);
}

SimulationHistory::SimulationHistory(Databuff *hitbuffer, int world_size){
	numFacet=0;
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	//normalFacets = std::vector<unsigned int>();

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

			if(!p->vipFacets.empty()){
				for(unsigned int i =0; i< p->vipFacets.size(); i++ ){
					if(int(numFacet)==p->vipFacets[i].first){
						f.sh.isVipFacet=true;
						break;
					}
				}
			}
			//if(!f.sh.isVipFacet)
			//	normalFacets.push_back(numFacet);
			numFacet+=1;
		}
	}
	coveringList.initList(numFacet);
	coveringList.appendCurrent(0);
	coveringList.initStatistics(numFacet);

	hitList.initList(numFacet);
	hitList.appendCurrent(0);

	desorbedList.initList(numFacet);
	desorbedList.appendCurrent(0);

	errorList_event.initList(numFacet);
	errorList_event.initCurrent(numFacet);
	errorList_event.appendCurrent(0.0);

	errorList_covering.initList(numFacet);
	errorList_covering.initCurrent(numFacet);
	errorList_covering.appendCurrent(0.0);

	particleDensityList.initList(numFacet);
	particleDensityList.initCurrent(numFacet);
	particleDensityList.appendCurrent(0.0);

	pressureList.initList(numFacet);
	pressureList.initCurrent(numFacet);
	pressureList.appendCurrent(0.0);

	currentStep=0;
	stepSize=0.0;
	stepSize_outgassing=0.0;

	numSubProcess=world_size-1;
	smallCoveringFactor=1;

	//std::cout<<"Normal facets: ";
	//for (unsigned int i =0; i< normalFacets.size(); i++){
	//	std::cout <<"\t" << normalFacets[i];
	//}
	//std::cout <<std::endl;

}

void SimulationHistory::updateHistory(){

	nParticles=0;
	flightTime=0.0;

	boost::multiprecision::uint128_t covering;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering = getCovering(&f);

			f.tmpCounter[0].hit.covering=covering.convert_to<llong>();
		}
	}
	coveringList.historyList.first.clear();
	coveringList.historyList.second.clear();
	coveringList.initList(numFacet);
	coveringList.appendCurrent(0);

	hitList.reset();
	hitList.initList(numFacet);
	hitList.initCurrent(numFacet);
	hitList.appendCurrent(0);

	desorbedList.reset();
	desorbedList.initList(numFacet);
	desorbedList.initCurrent(numFacet);
	desorbedList.appendCurrent(0);

	errorList_event.reset();
	errorList_event.initList(numFacet);
	errorList_event.initCurrent(numFacet);

	errorList_covering.reset();
	errorList_covering.initList(numFacet);
	errorList_covering.initCurrent(numFacet);

	particleDensityList.reset();
	particleDensityList.initList(numFacet);
	particleDensityList.initCurrent(numFacet);
	particleDensityList.appendCurrent(0.0);

	pressureList.reset();
	pressureList.initList(numFacet);
	pressureList.initCurrent(numFacet);
	pressureList.appendCurrent(0.0);

	updateStepSize();
	//lastTime=0.0;
}

void SimulationHistory::updateStepSize(){
	stepSize = getStepSize();

	if(lastTime+stepSize<=p->outgassingTimeWindow){
		stepSize_outgassing = stepSize;
	}
	else if(lastTime<p->outgassingTimeWindow){
		stepSize_outgassing = (p->outgassingTimeWindow-lastTime);
	}
	else{
		stepSize_outgassing = 0.0;
	}
}

void SimulationHistory::appendList(double time){

	if(time==-1.0) //One step
		time=coveringList.historyList.first.back()+1.0;

	coveringList.appendCurrent(time);
}

void SimulationHistory::erase(int idx){
	coveringList.erase(idx);
	errorList_event.erase(idx);
	errorList_covering.erase(idx);
	particleDensityList.erase(idx);
	pressureList.erase(idx);
}

void SimulationHistory::print(bool write){
	std::vector<double> errorPerIt_event;//This is the error_event (Desorb + Hit)
	std::vector<double> errorPerIt_covering;//This is the error_event (Desorb + Adsorb)
	std::vector<boost::multiprecision::uint128_t> covPerIt;
	std::tie(errorPerIt_event, errorPerIt_covering,covPerIt) = CalcPerIteration();

	coveringList.print(p->outFile,covPerIt, "Accumulative covering", p->histSize, true, write);
	//hitList.print(p->outFile, "Accumulative number hits", p->histSize, true,write);//Since we do not accumulate hits anymore over all iterations, we do not need this anymore.
	//desorbedList.print(p->outFile, "Accumulative number desorbed", p->histSize,true,write);//Since we do not accumulate desorbs anymore over all iterations, we do not need this anymore.
	errorList_event.print(p->outFile,errorPerIt_event, "Error Event (Desorb + Hit) per iteration", p->histSize, true,write);
	errorList_covering.print(p->outFile, errorPerIt_covering, "Error Covering (Desorb + Adsorb) per iteration", p->histSize,true,write);
	particleDensityList.print(p->outFile, "Particle density per iteration", p->histSize,true,write);
	pressureList.print(p->outFile, "Pressure per iteration", p->histSize,true,write);
}

void SimulationHistory::write(std::string path){
	coveringList.write(path+"/covering.txt", p->histSize);
	errorList_covering.write(path+"/errorCovering.txt", p->histSize);
	particleDensityList.write(path+"/particleDensity.txt", p->histSize);
	pressureList.write(path+"/pressure.txt", p->histSize);
}

//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
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
	for(int j=0; j<p->maxSimPerIt && !(j>0 && simHistory->nParticles>targetParticles && checkErrorSub(targetError, totalError, pow(simHistory->numSubProcess,0.5)))&& !eos; j++){

		for(i=0; i<(double)(simutime) && !eos;i+=realtimestep){

			if(i>=(double(simutime)*0.99)){break;}
			if(simHistory->coveringList.empty()){
				simHistory->appendList(hitbuffer,i); //append list with initial covering
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

		if(j%(int)(30000/simutime)==0 || (simHistory->nParticles>targetParticles && checkErrorSub(targetError, totalError, pow(simHistory->numSubProcess,0.5)))|| eos){
			// Print information every 30s or if target reached
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"  "<<rank<<": Step "<<std::setw(4)<<std::right <<j <<"    &    Total time " <<std::setw(10)<<std::right <<totalTime <<"ms    &    Desorbed particles "<<std::setw(10)<<std::right<<simHistory->nParticles <<"    &    Total error "  <<std::setw(10)<<std::left<<totalError<<std::endl;

			if(!checkErrorSub(targetError, totalError, pow(simHistory->numSubProcess,0.5))){
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

std::tuple<bool, llong > SimulationHistory::updateHistory(Databuff *hitbuffer){
	nParticles=0;
	flightTime=0.0;
	lastTime=0.0;

	startNewParticle=false;


	double numHit;
	llong numDes;
	llong smallCoveringFactor=0;
	boost::multiprecision::uint128_t covering;

	boost::multiprecision::uint128_t mincov = boost::multiprecision::uint128_t(p->coveringMinThresh);

	bool smallCovering=false;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering = boost::multiprecision::uint128_t(getCovering(&f, hitbuffer));
			numHit=getHits(&f, hitbuffer);
			numDes=getnbDesorbed(&f, hitbuffer);

			coveringList.setCurrent(&f, covering);
			hitList.setCurrent(&f, numHit);
			desorbedList.setCurrent(&f, numDes);

			f.tmpCounter[0].hit.covering=covering.convert_to<llong>();

			if(llong(covering)<p->coveringMinThresh && f.sh.desorption>0.0 && covering >0){
				smallCovering=true;
				if(covering<mincov){
					mincov=covering;
				}
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

	errorList_event.reset();
	errorList_event.initCurrent(numFacet);
	errorList_covering.reset();
	errorList_covering.initCurrent(numFacet);

	//stepSize=manageStepSize();

	if(smallCovering && llong(mincov) < p->coveringMinThresh){
		smallCoveringFactor=llong(1.0+1.1*double(p->coveringMinThresh)/(double(mincov)));

		//std::cout <<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for mincov "<< mincov <<std::endl;
		//p->outFile<<"Small covering: multiply covering and threshold by " <<smallCoveringFactor <<" for mincov "<< mincov <<std::endl;

		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				if(f.sh.desorption==0) continue;
				f.tmpCounter[0].hit.covering*=smallCoveringFactor;
				sHandle->coveringThreshold[getFacetIndex(&f)]*=smallCoveringFactor;
			}
		}
	}

	return {std::make_tuple(smallCovering,smallCoveringFactor)};
}
*/

//----------deprecated functions because of new K_real/virt approach
/*
void UndoSmallCovering(Databuff *hitbuffer_sum){
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitSum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				facetHitSum->hit.covering/=simHistory->smallCoveringFactor;
				simHistory->coveringList.setLast(&f, simHistory->coveringList.getLast(&f)/boost::multiprecision::uint128_t(simHistory->smallCoveringFactor));
			}
	}
}
*/

//----------Not used anymore
/*
void SimulationHistory::appendList(Databuff *hitbuffer, double time){

	std::vector<boost::multiprecision::uint128_t> currentCov;
	currentCov =std::vector<boost::multiprecision::uint128_t> ();

	if(time==-1.0) //One step
		time=coveringList.historyList.back().first+1.0;

	boost::multiprecision::uint128_t covering;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			covering=boost::multiprecision::uint128_t(getCovering(&f, hitbuffer));
			currentCov.push_back(covering);
		}
	}
	coveringList.appendList(currentCov, time);
}
*/
