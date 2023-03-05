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
 * Iterative simulation with contamination
 */

#include "Simulation.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/float128.hpp>

//#include <boost/accumulators/statistics/rolling_variance.hpp>
//#include <boost/accumulators/statistics/rolling_mean.hpp>

static std::string year[]={"Years","years","Year","year","Yr","yr","Y","y"};
static std::string month[]={"Months","months","Month","month","Mth","mth","mo","Mo"};
static std::string day[]={"day","Day","days","Days","d","D"};
static std::string hour[]={"hour","Hour","hours","Hours","h","H","hr","Hr"};
static std::string min[]={"Minutes","minutes","Minute","minute","min","Min","m","M"};
static std::string sec[]={"Seconds","seconds","Second","second","sec","Sec","s","S"};


const double diameterH2O = 2.76E-10;
const double kb = 1.38E-23;
const double h= 6.626E-34;

//const double tau = 1E-13;
//const double tuneE=2.64665;//tanh(2.64665)~0,99

//-----------------------------------------------------------

std::string home_to_tilde(std::string path);
int getFacetIndex(SubprocessFacet *iFacet);

//-----------------------------------------------------------

template <typename T> class HistoryList{
public:
	std::pair< std::vector<double>,std::vector<std::vector<T>> > historyList; // pair: list of times, list of facets
	std::vector<T> currentList; // list of facets
	std::vector<T> predictList; // Predicted per facet covering values at the end of the current time step.
	std::vector<std::pair<boost::multiprecision::float128,boost::multiprecision::float128>> statisticsList; // list of mean-std pair
	unsigned int currIt;


	HistoryList(){
		historyList.first = std::vector<double>();
		historyList.second = std::vector<std::vector<T>>();
		currentList=std::vector<T>();
		predictList=std::vector<T>(); // list containing facet values at the end of the predictor step
		statisticsList=std::vector<std::pair<boost::multiprecision::float128,boost::multiprecision::float128>> ();
		currIt=0;

	}

	void reset(){
		historyList.first.clear();
		historyList.second.clear();
		currentList.clear();
		predictList.clear();
		statisticsList.clear();
		currIt=0;
	}

	// initialize lists
	void initList(unsigned int numFacet){
		for(unsigned int i=0; i<numFacet; i++){
			historyList.second.push_back(std::vector<T>());
		}
	}

	void initStatistics(unsigned int numFacet){
		for(unsigned int i=0; i<numFacet; i++){
			statisticsList.push_back(std::make_pair(static_cast<boost::multiprecision::float128>(0),static_cast<boost::multiprecision::float128>(0)));
		}
	}

	void initCurrent(unsigned int numFacet){
		for(unsigned int i=0; i<numFacet; i++){
			currentList.push_back(static_cast<T>(0));
		}
		//initStatistics(numFacet);
	}

	void initPredict(unsigned int numFacet){
		for(unsigned int i=0; i<numFacet; i++){
			predictList.push_back(static_cast<T>(0));
		}
	}

	// append lists
	void appendCurrent(double time=-1){ //append currentList to historyList
		if(time==-1.0) //one step
				time=historyList.first.back()+1.0;
		historyList.first.push_back(time);
		for(unsigned int j=0; j < currentList.size(); j++){
			historyList.second[j].push_back(currentList[j]);
		}
		currIt+=1;
	}

	void appendPredict(double time=-1){ //append predictList to historyList
			if(time==-1.0) //one step
					time=historyList.first.back()+1.0;
			historyList.first.push_back(time);
			for(unsigned int j=0; j < predictList.size(); j++){
				historyList.second[j].push_back(predictList[j]);
			}
			currIt+=1;
		}

	// Statistics: Calculate mean/std per facet
	void updateStatistics(int rollingWindowSize, unsigned int offset=0){
		//std::cout<< historyList.second.back().size() <<std::endl;
		if(offset==0 && currIt>historyList.first.size()){
			offset=1;
		}

		for(uint j=0; j<statisticsList.size();j++){
			if(historyList.first.size()<rollingWindowSize+offset){ //Not enough samples for rolling window
				// mean and std = 0
				statisticsList[j].first=static_cast<boost::multiprecision::float128>(0);
				statisticsList[j].second=static_cast<boost::multiprecision::float128>(0);
			}
			else if(rollingWindowSize==1){
				// mean current value, std 0
				statisticsList[j].first =static_cast<boost::multiprecision::float128>(historyList.second[j].back());
				statisticsList[j].second=static_cast<boost::multiprecision::float128>(0);
			}
			else{ //Do statistics
				// mean = sum(point)/N
				T sum = std::accumulate(std::end(historyList.second[j])-rollingWindowSize, std::end(historyList.second[j]), static_cast<T>(0));
				statisticsList[j].first = static_cast<boost::multiprecision::float128>(sum) / static_cast<boost::multiprecision::float128>(rollingWindowSize);
				// std = sqrt (sum((point -mean)^2)/(N-1))
				boost::multiprecision::float128 accum = boost::multiprecision::float128(0.0);
				std::for_each (std::end(historyList.second[j])-rollingWindowSize, std::end(historyList.second[j]), [&](const T d) {
				    accum += (static_cast<boost::multiprecision::float128>(d) - statisticsList[j].first) * (static_cast<boost::multiprecision::float128>(d) - statisticsList[j].first);
				});
				statisticsList[j].second = boost::multiprecision::sqrt(accum/static_cast<boost::multiprecision::float128>(rollingWindowSize-1));
			}
		}
	}

	// Statistics: sum over ratio std/mean weighted with area
	boost::multiprecision::float128 getAverageStatistics(Simulation *sHandle,bool opacityCheck, bool doFocusOnly=false, std::vector<int> focusFacets=std::vector<int>()){
		double totalArea=0.0;
		bool meanZero=true;
		boost::multiprecision::float128 totalStatistics= boost::multiprecision::float128(0);
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				double idx=getFacetIndex(&f);
				// skip if not in focusGroup
				if(doFocusOnly && std::find(std::begin(focusFacets),std::end(focusFacets),idx)==std::end(focusFacets)) continue;

				if( (opacityCheck && f.sh.opacity!=0.0) || !opacityCheck){ // If opacity has to be checked (e.g. for covering): only consider facet if opacity is larger than 0
					totalArea+=f.sh.area;
					// area * std/mean
					totalStatistics+=boost::multiprecision::float128(f.sh.area)*statisticsList[idx].second/statisticsList[idx].first;
					if(statisticsList[idx].first!= boost::multiprecision::float128(0)){
						meanZero=false;
					}
				}
			}
		}
		// if mean is zero: return inf, otherwise return totalStatistics/totalArea
		return meanZero? std::numeric_limits<boost::multiprecision::float128>::infinity():totalStatistics/boost::multiprecision::float128(totalArea);
	}

	// Print lists
	std::string convertTime(double time){//convert seconds
		bool empty=true;
		std::string final="";
		std::div_t divresult;
		if(time/60.0>1){
			divresult=std::div((int)(time+0.5), 60); //min sec
			int seconds=divresult.rem;
			divresult=std::div(divresult.quot, 60); //hours min
			int minutes=divresult.rem;
			divresult=std::div(divresult.quot, 24); //days h
			int hours=divresult.rem;
			divresult=std::div(divresult.quot, 30.4375); //month days
			int days=divresult.rem;
			divresult=std::div(divresult.quot, 12); //years month
			int months=divresult.rem;
			int years=divresult.quot;

			if(years!=0) {final=final+std::to_string(years)+"y"; empty=false;}else{final=empty?final+"   ":final+"---";}//3
			if(months!=0) {final=(months>9)?final+std::to_string(months)+"mo":final+"0"+std::to_string(months)+"mo";empty=false;}else{final=empty?final+"    ":final+"----";}//4
			if(days!=0) {final=(days>9)?final+std::to_string(days)+"d":final+"0"+std::to_string(days)+"d";empty=false;}else{final=empty?final+"   ":final+"---";}//3
			if(hours!=0) {final=(hours>9)?final+std::to_string(hours)+"h":final+"0"+std::to_string(hours)+"h";empty=false;}else{final=empty?final+"   ":final+"---";}//3
			if(minutes!=0) {final=(minutes>9)?final+std::to_string(minutes)+"min":final+"0"+std::to_string(minutes)+"min";empty=false;}else{final=empty?final+"     ":final+"-----";}//5
			if(seconds!=0) {final=(seconds>9)?final+std::to_string(seconds)+"s":final+"0"+std::to_string(seconds)+"s";empty=false;}else{final=empty?final+"   ":final+"---";}//3
		}
		if(final==""){
			final=std::to_string((int)(time+0.5))+"s";
		}

		return final;
	}

	// Print historyList
	void print(std::ostream& outstream, std::string msg= "", int histSize = std::numeric_limits<int>::infinity(), std::vector<T> totalvec =std::vector<T>(), int entryWidth=7, bool convertToFloat=true){
		std::ostringstream out (std::ostringstream::app);

		uint offset_table=0;
		if(histSize != std::numeric_limits<int>::infinity()&& currIt > (unsigned int)histSize+1){
			offset_table=currIt - uint(histSize)-1;
		}
		if(msg!="")
			out<<std::endl <<msg <<std::endl;

		out <<std::setw(9)<<std::right<<"Iteration";
		out <<std::setw(14)<<std::right<<"Time[s]";
		out <<std::setw(22)<<std::right<<"Time";
		for(uint i=0;i<historyList.first.size();i++)
		{
			if(i==0){
				for(uint j=0; j<historyList.second.size();j++){
					out <<"\t" <<"Facet" <<std::setw(entryWidth)<<std::setfill('-')<<std::right <<j;
					}
				if(totalvec.size()==historyList.first.size())
					out <<"\t" <<std::setw(20)<<std::setfill(' ')<<std::right<<"Total";
			}
			out<<std::endl;

			out<<std::setw(9)<<std::setfill(' ')<<std::right <<i+offset_table*(i>0?1:0);
			out<<std::setw(14)<<std::right <<historyList.first[i] ;
			out<<std::setw(22)<<std::right <<convertTime(historyList.first[i]);

			for(uint j=0; j<historyList.second.size();j++)
			{
				if(convertToFloat)
					out <<"\t" <<std::setw(5+entryWidth)<<std::right <<boost::multiprecision::float128(historyList.second[j][i]);
				else
					out <<"\t" <<std::setw(5+entryWidth)<<std::right <<historyList.second[j][i];
			}
			if(totalvec.size()==historyList.first.size())
				out<<"\t"<<std::setw(20)<<std::right<<totalvec[i];

		}
		out<<std::endl<<std::endl;
		outstream <<out.str();
	}
	// Print currentList
	void printCurrent(std::ostream& outstream, std::string msg= ""){
		std::ostringstream tmpstream (std::ostringstream::app);

		tmpstream<<"    " <<std::setw(20)<<std::left<<msg;

		for(uint i=0;i<currentList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<boost::multiprecision::float128(currentList[i]);
		}
		tmpstream<<std::endl;
		outstream<<tmpstream.str();
	}

	void printPredict(std::ostream& outstream, std::string msg=""){
		std::ostringstream tmpstream (std::ostringstream::app);

		tmpstream<<"    " <<std::setw(20)<<std::left<<msg;

		for(uint i=0;i<predictList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<boost::multiprecision::float128(predictList[i]);
		}
		tmpstream<<std::endl;
		outstream<<tmpstream.str();
	}

	// Print statistictsList
	void printStatistics(std::ostream& outstream, std::string msg= "", int textwidth=45){
		std::ostringstream tmpstream (std::ostringstream::app);

		tmpstream<<std::left<<msg <<std::endl;
		tmpstream<<std::setw(textwidth)<<std::right <<"mean";
		for(uint i=0;i<statisticsList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<statisticsList[i].first;
		}
		tmpstream<<std::endl<<std::setw(textwidth)<<std::right <<"std";
		for(uint i=0;i<statisticsList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<statisticsList[i].second;
		}
		tmpstream<<std::endl;
		outstream<<tmpstream.str();
	}
	// Write historyList to file
	void write(std::string filename, int histSize = std::numeric_limits<int>::infinity()){
		std::ostringstream tmpstream (std::ostringstream::app);
		print(tmpstream,"",histSize,std::vector<T>(),15, true);

		std::ofstream out(filename,std::ofstream::out|std::ios::trunc);
		out<<tmpstream.str();
		out.close();
		std::cout <<"Results saved to " <<home_to_tilde(filename) <<std::endl;
	}
	// Erase entries at index
	void erase(int idx){
		historyList.first.erase(historyList.first.begin()+idx);
		for(unsigned int j=0; j<historyList.second.size();j++){
			historyList.second[j].erase(historyList.second[j].begin()+idx);
		}
	}

	bool empty(){return historyList.first.empty();}

	void setCurrent(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	currentList[covidx]=newValue;}
	void setCurrent(int idx, T newValue){currentList[idx]=newValue;}
	void setPredict(SubprocessFacet *iFacet, T newValue){int idx = getFacetIndex(iFacet);	predictList[idx]=newValue;}
	//T getLast(int idx){return historyList.second[idx].back();}
	T getLast(SubprocessFacet *iFacet){int covidx = getFacetIndex(iFacet);return historyList.second[covidx].back();}
	T getForelast(SubprocessFacet *iFacet){
		int covidx = getFacetIndex(iFacet);
		int index_fl = historyList.second[covidx].size() - 2;//index of the fore_last_vector_element
		//last element's = size -1, since first element's index = 0
		return historyList.second[covidx].at(index_fl);
	}
	T getCurrent(int idx){return currentList[idx];}
	T getCurrent(SubprocessFacet *iFacet){int covidx = getFacetIndex(iFacet);return currentList[covidx];}
	T getPredict(int idx){return predictList[idx];}
	T getPredict(SubprocessFacet *iFacet){int idx = getFacetIndex(iFacet);return predictList[idx];}
	void setLast(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	historyList.second[covidx].back()=newValue;}
	int getlastindex(){return historyList.first.size() - 1;}//returns the last index of the vector saving the points in time

	//---------------------------------------------------
	// Not used anymore
	/*
	void appendList(std::vector<T> List, double time=-1){

		if(time==-1.0) //One step
				time=historyList.back().first+1.0;
		historyList.first.push_back(time);
		for(int j=0; j < currentList.size(); j++){
			historyList.second[j].push_back(List[j]);
		}
		currIt+=1;
	}*/

};

class ProblemDef{
public:
	//ProblemDef(int argc, char *argv[]);
	ProblemDef();

	void readArg(int argc, char *argv[], int rank=1);
	bool readInputfile(std::string filename, int rank=1, int save=1);
	void writeInputfile(std::string filename, int rank=1);
	void printInputfile(std::ostream& out, bool printConversion=true);
	void SetFocusGroup(int facets);

	bool saveResults;
	bool saveConsole;

	std::string resultPath;
	std::string contaminationFlowPath;
	//std::string resultbufferPath;
	std::ofstream outFile;

	// These can be given as parameters directly
	std::string loadbufferPath;
	std::string hitbufferPath;
	std::string coveringPath;
	double simulationTime;
	std::string unit;

	// These can be given through input file only
	int iterationNumber; //number of iterations
	int usePCMethod; // 0: Do not use PC-Method, 1: Use PC-Method
	double maxTime;
	std::string maxUnit;

	double particleDia;

	double E_de;
	double sticking;

	double H_vap;
	//double W_tr;

	std::string errorMode;
	llong targetParticles;
	double targetError;
	llong targetParticles_input;
	double targetError_input;
	double noupdateError;

//	double hitRatioLimit;
	double t_min; // [s]

	double t_max; // [s]
	int maxTimePerIt; //[s]

	llong coveringMinThresh;

	int histSize;

	double outgassingTimeWindow; //[s] uniform distribution of outgassing over outgassingTimeWindow
	double counterWindowPercent; // [%]
//	double desWindowPercent; // [%]

	int rollingWindowSize;
	double convergenceTarget;
	double convergenceTime; //[s]
	bool stopConverged;

	std::vector< std::pair<int,double> > vipFacets;
	std::vector< std::vector<int> > facetGroups; // vector that includes vector of facet groups
	std::pair<std::vector<int>,std::vector<int>> focusGroup; // pair: facet group indices and corresponding facet indices
	bool doFocusGroupOnly; // determines whether facets in focus groups are used for error/convergence or all facets

	//These cannot be given, but are computed from other variables
	int simulationTimeMS;
	double maxTimeS;
	bool doCoveringFile;

private:
	void createOutput(int save);
};

class SimulationHistory{
public:
	SimulationHistory(int);
	SimulationHistory(Databuff *hitbuffer,int);

	HistoryList<boost::multiprecision::uint128_t> coveringList;
	HistoryList<double> hitList;
	HistoryList<llong> desorbedList;
	HistoryList<llong> outgassedList;
	HistoryList<llong> adsorbedList;
	HistoryList<double> errorList_event;
	HistoryList<double> errorList_covering;
	HistoryList<double> particleDensityList;
	HistoryList<double> pressureList;

	//std::vector<unsigned int> normalFacets;

	llong smallCoveringFactor;

	unsigned int numFacet;
	int numSubProcess;

	double flightTime; // [s]
	int nParticles;
	int nLeaks;

	double lastTime; // [s]
	int currentStep;
	int pcStep; // 0: In predictor step. 1: In corrector step
	double stepSize; // [s]
	double stepSize_outgassing; //[s]

	//void appendList(Databuff *hitbuffer, double time=-1.0);
	//void appendList(double time=-1.0);
	void erase(int idx);

	void print();
	void write(std::string path);
	//std::tuple<bool, llong > updateHistory(Databuff *hitbuffer);
	void updateHistory();

	void updateStepSize_outgassing();

};

//-----------------------------------------------------------
//SimulationLinux.cpp
std::tuple<bool, std::vector<int>>  simulateSub2(Databuff *hitbuffer, int rank, int simutime);

double convertunit(double simutime, std::string unit);
void printStream(std::string string, bool print=true);
void checkSmallCovering(int rank, Databuff *hitbuffer_sum);
bool readCovering(Databuff* hitbuffer, std::string coveringFile, int rank);

//ProblemDef
//SimulationHistory

//-----------------------------------------------------------
//UpdateSubProcess.cpp

void UpdateSticking();
bool UpdateDesorption();
void UpdateSojourn();

double CalcErrorSub(std::string mode);

void UpdateMCSubHits(Databuff *databuffer, int rank);
void initbufftozero(Databuff *databuffer);

//-----------------------------------------------------------
//UpdateMainProcess.cpp

double getStepSize();

void UpdateCovering(Databuff *hitbuffer_sum);
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer, bool step_size_change);
void UpdateErrorMain();
void UpdateParticleDensityAndPressure(Databuff *hitbuffer_sum);

std::tuple<std::vector<double>,std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration();

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank);

//-----------------------------------------------------------
//SimulationCalc.cpp
FacetHitBuffer* getFacetHitBuffer(SubprocessFacet *iFacet, Databuff *hitbuffer);
double calcNmono(SubprocessFacet *iFacet);

llong getnbDesorbed(Databuff *hitbuffer_sum);
llong getnbOutgassed(Databuff *hitbuffer_sum);
llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer);
llong getnbOutgassed(SubprocessFacet *iFacet, Databuff *hitbuffer);
llong getnbAdsorbed(SubprocessFacet *iFacet, Databuff *hitbuffer);//In the original Molflow, particles were absorbed not adsorbed. In ContaminationFlow we regard all old code parts
// which are called 'Absorb' actually as an 'Adsorb'. But we did not rename them.
llong getnbHits(Databuff *hitbuffer_sum);

boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer);
boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet);
boost::multiprecision::uint128_t getPredictedCovering(SubprocessFacet *iFacet);

double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer);
//std::tuple<double, double, double> getVelocities(SubprocessFacet *iFacet, Databuff *hitbuffer_sum);


//----------
boost::multiprecision::float128 calctotalDesorption();
double calcOutgassingFactor(SubprocessFacet *iFacet);

boost::multiprecision::float128 calcCoverage(SubprocessFacet *iFacet);
boost::multiprecision::float128 calcPredictedCoverage(SubprocessFacet *iFacet);
boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum);
void calcSticking(SubprocessFacet *iFacet);

boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet);
double calcParticleDensity(Databuff *hitbuffer_sum , SubprocessFacet *f);
double calcPressure(Databuff *hitbuffer_sum , SubprocessFacet *f);

double calcStartTime(SubprocessFacet *iFacet, bool desorbed_b, bool printWarning=false);

//-----------------------------------------------------------
//Iteration.cpp

std::tuple<double,double> CalcErrorAll(int it=-1);
void UpdateErrorList(Databuff *hitbuffer_sum);
bool checkError(double targetError, double currentError, double factor, std::string mode);
HistoryList<double>* getErrorList(std::string mode);

//void setCoveringThreshold(Databuff *hitbuffer, int size, int rank);
void setCoveringThreshold(int size, int rank);
void initCoveringThresh();

std::tuple<bool, double> TimestepControl(Databuff *hitbuffer_sum);

//double estimateTmin(Databuff *hitbuffer);
double estimateAverageFlightTime(Databuff *hitbuffer_sum);


//-----------------------------------------------------------
//worker.cpp
void CalcTotalOutgassingWorker();
