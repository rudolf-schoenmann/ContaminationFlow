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
//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/multiprecision/float128.hpp>

static const char *year[]={"Years","years","Year","year","Yr","yr","Y","y"};
static const char *month[]={"Months","months","Month","month","Mth","mth","mo","Mo"};
static const char *day[]={"day","Day","days","Days","d","D"};
static const char *hour[]={"hour","Hour","hours","Hours","h","H","hr","Hr"};
static const char *min[]={"Minutes","minutes","Minute","minute","min","Min","m","M"};

const double carbondiameter = 2 *76E-12;
const double kb = 1.38E-23;
//const double tau = 1E-13;
const double h= 6.626E-34;
const double tuneE=2.64665;//tanh(2.64665)~0,99

template <typename T> class HistoryList{
public:
	HistoryList(){
		pointintime_list = std::vector< std::pair<double,std::vector<T>> >();
		currentList=std::vector<T>();
		currIt=0;
	}

	std::vector< std::pair<double,std::vector<T>> > pointintime_list;
	std::vector<T> currentList;
	unsigned int currIt;

	void reset(unsigned int numFacet=0){
		pointintime_list.clear();
		currentList.clear();
		currIt=0;
	}

	void initCurrent(unsigned int numFacet){
		for(unsigned int i=0; i<numFacet; i++){
			currentList.push_back(0);
		}
	}

	void appendCurrent(double time=-1){
		if(time==-1.0) //one step
				time=pointintime_list.back().first+1.0;
		pointintime_list.push_back(std::make_pair(time,currentList));
		currIt+=1;

	}

	void appendList(std::vector<T> List, double time=-1){

		if(time==-1.0) //One step
				time=pointintime_list.back().first+1.0;
		pointintime_list.push_back(std::make_pair(time,List));
		currIt+=1;

	}
	std::string convertTime(double time){
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
	void print(std::ostream& out, std::string msg= "", int histSize = std::numeric_limits<int>::infinity()){

		uint offset_table=0;
		if(histSize != std::numeric_limits<int>::infinity() && currIt > histSize+1){
			offset_table=currIt - uint(histSize)-1;
		}

		out<<std::endl <<msg <<std::endl;

		out <<std::setw(9)<<std::right<<"Iteration\t";
		out <<std::setw(11)<<std::right<<"Time[s]";
		out <<std::setw(22)<<std::right<<"Time";
		for(uint i=0;i<pointintime_list.size();i++)
		{
			if(i==0){
				for(uint j=0; j<pointintime_list[i].second.size();j++)
						{
					out <<"\t" <<std::setw(6)<<std::right <<"Facet " <<std::setw(8)<<std::right <<j;
						}
			}
			out<<std::endl;
			out<<std::setw(9)<<std::right <<i+offset_table*(i>0?1:0)<<"\t";
			out<<std::setw(11)<<std::right <<pointintime_list[i].first ;
			out<<std::setw(22)<<std::right <<convertTime(pointintime_list[i].first);

			for(uint j=0; j<pointintime_list[i].second.size();j++)
			{
				out <<"\t" <<std::setw(14)<<std::right <<boost::multiprecision::float128(pointintime_list[i].second[j]);

			}

		}
		out<<std::endl<<std::endl;
	}

	void print(std::ostream& out, std::vector<T> totalvec, std::string msg= "", int histSize = std::numeric_limits<int>::infinity()){

		uint offset_table=0;
		if(histSize != std::numeric_limits<int>::infinity()&& currIt > histSize+1){
			offset_table=currIt - uint(histSize)-1;
		}
		out<<std::endl <<msg <<std::endl;

		out <<std::setw(9)<<std::right<<"Iteration\t";
		out <<std::setw(11)<<std::right<<"Time[s]";
		out <<std::setw(22)<<std::right<<"Time";
		for(uint i=0;i<pointintime_list.size();i++)
		{
			if(i==0){
				for(uint j=0; j<pointintime_list[i].second.size();j++){
					out <<"\t" <<std::setw(6)<<std::right <<"Facet-" <<std::setw(8)<<std::setfill('-')<<std::right <<j;
					}

				out <<"\t" <<std::setw(14)<<std::setfill(' ')<<std::right<<"Total";
			}
			out<<std::endl;

			out<<std::setw(9)<<std::right <<i+offset_table*(i>0?1:0)<<"\t";
			out<<std::setw(11)<<std::right <<pointintime_list[i].first ;
			out<<std::setw(22)<<std::right <<convertTime(pointintime_list[i].first);

			for(uint j=0; j<pointintime_list[i].second.size();j++)
			{
				if(j==pointintime_list[i].second.size()-1)
					out <<"\t" <<std::setw(14)<<std::right <<pointintime_list[i].second[j];
				else
					out <<"\t" <<std::setw(14)<<std::right <<boost::multiprecision::float128(pointintime_list[i].second[j]);

			}
			out<<"\t"<<std::setw(14)<<std::right<<totalvec[i];

		}
		out<<std::endl<<std::endl;
	}

	void printCurrent(std::ostream& out, std::string msg= ""){
		std::ostringstream tmpstream (std::ostringstream::app);

		tmpstream<<"    " <<std::setw(20)<<std::left<<msg;

		for(uint i=0;i<currentList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<boost::multiprecision::float128(currentList[i]);
		}
		tmpstream<<std::endl;
		out<<tmpstream.str();
	}

	void printCurrent(std::string msg= ""){
		std::ostringstream tmpstream (std::ostringstream::app);

		tmpstream<<"    " <<std::setw(20)<<std::left<<msg;

		for(uint i=0;i<currentList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<currentList[i];
		}
		tmpstream<<std::endl;
		std::cout<<tmpstream.str();
	}

	void write(std::string filename, int histSize = std::numeric_limits<int>::infinity()){

		uint offset_table=0;
		if(histSize != std::numeric_limits<int>::infinity()&& currIt > histSize+1){
			offset_table=currIt - uint(histSize)-1;
		}

		std::ofstream out(filename,std::ofstream::out|std::ios::trunc);

		out <<std::setw(9)<<std::right<<"Iteration\t";
		out <<std::setw(11)<<std::right<<"Time[s]";
		out <<std::setw(22)<<std::right<<"Time";
		for(uint i=0;i<pointintime_list.size();i++)
		{
			if(i==0){
				for(uint j=0; j<pointintime_list[i].second.size();j++)
						{
					out <<"\t" <<std::setw(6)<<std::right <<"Facet-" <<std::setw(8)<<std::setfill('-')<<std::right <<j;
						}
			}
			out<<std::endl;

			out<<std::setw(9)<<std::setfill(' ')<<std::right <<i+offset_table*(i>0?1:0)<<"\t";
			out<<std::setw(11)<<std::right <<pointintime_list[i].first ;
			out<<std::setw(22)<<std::right <<convertTime(pointintime_list[i].first);

			for(uint j=0; j<pointintime_list[i].second.size();j++)
			{
				out <<"\t" <<std::setw(14)<<std::right <<pointintime_list[i].second[j];

			}

		}
		out.close();
		std::cout <<"Results saved to " <<filename <<std::endl;
	}
	void read(std::string filename, Databuff *hitbuffer, Simulation *sHandle){//Rudi: Not ready yet.
		pointintime_list.clear();
		std::string line;

		std::ifstream input(filename,std::ifstream::in);
		std::cout <<"Reading in covering history from " <<filename <<std::endl;
		while(std::getline(input,line)){
			std::vector<llong> currentstep;
			currentstep =std::vector<llong> ();

			llong covering;
			double time;
			std::istringstream is( line );

			is >> time;
			while(!is.eof()){
				is >> covering;
				currentstep.push_back(covering);

			}
			pointintime_list.push_back(std::make_pair(time,currentstep));
		}
		input.close();

		int i=0;
		double num_mol=0.0;
		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
					num_mol=pointintime_list.back().second[i]; //Rudi: Maybe wrong, since we changed covering and introduced coverage.
					f.tmpCounter[0].hit.covering = pointintime_list.back().second[i];
					//calcStickingnew(&f, hitbuffer); // calculate new sticking for new covering value
					// 1) Update the hitbuffer with the last covering value
					// 2) calcStickingnew(&f, hitbuffer);
					std::cout <<"Facet "<<i <<"\t covering: " <<f.tmpCounter[0].hit.covering <<"\t Corresponding number of particles: " <<num_mol <<std::endl;
					i+=1;
			}
		}
		std::cout <<std::endl;

	}

	bool empty(){return pointintime_list.empty();}

	void setCurrent(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	currentList[covidx]=newValue;}
	T getLast(int idx){return pointintime_list.back().second[idx];}
	T getLast(SubprocessFacet *iFacet){int covidx = getFacetIndex(iFacet);return pointintime_list.back().second[covidx];}
	T getCurrent(int idx){return currentList[idx];}
	T getCurrent(SubprocessFacet *iFacet){int covidx = getFacetIndex(iFacet);return currentList[covidx];}
	std::vector<T> getCurrent(){return pointintime_list.back().second;};
	void setLast(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	pointintime_list.back().second[covidx]=newValue;}

};

class ProblemDef{
public:
	ProblemDef(int argc, char *argv[]);
	ProblemDef();

	void createOutput(int save);
	void readArg(int argc, char *argv[], int rank=1);
	void readInputfile(std::string filename, int rank=1, int save=1);
	void writeInputfile(std::string filename, int rank=1);
	void printInputfile(std::ostream& out);

	bool saveResults;

	std::string resultpath;
	std::string resultbufferPath;
	std::ofstream outFile;

	// These can be given as parameters directly
	std::string loadbufferPath;
	std::string hitbufferPath;
	double simulationTime;
	std::string unit;

	// These can be given through input file only
	int iterationNumber; //number of iterations, all of length simulationTimeMS so far
	double maxTime;
	std::string maxUnit;

	double particleDia;

	double E_de;
	double sticking;

	double H_vap;
	double W_tr;

	int targetParticles;
	double targetError;

	double hitRatioLimit;
	double t_min; // [s]

	double t_max; // [s]
	int maxSimPerIt;

	llong coveringMinThresh;

	int histSize;

	double counterWindowPercent; // [%]

	std::vector< std::pair<int,double> > vipFacets;

	//These cannot be given, but are computed from other variables
	int simulationTimeMS;
	double maxTimeS;
};

class SimulationHistory{
public:
	SimulationHistory(int);
	SimulationHistory(Databuff *hitbuffer,int);

	HistoryList<boost::multiprecision::uint128_t> coveringList;
	HistoryList<double> hitList;
	HistoryList<llong> desorbedList;
	HistoryList<double> errorList_event;
	HistoryList<double> errorList_covering;

	std::vector<unsigned int> normalFacets;

	bool startNewParticle;
	llong smallCoveringFactor;

	unsigned int numFacet;
	int numSubProcess;

	double flightTime; // [s]
	int nParticles;

	double lastTime; // [s]
	int currentStep;
	double stepSize; // [s]

	void appendList(Databuff *hitbuffer, double time=-1.0);
	void appendList(double time=-1.0);

	void print(bool write=false);
	void write(std::string path);
	//std::tuple<bool, llong > updateHistory(Databuff *hitbuffer);
	void updateHistory();

};

class NullStream : public std::ostream {
    class NullBuffer : public std::streambuf {
    public:
        int overflow( int c ) { return c; }
    } m_nb;
public:
    NullStream() : std::ostream( &m_nb ) {}
};

//-----------------------------------------------------------
//SimulationLinux.cpp
//std::tuple<bool, std::vector<int> >  simulateSub(Databuff *hitbuffer, int rank, int simutime);
std::tuple<bool, std::vector<int> >  simulateSub2(Databuff *hitbuffer, int rank, int simutime);

double convertunit(double simutime, std::string unit);

void printConsole(std::string str,std::ofstream outFile);
bool checkSmallCovering(int rank, Databuff *hitbuffer_sum);
void UndoSmallCovering(Databuff *hitbuffer_sum);
//ProblemDef
//SimulationHistory

//-----------------------------------------------------------
//UpdateSubProcess.cpp

//void UpdateSticking(Databuff *hitbuffer);
void UpdateSticking();

//bool UpdateDesorptionRate (Databuff *hitbuffer);
bool UpdateDesorptionRate();

//void UpdateSojourn(Databuff *hitbuffer);
void UpdateSojourn();
std::tuple<double,double> UpdateErrorAll();
double UpdateError(std::string mode);
void UpdateErrorSub();
bool checkErrorSub(double targetError, double currentError, double factor, std::string mode="covering");


void UpdateMCSubHits(Databuff *databuffer, int rank);

void initbufftozero(Databuff *databuffer);

//-----------------------------------------------------------
//UpdateMainProcess.cpp

double manageStepSize();
double getStepSize();

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank);

void UpdateCovering(Databuff *hitbuffer_sum);
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer);

void UpdateErrorMain(Databuff *hitbuffer_sum);
std::tuple<std::vector<double>,std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration();

void printVelocities(Databuff *hitbuffer);

//-----------------------------------------------------------
//SimulationCalc.cpp

llong getnbDesorbed(Databuff *hitbuffer_sum);
llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer);
llong getnbAdsorbed(SubprocessFacet *iFacet, Databuff *hitbuffer);//In the original Molflow, particles were absorbed not adsorbed. In ContaminationFlow we regard all old code parts
// which are called 'Absorb' actually as an 'Adsorb'. But we did not rename them.
llong getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer);
boost::multiprecision::uint128_t getCovering(SubprocessFacet *iFacet);
double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer);
std::tuple<double, double, double> getVelocities(SubprocessFacet *iFacet, Databuff *hitbuffer_sum);

double calcStep(long double variable, double start, double end, double inflection_point, double Wtr);
//double calcEnergy(SubprocessFacet *iFacet, Databuff *hitbuffer);
double calcEnergy(SubprocessFacet *iFacet);

boost::multiprecision::float128 calcCoverage(SubprocessFacet *iFacet);

boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum);
//void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer);
void calcStickingnew(SubprocessFacet *iFacet);
//boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer);
//boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet);
boost::multiprecision::float128 calcDesorption(SubprocessFacet *iFacet);

double calcStartTime(SubprocessFacet *iFacet);

//-----------------------------------------------------------
//Iteration.cpp
/*
double estimateTmin(Databuff *hitbuffer);
*/
double estimateAverageFlightTime();

//void allocateCovering(Databuff *hitbuffer, int size, int rank);
void setCoveringThreshold(Databuff *hitbuffer, int size, int rank);
void setCoveringThreshold(int size, int rank);
void initCoveringThresh();

//-----------------------------------------------------------
//worker.cpp
void CalcTotalOutgassingWorker();
