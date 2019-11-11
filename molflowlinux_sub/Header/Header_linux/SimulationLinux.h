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
const double tuneE=2.7;

template <typename T> class HistoryList{
public:
	HistoryList(){
		pointintime_list = std::vector< std::pair<double,std::vector<T>> >();
		currentList=std::vector<T>();
	}

	std::vector< std::pair<double,std::vector<T>> > pointintime_list;
	std::vector<T> currentList;

	void reset(unsigned int numFacet=0){
		pointintime_list.clear();
		currentList.clear();
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

		}

	void appendList(std::vector<T> List, double time=-1){

		if(time==-1.0) //One step
				time=pointintime_list.back().first+1.0;
		pointintime_list.push_back(std::make_pair(time,List));

	}
	std::string convertTime(double time){
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

			if(years!=0) {final=final+std::to_string(years)+"y";}else{final=final+"   ";}//3
			if(months!=0) {final=(months>9)?final+std::to_string(months)+"mo":final+" "+std::to_string(months)+"mo";}else{final=final+"    ";}//4
			if(days!=0) {final=(days>9)?final+std::to_string(days)+"d":final+" "+std::to_string(days)+"d";}else{final=final+"   ";}//3
			if(hours!=0) {final=(hours>9)?final+std::to_string(hours)+"h":final+" "+std::to_string(hours)+"h";}else{final=final+"   ";}//3
			if(minutes!=0) {final=(minutes>9)?final+std::to_string(minutes)+"min":final+" "+std::to_string(minutes)+"min";}else{final=final+"     ";}//5
			if(seconds!=0) {final=(seconds>9)?final+std::to_string(seconds)+"s":final+" "+std::to_string(seconds)+"s";}else{final=final+"   ";}//3
		}
		if(final==""){
			final=std::to_string((int)(time+0.5))+"s";
		}

		return final;
	}
	void print(std::ostream& out, std::string msg= ""){
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

			out<<std::setw(9)<<std::right <<i<<"\t";
			out<<std::setw(11)<<std::right <<pointintime_list[i].first ;
			out<<std::setw(22)<<std::right <<convertTime(pointintime_list[i].first);

			for(uint j=0; j<pointintime_list[i].second.size();j++)
			{
				out <<"\t" <<std::setw(14)<<std::right <<boost::multiprecision::float128(pointintime_list[i].second[j]);

			}

		}
		out<<std::endl<<std::endl;
	}

	void print(std::ostream& out, std::vector<T> totalvec, std::string msg= ""){
			out<<std::endl <<msg <<std::endl;

			out <<std::setw(9)<<std::right<<"Iteration\t";
			out <<std::setw(11)<<std::right<<"Time[s]";
			out <<std::setw(22)<<std::right<<"Time";
			for(uint i=0;i<pointintime_list.size();i++)
			{
				if(i==0){
					for(uint j=0; j<pointintime_list[i].second.size();j++){
						out <<"\t" <<std::setw(6)<<std::right <<"Facet " <<std::setw(8)<<std::right <<j;
						}

					out <<"\t" <<std::setw(14)<<std::right<<"Total";
				}
				out<<std::endl;

				out<<std::setw(9)<<std::right <<i<<"\t";
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

		tmpstream<<"    " <<std::setw(12)<<std::left<<msg;

		for(uint i=0;i<currentList.size();i++)
		{
			tmpstream <<"\t" <<std::setw(12)<<std::right <<boost::multiprecision::float128(currentList[i]);
		}
		tmpstream<<std::endl;
		out<<tmpstream.str();
	}

	void write(std::string filename){

		std::ofstream out(filename,std::ofstream::out|std::ios::trunc);

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

			out<<std::setw(9)<<std::right <<i<<"\t";
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

	/*
	void newLine(unsigned int numFacet, double time=-1){
		if(time==-1){
			time=pointintime_list.back().first+1.0;
		}
		std::vector<T> tempvec; tempvec = std::vector<T>();
		for(unsigned int i=0; i<numFacet; i++){
			tempvec.push_back(0.0);
		}
		pointintime_list.push_back(std::make_pair(time,tempvec));

	}*/

	//void setCurrent(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	pointintime_list.back().second[covidx]=newValue;}
	void setCurrentList(SubprocessFacet *iFacet, T newValue){int covidx = getFacetIndex(iFacet);	currentList[covidx]=newValue;}
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

	//double s1;
	//double s2;
	double E_de;
	//double E_ad;
	//double d;
	double sticking;

	double H_vap;
	double W_tr;

	int targetParticles;
	double targetError;

	double hitRatioLimit;
	double Tmin;

	//TODO testing
	int coveringLimit;
	//double coveringMinFactor;
	//double coveringMaxFactor;
	llong coveringMinThresh;

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
	HistoryList<double> errorList;


	unsigned int numFacet;
	int numSubProcess;

	llong nbDesorbed_old;
	double flightTime;
	int nParticles;
	double lastTime;
	int currentStep;
	double stepSize;
	//double currentStepSizeFactor;

	//double StepSizeComputationTimeFactor;

	void appendList(Databuff *hitbuffer, double time=-1.0);
	void print(bool write=false);
	void write(std::string path);
	std::tuple<bool, llong > updateHistory(Databuff *hitbuffer);


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
std::tuple<bool, std::vector<int> >  simulateSub(Databuff *hitbuffer, int rank, int simutime);
double convertunit(double simutime, std::string unit);

void printConsole(std::string str,std::ofstream outFile);
//ProblemDef
//SimulationHistory

//-----------------------------------------------------------
//UpdateSubProcess.cpp

void UpdateSticking(Databuff *hitbuffer);
bool UpdateDesorptionRate (Databuff *hitbuffer);
void UpdateSojourn(Databuff *hitbuffer);
double UpdateError();
void UpdateErrorSub();


void UpdateMCSubHits(Databuff *databuffer, int rank);

void initbufftozero(Databuff *databuffer);

//-----------------------------------------------------------
//UpdateMainProcess.cpp

double manageStepSize(bool updateCurrentStep=false);

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, Databuff *physbuffer,int rank);
void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank);

void UpdateCovering(Databuff *hitbuffer, Databuff *hitbuffer_original, double time_step);
void UpdateCovering(Databuff *hitbuffer_sum);
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer);

void UpdateErrorMain(Databuff *hitbuffer_sum);
std::tuple<std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration();

//-----------------------------------------------------------
//SimulationCalc.cpp

llong getnbDesorbed(Databuff *hitbuffer_sum);
llong getnbDesorbed(SubprocessFacet *iFacet, Databuff *hitbuffer);
llong getCovering(SubprocessFacet *iFacet, Databuff *hitbuffer);
double getHits(SubprocessFacet *iFacet, Databuff *hitbuffer);

double calcStep(long double var, double start, double end, double step, double Wtr);
double calcEnergy(SubprocessFacet *iFacet, Databuff *hitbuffer);

boost::multiprecision::float128 GetMoleculesPerTP(Databuff *hitbuffer_sum, llong nbDesorbed_old);
void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer);
boost::multiprecision::float128 calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer);


//-----------------------------------------------------------
//Iteration.cpp
/*
double estimateTmin();
double estimateTmin_RudiTest(Databuff *hitbuffer);
*/
double estimateAverageFlightTime();

//void allocateCovering(Databuff *hitbuffer, int size, int rank);
void setCoveringThreshold(Databuff *hitbuffer, int size, int rank);
void initCoveringThresh();

//-----------------------------------------------------------
//worker.cpp
void CalcTotalOutgassingWorker();
