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

static const char *day[]={"day","Day","days","Days","d","D"};
static const char *hour[]={"hour","Hour","hours","Hours","h","H","hr","Hr"};
static const char *min[]={"Minutes","minutes","Minute","minute","min","Min","m","M"};

const double carbondiameter = 2 *76E-12;
const double kb = 1.38E-23;
const double tau = 1E-13;

class CoveringHistory{
public:
	CoveringHistory();
	CoveringHistory(Databuff *hitbuffer);
	std::vector< std::pair<double,std::vector<llong>> > pointintime_list;
	//std::vector< std::pair<double,std::vector<double>> > pointintime_list_read;

	void appendList(Databuff *hitbuffer, double time=-1.0);
	void appendList(std::vector<llong> List);
	void print();
	void write(std::string filename);
	void read(std::string filename, Databuff *hitbuffer);

	llong getCurrentCovering(int idx);
	llong getCurrentCovering(SubprocessFacet *iFacet);
	void setCurrentCovering(SubprocessFacet *iFacet, llong newCovering);
};

class ProblemDef{
public:
	ProblemDef(int argc, char *argv[]);
	void readInputfile(std::string filename);
	void writeInputfile(std::string filename);
	double simulationTime;
	int simulationTimeMS;
	std::string unit;

	int iterationNumber;

	std::string loadbufferPath;
	std::string hitbufferPath;
	std::string resultbufferPath;

};


bool simulateSub(Databuff *hitbuffer, int rank, int simutime);
double convertunit(double simutime, std::string unit);

void initbufftozero(Databuff *databuffer);
void initcounterstozero(Databuff *databuffer);

void UpdateMCSubHits(Databuff *databuffer, int rank);
void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, Databuff *physbuffer,int rank);
void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, CoveringHistory *history,int rank);

int getFacetIndex(SubprocessFacet *iFacet);
llong getnbDesorbed(Databuff *hitbuffer_sum);
void CalcTotalOutgassingWorker();

llong calcCovering(SubprocessFacet *iFacet, Databuff *hitbuffer);
void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer);
double calcDesorptionRate(SubprocessFacet *iFacet, Databuff *hitbuffer);

double estimateTmin();
double estimateTmin_RudiTest(Databuff *hitbuffer);

void UpdateSticking(Databuff *hitbuffer);
void UpdateDesorptionRate (Databuff *hitbuffer);

void UpdateCovering(Databuff *hitbuffer, Databuff *hitbuffer_original, double time_step);
void UpdateCovering(CoveringHistory *history, Databuff *hitbuffer_sum, double time_step, llong *nbDesorbed_old);
void UpdateCoveringphys(CoveringHistory *history, Databuff *hitbuffer_sum, Databuff *hitbuffer);

void testprintcovering(Databuff *hitbuffer);
