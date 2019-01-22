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

class CoveringHistory{
public:
	CoveringHistory();
	CoveringHistory(Databuff *hitbuffer);
	std::vector< std::pair<double,std::vector<double>> > pointintime_list;
	//std::vector< std::pair<double,std::vector<double>> > pointintime_list_read;

	void appendList(double time);
	void print();
	void write(std::string filename);
	void read(std::string filename, Databuff *hitbuffer);
};


bool simulateSub(Databuff *hitbuffer, int rank, int simutime);
double convertunit(double simutime, std::string unit);
void UpdateSubHits(Databuff *databuffer, int rank);
void UpdateSubMCHits(Databuff *databuffer, int rank, size_t nbMoments);
//void initbufftozero(size_t nbMoments, Databuff *buffer);

void UpdateMCmainHits(Databuff *mainbuffer, Databuff *subbuffer,int rank, size_t nbMoments);
void UpdateMainHits(Databuff *databuffer,Databuff *subbuffer, int rank);

void UpdateSticking(Databuff *hitbuffer);

void calcStickingnew(SubprocessFacet *iFacet, Databuff *hitbuffer);
//double calcDesorption(SubprocessFacet *iFacet);
double calcNmono(SubprocessFacet *iFacet);
double calcRealCovering(SubprocessFacet *iFacet);
double calcKrealvirt(SubprocessFacet *iFacet, int m);


double estimateTmin();
double estimateTmin_RudiTest();

