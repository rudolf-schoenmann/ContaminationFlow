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
 * Reduced Worker class for iterativve algorithm
 */
#pragma once

#include <string>
#include <vector>
#include "GLApp/GLTypes.h"
//#include "Smp.h"
#include "Buffer_shared.h" //LEAK, HIT

#include "Parameter.h"
#include "Vector.h" //moving parts
#include "MolflowTypes.h"

#define CDF_SIZE 100 //points in a cumulative distribution function


class Worker
{

public:

  // Constructor
  Worker();
  ~Worker();

  // Return a handle to the currently loaded geometry

  //std::vector<std::pair<double, double>> Generate_ID(int paramId);
  //int GenerateNewID(int paramId);
  std::vector<std::pair<double, double>> Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size);
  int GenerateNewCDF(double temperature);
  void CalcTotalOutgassing();
  int GetCDFId(double temperature);
  //int GetIDId(int paramId);

  std::vector<Parameter> parameters;
  //int displayedMoment;

  std::vector<std::vector<std::pair<double, double>>> CDFs; //cumulative distribution function for each temperature
  //std::vector<std::vector<std::pair<double, double>>> IDs; //integrated distribution function for each time-dependent desorption type
  std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
  //std::vector<double> moments;             //moments when a time-dependent simulation state is recorded
  //std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
  //std::vector<std::string> userMoments;    //user-defined text values for defining time moments (can be time or time series)


};
