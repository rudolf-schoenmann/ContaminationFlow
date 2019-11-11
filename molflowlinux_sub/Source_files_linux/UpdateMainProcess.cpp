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
 * This file contains the functions that adds up Databuff structs
 */

#include "SimulationLinux.h"
#include "GLApp/MathTools.h"
#include <math.h>

extern Simulation *sHandle;
extern ProblemDef* p;
extern SimulationHistory* simHistory;

// Step size for intervals
double getStepSize(){
	double T_min = p->Tmin;//set minimal time resolution to 1E-4 seconds.

	//Dynamical calculation of min_time is not straight forward, since 'manageTimeStep()' can change it.
	//Dynamical calculation can be done later, if it is regarded as useful.
	double t_start = T_min*exp((double)simHistory->currentStep*(log(p->maxTimeS/(T_min))/(double)p->iterationNumber));
	double t_stop = T_min*exp((double)(simHistory->currentStep+1)*(log(p->maxTimeS/(T_min))/(double)p->iterationNumber));
	double Delta = t_stop - t_start;
	return Delta;
	/*if(simHistory->currentStep==0){
		double T_min = estimateAverageFlightTime();
		return T_min*exp((double)simHistory->currentStep*(log(p->maxTimeS/T_min)/(double)p->iterationNumber));
		}
	else{
		return exp((double)simHistory->currentStep*(log(p->maxTimeS/simHistory->coveringList.pointintime_list[1].first)/(double)p->iterationNumber));
		}*/

}

double manageStepSize(bool updateCurrentStep){
	double step_size = getStepSize();
	bool incrCurrentStep=true;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(f.sh.desorption==0||f.sh.temperature==0)continue;

			int sizeList=simHistory->coveringList.pointintime_list.size();
			boost::multiprecision::uint128_t covering_phys = simHistory->coveringList.getLast(&f);
			boost::multiprecision::uint128_t covering_phys_before = simHistory->coveringList.pointintime_list[sizeList-2].second[getFacetIndex(&f)];

			if(updateCurrentStep){
				std::cout <<"Facet "<<getFacetIndex(&f) <<": " <<(f.sh.desorption/boost::multiprecision::float128(kb* f.sh.temperature))*boost::multiprecision::float128(step_size) +boost::multiprecision::float128(0.5) <<" >? " <<boost::multiprecision::float128(covering_phys) <<std::endl;
				p->outFile<<"Facet "<<getFacetIndex(&f) <<": " <<(f.sh.desorption/boost::multiprecision::float128(kb* f.sh.temperature))*boost::multiprecision::float128(step_size) +boost::multiprecision::float128(0.5) <<" >? " <<boost::multiprecision::float128(covering_phys) <<std::endl;
			}

			if ((boost::multiprecision::uint128_t)((f.sh.desorption/boost::multiprecision::float128(kb* f.sh.temperature))*boost::multiprecision::float128(step_size) +boost::multiprecision::float128(0.5))>covering_phys){

				boost::multiprecision::float128 test_size=boost::multiprecision::float128(covering_phys)/(f.sh.desorption/boost::multiprecision::float128(kb* f.sh.temperature));
				step_size=0.9 *test_size.convert_to<double>();

				incrCurrentStep=false;
			}
			/*
			if(covering_phys < covering_phys_before && sizeList>1) //covering_pyhs-covering_phys_before < 0
			{
				boost::multiprecision::float128 CovDiff=static_cast<boost::multiprecision::float128>(covering_phys_before-covering_phys);
				boost::multiprecision::float128 CurrTime= static_cast<boost::multiprecision::float128>(simHistory->coveringList.pointintime_list[sizeList-1].first-simHistory->coveringList.pointintime_list[sizeList-2].first);

				std::cout <<"Facet" <<getFacetIndex(&f)<<":\t CovDiff = Cov_before - Cov_now = "<<covering_phys_before <<" - " <<covering_phys<<" = "<<CovDiff <<"\t time_now-time_before = " <<CurrTime<<"\t If clause : "<<CovDiff*boost::multiprecision::float128(step_size)/CurrTime<<" >? " <<covering_phys<<std::endl;
				if(CovDiff*boost::multiprecision::float128(step_size)/CurrTime > static_cast<boost::multiprecision::float128>(covering_phys)){
					//boost::multiprecision::float128 test_size = static_cast<boost::multiprecision::float128>(covering_phys) * CurrTime/CovDiff;
					//step_size= 0.9 * test_size.convert_to<double>();
					//incrCurrentStep=false;
				}

			}*/

		}
	}

	if(incrCurrentStep&&updateCurrentStep){//needed here?  => JEIN
		simHistory->currentStep+=1;
		std::cout<<"Increase simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
		p->outFile<<"Increase simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
	}

	return step_size;
}

// Function that adapts timestep if needed, to avoid negative covering
//This function is not used anymore
double manageTimeStep(Databuff *hitbuffer_sum, double Krealvirt){
	//double test_time_step = pow(10,-14);
	double test_time_step = estimateAverageFlightTime();
	double stepSize=getStepSize();
	bool incrCurrentStep=false;

	boost::multiprecision::uint128_t covering_phys;
	boost::multiprecision::uint128_t covering_sum;
	boost::multiprecision::float128 covering_check;

	//double	covering_diff, covering_diff_min = 0;
	//llong facet_count_1 = 0;
	//llong facet_count_2 = 0;
	//bool decreased_time_step, increased_test_step;
	//decreased_time_step = increased_test_step = false;
	//double minimum_time_step_increase = 0;
	
	/*
	//Man könnte sich auch den 'increase time step check' pro Facet sparen, indem man den time step sofort erhöht, wenn man sieht, dass ein virtuelles Testteilchen weniger als einem realen Teilchen entspricht:
	if (test_time_step*Krealvirt < 1){
		test_time_step = test_time_step * (1/(test_time_step*Krealvirt)) *20; // 20 is just a chosen random value to increase the output per iteration
		std::cout << "increased time step to " << test_time_step << " s." << std::endl;
		p->outFile << "increased time step to " << test_time_step << " s." << std::endl;
	}*/
	//Damit hat man aber den time step stärker nach unten hin begrenzt, als mit dem 'increase time step check' pro Facet. Mit dem 'increase time step check' pro Facet kann es auch ein Verhältnis von realen zu Testteilchen
	//kleiner eins geben. Damit hat man eine besserer Statistik.
	
	/*
	//Check, if the time step is too short for more than 2 facets. We avoid changing the covering counter by a value of +- 0 (for more than N-2 facets) by increasing the time step.
	if (test_time_step*Krealvirt < 1){
		for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					covering_phys = history->coveringList.getCurrent(&f);
					covering_sum = getCovering(&f, hitbuffer_sum);

					if (covering_sum < covering_phys){
						covering_diff = (covering_phys - covering_sum)*Krealvirt*test_time_step;
						if (covering_diff_min == 0 || covering_diff_min > covering_diff){
							covering_diff_min = covering_diff;
						}
					}
					if (covering_sum > covering_phys){
						covering_diff = (covering_sum - covering_phys)*Krealvirt*test_time_step;
						if (covering_diff_min == 0 || covering_diff_min > covering_diff){
							covering_diff_min = covering_diff;
						}
					}
					if (covering_sum == covering_phys){
						facet_count_2 += 1;
					}					
					if ((llong)covering_diff == 0){
							facet_count_1 += 1;
					}
			}
		}
		if (sHandle->sh.nbFacet - facet_count_1 < 2){
			if (covering_diff_min == 0){
				if (sHandle->sh.nbFacet == facet_count_2){
					std::cout << "Covering does not change in this iteration step!";
				}
				else{			
					std::cout << "Error: Overflow of double value 'covering_diff'!" << std::endl;
				}
			}
			else{
				minimum_time_step_increase = (1/covering_diff_min) * test_time_step;
				test_time_step = 1000 * (1/covering_diff_min) * test_time_step; // 1000 is just a chosen random value to increase the output per iteration
				increased_test_step = true;
				std::cout << "Increased time step."<< std::endl;
				std::cout << "Time step is now "<<test_time_step<< " s."<< std::endl;
			}
		}
	}
	*/

	// Select larger value between test_time_step and stepSize
	if(stepSize>test_time_step){
		test_time_step=stepSize;
		std::cout<<"Replace test_time_step with stepSize: "<<stepSize <<std::endl;
		p->outFile<<"Replace test_time_step with stepSize: "<<stepSize <<std::endl;
		incrCurrentStep=true;
	}
	else{
		simHistory->currentStep+=1;
	}

	//Check, if the time step is too long. We avoid the covering counter going negative (=overflow of llong) by decreasing the time step.
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
					covering_phys = simHistory->coveringList.getLast(&f);
					covering_sum = static_cast<boost::multiprecision::uint128_t>(getCovering(&f, hitbuffer_sum));

					if (covering_sum < covering_phys){

						covering_check = boost::multiprecision::float128(covering_phys) - boost::multiprecision::float128(covering_phys - covering_sum)*boost::multiprecision::float128(Krealvirt*test_time_step);

						if(covering_check<0.0){
							boost::multiprecision::float128 tmp_step=boost::multiprecision::float128(covering_phys)/(boost::multiprecision::float128(covering_phys - covering_sum)*boost::multiprecision::float128(Krealvirt));
							test_time_step=0.2*tmp_step.convert_to<double>();//0.05 (decreasing multiplier) is a trade off between fast convergence and small oscillations
							//decreased_time_step = true;
							std::cout<<"Decreased Tmin: "<<test_time_step <<std::endl;
							p->outFile<<"Decreased Tmin: "<<test_time_step <<std::endl;
							incrCurrentStep=false;
							//std::cout << covering_check << std::endl;
						}	
					}
					
			}
	}
	//increment currentStep if stepSize chosen and not decreased
	if(incrCurrentStep){
		simHistory->currentStep+=1;
		std::cout<<"Increase simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
		p->outFile<<"Increase simHistory->currentStep: "<<simHistory->currentStep <<std::endl;
	}
	
	/*
	if (increased_test_step && decreased_time_step){
		if (test_time_step < minimum_time_step_increase){
			std::cout << "Note: Some particles might be lost due to choice of time step." << std::endl;
		}
	}*/
	
	return test_time_step;

}
//-----------------------------------------------------------
//Update Covering
/*
//Buffer version
void UpdateCovering(Databuff *hitbuffer_phys, Databuff *hitbuffer_sum, double time_step){//Updates Covering after one Iteration using Krealvirt, resets other counters
	//If one wants to read out pressure and particle density, this must be done before calling UpdateCovering.
	//Calculates with the summed up counters of hitbuffer_sum how many test particles are equivalent to one physical particle.
	//Then the physical values are stored in the hitbuffer.
	double Krealvirt = GetMoleculesPerTP(hitbuffer_sum,0);
	std::cout <<"Krealvirt = " << Krealvirt << std::endl;
	llong covering_phys;
	llong covering_sum;
	double covering_check;
	BYTE *buffer_phys;
	buffer_phys = hitbuffer_phys->buff;
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;
	double test_time_step = pow(10,-14);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer_phys = (FacetHitBuffer *)(buffer_phys + f.sh.hitOffset);
				covering_phys = facetHitBuffer_phys->hit.covering;
				FacetHitBuffer *facetHitBuffer_sum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_sum = facetHitBuffer_sum->hit.covering;
				std::cout<<std::endl << "Facet " << &f << std::endl;
				std::cout << "covering_sum = " << covering_sum << std::endl;
				std::cout<< "covering_phys_before = " << covering_phys << std::endl;
				if (covering_sum > covering_phys){
					llong covering_delta = static_cast < llong > ((covering_sum - covering_phys)*Krealvirt*time_step);
					covering_phys += covering_delta;
					std::cout << "covering rises"<< std::endl;
				}
				else{
					covering_check = covering_phys + (covering_phys - covering_sum)*Krealvirt*(-1)*time_step;
					std::cout <<"covering_check = " << covering_check << std::endl;
					if(!(covering_check<0)){
						llong covering_delta = static_cast < llong > ((covering_phys - covering_sum)*Krealvirt*time_step);
						covering_phys -= covering_delta;
						std::cout << "covering decreases but remains positive" << std::endl;
					}
					else {
						std::cout<<"Upps! Covering darf nicht negativ sein. Iteration wird nicht upgedated."<<std::endl;
						std::cout<<"test: "<<(double)(covering_phys/(covering_phys - covering_sum)*Krealvirt) <<std::endl;
						//std::cout << covering_check << std::endl;
						//nichts updaten
						//iteration neu starten mit weniger nbSteps; Wie viel weniger? 1/10 der vorigen Anzahl?
					}
				}
				std::cout<< "covering_phys_after = " << covering_phys << std::endl;
				facetHitBuffer_phys->hit.covering = covering_phys;
				//Reset of Hitbuffer_phys for the next Iteration Step
				facetHitBuffer_phys->hit.nbAbsEquiv = 0;
				facetHitBuffer_phys->hit.nbDesorbed = 0;
				facetHitBuffer_phys->hit.nbMCHit = 0;
				facetHitBuffer_phys->hit.nbHitEquiv = 0;
				facetHitBuffer_phys->hit.sum_1_per_ort_velocity = 0;
				facetHitBuffer_phys->hit.sum_v_ort = 0;
				facetHitBuffer_phys->hit.sum_1_per_velocity = 0;
				//Reset of Hitbuffer_sum for the next Iteration Step
				facetHitBuffer_sum->hit.nbAbsEquiv = 0;
				facetHitBuffer_sum->hit.nbDesorbed = 0;
				facetHitBuffer_sum->hit.nbMCHit = 0;
				facetHitBuffer_sum->hit.nbHitEquiv = 0;
				facetHitBuffer_sum->hit.sum_1_per_ort_velocity = 0;
				facetHitBuffer_sum->hit.sum_v_ort = 0;
				facetHitBuffer_sum->hit.sum_1_per_velocity = 0;
		}
	}
	if(covering_check){ //TODO always true?
	//Reset GlobalHitBuffer
	GlobalHitBuffer *gHits_phys;
	gHits_phys = (GlobalHitBuffer *)buffer_phys;
	gHits_phys->globalHits.hit.nbMCHit = 0;
	gHits_phys->globalHits.hit.nbHitEquiv = 0;
	gHits_phys->globalHits.hit.nbAbsEquiv = 0;
	gHits_phys->globalHits.hit.nbDesorbed = 0;
	GlobalHitBuffer *gHits_sum;
	gHits_sum = (GlobalHitBuffer *)buffer_sum;
	gHits_sum->globalHits.hit.nbMCHit = 0;
	gHits_sum->globalHits.hit.nbHitEquiv = 0;
	gHits_sum->globalHits.hit.nbAbsEquiv = 0;
	gHits_sum->globalHits.hit.nbDesorbed = 0;
	}
}
*/
//Simhistory version
void UpdateCovering(Databuff *hitbuffer_sum){//Updates Covering after one Iteration using Krealvirt, resets other counters
	//If one wants to read out pressure and particle density, this must be done before calling UpdateCovering.
	//Calculates with the summed up counters of hitbuffer_sum how many test particles are equivalent to one physical particle.
	//Then the physical values are stored in the hitbuffer.
	//simTime in ms

	boost::multiprecision::float128 Krealvirt = GetMoleculesPerTP(hitbuffer_sum, simHistory->nbDesorbed_old);
	llong nbDesorbed = getnbDesorbed(hitbuffer_sum)-simHistory->nbDesorbed_old;
	//std::cout <<"nbDesorbed before and after:\t" << history->nbDesorbed_old <<'\t';
	simHistory->nbDesorbed_old = getnbDesorbed(hitbuffer_sum);
	//std::cout << history->nbDesorbed_old <<std::endl;

	boost::multiprecision::uint128_t covering_phys;
	boost::multiprecision::uint128_t covering_sum;
	boost::multiprecision::float128 covering_check;


	//std::cout << "Tmin = " << estimateAverageFlightTime() << " s."<< std::endl;
	//p->outFile << "Tmin = " << estimateAverageFlightTime() << " s."<< std::endl;

	double time_step;
	if(Krealvirt==boost::multiprecision::float128(0.0)){ //if no Krealvirt(no desorption), do not increase currentStep, time_step=0
		time_step=0;
	}
	else{
		double error=0.0;
		double area=0.0;

		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				if(simHistory->errorList.getCurrent(&f)== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
					continue;

				error+=simHistory->errorList.getCurrent(&f)*f.sh.area;
				area+=f.sh.area;
			}
		}
		// Print total error and error per facet of this iteration
		std::cout <<"Total Error "<<error/area <<std::endl;
		simHistory->errorList.printCurrent(std::cout);

		//if targetError not reached: do not update currentstep
		if(error/area < /*1.05**/p->targetError)
			{time_step = manageStepSize(true);}
		else
			{time_step = manageStepSize(false);}
	}

	std::cout <<"Krealvirt = " << Krealvirt << std::endl;
	std::cout << "Covering difference will be multiplied by Krealvirt*(time step): " << Krealvirt*boost::multiprecision::float128(time_step) << std::endl;

	p->outFile <<"Krealvirt = " << Krealvirt << std::endl;
	p->outFile << "Covering difference will be multiplied by Krealvirt*(time step): " << Krealvirt*boost::multiprecision::float128(time_step) << std::endl;
	//std::cout <<"testing timestep: " <<time_step <<'\t' <<estimateAverageFlightTime() <<std::endl;

	//double rounding=1/simHistory->numFacet;
	boost::multiprecision::float128 rounding(0.0);
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {

				covering_phys = simHistory->coveringList.getLast(&f);
				covering_sum = boost::multiprecision::uint128_t(getCovering(&f, hitbuffer_sum));

				std::cout<<std::endl << "Facet " << getFacetIndex(&f)<< std::endl;
				std::cout << "covering_sum = " << covering_sum  << " = "<< boost::multiprecision::float128(covering_sum) << std::endl;
				std::cout<< "covering_phys_before = " << covering_phys << " = "<< boost::multiprecision::float128(covering_phys) << std::endl;

				p->outFile<<std::endl << "Facet " << getFacetIndex(&f) << std::endl;
				p->outFile<< "covering_sum = " << covering_sum  << " = "<< boost::multiprecision::float128(covering_sum) << std::endl;
				p->outFile<< "covering_phys_before = " << covering_phys << " = "<< boost::multiprecision::float128(covering_phys) << std::endl;

				if (covering_sum > covering_phys){
					//+0.5 for rounding
					boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding + boost::multiprecision::float128(covering_sum - covering_phys)*Krealvirt*boost::multiprecision::float128(time_step));
					covering_phys += covering_delta;
					std::cout << "covering rises by " <<covering_delta << " = "<<boost::multiprecision::float128(covering_delta) << std::endl;
					p->outFile << "covering rises by " <<covering_delta << " = "<<boost::multiprecision::float128(covering_delta) << std::endl;
				}
				else{
					//+0.5 for rounding
					boost::multiprecision::uint128_t covering_delta = static_cast < boost::multiprecision::uint128_t > (rounding + boost::multiprecision::float128(covering_phys - covering_sum)*Krealvirt*boost::multiprecision::float128(time_step));
					if(covering_phys+1==covering_delta){
						std::cout<<"!!!-----Correct covering_delta by 1: "<<covering_phys<<" - "<<covering_delta <<" because of rounding-----!!!"<<std::endl;
						p->outFile<<"!!!-----Correct covering_delta by 1: "<<covering_phys<<" - "<<covering_delta <<" because of rounding-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}
					else if(covering_phys<covering_delta && covering_phys+10>=covering_delta){
						std::cout<<"!!!-----Correct covering_delta: "<<covering_phys<<" - "<<covering_delta <<" because of bad statistic-----!!!"<<std::endl;
						p->outFile<<"!!!-----Correct covering_delta: "<<covering_phys<<" - "<<covering_delta <<" because of bad statistic-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}
					else if(covering_phys<covering_delta){
						std::cout<<"!!!-----Covering gets negative: "<<covering_phys<<" - "<<covering_delta <<". Correct Covering=0-----!!!"<<std::endl;
						p->outFile<<"!!!-----Covering gets negative: "<<covering_phys<<" - "<<covering_delta <<". Correct Covering=0-----!!!"<<std::endl;
						covering_delta=covering_phys;
					}

					covering_phys -= covering_delta;
					std::cout << "covering decreases by "<<covering_delta << " = " << boost::multiprecision::float128(covering_delta) << std::endl;
					p->outFile << "covering decreases by "<<covering_delta << " = " << boost::multiprecision::float128(covering_delta) << std::endl;

				}
				std::cout<< "covering_phys_after = " << covering_phys << " = " << boost::multiprecision::float128(covering_phys) << std::endl;
				std::cout<< "coveringThreshhold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << " = " << boost::multiprecision::float128(sHandle->coveringThreshold[getFacetIndex(&f)]) << std::endl;
				p->outFile<< "covering_phys_after = " << covering_phys << " = " << boost::multiprecision::float128(covering_phys) << std::endl;
				p->outFile<< "coveringThreshhold = " << sHandle->coveringThreshold[getFacetIndex(&f)] << " = " << boost::multiprecision::float128(sHandle->coveringThreshold[getFacetIndex(&f)]) << std::endl;

				simHistory->coveringList.setCurrentList(&f, covering_phys);
		}
	}
	simHistory->coveringList.appendCurrent(simHistory->lastTime+time_step);
	simHistory->lastTime+=time_step;
	simHistory->hitList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->errorList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering
	simHistory->desorbedList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateErrorMain before UpdateCovering

	simHistory->stepSize=time_step;
}

void UpdateErrorMain(Databuff *hitbuffer_sum){

	double num_hit_it=0;

	//save current num total hits in currentList, add difference current-old to num_hit_it
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			num_hit_it+=f.sh.opacity * (getHits(&f,hitbuffer_sum)-simHistory->hitList.getLast(&f) + getnbDesorbed(&f, hitbuffer_sum) - simHistory->desorbedList.getLast(&f));
		}
	}
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			double num_hit_f=f.sh.opacity * ( getHits(&f,hitbuffer_sum)-simHistory->hitList.getLast(&f) + getnbDesorbed(&f, hitbuffer_sum) - simHistory->desorbedList.getLast(&f));

			if(num_hit_f/num_hit_it<p->hitRatioLimit){// threshold. If reached, small number of hits neglected
				num_hit_it-=num_hit_f; //TODO also adapt facet counters??
				num_hit_f=0;
			}

			if(f.sh.opacity==0){simHistory->errorList.setCurrentList(&f, 0.0);}
			else{
				double error=pow((1/num_hit_f)*(1-num_hit_f/num_hit_it),0.5);
				simHistory->errorList.setCurrentList(&f, error);
			}
			simHistory->hitList.setLast(&f,getHits(&f,hitbuffer_sum));
			simHistory->desorbedList.setLast(&f,getnbDesorbed(&f,hitbuffer_sum));
		}
	}

	simHistory->errorList.appendCurrent(simHistory->lastTime);
	//simHistory->hitList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateCovering before UpdateErrorMain
	//simHistory->desorbedList.pointintime_list.back().first=simHistory->lastTime; // Uncomment if UpdateCovering before UpdateErrorMain

}

std::tuple<std::vector<double>,std::vector<boost::multiprecision::uint128_t>>  CalcPerIteration(){
	std::vector<double> errorPerIt;
	errorPerIt =std::vector<double> ();

	std::vector<boost::multiprecision::uint128_t> covPerIt;
	covPerIt =std::vector<boost::multiprecision::uint128_t> ();

	for(int it=0; it<simHistory->errorList.pointintime_list.size();it++){
		// Total error/covering for each iteration
		double error=0.0;
		double area=0.0;
		boost::multiprecision::uint128_t covering=0;

		for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				int idx=getFacetIndex(&f);
				covering+=simHistory->coveringList.pointintime_list[it].second[idx];

				double err=simHistory->errorList.pointintime_list[it].second[idx];
				if(err== std::numeric_limits<double>::infinity()||f.sh.opacity==0)//ignore facet if no hits (=inf error)
					continue;

				error+=err*f.sh.area;
				area+=f.sh.area;
			}
		}
		errorPerIt.push_back(error/area);
		covPerIt.push_back(covering);
	}
	return std::make_tuple(errorPerIt,covPerIt);
}

// Copy covering to buffer
void UpdateCoveringphys(Databuff *hitbuffer_sum, Databuff *hitbuffer){
	boost::multiprecision::uint128_t covering_phys;
	BYTE *buffer_sum;
	buffer_sum = hitbuffer_sum->buff;

	BYTE *buffer;
	buffer = hitbuffer->buff;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				FacetHitBuffer *facetHitSum = (FacetHitBuffer *)(buffer_sum + f.sh.hitOffset);
				covering_phys = simHistory->coveringList.getLast(&f);
				facetHitBuffer->hit.covering=covering_phys.convert_to<llong>();
				facetHitSum->hit.covering=covering_phys.convert_to<llong>();
			}
	}

	simHistory->flightTime=0.0;
	simHistory->nParticles=0;
}

//-----------------------------------------------------------
// Update Buffer of main process using buffer of sub process
// Physbuffer version
void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, Databuff *physbuffer,int rank) {
	BYTE *buffer, *subbuff, *buffer_phys;
	GlobalHitBuffer *gHits, *subHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif


	buffer = mainbuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	//added subbuffer that contains simulation results from a subprocess, to be added to mainbuffer
	subbuff=subbuffer->buff;
	subHits=(GlobalHitBuffer *)subbuff;

	//buffer_phys holds the physical values of covering before the iteration step
	buffer_phys = physbuffer->buff;

	size_t nbMoments=(size_t)sHandle->moments.size();
/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl <<std::endl;

	std::cout <<subHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<subHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<subHits->globalHits.hit.covering <<std::endl;
	std::cout <<subHits->distTraveled_total  <<std::endl;
	std::cout <<subHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;*/

	// Global hits and leaks: adding local hits to shared memory
	sHandle->tmpGlobalResult.globalHits.hit.nbMCHit=gHits->globalHits.hit.nbMCHit += subHits->globalHits.hit.nbMCHit;
	sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv=gHits->globalHits.hit.nbHitEquiv += subHits->globalHits.hit.nbHitEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv=gHits->globalHits.hit.nbAbsEquiv += subHits->globalHits.hit.nbAbsEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed=gHits->globalHits.hit.nbDesorbed += subHits->globalHits.hit.nbDesorbed;
	sHandle->tmpGlobalResult.distTraveled_total=gHits->distTraveled_total += subHits->distTraveled_total;
	sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly=gHits->distTraveledTotal_fullHitsOnly += subHits->distTraveledTotal_fullHitsOnly;

	/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;
	std::cout <<gHits->hitCacheSize <<std::endl;*/

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	//sHandle->wp.sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < subHits->leakCacheSize; leakIndex++)
		sHandle->tmpGlobalResult.leakCache[leakIndex]=gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = subHits->leakCache[leakIndex];
	sHandle->tmpGlobalResult.nbLeakTotal=gHits->nbLeakTotal += subHits->nbLeakTotal;
	sHandle->tmpGlobalResult.lastLeakIndex=gHits->lastLeakIndex = (gHits->lastLeakIndex + subHits->leakCacheSize) % LEAKCACHESIZE;
	sHandle->tmpGlobalResult.leakCacheSize=gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + subHits->leakCacheSize);


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	if (rank == 0) {
		for (size_t hitIndex = 0; hitIndex < subHits->hitCacheSize; hitIndex++)
			sHandle->tmpGlobalResult.hitCache[hitIndex]=gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = subHits->hitCache[hitIndex];

		if (subHits->hitCacheSize > 0) {
			sHandle->tmpGlobalResult.lastHitIndex = gHits->lastHitIndex = (gHits->lastHitIndex + subHits->hitCacheSize) % HITCACHESIZE;
			sHandle->tmpGlobalResult.hitCache[gHits->lastHitIndex].type = gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			sHandle->tmpGlobalResult.hitCacheSize = gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + subHits->hitCacheSize);
		}
	}
	/*
	std::cout <<gHits->hitCacheSize <<std::endl;
	std::cout <<gHits->leakCacheSize <<std::endl;
	std::cout <<gHits->nbLeakTotal <<std::endl<<std::endl;*/


	//Global histograms

		for (unsigned int m = 0; m < (1 + nbMoments); m++) {
			BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
			BYTE *subhist = subbuff + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
			if (sHandle->wp.globalHistogramParams.recordBounce) {
				double* nbHitsHistogram = (double*)histCurrentMoment;
				double* nbHitsSub=(double*)subhist;
				for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
					sHandle->tmpGlobalHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
				}
			}

			if (sHandle->wp.globalHistogramParams.recordDistance) {
				double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
				double* distanceSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					sHandle->tmpGlobalHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
				}
			}
			if (sHandle->wp.globalHistogramParams.recordTime) {
				double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
				double* timeSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					sHandle->tmpGlobalHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
				}
			}

		}


	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	//int num=0;
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {

		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//if (f.hitted) {
				//std::cout <<"Facet " <<num <<std::endl; // Da wird immer "0" angezeigt. Was ist der Sinn?
				for (unsigned int m = 0; m < (1 + nbMoments); m++) { // Add hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitSub = (FacetHitBuffer *)(subbuff + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitphys = (FacetHitBuffer *)(buffer_phys +f.sh.hitOffset + m * sizeof(FacetHitBuffer));
/*
					std::cout <<sizeof(GlobalHitBuffer) <<std::endl;
					std::cout <<f.sh.hitOffset  <<std::endl;
					std::cout <<sizeof(FacetHitBuffer) <<std::endl;
					std::cout <<facetHitSub->hit.covering <<std::endl;*/

					/*if(f.globalId == 1){
					std::cout <<"facet number\t\t\t"<< f.globalId <<std::endl;
					std::cout <<"buffer before\t\t\t" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/

					f.tmpCounter[m].hit.nbAbsEquiv = facetHitBuffer->hit.nbAbsEquiv += facetHitSub->hit.nbAbsEquiv;
					f.tmpCounter[m].hit.nbDesorbed = facetHitBuffer->hit.nbDesorbed += facetHitSub->hit.nbDesorbed;
					f.tmpCounter[m].hit.nbMCHit = facetHitBuffer->hit.nbMCHit += facetHitSub->hit.nbMCHit;
					f.tmpCounter[m].hit.nbHitEquiv = facetHitBuffer->hit.nbHitEquiv += facetHitSub->hit.nbHitEquiv;;
					f.tmpCounter[m].hit.sum_1_per_ort_velocity = facetHitBuffer->hit.sum_1_per_ort_velocity += facetHitSub->hit.sum_1_per_ort_velocity;
					f.tmpCounter[m].hit.sum_v_ort = facetHitBuffer->hit.sum_v_ort += facetHitSub->hit.sum_v_ort;
					f.tmpCounter[m].hit.sum_1_per_velocity = facetHitBuffer->hit.sum_1_per_velocity += facetHitSub->hit.sum_1_per_velocity;

					//facetHitBuffer->hit.covering += facetHitSub->hit.covering; //We do that in another way.
					if (facetHitSub->hit.covering > facetHitphys->hit.covering){
					facetHitBuffer->hit.covering += (facetHitSub->hit.covering - facetHitphys->hit.covering);
					}
					else{
						if(facetHitBuffer->hit.covering > (facetHitphys->hit.covering - facetHitSub->hit.covering)){
							facetHitBuffer->hit.covering -= (facetHitphys->hit.covering - facetHitSub->hit.covering);
							}
						else facetHitBuffer->hit.covering = 0;//Counter cannot be negative! Maybe we could interrupt the iteration here?
					}
					f.tmpCounter[m].hit.covering = facetHitBuffer->hit.covering;
					/*
					if(f.globalId == 1){
					std::cout <<"buffer afterwards" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/
				}

				if (f.sh.isProfile) { //(MY) Add profiles
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						ProfileSlice *shProfileSub = (ProfileSlice *)(subbuff + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							f.profile[m][j] = shProfile[j] += shProfileSub[j];
						}
					}
				}

				if (f.sh.isTextured) {// Add texture
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						TextureCell *shTextureSub = (TextureCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)
						//std::cout <<"timecorrection: " <<timeCorrection <<std::endl<<std::endl;

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Add temporary hit counts
								f.texture[m][add] = shTexture[add] += shTextureSub[add];

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f.textureCellIncrements[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f.largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}
							}
						}
					}
				}

				if (f.sh.countDirection) {
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						DirectionCell *shDirSub = (DirectionCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;
								f.direction[m][add].dir.x = shDir[add].dir.x += shDirSub[add].dir.x;
								f.direction[m][add].dir.y = shDir[add].dir.y += shDirSub[add].dir.y;
								f.direction[m][add].dir.z = shDir[add].dir.z += shDirSub[add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								f.direction[m][add].count = shDir[add].count += shDirSub[add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					size_t *shAngleMapSub = (size_t *)(subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							f.angleMap.pdf[add] = shAngleMap[add] += shAngleMapSub[add];
						}
					}
				}

				//Facet histograms

					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
						BYTE *histSub = subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
						if (f.sh.facetHistogramParams.recordBounce) {
							double* nbHitsHistogram = (double*)histCurrentMoment;
							double* nbHitsSub = (double*)histSub;
							for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
								f.tmpHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
							}
						}
						if (f.sh.facetHistogramParams.recordDistance) {
							double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
							double* distanceSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
								f.tmpHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
							}
						}
						if (f.sh.facetHistogramParams.recordTime) {
							double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							double* timeSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
								f.tmpHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
							}
						}
					}

			//} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
		sHandle->tmpGlobalResult.texture_limits[v]=gHits->texture_limits[v];
	}


	//ResetTmpCounters();
	//extern char* GetSimuStatus();
	//SetState(NULL, GetSimuStatus(), false, true); // (Rudi) Don't need that.

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}

//Simhistory version

void UpdateMCMainHits(Databuff *mainbuffer, Databuff *subbuffer, SimulationHistory *history,int rank) {
	BYTE *buffer, *subbuff;
	GlobalHitBuffer *gHits, *subHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif


	size_t nbMoments = (size_t)sHandle->moments.size();

	buffer = mainbuffer->buff;
	gHits = (GlobalHitBuffer *)buffer;

	//added subbuffer that contains simulation results from a subprocess, to be added to mainbuffer
	subbuff=subbuffer->buff;
	subHits=(GlobalHitBuffer *)subbuff;

/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl <<std::endl;

	std::cout <<subHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<subHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<subHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<subHits->globalHits.hit.covering <<std::endl;
	std::cout <<subHits->distTraveled_total  <<std::endl;
	std::cout <<subHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;*/

	// Global hits and leaks: adding local hits to shared memory
	sHandle->tmpGlobalResult.globalHits.hit.nbMCHit= gHits->globalHits.hit.nbMCHit += subHits->globalHits.hit.nbMCHit;
	sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv=gHits->globalHits.hit.nbHitEquiv += subHits->globalHits.hit.nbHitEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv=gHits->globalHits.hit.nbAbsEquiv += subHits->globalHits.hit.nbAbsEquiv;
	sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed=gHits->globalHits.hit.nbDesorbed += subHits->globalHits.hit.nbDesorbed;
	sHandle->tmpGlobalResult.distTraveled_total=gHits->distTraveled_total += subHits->distTraveled_total;
	sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly=gHits->distTraveledTotal_fullHitsOnly += subHits->distTraveledTotal_fullHitsOnly;

	/*
	std::cout <<gHits->globalHits.hit.nbMCHit  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbHitEquiv   <<std::endl;
	std::cout <<gHits->globalHits.hit.nbAbsEquiv  <<std::endl;
	std::cout <<gHits->globalHits.hit.nbDesorbed <<std::endl;
	std::cout <<gHits->globalHits.hit.covering <<std::endl;
	std::cout <<gHits->distTraveled_total  <<std::endl;
	std::cout <<gHits->distTraveledTotal_fullHitsOnly <<std::endl<<std::endl;
	std::cout <<gHits->hitCacheSize <<std::endl;*/

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	//sHandle->wp.sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < subHits->leakCacheSize; leakIndex++)
		sHandle->tmpGlobalResult.leakCache[leakIndex]=gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = subHits->leakCache[leakIndex];
	sHandle->tmpGlobalResult.nbLeakTotal=gHits->nbLeakTotal += subHits->nbLeakTotal;
	sHandle->tmpGlobalResult.lastLeakIndex=gHits->lastLeakIndex = (gHits->lastLeakIndex + subHits->leakCacheSize) % LEAKCACHESIZE;
	sHandle->tmpGlobalResult.leakCacheSize=gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + subHits->leakCacheSize);


	// HHit (Only prIdx 0) //Rudi: I think that's some Hit-History stuff. Not necessary to comment out (presumably).
	if (rank == 0) {
			for (size_t hitIndex = 0; hitIndex < subHits->hitCacheSize; hitIndex++)
				sHandle->tmpGlobalResult.hitCache[hitIndex]=gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = subHits->hitCache[hitIndex];

			if (subHits->hitCacheSize > 0) {
				sHandle->tmpGlobalResult.lastHitIndex = gHits->lastHitIndex = (gHits->lastHitIndex + subHits->hitCacheSize) % HITCACHESIZE;
				sHandle->tmpGlobalResult.hitCache[gHits->lastHitIndex].type = gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
				sHandle->tmpGlobalResult.hitCacheSize = gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + subHits->hitCacheSize);
			}
		}
	/*
	std::cout <<gHits->hitCacheSize <<std::endl;
	std::cout <<gHits->leakCacheSize <<std::endl;
	std::cout <<gHits->nbLeakTotal <<std::endl<<std::endl;*/


	//Global histograms
	for (unsigned int m = 0; m < (1 + nbMoments); m++) {
		BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
		BYTE *subhist = subbuff + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
		if (sHandle->wp.globalHistogramParams.recordBounce) {
			double* nbHitsHistogram = (double*)histCurrentMoment;
			double* nbHitsSub=(double*)subhist;
			for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
				sHandle->tmpGlobalHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
			}
		}

		if (sHandle->wp.globalHistogramParams.recordDistance) {
			double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			double* distanceSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
			for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
				sHandle->tmpGlobalHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
			}
		}
		if (sHandle->wp.globalHistogramParams.recordTime) {
			double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
			double* timeSub = (double*)(subhist + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
			for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
				sHandle->tmpGlobalHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
			}
		}

	}



	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	//int num=0;
	for (s = 0; s < (int)sHandle->sh.nbSuper; s++) {

		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			//if (f.hitted) {
				//std::cout <<"Facet " <<num <<std::endl; // Da wird immer "0" angezeigt. Was ist der Sinn?
				for (unsigned int m = 0; m < (1 + nbMoments); m++) { // Add hits
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					FacetHitBuffer *facetHitSub = (FacetHitBuffer *)(subbuff + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					llong covering_phys= history->coveringList.getLast(&f).convert_to<llong>();
					llong covering_sum = facetHitSub->hit.covering;
/*
					std::cout <<sizeof(GlobalHitBuffer) <<std::endl;
					std::cout <<f.sh.hitOffset  <<std::endl;
					std::cout <<sizeof(FacetHitBuffer) <<std::endl;
					std::cout <<facetHitSub->hit.covering <<std::endl;*/

					/*if(f.globalId == 1){
					std::cout <<"facet number\t\t\t"<< f.globalId <<std::endl;
					std::cout <<"buffer before\t\t\t" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/

					f.tmpCounter[m].hit.nbAbsEquiv = facetHitBuffer->hit.nbAbsEquiv += facetHitSub->hit.nbAbsEquiv;
					f.tmpCounter[m].hit.nbDesorbed = facetHitBuffer->hit.nbDesorbed += facetHitSub->hit.nbDesorbed;
					f.tmpCounter[m].hit.nbMCHit = facetHitBuffer->hit.nbMCHit += facetHitSub->hit.nbMCHit;
					f.tmpCounter[m].hit.nbHitEquiv = facetHitBuffer->hit.nbHitEquiv += facetHitSub->hit.nbHitEquiv;;
					f.tmpCounter[m].hit.sum_1_per_ort_velocity = facetHitBuffer->hit.sum_1_per_ort_velocity += facetHitSub->hit.sum_1_per_ort_velocity;
					f.tmpCounter[m].hit.sum_v_ort = facetHitBuffer->hit.sum_v_ort += facetHitSub->hit.sum_v_ort;
					f.tmpCounter[m].hit.sum_1_per_velocity = facetHitBuffer->hit.sum_1_per_velocity += facetHitSub->hit.sum_1_per_velocity;
					//facetHitBuffer->hit.covering += facetHitSub->hit.covering; //We do that in another way.

					if (facetHitSub->hit.covering > covering_phys){
					facetHitBuffer->hit.covering += (covering_sum - covering_phys);
					}
					else{
						if(facetHitBuffer->hit.covering > (covering_phys - covering_sum)){
							facetHitBuffer->hit.covering -= (covering_phys- covering_sum);
							}
						else{
							facetHitBuffer->hit.covering = 0;//Counter cannot be negative! Maybe we could interrupt the iteration here?
						}
					}
					f.tmpCounter[m].hit.covering = facetHitBuffer->hit.covering;

					/*
					if(f.globalId == 1){
					std::cout <<"buffer afterwards" <<std::endl;
					std::cout <<"hit.nbAbsEquiv\t\t\t"<<facetHitBuffer->hit.nbAbsEquiv <<"\t (as nbMCHitEquiv)"<<std::endl;
					std::cout <<"hit.nbDesorbed\t\t\t"<<facetHitBuffer->hit.nbDesorbed <<std::endl;
					std::cout <<"hit.nbMCHit\t\t\t"<<facetHitBuffer->hit.nbMCHit <<std::endl;
					std::cout <<"hit.nbHitEquiv\t\t\t"<<facetHitBuffer->hit.nbHitEquiv <<"\t gives the number of equivalent MC hits in low flux mode"<<std::endl;
					std::cout <<"hit.sum_1_per_ort_velocity\t"<<facetHitBuffer->hit.sum_1_per_ort_velocity <<std::endl;
					std::cout <<"hit.sum_v_ort\t\t\t"<<facetHitBuffer->hit.sum_v_ort <<std::endl;
					std::cout <<"hit.sum_1_per_velocity\t\t"<<facetHitBuffer->hit.sum_1_per_velocity <<std::endl;
					std::cout <<"hit.covering\t\t\t"<<facetHitBuffer->hit.covering <<std::endl;
					std::cout <<"hit.covering [Number of particles]\t" << facetHitBuffer->hit.covering * calcNmono(&f) << std::endl;
					}*/
				}

				if (f.sh.isProfile) { //(MY) Add profiles
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						ProfileSlice *shProfileSub = (ProfileSlice *)(subbuff + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < (int)PROFILE_SIZE; j++) {
							f.profile[m][j] = shProfile[j] += shProfileSub[j];
						}
					}
				}

				if (f.sh.isTextured) {// Add texture
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						TextureCell *shTextureSub = (TextureCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)
						//std::cout <<"timecorrection: " <<timeCorrection <<std::endl<<std::endl;

						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Add temporary hit counts
								f.texture[m][add] = shTexture[add] += shTextureSub[add];

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f.textureCellIncrements[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f.largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}
							}
						}
					}
				}

				if (f.sh.countDirection) {
					for (unsigned int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						DirectionCell *shDirSub = (DirectionCell *)(subbuff + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						for (y = 0; y < (int)f.sh.texHeight; y++) {
							for (x = 0; x < (int)f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;
								f.direction[m][add].dir.x = shDir[add].dir.x += shDirSub[add].dir.x;
								f.direction[m][add].dir.y = shDir[add].dir.y += shDirSub[add].dir.y;
								f.direction[m][add].dir.z = shDir[add].dir.z += shDirSub[add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								f.direction[m][add].count = shDir[add].count += shDirSub[add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					size_t *shAngleMapSub = (size_t *)(subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (int)(f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < (int)f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							f.angleMap.pdf[add] = shAngleMap[add] += shAngleMapSub[add];
						}
					}
				}

				//Facet histograms

				for (unsigned int m = 0; m < (1 + nbMoments); m++) {
					BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
					BYTE *histSub = subbuff + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
					if (f.sh.facetHistogramParams.recordBounce) {
						double* nbHitsHistogram = (double*)histCurrentMoment;
						double* nbHitsSub = (double*)histSub;
						for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
							f.tmpHistograms[m].nbHitsHistogram[i] = nbHitsHistogram[i] += nbHitsSub[i];
						}
					}
					if (f.sh.facetHistogramParams.recordDistance) {
						double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
						double* distanceSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize());
						for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
							f.tmpHistograms[m].distanceHistogram[i] = distanceHistogram[i] += distanceSub[i];
						}
					}
					if (f.sh.facetHistogramParams.recordTime) {
						double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
						double* timeSub = (double*)(histSub + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
						for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
							f.tmpHistograms[m].timeHistogram[i] = timeHistogram[i] += timeSub[i];
						}
					}
				}

			//} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
		sHandle->tmpGlobalResult.texture_limits[v]=gHits->texture_limits[v];
	}

	//ResetTmpCounters(); //This will reset counters -> but needed for Tmin

	//extern char* GetSimuStatus();
	//SetState(NULL, GetSimuStatus(), false, true); // (Rudi) Don't need that.

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}
