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
 * This file contains calculations for the iterative algorithm
 */

#include "SimulationContaminationFlow.h"
#include <math.h>
//#include "levmar.h"
extern SimulationHistory* simHistory;
extern Simulation *sHandle;
extern ProblemDef *p;


//Update Error lists
std::tuple<double,double,double,double> getErrorVariables(SubprocessFacet* f, Databuff *hitbuffer_sum){
	if(hitbuffer_sum==NULL){
		return std::make_tuple(f->tmpCounter[0].nbHitEquiv, (double)f->tmpCounter[0].nbDesorbed, (double)f->tmpCounter[0].nbOutgassed, f->tmpCounter[0].nbAbsEquiv);
	}
	else{
		return std::make_tuple(getHits(f,hitbuffer_sum),(double)getnbDesorbed(f, hitbuffer_sum),(double)getnbOutgassed(f, hitbuffer_sum),getnbAdsorbed(f,hitbuffer_sum));
	}
}

void UpdateErrorList(Databuff *hitbuffer_sum){ // hitbuffer_sum==NULL: subprocess error, else: main error
	double nbhits_all=0; double nbdes_all=0;
	double nbout_all=0; double nbads_all=0;

	double nbhits=0.0; double nbdes=0.0; double nbout=0.0;double nbads=0.0;

	//double factor=hitbuffer_sum==NULL?pow(simHistory->numSubProcess,0.5):1.0; // To be consistent with the ignored facets for calculating the error after summation over all subprocesses, here the hitRationLimit must be reduced with the correction factor due to multiple subprocesses.

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			std::tie(nbhits,nbdes,nbout,nbads)=getErrorVariables(&f, hitbuffer_sum);
			nbhits_all+=nbhits;
			nbdes_all += nbdes;
			nbout_all += nbout;
			nbads_all +=nbads;
		}
	}

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			std::tie(nbhits,nbdes,nbout,nbads)=getErrorVariables(&f, hitbuffer_sum);
			if(nbdes_all==0){
				nbdes_all = 1;// prevent error being 'nan'
			}
			if(nbads_all==0){
				nbads_all = 1;// prevent error being 'nan'
			}
			if(nbout_all==0){
				nbout_all = 1;// prevent error being 'nan'
			}
			if(nbhits_all==0){
				nbhits_all = 1;// prevent error being 'nan'
			}
			double error_event=(pow(((nbdes*(1-nbdes/nbdes_all))+(nbout*(1-nbout/nbout_all))+(nbhits*(1-nbhits/nbhits_all))),0.5)/(nbdes_all+nbout_all+nbhits_all));
			double error_covering=(pow(((nbdes*(1-nbdes/nbdes_all))+(nbads*(1-nbads/nbads_all))),0.5)/(nbdes_all+nbads_all));
			//if counters are zero, the error formula yields zero. This has to be excluded.
			if(error_event == 0){
				if((nbdes == nbdes_all && nbdes_all != 0) || (nbout == nbout_all && nbout_all != 0) || (nbhits == nbhits_all&& nbhits_all != 0)){
					error_event=0; //error would be zero due to formula. => will not be set to 1, since there is some kind of statistics
				}
				else{
					error_event=1;
				}
			}
			if(error_covering == 0){
				if((nbdes == nbdes_all && nbdes_all != 0) || (nbads == nbads_all && nbads_all != 0)){
					error_covering=0; //error would be zero due to formula. => will not be set to 1, since there is some kind of statistics
				}
				else{
					error_covering=1;
				}
			}
			//scale the error with the number of MC test-particles => error goes to zero with
			//MC test-particles -> inf
			simHistory->errorList_event.setCurrent(&f, error_event);
			simHistory->errorList_covering.setCurrent(&f, error_covering);

			//save current hits in currentList,
			if(hitbuffer_sum==NULL){ // Sub process: set currentList values
				simHistory->hitList.setCurrent(&f,nbhits);
				simHistory->desorbedList.setCurrent(&f,nbdes);
				simHistory->outgassedList.setCurrent(&f,nbout);
				simHistory->adsorbedList.setCurrent(&f,nbads);
			}
			else{ // Main process: set historyList values
				simHistory->hitList.setLast(&f,nbhits);
				simHistory->desorbedList.setLast(&f,nbdes);
				simHistory->outgassedList.setLast(&f,nbout);
				simHistory->adsorbedList.setLast(&f,nbads);
			}
		}
	}
}

// Calculate total error
std::tuple<double,double> CalcErrorAll(int it){//calculates the averaged total error weighted with the facets area to decide, if the desired uncertainty level is reached
	double error_event=0.0; double error_covering=0.0;
	//double area_event=0.0; double area_covering=0.0;

	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			int idx=getFacetIndex(&f);
			if(f.sh.opacity==0||(p->doFocusGroupOnly && std::find(std::begin(p->focusGroup.second),std::end(p->focusGroup.second),idx)==std::end(p->focusGroup.second))) {continue;}
			//ignore facet if no opacity or if facet is not in focus group
			//if no focus group is given, than all facets are the focus group
			double error_event_facet=(it==-1)?simHistory->errorList_event.getCurrent(idx):simHistory->errorList_event.historyList.second[idx][it];
			if(!std::isinf(error_event_facet)&&!std::isnan(error_event_facet)){//ignore facet if no hits (=inf error)
				//error_event+=error_event_facet*pow(f.sh.area,1/2);
				//area_event+=f.sh.area;
				error_event+=error_event_facet;
			}

			double error_covering_facet=(it==-1)?simHistory->errorList_covering.getCurrent(idx):simHistory->errorList_covering.historyList.second[idx][it];
			if(!std::isinf(error_covering_facet)&&!std::isnan(error_covering_facet)){//ignore facet if no hits (=inf error)
				//error_covering+=error_covering_facet*pow(f.sh.area,1/2);
				//area_covering+=f.sh.area;
				error_covering+=error_covering_facet;
			}
		}
	}
	//return std::make_tuple(error_event/(pow(area_event,0.5)), error_covering/(pow(area_covering,0.5)));
	return std::make_tuple(error_event/(pow(simHistory->numFacet,0.5)), error_covering/(pow(simHistory->numFacet,0.5)));
}

// Check if error reached targets
bool checkError(double targetError, double currentError, double factor, std::string mode){
	bool vipCheck = currentError<=targetError;
	if(!p->vipFacets.empty()){
		HistoryList<double> *listptr;
		listptr = getErrorList(mode);
		if(listptr==NULL){
			return true;
		}

		for(unsigned int i = 0; i < p->vipFacets.size(); i++){
			//vipCheck = vipCheck && (listptr->getCurrent(p->vipFacets[i].first)== std::numeric_limits<double>::infinity() || listptr->getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
			vipCheck = vipCheck && (listptr->getCurrent(p->vipFacets[i].first) <= p->vipFacets[i].second * factor);
		}
	}
	return vipCheck;
}

HistoryList<double>* getErrorList(std::string mode){
	if(mode=="covering"){
		return &simHistory->errorList_covering;
	}
	else if(mode=="event"){
		return &simHistory->errorList_event;
	}
	else{
		std::cout<<"------------! Error mode '"<<mode <<"' not implemented !------------\n";
		return NULL;
	}
}

//-----------------------------------------------------------
//Functions for covering threshold
void initCoveringThresh(){
	// Initializes size of cveringThreshold
	sHandle->coveringThreshold=std::vector<llong> ();
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for(uint i=0; i<sHandle->structures[j].facets.size(); i++){
			sHandle->coveringThreshold.push_back(0);
		}
	}
}

void setCoveringThreshold(int size, int rank){
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet per subprocess
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if(rank==1)
				{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim + simHistory->coveringList.getCurrent(&f).convert_to<llong>()%num_sim;}
			else
				{cov_sim=simHistory->coveringList.getCurrent(&f).convert_to<llong>()/num_sim;}

			// Covering threshold: current covering of face - particles that can be desorbed
			sHandle->coveringThreshold[getFacetIndex(&f)]=simHistory->coveringList.getCurrent(&f).convert_to<llong>()-cov_sim;
		}
	}
}

//-----------------------------------------------------------

std::tuple<bool, double> TimestepControl(Databuff *hitbuffer_sum){//Return value => bool: is a change of the time step necessary; double: new time step
	std::ostringstream tmpstream (std::ostringstream::app);
	tmpstream << "--------------------------------------" << std::endl;
	tmpstream << "TimestepControl function is executed!: " << std::endl << std::endl;
	bool repetition = 0;
	double stepSize_recom = simHistory->stepSize;//recom (recommendation of step size)
	std::vector<bool> cases = {0, 0, 0, 0, 0};
	double Krealvirt = GetMoleculesPerTP(hitbuffer_sum).convert_to<double>();
	//Check possibilities of different error cases

	//------------------------ CASE I --------------------
	double decreasethreshold = 0.10; //if particles of a facet decrease stronger than implemented by the model,
	//the iteration step has to be repeated with better statistics.
	//decreasethreshold describes the uncertainty.
	double desorption_f = 0;
	//double cov_f_b4 = 0;
	//double cov_f_aft = 0;
	BYTE *buffer = hitbuffer_sum->buff;


	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
			desorption_f = double(f.sh.desorption);
			//cov_f_b4 = double(simHistory->coveringList.getCurrent(getFacetIndex(&f)));
			//cov_f_aft = double(simHistory->coveringList.getPredict(getFacetIndex(&f)));
			if((facetHitBuffer->nbDesorbed == 0)&&(desorption_f > 0)){
				tmpstream << "Facet "<< getFacetIndex(&f) << " did not desorb test-particles!" <<std::endl;
				tmpstream << "Calculated desorption (physical particles/Krealvirt) of Facet "<< getFacetIndex(&f) << " would be " << (desorption_f/Krealvirt) << " particles."<<std::endl;
				tmpstream << "CASE I is not triggered!" <<std::endl;
			}
			else if((facetHitBuffer->nbDesorbed == 0)&&(desorption_f == 0)){
				continue;
				}
			else if((facetHitBuffer->nbDesorbed*Krealvirt)<(desorption_f*(1-decreasethreshold))){
			tmpstream << "CASE I is detected for facet " <<getFacetIndex(&f) <<"."<<std::endl;
			cases.at(0) = 1;
			repetition = true;
			}
		}
	}
	if(repetition){
		p->targetError *= 0.5;
		p->targetParticles *= 4;
		tmpstream << "CASE II will not be checked." <<std::endl;
	}
	else{//if(!repetition){
		//------------------------ CASE II -------------------
			double growththreshold = 0.05; //if total amount of particles is not conserved within the growth threshold,
			//the iteration step has to be repeated with better statistics.
			//total number of particles of last simulation:
			boost::multiprecision::uint128_t total_covering_before = 0;
			for(unsigned int j=0; j<simHistory->coveringList.currentList.size();j++){
				total_covering_before +=simHistory->coveringList.getCurrent(j);
			}
			double tot_cov_b4 = double(total_covering_before);
			//total number of particles of last simulation:
			boost::multiprecision::uint128_t total_covering_after = 0;
			for(unsigned int j=0; j<simHistory->coveringList.predictList.size();j++){
				total_covering_after += simHistory->coveringList.getPredict(j);
			}
			double tot_cov_aft = double(total_covering_after);
			if((tot_cov_b4+getnbOutgassed(hitbuffer_sum))*(1+growththreshold)<tot_cov_aft+simHistory->nLeaks*Krealvirt){
				p->targetError *= 0.5;
				p->targetParticles *= 4;
				repetition = true;
				tmpstream << "CASE II detected." <<std::endl;
				cases.at(1) = 1;
			}
			else if((tot_cov_b4+getnbOutgassed(hitbuffer_sum))*(1-growththreshold)>tot_cov_aft+simHistory->nLeaks*Krealvirt){
				p->targetError *= 0.5;
				p->targetParticles *= 4;
				repetition = true;
				tmpstream << "CASE II detected." <<std::endl;
				cases.at(1) = 1;
						}
			else{//Cases I and II are not triggered
				if(p->targetError<p->targetError_input){
					p->targetError *=2;
				}
				if(p->targetParticles>p->targetParticles_input){
					p->targetParticles *=0.25;
				}
			}

	}
	//------------------------ CASE IV -------------------
	double transitionthreshold = 0.05;//If change of coverage is larger than transitionthreshold in case of a
	// transition [=> if(coverage_f_b4 > 0.01 && coverage_f_b4 < 1.1)], reduce the time step
	double coverage_f_b4 = 0;
	double coverage_f_aft = 0;
	/*double tau_0= 0;
	double residence_energy = 0;
	double tau = 0;
	*/
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				coverage_f_b4 = double(calcCoverage(&f));
				coverage_f_aft = double(calcPredictedCoverage(&f));
				//std::cout << "Facet " << getFacetIndex(&f) << std::endl;
				//std::cout << "coverage_f_b4 = " << coverage_f_b4 << std::endl;
				//std::cout << "coverage_f_aft = " << coverage_f_aft << std::endl;
				//std::cout << "std::abs(cov_f_aft-cov_f_b4) = " << std::abs(coverage_f_aft-coverage_f_b4) << std::endl;
				if((coverage_f_b4 > 0.01 && coverage_f_b4 < 1.1)||(coverage_f_aft > 0.01 && coverage_f_aft < 1.1)){
					if(std::abs(coverage_f_aft-coverage_f_b4)>transitionthreshold){
						/*if(coverage_f_b4 >= 1){
							residence_energy = p->H_vap;
						}
						else {
							residence_energy = p->E_de;
						}
						tau_0 = double(h/(kb*f.sh.temperature));
						tau = tau_0 *exp(residence_energy / (kb*f.sh.temperature));
						if (tau < stepSize_recom){
							repetition = true;
							stepSize_recom = std::max(tau,p->t_min);
						}
						*/
						if(stepSize_recom>p->t_min){
							repetition = true;
							tmpstream << "CASE IV is detected for facet " <<getFacetIndex(&f) <<"."<<std::endl;
							tmpstream << "Other facets will not be checked!"<<std::endl;
							cases.at(3) = 1;
							stepSize_recom *= 0.1;
							if(stepSize_recom<p->t_min){
								stepSize_recom = p->t_min;
								tmpstream << "Step size is limited by t_min."<<std::endl;
							}
							break;
						}

					}
				}
			}
	}

	//------------------------ CASE V --------------------

	double f_HitEquiv = 0;
	llong f_Des = 0;
	double f_AbsEquiv = 0;
	bool too_large_reemission = false;
	bool too_large_desorption = false;
	double tau_0= 0;
	double residence_energy = 0;
	double tau = 0;
	coverage_f_b4 = 0;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			f_HitEquiv = getFacetHitBuffer(&f,hitbuffer_sum)->nbHitEquiv;
			f_Des = getFacetHitBuffer(&f,hitbuffer_sum)->nbDesorbed;
			f_AbsEquiv = getFacetHitBuffer(&f,hitbuffer_sum)->nbMCHit;
			coverage_f_b4 = double(calcCoverage(&f));
			tau_0 = double(h/(kb*f.sh.temperature));
			if(coverage_f_b4 >= 1){
				residence_energy = p->H_vap;
				tau = tau_0 *exp(residence_energy / (kb*f.sh.temperature));
			}
			else {
				tau = tau_0 *(coverage_f_b4*exp(p->H_vap / (kb*f.sh.temperature))+(1-coverage_f_b4)*exp(p->E_de / (kb*f.sh.temperature)));
			}
			too_large_reemission = (((f_HitEquiv-f_AbsEquiv)*(Krealvirt/calcNmono(&f))*(tau/simHistory->stepSize))>1);
			too_large_desorption = (((f_Des*Krealvirt/calcNmono(&f))>1)&&(f_HitEquiv*(Krealvirt/calcNmono(&f))*(tau/simHistory->stepSize)>1));
			if (too_large_reemission || too_large_desorption){
				if(stepSize_recom>p->t_min){
					repetition = true;
					tmpstream << "CASE V is detected for facet " <<getFacetIndex(&f) <<"."<<std::endl;
					tmpstream << "Other facets will not be checked!"<<std::endl;
					cases.at(4) = 1;
					stepSize_recom *= 0.1;
					if(stepSize_recom<p->t_min){
						stepSize_recom = p->t_min;
						tmpstream << "Step size is limited by t_min."<<std::endl;
					}
					break;
				}

			}
		}
	}

	//------------------------ CASE III ------------------
	//has to be at the end, since it overwrites stepSize_recom with p->t_min
	double avg_flighttime = estimateAverageFlightTime(hitbuffer_sum);
	std::cout << "avg_flighttime = " << avg_flighttime << std::endl;
	if (avg_flighttime > p->t_min){
		p->t_min = 10*avg_flighttime;// factor 10 (randomly chosen) to have more of the flight time distribution covered
		repetition = true;
		tmpstream << "CASE III detected." <<std::endl;
		stepSize_recom = p->t_min;
		cases.at(2) = 1;
	}
	simHistory->flightTime=0.0;
	//------------------------ End of CASES I to V -------

	if(repetition){//revert the value of lastTime in case of repetition
		simHistory->lastTime-=simHistory->stepSize;
	}
	tmpstream << "Cases: ";
	for (int elem : cases){
		tmpstream << " " << elem;
	}
		tmpstream << std::endl;
	printStream(tmpstream.str());
	return std::make_tuple(repetition,stepSize_recom);
}


/*not needed anymore
//Estimation of Tmin
double estimateTmin(Databuff *hitbuffer){ //not ready yet => finish //TODO
//Ich muss der Funktion noch einen Hitbuffer übergeben. Ich brauche ja 'covering'.
	double tmin=0;
	double sum_v_avg = 0;
	double normalization_factor_v = 0;
	//llong covering;
	//avarage of <v> (<v> is the average velocity on a single facet depending on temperature)
    // over all facets weighted with the rate of outgoing particles (outgassing + desorption) per facet
	double tmin_particles_out = 0;
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
			for (SubprocessFacet& f : sHandle->structures[s].facets) {
				boost::multiprecision::uint128_t covering = getFacetHitBuffer(&f,hitbuffer)->covering;
				double v_avg_therm = 3*pow((3.14159265359/8),0.5)* pow((8.314472* f.sh.temperature/(sHandle->wp.gasMass*0.001)),0.5); //0.001 to convert MolarMass form g to kg
				sum_v_avg += v_avg_therm * (f.sh.outgassing + f.sh.desorption.convert_to<double>());
				normalization_factor_v += f.sh.outgassing + f.sh.desorption.convert_to<double>();
				if ((f.sh.outgassing + f.sh.desorption.convert_to<double>()) > 0){ //avoid division by 0
					double ttemp= (double)((double) covering/((f.sh.outgassing + f.sh.desorption.convert_to<double>())/ (kb*f.sh.temperature)));
					if (!tmin_particles_out){
						tmin_particles_out = (ttemp);
					}
					if (tmin_particles_out > ttemp){
						tmin_particles_out = ttemp;}
				}				
			 }
	}
	//if (normalization_factor_v == 0) std::cout << "divide by zero!" << std::endl;
	double v_avg = sum_v_avg/normalization_factor_v;
	std::cout << "v_avg = sum_v_avg/normalization_factor_v = " << sum_v_avg << "/" << normalization_factor_v << std::endl;
	double av_path_length = sHandle->tmpGlobalResult.distTraveled_total*0.01/sHandle->tmpGlobalResult.globalHits.nbMCHit;// The 0.01 is because the distance should be saved in cm. ToDo: Check up, if that's right.
	std::cout <<  "av_path_length = sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.nbMCHit = " << sHandle->tmpGlobalResult.distTraveled_total << "/" << sHandle->tmpGlobalResult.globalHits.nbMCHit << std::endl;
	std::cout << "tmin = av_path_length /v_avg = " << av_path_length << " / " << v_avg << std::endl;
	tmin = av_path_length /v_avg;

	//Time step information for developing

	std::cout << "Time step information for developing" << std::endl;
	std::cout << "estimateTmin: tmin = " <<tmin<< "ms"<< std::endl;
	std::cout << "estimateTmin: tmin_particles_out = " <<tmin_particles_out<< "ms"<< std::endl;
	return tmin;
}
*/

double estimateAverageFlightTime(Databuff *hitbuffer_sum){
	//return simHistory->flightTime/simHistory->nParticles;
	llong nbtotalHits = getnbHits(hitbuffer_sum);
	return simHistory->flightTime/nbtotalHits;
}


//----------deprecated functions because hitbuffer not sent to sub processes anymore
/*
void setCoveringThreshold(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong num_sim=size-1;
	llong cov_sim; //number of particles that can be desorbed from facet
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{cov_sim=facetHitBuffer->covering/num_sim + facetHitBuffer->covering%num_sim;}
				else
					{cov_sim=facetHitBuffer->covering/num_sim;}

				sHandle->coveringThreshold[getFacetIndex(&f)]=facetHitBuffer->covering-cov_sim;
			}
	}
}
*/

//----------Not used anymore
/*
void allocateCovering(Databuff *hitbuffer, int size, int rank){
	BYTE *buffer;
	buffer = hitbuffer->buff;
	llong temp=size+1;
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
			for (SubprocessFacet& f : sHandle->structures[j].facets) {
				FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset);
				if(rank==1)
					{facetHitBuffer->covering=facetHitBuffer->covering/temp + facetHitBuffer->covering%temp;}
				else
					{facetHitBuffer->covering=facetHitBuffer->covering/temp;}
			}
	}
}
*/
