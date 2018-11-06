#include "worker.h"
#include "Simulation.h"
#include <vector>

//MolFlow *mApp;
extern Simulation *sHandle; //delcared in molflowSub.cpp

Worker::Worker() {

	//Molflow specific
	temperatures = std::vector<double>();
	//moments = std::vector<double>();
	//desorptionParameterIDs = std::vector<size_t>();
	//userMoments = std::vector<std::string>(); //strings describing moments, to be parsed
	CDFs = std::vector<std::vector<std::pair<double, double>>>();
	//IDs = std::vector<std::vector<std::pair<double, double>>>();
	parameters = std::vector<Parameter>();
	//displayedMoment = 0; //By default, steady-state is displayed
	/*
	sHandle->wp.timeWindowSize = 1E-10; //Dirac-delta desorption pulse at t=0
	sHandle->wp.useMaxwellDistribution = true;
	sHandle->wp.calcConstantFlow = true;
	sHandle->wp.gasMass = 28.0;
	sHandle->wp.enableDecay = false;
	sHandle->wp.halfLife = 1;
	sHandle->wp.finalOutgassingRate = sHandle->wp.finalOutgassingRate_Pa_m3_sec = sHandle->wp.totalDesorbedMolecules = 0.0;
	sHandle->wp.motionType = 0;
	sHandle->wp.sMode = MC_MODE;*/
}


std::vector<std::pair<double, double>> Worker::Generate_CDF(double gasTempKelvins, double gasMassGramsPerMol, size_t size){
	std::vector<std::pair<double, double>> cdf; cdf.reserve(size);
	double Kb = 1.38E-23;
	double R = 8.3144621;
	double a = sqrt(Kb*gasTempKelvins / (gasMassGramsPerMol*1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

	//Generate cumulative distribution function
	double mostProbableSpeed = sqrt(2 * R*gasTempKelvins / (gasMassGramsPerMol / 1000.0));
	double binSize = 4.0*mostProbableSpeed / (double)size; //distribution generated between 0 and 4*V_prob

	for (size_t i = 0; i < size; i++) {
		double x = (double)i*binSize;
		double x_square_per_2_a_square = pow(x, 2) / (2 * pow(a, 2));
		cdf.push_back(std::make_pair(x, 1 - exp(-x_square_per_2_a_square)*(x_square_per_2_a_square + 1)));

	}


	return cdf;
}

int Worker::GetCDFId(double temperature) {

	int i;
	for (i = 0; i<(int)temperatures.size() && (abs(temperature - (double)(temperatures[i]))>1E-5); i++); //check if we already had this temperature
	if (i >= (int)temperatures.size()) i = -1; //not found
	return i;
}

int Worker::GenerateNewCDF(double temperature){
	size_t i = temperatures.size();
	temperatures.push_back(temperature);
	CDFs.push_back(Generate_CDF(temperature, sHandle->wp.gasMass, CDF_SIZE));
	return (int)i;
}

void Worker::CalcTotalOutgassing() {
	// Compute the outgassing of all source facet
	sHandle->wp.totalDesorbedMolecules = sHandle->wp.finalOutgassingRate_Pa_m3_sec = sHandle->wp.finalOutgassingRate = 0.0;

	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if (f.sh.desorbType != DES_NONE) { //there is a kind of desorption
				if (f.sh.useOutgassingFile) { //outgassing file
					for (unsigned int l = 0; l < (f.sh.outgassingMapWidth*f.sh.outgassingMapHeight); l++) {
						sHandle->wp.totalDesorbedMolecules += sHandle->wp.latestMoment * f.outgassingMap[l] / (1.38E-23*f.sh.temperature);
						sHandle->wp.finalOutgassingRate += f.outgassingMap[l] / (1.38E-23*f.sh.temperature);
						sHandle->wp.finalOutgassingRate_Pa_m3_sec += f.outgassingMap[l];
					}
				}
				else { //regular outgassing
					if (f.sh.outgassing_paramId == -1) { //constant outgassing
						sHandle->wp.totalDesorbedMolecules += sHandle->wp.latestMoment * f.sh.outgassing / (1.38E-23*f.sh.temperature);
						sHandle->wp.finalOutgassingRate += f.sh.outgassing / (1.38E-23*f.sh.temperature);  //Outgassing molecules/sec
						sHandle->wp.finalOutgassingRate_Pa_m3_sec += f.sh.outgassing;
					}
					else { //time-dependent outgassing
						//sHandle->wp.totalDesorbedMolecules += IDs[f.sh.IDid].back().second / (1.38E-23*f.sh.temperature);
						size_t lastIndex = parameters[f.sh.outgassing_paramId].GetSize() - 1;
						double finalRate_mbar_l_s = parameters[f.sh.outgassing_paramId].GetY(lastIndex);
						sHandle->wp.finalOutgassingRate += finalRate_mbar_l_s *0.100 / (1.38E-23*f.sh.temperature); //0.1: mbar*l/s->Pa*m3/s
						sHandle->wp.finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s *0.100;
					}
				}
			}
		}
	}
	//if (mApp->globalSettings) mApp->globalSettings->UpdateOutgassing();

}



