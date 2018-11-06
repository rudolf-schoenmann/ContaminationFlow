#include "SimulationLinux.h"

extern Simulation *sHandle;

double estimateTmin(){//TODO something is wrong here
	double tmin=1;
	int facetcounter=0;
	double sum_1_v_ort=0.0;
	double sum_abs=0.0;
	//double sum_hits=0.0;

	//mittelwert 1/v_ort
	size_t nbMoments = sHandle->moments.size();
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			for (size_t m = 0; m <= nbMoments; m++) {
				//sum_1_v_orth
				sum_1_v_ort+=f.tmpCounter[m].hit.sum_1_per_ort_velocity*f.sh.area; //(s/m)*m^2
				facetcounter++;

				//covering
				double coveringtemp=f.tmpCounter[m].hit.covering;
				if(coveringtemp>0)
					sum_abs+=calcNmono(&f)*coveringtemp; //mass
				//sum_hits+=f.tmpCounter[m].hit.nbHitEquiv;//or nbMChit
			}
		}
	}
	sum_1_v_ort /=facetcounter;//(s*m)
	//double Nmean=sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv/sum_abs;//(1/mass) //nbHitEquiv or nbMCHit
	double Nmean=1/sum_abs;//(1/mass) //nbHitEquiv or nbMCHit
	std::cout <<"test " <<sum_1_v_ort <<std::endl;
	std::cout <<"test " <<Nmean <<std::endl;

	//TODO durchschnittstrecke
	//double dmean=sHandle->tmpGlobalResult.distTraveled_total/sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;//m
	double dmean=sHandle->tmpGlobalResult.distTraveled_total;//m
	std::cout <<"test " <<dmean <<std::endl;

	tmin=dmean*sum_1_v_ort*Nmean;//s
	return tmin*1000;//ms
}
