#include "Simulation.h"
#include "GLApp/MathTools.h"
#include <math.h>

// Global handles
extern Simulation* sHandle; //Declared at molflowSub.cpp

double updatecovering(SubprocessFacet *iFacet)
{
	//TODO: adapt units, this may not yet be the correct result
	double N_mono= iFacet->sh.area/Sqr(76*pow(10.0,-12.0));
	double dN_surf=sHandle->wp.gasMass/12.011;

	return dN_surf/N_mono;

}
