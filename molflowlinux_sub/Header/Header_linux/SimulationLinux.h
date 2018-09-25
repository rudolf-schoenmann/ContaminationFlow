#include "Simulation.h"

static const char *day[]={"day","Day","days","Days","d","D"};
static const char *hour[]={"hour","Hour","hours","Hours","h","H","hr","Hr"};
static const char *min[]={"Minutes","minutes","Minute","minute","min","Min","m","M"};

bool simulateSub(Databuff *hitbuffer, int rank, int simutime);
double convertunit(double simutime, std::string unit);
void UpdateSubHits(Databuff *databuffer, int rank);
void UpdateSubMCHits(Databuff *databuffer, int rank, size_t nbMoments);
void initbufftozero(size_t nbMoments, BYTE *buffer,SubprocessFacet *f);
