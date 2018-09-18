#include "Simulation.h"

static const char *day[]={"day","Day","days","Days","d","D"};
static const char *hour[]={"hour","Hour","hours","Hours","h","H","hr","Hr"};
static const char *min[]={"Minutes","minutes","Minute","minute","min","Min","m","M"};

bool simulateSub(Databuff *hitbuffer, int rank, double simutime, std::string unit);
double convertunit(double simutime, std::string unit);
