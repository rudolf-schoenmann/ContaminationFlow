# Quick start ContaminationFlow

## Copy git reposity

```
$ git clone https://gitlab.lrz.de/rudolf_schoenmann/ContaminationFlow.git

```

## Dependencies
Install dependencies

```
# Examples are given for RockyLinux and the dnf package manager:

# Running the ContaminationFlow:
$ dnf install openmpi

#Building ContaminationFlow
$ dnf install gsl
$ dnf install gcc
$ dnf install gcc-c++
$ dnf install make

#Version control with git and SSH
$ dnf install git
$ dnf install openssh

#Integrated Development Environment - Eclipse
$ dnf install java-17-openjdk



```

## Build

We use [Eclipse IDE for C/C++ Developers](https://www.eclipse.org/downloads/packages/release/photon/r/eclipse-ide-cc-developers) version photon release (4.8.0).
Open Eclipse via terminal:

```
# Load MPI -- to find mpic++
module load mpi


# Optional: check if loaded correctly
which mpic++ #e.g., /usr/lib64/openmpi/bin/mpic++

# Open eclipse
cd /path-to-eclipse-executable/
./eclipse
```
### Create ContaminationFlow project
* Create new workspace (**File->Switch Workspace->Other...**).
* Add project (**File->Open Projects from File System...** and select ContaminationFlow).


### Build Execulatbe
The ContaminationFlow executable can be build via terminal.

```
# Load MPI
module load mpi


# Build ContaminationFlow executable
cd /path-to-ContaminationFlow/Debug/
make all
```

The ContaminationFlow executable can be build via Ecplise specifying the 'Build directory' as '${workspace_loc:/ContaminationFlow}/Debug' in the Eclipse project properties in the 'C/C++ Build' settings.


## Run

```
# Load mpi
$ module load mpi
source /path-to-ContaminationFlow/StartContaminationFlow.tcsh

# Or run ContaminationFlow executable directly
mpirun -n N /path-to-ContaminationFlow/Debug/ContaminationFlow /path-to-an-Inputfile/InputFile.txt

```

with

|Parameter| Valid Input | Description|
|------------ | -------------| -------------|
|N| int value, N > 1|Number of processes. 1 main process, N-1 simulation processes|
|path-to-ContaminationFlow| Directory | path to github repository|
|/path-to-an-Inputfile/InputFile.txt| Readable text file | path to an input file that defines simulation parameters|
|save| int value | 1 (default) simulation results saved, 0 not saved|

Important input file parameters:

|Parameter Name| Valid Input | Description|
|------------ | -------------| -------------|
|loadbufferPath|Readable file|Loadbuffer file|
|hitbufferPath|Readble file | Hitbuffer file|
|coveringPath|Readble file | Covering file, sets covering per facet. Hitbuffer not imported if given|
|simulationTime|Double value | "Number" of simulation time (=computation time before targets are reached)|
|unit|String|"Unit" of simulation time (=computation time before targets are reached)|
|iterationNumber|Int value|Number of desired iterations|
|usePCMethod|Int value|0: Do not use PC method, 1: Use PC Method v1, 2: Use PC Method v2|
|maxTime| Double value| "Number" of maximum simulated time (=total time in simulated system)|
|maxUnit| String|"unit" of maximum simulated time (=total time in simulated system)|
|errorMode|String: covering,event| Desired error that is tracked|
|targetParticles|Int value| Number of target particles per iteration|
|targetError|Double value| Target error per iteration|
|t_min|Double value| Minimum time step (=simulated time between desorb and adsorb)|
|rollingWindowSize|Int value| Number of precedings iterations used to calculate statistics|
|vipFacets|Alterning sequence: int double| Alterning sequence of facet number and respective target errors, seperated by tabs|

For more information on input file see **ContaminationFlow_doc/main.pdf**.

Notes:
* The simulation will not start in the following cases
  * N <= 1
  * More than zero moments defined in geometry
  * Invalid parameter names in input file
  * Invalid values for certain parameters in input file
  * Invalid hitbuffer/loadbuffer
  * Two-sided facet with opacity exists
  * Invalid covering file
* Input file
  * One parameter per line
  * Parameter name and value are seperated by a tab (no =)
  * Lines can be commented in input file by adding \# at the start of the line
* Covering file
  * Seperate text file
  * Two options:
    * Set covering: covering followed by covering values per facet. Example for 7 facet geometry:
      * covering 0 1312749422390254 1312749422390254 1312749422390254 1312749422390254 1312749422390254 1312749422390254
    * Set coverage: coverage followed by coverage values per facet. Example for 7 facet geometry:
      * coverage 0 1 1 1 1 1 1
* Export from ContaminationFlow Windows
  * Hitbuffer and Loadbuffer files
  * Covering file (directly copy from GUI or save to file)
  * Facet groups for input file (directly copy from GUI or save to file)
* In case an iteration takes too long (e.g., small residence time -> long time until adsorbtion), the following can be done
  * Decrease *targetError* and/or *targetParticles* in input file (you should probably not go below 0.2/50 though)
  * Decrease *maxTimePerIt* in input file
  * Increase *hitRatioLimit* in input file (you should probably not go above 0.01 though)
* An iteration ends in the following cases
  * Last event was not a hit and targets are reached
  * covering threshold is reached -> can lead to negative covering if not stopped
* To stop the simulation if converged, the following conditions have to hold
  * *stopConverged* is set to 1 in input file (true)
  * Minimum of *rollingWindowSize* iterations have elapsed
  * *convergenceTarget*, i.e., the ratio of std/mean per facet, is reached
  * *convergenceTime*, i.e., the simulated time since the start, is reached
* The simulation ends in the following cases
  * Maximum simulated time is reached
  * Virtually no desorption/outgassing
  
