# Quick start ContaminationFlow

## Copy github reposity

```
git clone https://gitlab.lrz.de/rudolf_schoenmann/MolflowLinux.git

```

## Dependencies
Install dependencies (_gsl-devel_, _openmpi-devel_)

```
# Fedora & dnf package manager
dnf install gsl-devel
dnf install openmpi-devel

```

## Build

We use [Eclipse IDE for C/C++ Developers](https://www.eclipse.org/downloads/packages/release/photon/r/eclipse-ide-cc-developers) version photon release (4.8.0).
Open Eclipse via terminal:

```
# Optional: Check possible modules
module avail #list should include mpi/openmpi-x86_64

# Load mpi -- to find mpic++
module load mpi/openmpi-x86_64

# Optional: check if loaded correctly
which mpic++ #e.g., /usr/lib64/openmpi/bin/mpic++

# Open eclipse
cd /path-to-eclipse-executable/
./eclipse
```
### Create Molflow project
* Create new workspace (**File->Switch Workspace->Other...** and insert the parent directory of MolflowLinux).
* Add project (**File->Open Projects from File System...** and select MolflowLinux).
* Build


### Build Execulatbe
The Molflow executable can be build via Ecplise (Note: always follow the steps above to open Eclipse).
Alternatively, after the first build, it can be build via terminal:

```
# Load mpi
module load mpi/openmpi-x86_64

# Build MolflowLinux executable
cd /path-to-MolflowLinux/Debug/
make MolflowLinux
```


## Run

```
# Load mpi
module load mpi/openmpi-x86_64

# Run MolflowLinux executable
mpirun -n N /path-to-MolflowLinux/Debug/MolflowLinux /path-to-an-Inputfile/InputFile.txt save

```

with

|Parameter| Valid Input | Description|
|------------ | -------------| -------------|
|N| int value, N > 1|Number of processes. 1 main process, N-1 simulation processes|
|path-to-MolflowLinux| Directory | path to github repository|
|/path-to-an-Inputfile/InputFile.txt| Readable text file | path to an input file that defines simulation parameters|
|save| int value | 1 (default) simulation results saved, 0 not saved|

Important input file parameters:

|Parameter Name| Valid Input | Description|
|------------ | -------------| -------------|
|loadbufferPath|Readable file|Loadbuffer file|
|hitbufferPath|Readble file | Hitbuffer file|
|simulationTime|Double value | "Number" of simulation time (=computation time before targets are reached)|
|unit|String|"Unit" of simulation time (=computation time before targets are reached)|
|iterationNumber|Int value|Number of desired iterations|
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
* Input file
  * One parameter per line
  * Parameter name and value are seperated by a tab (no =)
  * Lines can be commented in input file by adding \# at the start of the line
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
  