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


#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include "Buffer.h"
#include <unistd.h>
#include "Simulation.h"
#include "SimulationLinux.h"


typedef void *HANDLE;

// Global process variables
Simulation* sHandle; //Global handle to simulation, one per subprocess
SimulationHistory* simHistory;
ProblemDef* p;

//This function checks if the correct number of arguments has been passed
//does not check their validity, e.g. right type such as double/string, correct filename, etc

bool parametercheck(int argc, char *argv[], ProblemDef *p, int rank) {
	int i;
	if(rank==0){
		printf("argc: %d\n", argc);
		for (i = 0; i < argc; i++) {
			printf("argv[%d]: %s\n", i, argv[i]);
		}
	}
	if(argc>3){ // list of parameters given
		std::cout <<std::endl;
		if (argc < 5 || argc > 6) { //Check correct number of arguments
			if(rank==0){
				std::cout << "MolflowLinux requires 4 mandatory arguments and 1 optional argument."<< std::endl;
				std::cout << "Please pass these arguments to MolflowLinux:"<< std::endl;
				std::cout << "1. Name of load-buffer file to read in (e.g. loadbuffer)."<< std::endl;
				std::cout << "2. Name of hit-buffer file to read in (e.g. hitbuffer)."<< std::endl;
				std::cout << "3. Save results: 0=false, 1=true."<< std::endl;
				std::cout << "4. The total simulation time (e.g 2.5)." << std::endl;
				std::cout << "5. [OPTIONAL] Simulation time unit (e.g. seconds, minutes, hours, days). Default set to seconds." << std::endl;
				std::cout << "MolflowLinux is terminated now." << std::endl;
			}
		}
		else{ //Read input arguments
			if(rank==0){std::cout<<"Read arguments" <<std::endl;}

			if(!checkReadable(argv[1])||!checkReadable(argv[2])||std::atof(argv[4])<0.0) // check if parameters are feasible
				{return false;}
			p->readArg(argc, argv, rank);
			return true;
		}
	}
	else if(argc<4 && argc>1){ // Read input file
		if(checkReadable(argv[1])){
			p->readInputfile(argv[1],rank, argc==3?(int)std::atof(argv[2]):1);
			if(!checkReadable(p->hitbufferPath)||!checkReadable(p->loadbufferPath)){return false;}
			return true;}
		}
	return false;
	}

//-----------------------------------------------------------
//Main Function
int main(int argc, char *argv[]) {

	// Initialise data buffers
	Databuff hitbuffer; //Hitbuffer for the data of the subprocesses
	hitbuffer.buff=NULL;

	Databuff hitbuffer_sum; //Hitbuffer to sum up all of the subprocesses' data
	hitbuffer_sum.buff=NULL;

	Databuff loadbuffer; //Loadbuffer to read in data of geometry and physical parameters
	loadbuffer.buff=NULL;

	double t0,t1;
	double computationTime=0.0; // total time (measured by system clock) used for simulation


	/* Create child processes, each of which has its own variables.
	 * From this point on, every process executes a separate copy
	 * of this program.  Each process has a different process ID,
	 * ranging from 0 to num_procs minus 1, and COPIES of all
	 * variables defined in the program. No variables are shared.
	 **/
	MPI_Init(NULL, NULL);
	/* find out MY process ID, and how many processes were started. */

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(world_size<=1){
		if (rank == 0){
			std::cout <<"Minimum number of 2 processes needed. Currently "<<world_size <<std::endl;
		}
		MPI_Finalize();
		return 0;
	}

//---- Initialize ProblemDef that defines relevant simulation parameters for MolflowLinux
	p = new ProblemDef();
	// Check parameters (defined by command line arguments or input file) to be written to ProblemDef p
	if (!parametercheck(argc, argv,p,rank)) {
		if (rank == 0){
			std::cout <<"Required parameters cannot be read." <<std::endl;
			std::cout <<"Check command line arguments and/or parameters in input file." <<std::endl;
			std::cout <<"Ending Simulation." <<std::endl;
		}

		MPI_Finalize();
		return 0;
	}

	if (rank == 0) {
		p->printInputfile(p->outFile);

		//Read in buffer file (exported by Windows-Molflow). File given as first argument to main().
		importBuff(p->loadbufferPath,&loadbuffer);
		importBuff(p->hitbufferPath,&hitbuffer);

		std::cout << "Buffers sent. Wait for a few seconds. " << std::endl;
	}

//---- Send load-buffer to all other processes
	// Send size of loadbuffer from main process to subprocesses
	MPI_Bcast(&loadbuffer.size, sizeof(loadbuffer.size), MPI::BYTE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// Allocate memory for loadbuffer in subprocesses
	if (rank != 0) {
		loadbuffer.buff = new BYTE[loadbuffer.size];
	}

	// Send laodbuffer content to all subprocesses
	MPI_Bcast(loadbuffer.buff, loadbuffer.size, MPI::BYTE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

//----Send hit-buffer size to all other processes.
	// Send size of hitbuffer from main process to subprocesses in subprocesses
	MPI_Bcast(&hitbuffer.size, sizeof(hitbuffer.size), MPI::BYTE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// Allocate memory for hitbuffer
	if (rank != 0) {
		hitbuffer.buff = new BYTE[hitbuffer.size];
	}

//----Create Simulation handle and preprocess hitbuffer
	if (p->simulationTimeMS != 0) {
		//Creates sHandle instance for process 0 and all subprocesses (before the first iteration step starts)
		InitSimulation();
		// Load geometry from buffer to sHandle
		if (!LoadSimulation(&loadbuffer)) {
			std::cout << "Geometry not loaded." << std::endl;
			std::cout << "MolflowLinux is terminated now." << std::endl;
			MPI_Finalize();
			return 0;
		}
		initCoveringThresh();
		UpdateSojourn();

		if(rank==0){
			//Initialize hitbuffer_sum as copy of hitbuffer in main process
			//The hitbuffers of all subprocesses will be added up and written in the hitbuffer_sum
			hitbuffer_sum.buff = new BYTE[hitbuffer.size];
			memcpy(hitbuffer_sum.buff,hitbuffer.buff,hitbuffer.size);
			hitbuffer_sum.size =hitbuffer.size;

			// Initialize simHistory for main process from hitbuffer
			//simHistory contains the relevant results/quantities of the simulation. E.g., covering history, total simulated time, etc.
			simHistory = new SimulationHistory (&hitbuffer, world_size);
			//TODO: maybe add possibility of covering.txt file input
		}
		else{
			//Initialize simHistory for subprocesses
			simHistory = new SimulationHistory(world_size);
		}
		//Send hitbuffer to all subprocesses
		MPI_Bcast(hitbuffer.buff, hitbuffer.size, MPI::BYTE, 0, MPI_COMM_WORLD);
	}

//----Simulation
	int it = 0;
	double currentRatio=0.0;
	while(true){
		it++;
		// Start of Simulation
		if (p->simulationTimeMS != 0) {
			usleep(100);
			MPI_Barrier(MPI_COMM_WORLD);

			if(rank == 0){
				std::ostringstream tmpstream (std::ostringstream::app);
				tmpstream <<std::endl <<"----------------Starting iteration " <<it <<"----------------"<<std::endl;
				printStream(tmpstream.str());
			}

			//----Send coveringList content to all subprocesses
			// reset hitbuffer_sum (except covering) before sending to sub processes
			initbufftozero(&hitbuffer);
			if(rank==0){
				initbufftozero(&hitbuffer_sum);
			}

			// Send each coveringList entry one at a time
			for(unsigned int i=0; i<simHistory->numFacet;i++){
				MPI_Bcast(&simHistory->coveringList.currentList[i], 16, MPI::BYTE,0,MPI_COMM_WORLD);
			}

			/*
			for(int i=0; i<world_size;i++){
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank==i)
					simHistory->coveringList.printCurrent(std::cout,std::to_string(rank)+": coveringList at beginning of iteration");}
			*/
			MPI_Bcast(&simHistory->currentStep, 1, MPI::INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			// Set covering threshold (covering -covering/(size-1)). Iteration in subprocess is ended if this threshold is reached
			setCoveringThreshold(world_size, rank);

			UpdateSticking(); // Write sticking factor into sHandle for all subprocesses
			if(rank!=0){
				simHistory->updateHistory();// Write the current covering values from the simHistory to the sHandle and calculate stepSize (normal and outgassing).
			}
			else{
				simHistory->updateStepSize(); // Calculate stepSize (normal and outgassing) for this iteration
			}

			CalcTotalOutgassingWorker();// Calculate outgassing values for this iteration
			if(!UpdateDesorption()){// Write desorption into sHandle for all subprocesses
				// End simulation for very small desorption + outgassing particles
				if(rank==0) {
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream <<"Desorption smaller than 1E-50. Ending Simulation." <<std::endl;
					tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
					printStream(tmpstream.str());
				}
				break;
			}


			//----Simulation on subprocesses
			if (rank != 0) {
				sHandle->posCovering=true;//assumption that negative covering has been resolved before: implemented through covering threshold in sub processes and smallCovering workaoround
				checkSmallCovering(rank, &hitbuffer); // Calculate smallCoveringFactor for this iteration

				//Do the simulation
				bool eos; std::vector<int> facetNum;
				std::tie(eos, facetNum) = simulateSub2(&hitbuffer, rank, p->simulationTimeMS);

				std::ostringstream tmpstream (std::ostringstream::app);
				if (eos) {
					if(sHandle->posCovering)
						{tmpstream << "Maximum desorption reached." << std::endl;}
					else{
						for (uint facets=0; facets < facetNum.size(); facets++){
							tmpstream <<"Facet " <<facetNum[facets] <<" reached threshold " <<sHandle->coveringThreshold[facetNum[facets]] <<" for process " <<rank <<std::endl;
						}
					}
				} else {
					tmpstream << "Simulation for process " << rank << " for iteration " << it << " finished."<< std::endl;
				}
				printStream(tmpstream.str());
				MPI_Barrier(MPI_COMM_WORLD);
			}
			else{
				t0 = GetTick();
				checkSmallCovering(rank, &hitbuffer_sum); // Calculate smallCoveringFactor for this iteration
				MPI_Barrier(MPI_COMM_WORLD);
				t1 = GetTick();
				computationTime+=t1-t0; // Add calculation time of current iteration

			}

			//----iteratively add hitbuffer from subprocesses to main buffer
			for (int i = 1; i < world_size; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				if (rank == i) {
					//Process i sends hitbuffer to main process 0
					MPI_Send(hitbuffer.buff, hitbuffer.size, MPI::BYTE, 0, 0,MPI_COMM_WORLD);

					//Process i send flightTIme and number of particles to main process 0
					MPI_Send(&simHistory->flightTime, 1, MPI::DOUBLE, 0, 0,MPI_COMM_WORLD);
					MPI_Send(&simHistory->nParticles, 1, MPI::INT, 0, 0,MPI_COMM_WORLD);

				} else if (rank == 0) {
					//Main process 0 receives hitbuffer from Process i
					MPI_Recv(hitbuffer.buff, hitbuffer.size, MPI::BYTE, i, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					UpdateMCMainHits(&hitbuffer_sum, &hitbuffer, simHistory ,0);
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream << "Updated hitbuffer with process " << i <<std::endl;

					// Calculate flightTime and nParticles over all subprocesses -> These values are currently not used
					double old_flightTime=simHistory->flightTime;
					int old_nParticles = simHistory->nParticles; //These values are reset in UpdateCoveringPhys()
					MPI_Recv(&simHistory->flightTime, 1, MPI::DOUBLE, i, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(&simHistory->nParticles, 1, MPI::INT, i, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					simHistory->flightTime += old_flightTime;
					simHistory->nParticles+=old_nParticles;

					if(i==world_size-1){
						tmpstream <<std::endl << "flightTime " << simHistory->flightTime << std::endl;
						tmpstream << "nParticles " << simHistory->nParticles << std::endl <<std::endl;
					}
					printStream(tmpstream.str());
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);


			//----Update History: particle density, pressure, error, covering
			if (rank == 0) {
				UpdateParticleDensityAndPressure(&hitbuffer_sum);
				UpdateErrorMain(&hitbuffer_sum); // !! If order changes, adapt "time" entry in errorList !!
				UpdateCovering(&hitbuffer_sum);

				UpdateCoveringphys(&hitbuffer_sum, &hitbuffer);

				// Adapt size of history lists if p->histSize is exceeded
				if(p->histSize != std::numeric_limits<int>::infinity() && simHistory->coveringList.historyList.first.size() > uint(p->histSize+1)){
						simHistory->erase(1);
				}
				// print current coveringList
				simHistory->coveringList.print(p->outFile,"Accumulative covering after iteration "+std::to_string(it),p->histSize,true);

				simHistory->coveringList.updateStatistics(p->rollingWindowSize);

				currentRatio=double(simHistory->coveringList.getAverageStatistics(sHandle,true));
				simHistory->coveringList.printStatistics(p->outFile, "Rolling time window statistics over last "+std::to_string(p->rollingWindowSize)+" iterations. Mean ratio std/mean = "+std::to_string(currentRatio)+" with target ratio for convergence "+std::to_string(p->convergenceTarget), true);
			}

			if (rank == 0) {std::cout << "ending iteration " << it <<std::endl;}

			// Check if maximum simulation time is reached
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&simHistory->lastTime, 1, MPI::DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&currentRatio, 1, MPI::DOUBLE, 0, MPI_COMM_WORLD);
			if((int)(simHistory->lastTime+0.5) >= p->maxTimeS){
				if(rank==0) {
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream <<"Maximum simulation time reached: " <<simHistory->lastTime  <<" >= " <<p->maxTimeS <<std::endl;
					tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
					printStream(tmpstream.str());
				}
				break;
			} else if (currentRatio<=p->convergenceTarget){
				if(rank==0) {
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream <<"Simulation converged. Average ratio std/mean target reached: " <<currentRatio <<" <= "<<p->convergenceTarget <<std::endl;
					if(p->stopConverged)
						tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
					else
						tmpstream <<"Continue Simulation."<<std::endl;
					printStream(tmpstream.str());
				}
				if(p->stopConverged)
					break;
			}

		} else {
			std::cout << "Simulation time = 0.0 seconds. Nothing to do." << std::endl;
			break;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		//----Write simulation results to new buffer file. This has to be read in by Windows-Molflow.
		// Print history lists
		simHistory->print(true);

		if(p->saveResults){
			// Write result files
			simHistory->write(p->resultpath);

			// Export hitbuffer
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "Process 0 exporting final hitbuffer" << std::endl <<std::endl;
			printStream(tmpstream.str());
			exportBuff(p->resultbufferPath,&hitbuffer_sum);//export hitbuffer_sum
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//--delete buffers
	if (hitbuffer.buff != NULL) {
		delete[] hitbuffer.buff;
		hitbuffer.buff = NULL;
	}
	if (hitbuffer_sum.buff != NULL) {
		delete[] hitbuffer_sum.buff;
		hitbuffer_sum.buff = NULL;
	}
	if (loadbuffer.buff != NULL) {
		delete[] loadbuffer.buff;
		loadbuffer.buff = NULL;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		std::cout <<std::endl << "Closing MPI now." << std::endl;
	MPI_Finalize();
	if (rank == 0)
		std::cout << "Program finished." << std::endl;

	return 0;
}
