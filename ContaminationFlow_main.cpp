/*
Program:     ContaminationFlow
Description: Monte Carlo simulator for satellite contanimation studies
Authors:     Rudolf Schönmann / Hoai My Van
Copyright:   TU Munich
Forked from: Molflow (CERN) (https://cern.ch/molflow)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.B

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/


#include <iostream>
#include <string>
#include <mpi.h>
#include "Buffer.h"
#include <unistd.h>
#include "Simulation.h"
#include "SimulationLinux.h"

//typedef void *HANDLE;

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
				std::cout << "ContaminationFlowLinux requires 4 mandatory arguments and 1 optional argument."<< std::endl;
				std::cout << "Please pass these arguments to ContaminationFlowLinux:"<< std::endl;
				std::cout << "1. Name of load-buffer file to read in (e.g. loadbuffer)."<< std::endl;
				std::cout << "2. Name of hit-buffer file to read in (e.g. hitbuffer)."<< std::endl;
				std::cout << "3. Save results: 0=false, 1=true."<< std::endl;
				std::cout << "4. The total simulation time (e.g 2.5)." << std::endl;
				std::cout << "5. [OPTIONAL] Simulation time unit (e.g. seconds, minutes, hours, days). Default set to seconds." << std::endl;
				std::cout << "ContaminationFlowLinux is terminated now." << std::endl;
			}
		}
		else{ //Read input arguments
			if(rank==0){std::cout<<"Read arguments" <<std::endl;}

			if(!checkReadable(argv[1],rank)||!checkReadable(argv[2],rank)||std::atof(argv[4])<0.0) // check if parameters are feasible
				{return false;}
			p->readArg(argc, argv, rank);
			return true;
		}
	}
	else if(argc<4 && argc>1){ // Read input file
		if(checkReadable(argv[1],rank)){
			bool valid = p->readInputfile(argv[1],rank, argc==3?(int)std::atof(argv[2]):1);
			if((!p->doCoveringFile&&!checkReadable(p->hitbufferPath,rank))||!checkReadable(p->loadbufferPath,rank)||(p->doCoveringFile&&!checkReadable(p->coveringPath,rank))){return false;}
			return valid;}
		}
	return false;
	}

bool loadAndCheckSHandle(int rank, Databuff* hitbuffer, Databuff* loadbuffer){
	bool valid=true;

	// Load geometry from buffer to sHandle
	if (!LoadSimulation(loadbuffer)) { // Check if geometry in loadbuffer can be loaded
		if(rank==0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "Geometry not loaded." << std::endl;
			tmpstream << "ContaminationFlowLinux is terminated now." << std::endl;
			printStream(tmpstream.str());
		}
		valid=false;
	}
	else{
		p->particleDia = sHandle->wp.gasDiameter;
		if(rank == 0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"Particle diameter from input file is overwritten by value imported from loadbuffer! particleDia = "<< p->particleDia <<std::endl;
			printStream(tmpstream.str());
		}
	}

	size_t hitsize=sHandle->GetHitsSize();
	if(p->doCoveringFile){
		initBuffSize(hitbuffer,hitsize);
		initbufftozero(hitbuffer);
	}

	// Check for inconsistent hitbuffer size
	if(hitsize!=(unsigned int)hitbuffer->size){
		if(rank==0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "Hitbuffer size not correctly calculated." << std::endl;
			tmpstream << "ContaminationFlowLinux is terminated now." << std::endl;
			printStream(tmpstream.str());
		}
		valid=false;
	}

	// Check if there are zero moments
	if(sHandle->moments.size()){
		if(rank==0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "Number of moments " << sHandle->moments.size()<<" > 0. ContaminationFlow only works with 0 moments." << std::endl;
			printStream(tmpstream.str());
		}
		valid=false;
	}
	//Check for activation of 'Calculate constant flow':
	/* This has to activated otherwise the 'latest moment' can somehow be earlier than the end of the iteration
	 * time step. Then, the test-particle would be traced until the end of the iteration but to the 'latest moment'.
	 */
	if(!sHandle->wp.calcConstantFlow){
		sHandle->wp.calcConstantFlow = true;
		if(rank == 0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "'calcConstantFlow' is automatically enabled." << std::endl;
			tmpstream << "It was not checked when setting up the simulation." << std::endl;
			printStream(tmpstream.str());
		}
	}


	// Check for two sided facet with opacity
	for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if(f.sh.is2sided && f.sh.opacity>0.0){
				if(rank==0){
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream << "There is a two sided facet with opacity. This might cause undefined behavior." << std::endl;
					tmpstream << "ContaminationFlowLinux is terminated now." << std::endl;
					printStream(tmpstream.str());
				}
				valid=false;
			}
		}
	}

	//Read coveringFile
	if(p->doCoveringFile){
		if(!readCovering(hitbuffer,p->coveringPath,rank)){
			if(rank==0){
				std::ostringstream tmpstream (std::ostringstream::app);
				tmpstream << "Invalid covering file "<<home_to_tilde(p->coveringPath) << std::endl;
				tmpstream << "ContaminationFlowLinux is terminated now." << std::endl;
				printStream(tmpstream.str());
			}
			valid=false;
		}
	}
	return valid;
}

void clearAll(Databuff* hitbuffer, Databuff* hitbuffer_sum, Databuff* loadbuffer){
	//--delete buffers
	if (hitbuffer->buff != NULL) {
		delete[] hitbuffer->buff;
		hitbuffer->buff = NULL;
	}
	if (hitbuffer_sum->buff != NULL) {
		delete[] hitbuffer_sum->buff;
		hitbuffer_sum->buff = NULL;
	}
	if (loadbuffer->buff != NULL) {
		delete[] loadbuffer->buff;
		loadbuffer->buff = NULL;
	}
	if (simHistory!=NULL){
		delete simHistory;
		simHistory=NULL;
	}
	if(p!=NULL){
		delete p;
		p=NULL;
	}
	if(sHandle!=NULL){
		delete sHandle;
		sHandle=NULL;
	}
}

//-----------------------------------------------------------
//Main Function
int main(int argc, char *argv[]) {

	std::cout << "Hello!" << std::endl;
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
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"Minimum number of 2 processes needed. Currently "<<world_size <<std::endl;
			printStream(tmpstream.str());
		}
		clearAll(&hitbuffer, &hitbuffer_sum,&loadbuffer);
		MPI_Finalize();
		return 0;
	}

//---- Initialize ProblemDef that defines relevant simulation parameters for ContaminationFlowLinux
	p = new ProblemDef();
	// Check parameters (defined by command line arguments or input file) to be written to ProblemDef p

	if (!parametercheck(argc, argv,p,rank)) {
		if (rank == 0){
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream <<"Required parameters cannot be read." <<std::endl;
			tmpstream <<"Check command line arguments and/or parameters in input file." <<std::endl;
			tmpstream <<"Ending Simulation." <<std::endl;
			printStream(tmpstream.str());
		}
		clearAll(&hitbuffer, &hitbuffer_sum,&loadbuffer);
		MPI_Finalize();
		return 0;
	}

	if (rank == 0) {
		std::ostringstream tmpstream (std::ostringstream::app);
		p->printInputfile(tmpstream);
		printStream(tmpstream.str(),false);

		//Read in buffer files
		importBuff(p->loadbufferPath,&loadbuffer);
		if(!p->doCoveringFile){
			importBuff(p->hitbufferPath,&hitbuffer);
		}

		std::cout << "Sending buffers. Wait for a few seconds. " << std::endl;
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
	if(!p->doCoveringFile){
		// Send size of hitbuffer from main process to subprocesses in subprocesses
		MPI_Bcast(&hitbuffer.size, sizeof(hitbuffer.size), MPI::BYTE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		// Allocate memory for hitbuffer
		if (rank != 0) {
			hitbuffer.buff = new BYTE[hitbuffer.size];
		}
	}

//----Create Simulation handle and preprocess hitbuffer
	if (p->simulationTimeMS != 0) {
		//Creates sHandle instance for process 0 and all subprocesses (before the first iteration step starts)
		InitSimulation();

		// Load geometry from buffer to sHandle and check for invalid values
		if(!loadAndCheckSHandle(rank,&hitbuffer,&loadbuffer)){
			MPI_Barrier(MPI_COMM_WORLD);
			clearAll(&hitbuffer, &hitbuffer_sum,&loadbuffer);
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
			hitbuffer_sum.size = hitbuffer.size;
			
			// Initialize simHistory for main process from hitbuffer
			//simHistory contains the relevant results/quantities of the simulation. E.g., covering history, total simulated time, etc.
			simHistory = new SimulationHistory (&hitbuffer, world_size);
			std::cout <<p->focusGroup.second.size() <<" facet(s) in focusGroup"<<std::endl;
			for(int grpidx:p->focusGroup.second){
				std::cout <<"\t"<<grpidx;
			}
			std::cout <<std::endl;
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
	std::string monitoredFacets=(!p->doFocusGroupOnly||p->focusGroup.second.size()==simHistory->numFacet)?"all":"selected";
	while(true){
		it++;
		// Start of Simulation
		MPI_Barrier(MPI_COMM_WORLD);
		if (p->simulationTimeMS != 0) {
			for (;simHistory->pcStep <= (p->usePCMethod?1:0); simHistory->pcStep += 1) { // Predictor-corrector loop
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 0){
					std::ostringstream tmpstream (std::ostringstream::app);
					if (!p->usePCMethod)
						tmpstream <<std::endl <<"----------------Starting iteration " <<it <<"----------------"<<std::endl;
					else if (simHistory->pcStep == 0)
						tmpstream <<std::endl <<"----------------Starting predictor step of iteration " <<it <<"----------------"<<std::endl;
					else
						tmpstream <<std::endl <<"----------------Starting corrector step of iteration " <<it <<"----------------"<<std::endl;
					simHistory->coveringList.printCurrent(tmpstream, "coveringList.currentList: "); // (Berke): Will be removed later on
					simHistory->coveringList.printPredict(tmpstream, "coveringList.predictList: "); // (Berke): Will be removed later on
					printStream(tmpstream.str());
				}
				//---- Reset buffers and send coveringList content to all subprocesses
				// reset hitbuffer_sum (except covering) before sending to sub processes
				initbufftozero(&hitbuffer);
				if(rank==0){
					initbufftozero(&hitbuffer_sum);
				}
				// Send coveringList to subprocesses
				MPI_Bcast(&(simHistory->coveringList.currentList.front()), simHistory->coveringList.currentList.size()*16, MPI::BYTE,0,MPI_COMM_WORLD);
				// Send predictList to subprocess
				MPI_Bcast(&(simHistory->coveringList.predictList.front()), simHistory->coveringList.predictList.size()*16, MPI::BYTE,0,MPI_COMM_WORLD);
				// Send currentStep -> used to calculate stepSize
				MPI_Bcast(&simHistory->currentStep, 1, MPI::INT, 0, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				// Set covering threshold (covering -covering/(size-1)). Iteration in subprocess is ended if this threshold is reached
				setCoveringThreshold(world_size, rank);

				if(rank!=0){
					simHistory->updateHistory();// Write the current covering values from the simHistory to the sHandle and calculate stepSize (normal and outgassing).
				}
				else{
					if (simHistory->pcStep == 0)
						simHistory->updateStepSize(); // Calculate stepSize (normal and outgassing) for this iteration
				}

				if (simHistory->pcStep == 0) {
					// Not calculating these values again in corrector step
					UpdateSticking(); // Write sticking factor into sHandle for all subprocesses
					CalcTotalOutgassingWorker();// Calculate outgassing values for this iteration
					if(!UpdateDesorption()){// Write desorption into sHandle for all subprocesses
						// End simulation for very small desorbed + outgassing particles
						if(rank==0) {
							std::ostringstream tmpstream (std::ostringstream::app);
							tmpstream <<"Desorption smaller than 1E-50. Ending Simulation (is deactivated)." <<std::endl;
							tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
							printStream(tmpstream.str());
						}
						//break; Do not end the simulation
					}
				}

				//----Simulation on subprocesses
				if (rank != 0) {
					sHandle->posCovering=true;//assumption that negative covering has been resolved before: implemented through covering threshold in sub processes and smallCovering workaoround
					checkSmallCovering(rank, &hitbuffer); // Calculate smallCoveringFactor for this iteration

					MPI_Barrier(MPI_COMM_WORLD);
					//Do the simulation
					bool eos; std::vector<int> facetNum;
					std::tie(eos, facetNum) = simulateSub2(&hitbuffer, rank, p->simulationTimeMS); //coveringList, hitBuffer covering aynı değerde

					std::ostringstream tmpstream (std::ostringstream::app);
					if (eos) {
						tmpstream << "Iteration ended early for process "<<rank <<"."<<std::endl;
						if(sHandle->posCovering)
							{tmpstream << "Maximum desorption reached for process "<<rank << "."<< std::endl;}
						else{
							for (uint facets=0; facets < facetNum.size(); facets++){
								tmpstream <<"Facet " <<facetNum[facets] <<" reached threshold " <<sHandle->coveringThreshold[facetNum[facets]] <<" for process " <<rank <<"."<<std::endl;
							}
						}
					} else {
						if (!p->usePCMethod)
							tmpstream << "Simulation for process " << rank << " for iteration " << it << " finished."<< std::endl;
						else if (simHistory->pcStep == 0)
							tmpstream << "Simulation for process " << rank << " for predictor step of iteration " << it << " finished."<< std::endl;
						else
							tmpstream << "Simulation for process " << rank << " for corrector step " << simHistory->pcStep << " of iteration " << it << " finished."<< std::endl;
					}
					tmpstream <<std::endl;
					printStream(tmpstream.str());
					MPI_Barrier(MPI_COMM_WORLD);
				}
				else{
					t0 = GetTick();
					checkSmallCovering(rank, &hitbuffer_sum); // Calculate smallCoveringFactor for this iteration
					for (int s = 0; s < (int)sHandle->sh.nbSuper; s++) {
						for (SubprocessFacet& f : sHandle->structures[s].facets) {
							calcStartTime(&f,true,true);
						}
					}
					MPI_Barrier(MPI_COMM_WORLD);
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
						/*if(i==world_size-1){
							tmpstream <<std::endl << "flightTime " << simHistory->flightTime << std::endl;
							tmpstream << "nParticles " << simHistory->nParticles << std::endl <<std::endl;
						}*/
						printStream(tmpstream.str());
					}
				}

				MPI_Barrier(MPI_COMM_WORLD);

				//----Update History: particle density, pressure, error, covering
				if (rank == 0) {
					if (simHistory->pcStep == (p->usePCMethod?1:0)) {
						UpdateParticleDensityAndPressure(&hitbuffer_sum); // !! If order changes, adapt "time" entry in pressure/density lists !!
						UpdateErrorMain(&hitbuffer_sum); // !! If order changes, adapt "time" entry in errorLists !!
					}
					UpdateCovering(&hitbuffer_sum); // Calculate real covering after iteration
					UpdateCoveringphys(&hitbuffer_sum, &hitbuffer); // Update real covering in buffers
					if (simHistory->pcStep == (p->usePCMethod?1:0)) {
						// Adapt size of history lists if p->histSize is exceeded
						if(p->histSize != std::numeric_limits<int>::infinity() && simHistory->coveringList.historyList.first.size() > uint(p->histSize+1)){
								simHistory->erase(1);
						}
						// print current coveringList
						std::ostringstream tmpstream (std::ostringstream::app);
						simHistory->coveringList.print(tmpstream,"Accumulative covering after iteration "+std::to_string(it),p->histSize);

						// Calculate and print statistics
						simHistory->coveringList.updateStatistics(p->rollingWindowSize);

						currentRatio=double(simHistory->coveringList.getAverageStatistics(sHandle,true, p->doFocusGroupOnly,p->focusGroup.second));
						//tmpstream <<"Rolling time window statistics over last "+std::to_string(p->rollingWindowSize)+" iterations. Mean ratio std/mean = "+std::to_string(double(simHistory->coveringList.getAverageStatistics(sHandle,true, !p->doFocusGroupOnly,p->focusGroup.second)))+" with target ratio for convergence "+std::to_string(p->convergenceTarget)<<std::endl;
						simHistory->coveringList.printStatistics(tmpstream, "Rolling time window statistics over last "+std::to_string(p->rollingWindowSize)+" iterations for "+monitoredFacets+" facets. Mean ratio std/mean = "+std::to_string(currentRatio)+" with target ratio for convergence "+std::to_string(p->convergenceTarget));
						printStream(tmpstream.str());
					}
				}
				if (rank == 0 && simHistory->pcStep == 0 && p->usePCMethod) {std::cout << "ending prediction step " <<std::endl;}
				else if (rank == 0 && simHistory->pcStep == 1 && p->usePCMethod) {std::cout << "ending correction step " <<std::endl;} 
			} //End of predictor-corrector loop

			if (rank == 0) {
				std::cout << "ending iteration "  << it <<std::endl;
				simHistory->coveringList.predictList.clear(); // (Berke): Will be removed later on
				simHistory->coveringList.initPredict(simHistory->numFacet); // (Berke): Will be removed later on
			}
			simHistory->pcStep = 0;

			// Check if simulation has ended
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&simHistory->lastTime, 1, MPI::DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&currentRatio, 1, MPI::DOUBLE, 0, MPI_COMM_WORLD);
			//if((int)(simHistory->lastTime+0.5) >= p->maxTimeS){ // Maximum simulated time reached
			if(simHistory->lastTime >= 1.000001*p->maxTimeS){ // Maximum simulated time reached (at least 99% of it)
				if(rank==0) {
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream <<"Maximum simulated time reached: " <<simHistory->lastTime  <<" >= " <<p->maxTimeS <<std::endl;
					tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
					printStream(tmpstream.str());
				}
				break;
			} else if (currentRatio<=p->convergenceTarget){ // Simulation has converged
				if(rank==0) {
					std::ostringstream tmpstream (std::ostringstream::app);
					tmpstream <<"Simulation converged. Average ratio std/mean target reached: " <<currentRatio <<" <= "<<p->convergenceTarget <<std::endl;
					if(p->stopConverged && simHistory->lastTime>=p->convergenceTime)
						tmpstream <<"Computation Time (Simulation only): " <<computationTime/1000.0<<"s = "<<simHistory->coveringList.convertTime(computationTime/1000.0) <<std::endl;
					else if(!p->stopConverged)
						tmpstream <<"p->stopConverged=false. Continue Simulation."<<std::endl;
					else
						tmpstream <<"Convergence time not reached: "<<simHistory->lastTime <<"s < " <<p->convergenceTime <<"s. Continue Simulation."<<std::endl;
					printStream(tmpstream.str());
				}
				if(p->stopConverged && simHistory->lastTime>=p->convergenceTime) // End simulation only if "allowed" and convergenceTime reached
					break;
			}
			MPI_Barrier(MPI_COMM_WORLD); // 
		} else {
			std::cout << "Simulation time = 0.0 seconds. Nothing to do." << std::endl;
			break;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		//----Write simulation results to new buffer file. This has to be read in by Windows-Molflow.
		// Print history lists
		simHistory->print();

		if(p->saveResults){
			// Write result files
			std::ostringstream tmpstream (std::ostringstream::app);
			tmpstream << "Process 0 exporting simHistory" << std::endl <<std::endl;
			printStream(tmpstream.str());

			simHistory->write(p->resultPath);
			//exportBuff(p->resultbufferPath,&hitbuffer_sum);//export hitbuffer_sum
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// Clear memory
	clearAll(&hitbuffer, &hitbuffer_sum,&loadbuffer);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		std::cout <<std::endl << "Closing MPI now." << std::endl;
	MPI_Finalize();
	if (rank == 0)
		std::cout << "Program finished." << std::endl;

	return 0;
}
