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
 * Main function
 */
#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include "Buffer.h"
#include <unistd.h>
#include "Simulation.h"
#include "SimulationLinux.h"

/*#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <cstring>
//#include <tuple>
#include "Types.h"
//#include "tuple.hpp"
#include "File.h"
//#includ <mpi.h>
#include <math.h>
//#include <malloc.h>
#include "MolFlow.h"
#include "Facet_shared.h"
#include "MolflowGeometry.h"
*/

/*
#include "GLApp/GLFileBox.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/MathTools.h"
#include "GLApp/GLMenuBar.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTextField.h"

#include "RecoveryDialog.h"
#include "direct.h"*/

/*
#include <vector>
#include <string>
//#include <io.h>
#include <thread>
#include <numeric> //std::iota

#include "Interface.h"
#include "Worker.h"
//#include "ImportDesorption.h"
//#include "TimeSettings.h"
//#include "Movement.h"
#include "FacetAdvParams.h"
#include "FacetDetails.h"
#include "Viewer3DSettings.h"
#include "TextureScaling.h"
//#include "GlobalSettings.h"
#include "ProfilePlotter.h"
*/
typedef void *HANDLE;

// Global process variables
Simulation* sHandle; //Global handle to simulation, one per subprocess
CoveringHistory* covhistory;

//This function checks if the correct number of arguments has been passed
//does not check their validity, e.g. right type such as double/string, correct filename, etc
bool parametercheck(int argc, char *argv[]) {
	int i;
	printf("argc: %d\n", argc);
	for (i = 0; i < argc; i++) {
		printf("argv[%d]: %s\n", i, argv[i]);
	}
	std::cout <<std::endl;
	if (argc < 5 || argc > 6) {
		std::cout
				<< "MolflowLinux requires 4 mandatory arguments and 1 optional argument."
				<< std::endl;
		std::cout << "Please pass these arguments to MolflowLinux:"
				<< std::endl;
		std::cout << "1. Name of load-buffer file to read in (e.g. loadbuffer)."
				<< std::endl;
		std::cout << "2. Name of hit-buffer file to read in (e.g. hitbuffer)."
				<< std::endl;
		std::cout
				<< "3. Choose a name for the buffer file to export the simulation results (e.g. resultbuffer)."
				<< std::endl;
		std::cout << "4. The total simulation time (e.g 2.5)." << std::endl;
		std::cout
				<< "5. [OPTIONAL] Simulation time unit (e.g. seconds, minutes, hours, days). Default set to seconds."
				<< std::endl;
		std::cout << "MolflowLinux is terminated now." << std::endl;
		return false;
	}
	return true;
}

int main(int argc, char *argv[]) {

	// Initialise buffer
	Databuff hitbuffer;
	hitbuffer.buff=NULL;

	Databuff hitbuffer_original;
	hitbuffer_original.buff=NULL;

	Databuff loadbuffer;
	loadbuffer.buff=NULL;

	// Init simulation time and unit
	double SimulationTime;
	int newsimutime;
	std::string unit;

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

	if (rank == 0) {
		// Parameter check for MolflowLinux
		if (!parametercheck(argc, argv)) {
			MPI_Finalize();
			return 0;
		}

		int a = 1;
		long b = 1;
		llong c = 1;
		float d = 1;
		double e = 1;
		long double f = 1;
		char g = 1;
		long long h = 1;
		std::cout<<"size of int a = "<< sizeof(a)<<std::endl;
		std::cout<<"size of long b = "<< sizeof(b)<<std::endl;
		std::cout<<"size of uint_64 c = "<< sizeof(c)<<std::endl;
		std::cout<<"size of float d = "<< sizeof(d)<<std::endl;
		std::cout<<"size of double e = "<< sizeof(e)<<std::endl;
		std::cout<<"size of long double f = "<< sizeof(f)<<std::endl;
		std::cout<<"size of char g = "<< sizeof(g)<<std::endl;
		std::cout<<"size of long long h = "<< sizeof(h)<<std::endl;



		//Read in buffer file (exported by Windows-Molflow). File given as first argument to main().
		importBuff(argv[1],&loadbuffer);
		importBuff(argv[2],&hitbuffer);

		//loadbuffer.importBuff(argv[1]);
		//hitbuffer.importBuff(argv[2]);
		//importBuff(argv[2], &hitbuffer_original); //TODO: copy content from pointer rather than read again?

		//TODO Copy constructor?
		hitbuffer_original.buff = new BYTE[hitbuffer.size];
		memcpy(hitbuffer_original.buff,hitbuffer.buff,hitbuffer.size);
		hitbuffer_original.size =hitbuffer.size;


		//hitbuffer_original=hitbuffer;
		exportBuff(argv[3],&hitbuffer_original);
		//hitbuffer_original.importBuff("/home/van/hitbuffertest");

		/*
		 * show informations about the loading
		 * just interesting for debugging => build some conditional (if debug, then show)?
		 * leave out or put in function?*/
		std::cout << "size of " << argv[2] << " = " << hitbuffer.size
				<< std::endl;
		std::cout << "size of " << argv[1] << " = " << loadbuffer.size
				<< std::endl;

		std::cout << "Buffers sent. Wait for a few seconds. " << std::endl<< std::endl;
	}

	// Send load-buffer to all other processes
	// Send size of buffer
	MPI_Bcast(&loadbuffer.size, sizeof(loadbuffer.size), MPI::BYTE, 0,
			MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// Allocate memory for buffer
	if (rank != 0) { /* do work in any remaining processes */
		loadbuffer.buff = new BYTE[loadbuffer.size];
	}

	// Send buffer content
	MPI_Bcast(loadbuffer.buff, loadbuffer.size, MPI::BYTE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//Send hit-buffer to all other processes.
	// Send size of buffer
	MPI_Bcast(&hitbuffer.size, sizeof(hitbuffer.size), MPI::BYTE, 0,
			MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// Allocate memory for buffer
	if (rank != 0) { /* do work in any remaining processes */
		hitbuffer.buff = new BYTE[hitbuffer.size];
	}

	for(int it=0;it<1;it++){ //TODO parameterübergabe, simulationszeit anpassen

		if(rank == 0){
		std::cout <<std::endl <<"Starting iteration " <<it <<std::endl;
		}

		// Send buffer content
		MPI_Bcast(hitbuffer.buff, hitbuffer.size, MPI::BYTE, 0, MPI_COMM_WORLD);

		// extract Simulation time and unit
		if(it==0){
			if (argc == 5)
				unit = "s";
			else
				unit = argv[5];
			SimulationTime = std::atof(argv[4]);

			//compute simulation time in seconds
			newsimutime = (int) (convertunit(SimulationTime, unit) + 0.5);
			if (rank == 0)
				std::cout << "Simulation time " << SimulationTime << unit
						<< " converted to " << newsimutime << "ms" << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Start of Simulation
		if (newsimutime != 0) {
			if(it==0){
				//Creates sHandle instance
				InitSimulation();

				//SetReady();

				// Load geometry from buffer to sHandle
				// nur wenn it!=0 !?!
				if (!LoadSimulation(&loadbuffer)) {
					std::cout << "Geometry not loaded." << std::endl;
					std::cout << "MolflowLinux is terminated now." << std::endl;
					MPI_Finalize();
					return 0;
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
			//Simulation on subprocesses
			if (rank != 0) {
				/* do work in any remaining processes */
				std::cout <<std::endl << "Process " << rank << " starting iteration "<< it <<" now."<< std::endl;

				//Do the simulation
				if (!simulateSub(&hitbuffer, rank, newsimutime)) {
					std::cout << "Maximum desorption reached." << std::endl;
				} else {
					std::cout << "Simulation for process " << rank << " finished."<< std::endl;
				}
			}

			//iteratively add hitbuffer from subprocesses
			for (int i = 1; i < world_size; i++) {
				MPI_Barrier(MPI_COMM_WORLD);
				if (rank == i) {
					//Process i sends hitbuffer to Process 0
					MPI_Send(hitbuffer.buff, hitbuffer.size, MPI::BYTE, 0, 0,MPI_COMM_WORLD);

				} else if (rank == 0) {
					// Process 0 receives hitbuffer from Process i
					//delete[] hitbuffer.buff; hitbuffer.buff = new BYTE[hitbuffer.size];
					MPI_Recv(hitbuffer.buff, hitbuffer.size, MPI::BYTE, i, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					//std::string exporthit = "/home/van/hitbuffer0" + std::to_string(i);
					//exportBuff(exporthit, &hitbuffer);
					//std::cout << "Received hitbuffer from process" <<i << std::endl;

					//sleep(1);
					UpdateMainHits(&hitbuffer_original, &hitbuffer, 0);
					std::cout << "Updated hitbuffer with process " << i <<std::endl
							<< std::endl;

					//std::string exporthit = "~/resultbuffer0" + std::to_string(i);
					//exportBuff(exporthit, &hitbuffer_original);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			if (rank == 0) {
				//hitbuffer.buff = new BYTE[hitbuffer_original.size];
				memcpy(hitbuffer.buff,hitbuffer_original.buff,hitbuffer_original.size);
				std::cout << "ending iteration " << it <<std::endl;
				//________________________________________________________________________

			}
		} else {
			std::cout << "Simulation time = 0.0 seconds. Nothing to do."
					<< std::endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		//Write simulation results to new buffer file. This has to be read in  by Windows-Molflow.
		std::cout << "Process 0 exporting final hitbuffer" << std::endl <<std::endl;
		exportBuff(argv[3],&hitbuffer_original);
	}

	if (hitbuffer.buff != NULL) {
		delete[] hitbuffer.buff;
		hitbuffer.buff = NULL;
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
