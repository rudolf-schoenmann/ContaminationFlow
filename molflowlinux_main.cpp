#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include "Buffer.h"
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

bool parametercheck(int argc, char *argv[])
      {
    	    int i;
    	  	printf("argc: %d\n", argc);
    	  	for(i=0; i < argc; i++) {
    	  		printf("argv[%d]: %s\n", i, argv[i]);
    	  		}
    	  	if(argc < 3){
    	  		std::cout << "Please pass 2 arguments to MolflowLinux:"<< std::endl;
    	  		std::cout << "1. Name of buffer file to read in." << std::endl;
    	  		std::cout << "2. Choose a name for the buffer file to export the simulation results." << std::endl;
    	  		std::cout << "MolflowLinux is terminated now." << std::endl;
    	  		return false;
    		 	}
    	  	return true;
      }



int main(int argc, char *argv[]) {

	  /*Parameters passed to main:
	   *
	   * 1. Name of buffer file to read in.
	   * 2. Choose a name for the buffer file to export the simulation results.
	   * 3. Simulation time.
	   * 4.
	   * 5.
	   *
       */

      /* Create child processes, each of which has its own variables.
       * From this point on, every process executes a separate copy
       * of this program.  Each process has a different process ID,
       * ranging from 0 to num_procs minus 1, and COPIES of all
       * variables defined in the program. No variables are shared.
       **/

	  // Initialise buffer
	  Databuff buff;
	  buff.buffer = NULL;



      MPI_Init(NULL, NULL);
      /* find out MY process ID, and how many processes were started. */

      // Get the number of processes
      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);

      // Get the rank of the process
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      if(rank == 0) {
    	  /* do some work as process 0 */
    	  std::cout << "Hello! I'm head process "<< rank << std::endl;
    	  //std::cout << "Number of started processes: "<< world_size << std::endl;

    	  // Parameter check for MolflowLinux
    	  if (!parametercheck(argc, argv))
    	  {		MPI_Finalize();
      	  	return 0;
    	  }

    	  //Read in buffer file (exported by Windows-Molflow). File given as first argument to main().
    	  importBuff(argv[1], &buff);

    	  /*
    	   * show informations about the loading
    	   * just interesting for debugging => build some conditional (if debug, then show)?
    	   * leave out or put in function?
    	  std::cout << "size of " << argv[1] << " = " << buff.size <<std::endl;
    	  //std::cout << "size of buff = " << sizeof(buff) <<std::endl;
    	  std::cout << argv[1] << ": ";
    	  int i;
    	  for(i=0; i < buff.size; i++) {
    		  	  	  	if (i != (buff.size -1)){
    	      	  		std::cout<<buff.buffer[i];}
    		  	  	  	else {std::cout<<buff.buffer[i]<<std::endl;}
    	      	  		}
    	  //std::cout << argv[1] << ": " << buff.buffer << std::endl;
    	  */
      	  }


      //Send buffer to all other processes.
      //MPI_Bcast(&buff, sizeof(buff), MPI::BYTE, 0, MPI_COMM_WORLD);



          /*Sharing buffer (Geometry and Parameters) with the other processes
             Send Buffer
              * Send SimulationTime
              * Tell the other processes via MPI, that they should execute COMMAND_LOAD (???)
              * If all Processes successfully executed COMMAND_LOAD, then tell them to execute COMMAND_START
          The other processes are simulating. I do nothing.
          After PROCESS_DONE Messages from the other processes, Receive Data from the other processes
          Sum up Data from the other processes
             First Step: Stationary Simulation => do nothing here
             Second Step: Time dependend Mode (Maybe copy Algorithm from Marton,  when ready)
                 Iterative Algorithm: Update Parameters (Is parallisation necessary/desirable for updating?) + Sharing new Parameters with the other processes
          write Result a bufferfile (or maybe??? in a file or .zip archive)*/


      if (rank != 0){
    	 /* do work in any remaining processes */
    	 std::cout << "Hello! I'm worker process "<< rank << std::endl;
    	 //std::cout << "size of " << argv[1] << " = " << buff.size <<std::endl;
    	 //std::cout << "Received buffer: " << buff.buffer << std::endl;

    	 //delete[] buff.buffer;
         //(Receive the Buffer)

    	 // Sketch of the worker algorithm:

    	  /* If COMMAND_LOAD is executed successfully, tell process 0 (via MPI)
          * 
          */
         /*Do the simulation*/
            /* int m = 0 
             * case PROCESS_RUN:
             *      //SetStatus(GetSimuStatus) deactivate HitUpdate function
             *      if(m == SimulationTime) {
             *          Tell process 0 via MPI PROCESS_DONE, MPI update Hits an den Hauptprozess
                        };
             *      m++
             */

    	 // Here is a simplified copy of the "Main loop" of the Molflow+ subprocess (in molflowSub.cpp) of the Windows application.

    	  /*
    	   *
    	   * Main loop

  	  	  	 while( !end ) {
	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("COMMAND: LOAD (%zd,%llu)\n",prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          //sHandle->desorptionLimit = prParam2; // 0 for endless
          SetReady();
        }
        break;

      case COMMAND_START:
        printf("COMMAND: START (%zd,%llu)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(prParam) )
            SetState(PROCESS_RUN,GetSimuStatus());
          else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;

      case COMMAND_PAUSE:
        printf("COMMAND: PAUSE (%zd,%llu)\n",prParam,prParam2);
        if( !sHandle->lastHitUpdateOK ) {
          // Last update not successful, retry with a longer timeout
			if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,60000);
        }
        SetReady();
        break;

      case COMMAND_EXIT:
        printf("COMMAND: EXIT (%zd,%llu)\n",prParam,prParam2);
        end = true;
        break;

      case PROCESS_RUN:
        SetStatus(GetSimuStatus()); //update hits only
        eos = SimulationRun();      // Run during 1 sec
		if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,20); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
        if(eos) {
          if( GetLocalState()!=PROCESS_ERROR ) {
            // Max desorption reached
            SetState(PROCESS_DONE,GetSimuStatus());
            printf("COMMAND: PROCESS_DONE (Max reached)\n");
          }
        }
        break;

      default:
        Sleep(WAITTIME);
        break;
    }
  }



    	   * */

      }
      /* Stop this process */

      if(rank == 0) {
          	 //___________________________________________________________________________________________________________________________________________________________
          	 //Write simulation results to new buffer file. This has to be read in  by Windows-Molflow.

          	 //char fileexport[] = "/smbhome/schoenmann/buffertest";
          	 exportBuff(argv[2], &buff);
          	 // Build in safety check to not loosing simulation results, if the buffer export does not work?
          	 delete[] buff.buffer;
          	std::cout <<"____________________________________________________________________________________________________" << std::endl;
            }


    MPI_Finalize();
    

    return 0;
}
