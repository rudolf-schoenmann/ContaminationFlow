#include <iostream>
#include <mpi.h>
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


int main(int argc, char **argv) {                           // Parameters for Main: File to load + How long should the simulation run (# of Hits)


    //char filename = argv[1];
    //int processNumber = static_cast<int>(argv[2]);
    //int SimulationTime = static_cast<int>(argv[3];
    //argv[4] should be used to load the buffer file.
  
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
      int world_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

      if( world_rank == 0 ) {
    	  std::cout << "Hello, world! Ich bin das molflowlinux_project!" << std::endl;
    	  std::cout << "Hello! I'm head process "<< world_rank << std::endl;
    	  std::cout << "Hello! Number of processes: "<< world_size << std::endl;
          /* do some work as process 0 */
          /*load the buffer (Buffer has to be exported in the WindowsMolflow first)*/
          /*Sharing Geometry and Parameters with the other processes*/
             /* Send Buffer
              * Send SimulationTime
              * Tell the other processes via MPI, that they should execute COMMAND_LOAD
              * If all Processes successfully executed COMMAND_LOAD, then tell them to execute COMMAND_START */
          /*The other processes are simulating. I do nothing.*/
          /*After PROCESS_DONE Messages from the other processes, Receive Data from the other processes*/
          /*Sum up Data from the other processes*/
             /*First Step: Stationary Simulation => do nothing here*/
             /*Second Step: Time dependend Mode (Maybe copy Algorithm from Marton,  when ready)*/
                 /*Iterative Algorithm: Update Parameters (Is parallisation necessary/desirable for updating?) + Sharing new Parameters with the other processes*/
          /*write Result a bufferfile (or maybe??? in a file or .zip archive)*/
    	  }
    
      else {
    	  std::cout << "Hello! I'm worker process "<< world_rank << std::endl;
         /* do work in any remaining processes */
         /* execute int main(...) from molflowSub.cpp
         (Receive the Buffer)*/

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

    MPI_Finalize();
    
    
    return 0;
}
