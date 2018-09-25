/* Linux-Molflow is based on Molflow+ 2.6.70*/


#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include "Buffer.h"
#include <unistd.h>

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
#include "direct.h"

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
    	  	f(argc < 5 || argc > 6){
    	  	    	  		std::cout << "MolflowLinux requires 4 mandatory arguments and 1 optional argument."<< std::endl;
    	  	    	  		std::cout << "Please pass these arguments to MolflowLinux:"<< std::endl;
    	  	    	  		std::cout << "1. Name of load-buffer file to read in (e.g. loadbuffer)." << std::endl;
    	  	    	  		std::cout << "2. Name of hit-buffer file to read in (e.g. hitbuffer)." << std::endl;
    	  	    	  		std::cout << "3. Choose a name for the buffer file to export the simulation results (e.g. resultbuffer)." << std::endl;
    	  	    	  		std::cout << "4. The total simulation time (e.g 2.5)." << std::endl;
    	  	    	  		std::cout << "5. [OPTIONAL] Simulation time unit (e.g. seconds, minutes, hours, days). Default set to seconds." << std::endl;
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

	  // Initialise 'load' buffer and 'hit' buffer
		Databuff databuffer;
		databuffer.buff = NULL;




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
    	  importBuff(argv[1], &databuffer);

    	  /*
    	   * show informations about the loading
    	   * just interesting for debugging => build some conditional (if debug, then show)?
    	   * leave out or put in function?*/
    	  std::cout << "size of " << argv[1] << " = " << databuffer.size <<std::endl;
    	  //std::cout << "size of databuffer = " << sizeof(databuffer) <<std::endl;
    	  /*std::cout << argv[1] << ": ";
    	  int i;
    	  for(i=0; i < databuffer.size; i++) {
    		  	  	  	if (i != (databuffer.size -1)){
    	      	  		std::cout<<databuffer.buff[i];}
    		  	  	  	else {std::cout<<databuffer.buff[i]<<std::endl;}
    	      	  		}*/
    	  //std::cout << argv[1] << ": " << databuffer.buff << std::endl;
    	  std::cout << "Buffer sent. Wait for 1 second. " <<std::endl;
    	  }

          //Send buffer to all other processes.
          MPI_Bcast(&databuffer.size, sizeof(databuffer.size), MPI::BYTE, 0, MPI_COMM_WORLD);
          //std::cout << "size of " << argv[1] << " = " << databuffer.size <<std::endl;
          sleep(1);

          if (rank !=0){ /* do work in any remaining processes */
              	  databuffer.buff = new BYTE[databuffer.size];
                }

          //MPI_Barrier(MPI_COMM_WORLD);
          MPI_Bcast(databuffer.buff, databuffer.size, MPI::BYTE, 0, MPI_COMM_WORLD);

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
    	 std::cout << "size of " << argv[1] << " = " << databuffer.size <<std::endl;
    	 //std::cout << argv[1] << std::endl;
    	 //char fileexport [] = "/smbhome/schoenmann/buffertest";
    	 //exportBuff(fileexport, &databuffer);
    	 /*int i;
    	 for(i=0; i < 10; i++) {
    	    	  	  	if (i != (databuffer.size -1)){
    	    	     	std::cout<<databuffer.buff[i];}
    	    	     	else {std::cout<<databuffer.buff[i]<<std::endl;}
    	 	 	 	 	}*/

    	  //InitSimulation(); //Creates sHandle instance, commented as error otherwise
    	  //SetReady(); // Rudi: Soll ich das übernehmen?

         /*Do the simulation*/
            /* int m = 0 
             * case PROCESS_RUN:
             *      //SetStatus(GetSimuStatus) deactivate HitUpdate function
             *      if(m == SimulationTime) {
             *          Tell process 0 via MPI PROCESS_DONE, MPI update Hits an den Hauptprozess
                        };
             *      m++
             */

      }


      sleep(1);
      //MPI_Barrier(MPI_COMM_WORLD);


      if(rank == 0) {
          	 //Write simulation results to new buffer file. This has to be read in  by Windows-Molflow.
    	     std::cout << "Hello! I'm head process "<< rank << std::endl;
          	 exportBuff(argv[2], &databuffer);
          	 // Build in safety check to not loosing simulation results, if the buffer export does not work?
          	 //delete[] databuffer.buff;
          	//std::cout <<"____________________________________________________________________________________________________" << std::endl;
            }

             MPI_Barrier(MPI_COMM_WORLD);
             delete[] databuffer.buff;




    MPI_Finalize();
    

    return 0;
}
