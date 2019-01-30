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
 * This file does nothing. molflowlinux_main.cpp was derived from this.
 */
#define NOMINMAX
//#include <windows.h>
//#include <tlhelp32.h>

//#include <iostream>

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

//#include "Simulation.h" // already included in SimulationLinux.h
#include "SimulationLinux.h"

#ifdef WIN
//#include <Process.h> // For _getpid()
#endif

typedef void *HANDLE;

// Global process variables
extern Simulation* sHandle; //Global handle to simulation, one per subprocess

#define WAITTIME    100  // Answer in STOP mode
//#define TIMEOUT     300  // Process kills itself after no heartbeat (seconds)

//static Dataport *dpControl=NULL;
//static Dataport *dpHit=NULL;
//static Dataport *dpLog = NULL;

static Databuff *hitbuffer=NULL;
static Databuff *controlbuffer=NULL;
static Databuff *loadbuffer=NULL;

//static int       noHeartBeatSince;
static int       prIdx;
static size_t       prState;
static size_t       prParam;
static llong     prParam2;
static DWORD     hostProcessId;
//static float       heartBeat;
//static HANDLE    masterHandle;
static char      ctrlDpName[32];
static char      loadDpName[32];
static char		 logDpName[32];
static char      hitsDpName[32];

bool end = false;
bool IsProcessRunning(DWORD pid);

void GetState() {
  prState = PROCESS_READY;
  prParam = 0;

  if( controlbuffer->buff!=NULL ) {
    SHCONTROL *master = (SHCONTROL *)controlbuffer->buff;
    prState = master->states[prIdx];
    prParam = master->cmdParam[prIdx];
    prParam2 = master->cmdParam2[prIdx];
    master->cmdParam[prIdx] = 0;
    master->cmdParam2[prIdx] = 0;

    //ReleaseDataport(dpControl);
	
	if (!IsProcessRunning(hostProcessId)) {
		printf("Host synrad.exe (process id %d) not running. Closing.",(int)hostProcessId);
		SetErrorSub("Host synrad.exe not running. Closing subprocess.");
		end = true;
	}
  } else {
	  printf("Subprocess couldn't connect to Molflow.\n");
	  SetErrorSub("No connection to main program. Closing subprocess.");
	  sleep(5000);
	  end = true;
  }
}

size_t GetLocalState() {
  return prState;
}

void SetState(size_t state,const char *status,bool changeState, bool changeStatus) {

	prState = state;
	if (changeState) printf("\n setstate %zd \n",state);
	if( controlbuffer->buff!=NULL) {
		SHCONTROL *master = (SHCONTROL *)controlbuffer->buff;
		if (changeState) master->states[prIdx] = state;
		if (changeStatus) {
			strncpy(master->statusStr[prIdx], status, 127);
			master->statusStr[prIdx][127] = 0;
		}
		/*if( state==PROCESS_RUNAC ) {
			master->cmdParam[prIdx] = sHandle->prgAC;
		}*/
		//ReleaseDataport(dpControl);
	}

}

void SetErrorSub(const char *message) {

  printf("Error: %s\n",message);
  SetState(PROCESS_ERROR,message);

}

char *GetSimuStatus() {

  static char ret[128];
  llong count = sHandle->totalDesorbed;
  llong max   = sHandle->ontheflyParams.desorptionLimit/sHandle->ontheflyParams.nbProcess;

  /* (Rudi) No AC mode in MolflowLinux
  if( GetLocalState()==PROCESS_RUNAC ) sHandle->wp.sMode = AC_MODE;
  */

  switch(sHandle->wp.sMode) {

    case MC_MODE:
      if( max!=0 ) {
        double percent = (double)(count)*100.0 / (double)(max);
        sprintf(ret,"(%s) MC %I64d/%I64d (%.1f%%)",sHandle->sh.name.c_str(),(int)count,(int)max,percent);
      } else {
        sprintf(ret,"(%s) MC %I64d",sHandle->sh.name.c_str(),(int)count);
      }
      break;

    /*case AC_MODE: //(my)no AC_MODE
      if( sHandle->prgAC<100 ) {
          sprintf(ret,"(%s) AC (%zdx%zd) (%zd%%)",sHandle->sh.name.c_str(),
                      sHandle->nbAC,sHandle->nbAC,sHandle->prgAC);
      } else {
        if( max!=0 ) {
          double percent = (double)(count)*100.0 / (double)(max);
          sprintf(ret,"(%s) AC (%zdx%zd) %I64d/%I64d (%.1f%%)",sHandle->sh.name.c_str(),
                      sHandle->nbAC,sHandle->nbAC,count,max,percent);
        } else {
          sprintf(ret,"(%s) AC (%zdx%zd) %I64d",sHandle->sh.name.c_str(),sHandle->nbAC,
                      sHandle->nbAC,count);
        }
      }
      break;*/

    default:
    	std::cout << "Simulation mode is not valid" <<std::endl;
    	break;

  }

  return ret;

}

void SetReady() {

  if(sHandle->loadOK)
    SetState(PROCESS_READY,GetSimuStatus());
  else
    SetState(PROCESS_READY,"(No geometry)");

}

void SetStatus(char *status) {

  if( controlbuffer->buff !=NULL) {
    SHCONTROL *master = (SHCONTROL *)controlbuffer->buff;
	strncpy(master->statusStr[prIdx], status, 127);
	master->statusStr[prIdx][127] = 0;
    //ReleaseDataport(dpControl);
  }

}
/*(My) not needed -> case command_loadac
void LoadAC() {

  Dataport *loader;
  SHELEM_OLD *map;

  if( !sHandle->loadOK ) {
    SetErrorSub("No geometry loaded");
    return;
  }

  ClearACMatrix();

  // Load mesh
  loader = OpenDataport(loadDpName,prParam);
  if( !loader ) {
    char err[512];
    sprintf(err,"Failed to open 'loader' dataport %s (%zd Bytes)",loadDpName, prParam);
    SetErrorSub(err);
    return;
  }

  // Connect the dataport
  if( !AccessDataport(loader) ) {
    SetErrorSub("Failed to connect to DP");
    return;
  }
  printf("Connected to %s\n",loadDpName);

  map = (SHELEM_OLD *)malloc( prParam );
  memcpy(map,loader->buff,prParam);
  ReleaseDataport(loader);
  CLOSEDP(loader);

  SetState(PROCESS_RUNAC,GetSimuStatus());
  ComputeACMatrix(map);
  if( GetLocalState()!=PROCESS_ERROR ) {
    char s[128];
    if( sHandle->prgAC==100 ) {
      sprintf(s,"AC matrix calcultion ok (%.3f s)",sHandle->calcACTime);
    } else {
      sprintf(s,"AC matrix calcultion interrupted");
    }
    SetState(PROCESS_DONE,s);
  }
  free(map);

}
*/
void Load() {

  //Databuff *databuffer; //loadbuffer already initialized/imported in main
  //databuffer=NULL; //(my) change to import databuffer
  //size_t hSize;

  // Load geometry
  //databuffer = OpenDataport(loadDpName,prParam); //Eigentlich brauch ich diese Funktion nicht. Vielleicht ne Sicherungskopie vom Buffer anlegen?
  /*// (Rudi) Don't need that.
  if( !loader ) {
    char err[512];
    sprintf(err,"Failed to connect to 'loader' dataport %s (%zd Bytes)",loadDpName, prParam);
    SetErrorSub(err);
    return;
  }
  
  //TODO Load geometry from buffer here?

  printf("Connected to %s\n",loadDpName);
*/
  if( !LoadSimulation(loadbuffer) ) {
    //CLOSEDP(loader); Rudi) Don't need that.
    return;
  }
  //CLOSEDP(loader); Rudi) Don't need that.

  /* (Rudi) Don't need that now. Maybe some other kind of log later.
  //Connect to log dataport
  if (sHandle->ontheflyParams.enableLogging) {
	  dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
	  if (!dpLog) {
		  char err[512];
		  sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		  SetErrorSub(err);
		  sHandle->loadOK = false;
		  return;
	  }
	  // *((size_t*)dpLog->buff) = 0; //Autofill with 0. Besides, we don't write without access!
  }
  */

  /* (Rudi) Ich muss eigentlich keine Hits laden, wenn keine da.
  // Connect to hit dataport
  hSize = GetHitsSize();
  dpHit = OpenDataport(hitsDpName,hSize);
  if( !dpHit ) {
	  char err[512];
	  sprintf(err, "Failed to connect to 'hits' dataport (%zd Bytes)", hSize);
	  SetErrorSub(err);
	sHandle->loadOK = false;
    return;
  }

  printf("Connected to %s (%zd bytes)\n",hitsDpName,hSize);
  */
}
/*(my) not needed -> case command_updateparams
bool UpdateParams() {

	// Load geometry
	Dataport *loader = OpenDataport(loadDpName, prParam);
	if (!loader) {
		char err[512];
		sprintf(err, "Failed to connect to 'loader' dataport %s (%zd Bytes)", loadDpName, prParam);
		SetErrorSub(err);
		return false;
	}
	printf("Connected to %s\n", loadDpName);

	bool result = UpdateOntheflySimuParams(loader);
	CLOSEDP(loader);

	if (sHandle->ontheflyParams.enableLogging) {
		dpLog = OpenDataport(logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
		if (!dpLog) {
			char err[512];
			sprintf(err, "Failed to connect to 'dpLog' dataport %s (%zd Bytes)", logDpName, sizeof(size_t) + sHandle->ontheflyParams.logLimit * sizeof(ParticleLoggerItem));
			SetErrorSub(err);
			return false;
		}
		// *((size_t*)dpLog->buff) = 0; //Autofill with 0, besides we would need access first
	}
	sHandle->tmpParticleLog.clear();
	sHandle->tmpParticleLog.shrink_to_fit();
	if (sHandle->ontheflyParams.enableLogging) sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);

	return result;
}*/


int submain(int argc,char* argv[])
{
  bool eos = false;

  if(argc!=3) {
    printf("Usage: molflowSub peerId index\n");
    return 1;
  }

  hostProcessId=atoi(argv[1]);
  prIdx = atoi(argv[2]);

  sprintf(ctrlDpName,"MFLWCTRL%s",argv[1]);
  sprintf(loadDpName,"MFLWLOAD%s",argv[1]);
  sprintf(hitsDpName,"MFLWHITS%s",argv[1]);
  sprintf(logDpName, "MFLWLOG%s", argv[1]);

  hitbuffer = new Databuff; hitbuffer->buff=NULL;
  controlbuffer = new Databuff;controlbuffer->buff=NULL;
  loadbuffer = new Databuff;loadbuffer->buff=NULL;

  //dpControl = OpenDataport(ctrlDpName/*,sizeof(SHCONTROL)*/);
  //TODO Import control,hit,load buffer

  if(controlbuffer->buff==NULL ) {
    printf("Usage: Cannot connect to MFLWCTRL%s\n",argv[1]);
    return 1;
  }

  printf("Connected to %s (%zd bytes), molflowSub.exe #%d\n",ctrlDpName,sizeof(SHCONTROL),prIdx);

  InitSimulation(); //Creates sHandle instance

  // Sub process ready
  SetReady();

  // Main loop
  while( !end ) {
	GetState();
    switch(prState) {

      case COMMAND_LOAD:
        printf("COMMAND: LOAD (%zd,%lu)\n",prParam,prParam2);
        Load();
        if( sHandle->loadOK ) {
          //sHandle->desorptionLimit = prParam2; // 0 for endless
          SetReady();
        }
        break;
/*(my) case not needed
      case COMMAND_LOADAC:
        printf("COMMAND: LOADAC (%zd)\n",prParam);
        LoadAC();
        break;*/
/*(my) case not needed
	  case COMMAND_UPDATEPARAMS:
		  printf("COMMAND: UPDATEPARAMS (%zd,%I64d)\n", prParam, prParam2);
		  if (UpdateParams()) {
			  SetState(prParam, GetSimuStatus());
		  }
		  break;*/
/*(my) case not needed
	  case COMMAND_RELEASEDPLOG:
		  printf("COMMAND: UPDATEPARAMS (%zd,%I64d)\n", prParam, prParam2);
		  CLOSEDP(dpLog);
		  SetState(prParam, GetSimuStatus());
		  break;*/

      case COMMAND_START:
        printf("COMMAND: START (%zd,%lu)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(hitbuffer))
            SetState(PROCESS_RUN,GetSimuStatus());
          else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;

      case COMMAND_PAUSE:
        printf("COMMAND: PAUSE (%zd,%lu)\n",prParam,prParam2);
        /* (Rudi) Don't need that. Replace by new code.
        if( !sHandle->lastHitUpdateOK ) {
          // Last update not successful, retry with a longer timeout
			if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,60000);
        }
        */
        // UpdateHits(dpHit,dpLog,prIdx,60000); // (Rudi) Write hits to buffer (once).
        SetReady(); // (Rudi) Maybe that's not necessary either.
        break;

      case COMMAND_RESET:
        printf("COMMAND: RESET (%zd,%lu)\n",prParam,prParam2);
        ResetSimulation();
        SetReady();
        break;

      case COMMAND_EXIT:
        printf("COMMAND: EXIT (%zd,%lu)\n",prParam,prParam2);
        end = true;
        break;

      case COMMAND_CLOSE:
        printf("COMMAND: CLOSE (%zd,%lu)\n",prParam,prParam2);
        ClearSimulation();
        /* (Rudi) Don't need that.
        CLOSEDP(dpHit);
		CLOSEDP(dpLog);
		*/
        SetReady();
        break;
/*(my) case not needed
      case COMMAND_STEPAC:
        // Debug command
        printf("COMMAND: STEPAC (%zd,%llu)\n",prParam,prParam2);
        if( sHandle->loadOK ) {
          if( StartSimulation(prParam) ) {
            SetState(PROCESS_RUN,GetSimuStatus());
            SimulationACStep(1);
            if(dpHit) UpdateHits(dpHit,dpLog,prIdx,20);
            SetReady();
          } else {
            if( GetLocalState()!=PROCESS_ERROR )
              SetState(PROCESS_DONE,GetSimuStatus());
          }
        } else
          SetErrorSub("No geometry loaded");
        break;*/

      case PROCESS_RUN:
        SetStatus(GetSimuStatus()); //update hits only
        std::tie(eos,std::ignore) = SimulationRun();      // Run during 1 sec
        /*Don't need that.
		if (dpHit && (GetLocalState() != PROCESS_ERROR)) UpdateHits(dpHit,dpLog,prIdx,20); // Update hit with 20ms timeout. If fails, probably an other subprocess is updating, so we'll keep calculating and try it later (latest when the simulation is stopped).
		*/
        if(eos) {
          if( GetLocalState()!=PROCESS_ERROR ) {
            // Max desorption reached
        	// (Rudi) Think of, what to do in that case (in respect of make time dependend iterative algorithm).
            SetState(PROCESS_DONE,GetSimuStatus());
            printf("COMMAND: PROCESS_DONE (Max reached)\n");
          }
        }
        break;

      default:
        sleep(WAITTIME);
        break;
    }
  }

  // Release
  SetState(PROCESS_KILLED,"");
  //CLOSEDP(dpControl);
  //CLOSEDP(dpHit);
  //Why bother closing dataports? Windows will release handles automatically.
  return 0;

}

bool IsProcessRunning(DWORD pid)
{
	/* (my) windows stuff. rewrite
	HANDLE process = OpenProcess(SYNCHRONIZE, false, pid);
	DWORD ret = WaitForSingleObject(process, 0);
	CloseHandle(process);
	return ret == WAIT_TIMEOUT;*/
	return true;
}
