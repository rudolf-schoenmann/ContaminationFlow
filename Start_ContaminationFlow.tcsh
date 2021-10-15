#! /bin/tcsh

module load mpi

#if (! $?LD_LIBRARY_PATH) then       
#  setenv LD_LIBRARY_PATH "$PWD"/include/lib
#else
#  if ("$LD_LIBRARY_PATH" == "")  then
#      setenv LD_LIBRARY_PATH "$PWD"/include/lib
#  else 
#      setenv LD_LIBRARY_PATH "$PWD"/include/lib\:$LD_LIBRARY_PATH
#  endif
#endif

mpirun -n $1 Debug/MolflowLinux $2

exit 0


#start simulation with the following command: 
# source Start_ContaminationFlow.tcsh parameter parameter
# e.g. 
# source Start_ContaminationFlow.tcsh 6 /scratch/schoenmann/MolflowLinuxInput/InputFileCF3.txt
#
# 1st parameter ($1): number of processes
# 2nd parameter ($2): path of input file
