# load gsl library
if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib
else
  if ("$LD_LIBRARY_PATH" == "")  then
      setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib
  else 
      setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib\:$LD_LIBRARY_PATH
  endif
endif


#allocate compute nodes.
salloc -N $1
# The varaible var1 caluclates the number of processes as a function of compute nodes.
@ var1 = 16 * $1
echo " Start ContaminationFlow with" $var1 "processes." 
# call application.
srun -n $var1 Debug/MolflowLinux $2
#deallocate compute nodes.
exit

exit 0


#start simulation with the following command: 
# source Start_ContaminationFlow.tcsh parameter parameter
# e.g. 
# source Start_ContaminationFlow.tcsh 6 /scratch/schoenmann/MolflowLinuxInput/InputFileCF3.txt
#
# 1st parameter ($1): number of compute nodes
# 2nd parameter ($2): path of input file
