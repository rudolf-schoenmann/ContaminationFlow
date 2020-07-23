module load mpi

if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib
else
  if ("$LD_LIBRARY_PATH" == "")  then
      setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib
  else 
      setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib\:$LD_LIBRARY_PATH
  endif
endif
