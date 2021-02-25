setenv PATH $HOME/mpi/bin\:$HOME/gcc4/bin\:$PATH

if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib\:$HOME/mpi/lib\:$HOME/gcc4/lib\:$HOME/gcc4/lib64
else
  if ("$LD_LIBRARY_PATH" == "")  then
      setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib\:$HOME/mpi/lib\:$HOME/gcc4/lib\:$HOME/gcc4/lib64
  else 
      setenv LD_LIBRARY_PATH $HOME/MolflowLinux/include/lib\:$HOME/mpi/lib\:$HOME/gcc4/lib\:$HOME/gcc4/lib64\:$LD_LIBRARY_PATH
  endif
endif
