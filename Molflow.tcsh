module load mpi

if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH "$PWD"/include/lib
else
  if ("$LD_LIBRARY_PATH" == "")  then
      setenv LD_LIBRARY_PATH "$PWD"/include/lib
  else 
      setenv LD_LIBRARY_PATH "$PWD"/include/lib\:$LD_LIBRARY_PATH
  endif
endif
