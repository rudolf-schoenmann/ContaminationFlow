#! /bin/tcsh

if (! $?LD_LIBRARY_PATH) then       
  setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib
else
  if ("$LD_LIBRARY_PATH" == "")  then
      setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib
  else 
      setenv LD_LIBRARY_PATH /opt/ohpc/pub/libs/gnu7/gsl/2.4/lib\:$LD_LIBRARY_PATH
  endif
endif

exit 0
