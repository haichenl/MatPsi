#!/bin/sh 

## IMPORTANT ##
## All $psi4 mentioned below better be /home/haichen/psi4wpgcc44changecode/ 
## as the original psi4 source code has been changed a little bit 
## IMPORTANT ##

# Compilation notice: 
# Please make sure: 
# (1) ./ contains read_options.o from $psi4/objdir/src/bin/psi4/ 
# (2) ./lib contains all the libXXX folders from $psi4/src/lib 
#       and all the libPSI_XXX.a static libraries from $psi4/objdir/lib/ 
# (3) ./include contains folders libderiv and libint from $psi4/objdir/include
#       and all the .h files from $psi4/include 

# Simply execute this file to compile and link. 
# Refer to original psi4 documents if met blas/lapack/boost issues. 

mex -DLinux -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE CFLAGS='$CFLAGS -O2 -fopenmp -fPIC -pthread' MatPsi_mex.cpp -I./lib -I/usr/include/python2.7 -I./include -I/usr/include LDFLAGS='$LDFLAGS -fopenmp -Wl,-export-dynamic -Wl,--whole-archive read_options.o ./lib/libPSI_mints.a ./lib/libPSI_trans.a ./lib/libPSI_dpd.a ./lib/libPSI_chkpt.a ./lib/libPSI_iwl.a ./lib/libPSI_psio.a ./lib/libPSI_qt.a ./lib/libPSI_ciomr.a ./lib/libPSI_options.a ./lib/libPSI_util.a ./lib/libPSI_deriv.a ./lib/libPSI_int.a ./lib/libPSI_parallel.a -Wl,--no-whole-archive -lm -lbsd -llapack -lcblas -lf77blas -latlas -L/usr/local/atlasfpic/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7 -L/usr/lib/x86_64-linux-gnu -L/usr/lib -L/lib/x86_64-linux-gnu -L/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib -L/usr/lib -lgfortran -lm -lf77blas -latlas -L/usr/local/atlasfpic/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7 -L/usr/lib/x86_64-linux-gnu -L/usr/lib -L/lib/x86_64-linux-gnu -L/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib -lgfortran -lm -lboost_filesystem -lboost_python -lboost_regex -lboost_serialization -lboost_system -lboost_thread -lrt -lpthread -L/usr/lib/python2.7/config -lpthread -ldl -lutil -lm -lpython2.7 -Xlinker -export-dynamic -Wl,-O2 -Wl,-Bsymbolic-functions' -cxx 


