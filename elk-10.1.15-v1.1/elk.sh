#!/bin/bash

# set the number of OpenMP threads per node equal to the number of cores
# (this environment variable does not normally need to be set)
#export OMP_NUM_THREADS=

# no binding of threads to cores
export OMP_PROC_BIND=false

# increase the OpenMP stack size
export OMP_STACKSIZE=512M

# increase the regular stack size
ulimit -Ss 524288

# Elk executable file
~/elk/src/elk
