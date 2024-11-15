MAKE = make
AR = ar
F90 = gfortran
MKL_PATH= /opt/intel/oneapi/mkl/2024.1/lib
#SRC_MKL = /opt/intel/oneapi/mkl/2024.1/ # 请确保这是 MKL 库的正确路径
F90_OPTS = -O3 -Wall -Wunused -Warray-temporaries -Waliasing -Wampersand -Wconversion -Winteger-division -Wintrinsics-std -Wsurprising -Wtabs -Wunused-parameter -Wfunction-elimination -Wrealloc-lhs-all -Wcompare-reals -fopenmp -fbounds-check -pg
F90_LIB = -L$(MKL_PATH) -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread
SRC_MPI = mpi_stub.f90
SRC_FFT = mkl_dfti.f90 zfftifc_mkl.f90 cfftifc_mkl.f90
#
SRC_LIBXC = libxcifc_stub.f90
SRC_W90S = w90_stub.f90
