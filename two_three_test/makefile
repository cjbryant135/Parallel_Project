#-------------------------------------------------------
COMP = mpif90
FLAGS = -fcheck=all
LAPACKFILES = -L/usr/lib/lapack -L/usr/lib/libblas \
  -L/usr/lib/libpblas \
  -l lapack -l blas -llapack #-l scalapack -l pblas -lscalapack
#------------------------------------------------------------------------------
main_files = scala_test.f90
#mod_files =  sort_mod.f90 oned_parameters.f90 oned_module.f90

all: test_exe

test_exe: $(main_files)
	$(COMP)	$(FLAGS)	$(mod_files)	$(main_files)	$(LAPACKFILES) -o $@

clean: 
	rm *_exe *.o *.mod
