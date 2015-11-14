#-------------------------------------------------------
COMP = mpif90
FLAGS = -fcheck=all
LAPACKFILES = -L/usr/lib/lapack -L/usr/lib/libblas -l lapack -l blas -llapack
#------------------------------------------------------------------------------
test_files = eig_test.f90
#mod_files = oned_module.f90 sort_mod.f90

all: test_exe

test_exe:	$(test_files)
	$(COMP)	$(FLAGS)	$(test_files)	$(LAPACKFILES) -o $@

clean: 
	rm *_exe *.o *.mod
