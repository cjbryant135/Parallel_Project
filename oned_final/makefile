#-------------------------------------------------------
COMP = mpif90
FLAGS = -fcheck=all
LAPACKFILES = -L/usr/lib/lapack -L/usr/lib/libblas -l lapack -l blas -llapack
#------------------------------------------------------------------------------
main_files = oned_main.f90
mod_files =  sort_mod.f90 oned_parameters.f90 oned_module.f90

all: main_exe

main_exe: $(mod_files)	$(main_files)
	$(COMP)	$(FLAGS)	$(mod_files)	$(main_files)	$(LAPACKFILES) -o $@

clean: 
	rm *_exe *.o *.mod
