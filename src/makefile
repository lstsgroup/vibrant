FC       = gfortran
#FCFLAGS  = -Og -g  -fcheck=bounds,do,recursion,pointer -I/usr/local/include
FCFLAGS  = -Og -g -fopenmp -fcheck=bounds,do,recursion,pointer -I/scratch/fftw-3.3.10/build/include -I/scratch/greenX/build/modules
#FCLIBS   = -L/usr/local/lib -lfftw3 
FCLIBS   = -L/scratch/fftw-3.3.10/build/lib -lfftw3 -L/scratch/greenX/build/lib -lgx_ac

PROGRAM  = vibspec

# "make" builds all
all: $(PROGRAM) realclean

#vibspec: setup.o read_traj.o dipole_calc.o vel_cor.o fin_diff.o spectra.o main.o
vibspec: setup.o pade_interpolation.o read_traj.o dipole_calc.o vel_cor.o fin_diff.o spectra.o main.o

	$(FC) $(FCFLAGS) -o $@ $^ $(FCLIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(FCLIBS)

# Utility targets
.PHONY: realclean distclean 

realclean:
	rm -rf *.o *.mod

distclean: realclean
	rm -f *~ $(PROGRAM) 
