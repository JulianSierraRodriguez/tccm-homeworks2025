
F90     = gfortran
# FFLAGS  = -O3 -fopenmp -I/usr/lib/x86_64-linux-gnu/openmpi/include
# LIBS    = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lmpi_mpifh
LIBS = 
FFLAGS  = -O3

# Defaults
SRCS = \
	acceleration.f90 \
	kinetic_energy.f90 \
	potential_energy.f90 \
	read_input.f90 \
	main.f90 
OBJS = $(SRCS:.f90=.o)
PROGRAM = program.exe
# NP ?= 1
# OMP ?= 1

.PHONY: all clean

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^ $(LIBS)

#Compile step: .f90 -> .0  /// and also the .mod for modules
%.o: %.f90
	$(F90) $(FFLAGS) -c $<

# run: $(PROGRAM)
# 	OMP_NUM_THREADS=$(OMP) mpirun -np $(NP) ./$(PROGRAM)

clean:
	rm -f $(OBJS) $(PROGRAM)

help: 
	@echo "To compile a file    ===> make SRC=file.f90"
	@echo "To run an executable ===> make run SRC=file.f90 MPI=X OMP=Y"