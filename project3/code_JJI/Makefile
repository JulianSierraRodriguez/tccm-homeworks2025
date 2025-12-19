
F90     = gfortran
FFLAGS  = -O3 -fopenmp -I/usr/lib/x86_64-linux-gnu/openmpi/include
LIBS    = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lmpi_mpifh

# Defaults
SRC ?= dot_product.f90
OBJ = $(SRC:.f90=.o)
PROGRAM = $(SRC:.f90=.x)
NP ?= 1
OMP ?= 1

.PHONY: all run clean

all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(F90) $(FFLAGS) $(OBJ) $(LIBS) -o $@

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

run: $(PROGRAM)
	OMP_NUM_THREADS=$(OMP) mpirun -np $(NP) ./$(PROGRAM)

clean:
	rm -f $(OBJ) $(PROGRAM)

help: 
	@echo "To compile a file    ===> make SRC=file.f90"
	@echo "To run an executable ===> make run SRC=file.f90 MPI=X OMP=Y"