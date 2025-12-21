
F90     = gfortran
# FFLAGS  = -O3 -fopenmp -I/usr/lib/x86_64-linux-gnu/openmpi/include
# LIBS    = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lmpi_mpifh
LIBS = 
FFLAGS  = -O3

OBJ_DIR = build/obj
MOD_DIR = build/mod


# Defaults
SRCS = \
	acceleration_mod.f90 \
	energies.f90 \
	error_treat.f90 \
	treat_input.f90 \
  treat_output.f90 \
  verlet_al.f90 \
	main.f90 

# OBJS = $(SRCS:.f90=.o)
OBJS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SRCS))
# MODS = $(SRCS:.f90=.mod)
# MODS = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(SRCS))
PROGRAM = program.exe



.PHONY: all clean help

all: dirs $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F90) $(FFLAGS) -I$(MOD_DIR) -o $@ $^ $(LIBS)

#Compile step: .f90 -> .0  /// and also the .mod for modules
$(OBJ_DIR)/%.o: %.f90 | dirs
	$(F90) $(FFLAGS) -c -J$(MOD_DIR) -I$(MOD_DIR) -o $@ $<

dirs:
	mkdir -p $(OBJ_DIR) $(MOD_DIR)

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(PROGRAM) 

help: 
	@echo "To compile a file    ===> make SRC=file.f90"
	@echo "clean the files after compilation ===> make clean"
