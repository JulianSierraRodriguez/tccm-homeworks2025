
F90     = gfortran
LIBS = 
FFLAGS  = -O3

SRC_DIR = src
OBJ_DIR = build/obj
MOD_DIR = build/mod


# Where we have our code
SRCS = \
	$(SRC_DIR)/acceleration_mod.f90 \
	$(SRC_DIR)/energies.f90 \
	$(SRC_DIR)/error_treat.f90 \
	$(SRC_DIR)/treat_input.f90 \
  $(SRC_DIR)/treat_output.f90 \
  $(SRC_DIR)/verlet_al.f90 \
	$(SRC_DIR)/main.f90 

# We make a substitution of the f90 files path and extension into the one of the objects.
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

# Program name
PROGRAM = program.exe

# Tells makefile that these are not files, only commands
.PHONY: all clean dirs help

# first checks that the folders exists and then builds the program.
all: dirs $(PROGRAM)

# Link step, it needs the program and objects names and paths, then it links them, $@ is the target (program name) and $^ is the prerequisites (all the .o files)
$(PROGRAM): $(OBJS)
	$(F90) $(FFLAGS) -I$(MOD_DIR) -o $@ $^ $(LIBS)


#Compile step: .f90 -> .0  /// and also the .mod for modules
# -J tells where to compile the modules and -I tells where to find the modules  /// $@ the target (here the .o files) and $< the source files (the .f90 files)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | dirs
	$(F90) $(FFLAGS) -c -J$(MOD_DIR) -I$(MOD_DIR) -o $@ $<

# creates the directories if they dont exist
dirs:
	mkdir -p $(OBJ_DIR) $(MOD_DIR)

# cleans the compilation files and the executable
clean:
	rm -rf build $(PROGRAM) 

help: 
	@echo "To compile a file    ===> make"
	@echo "To create needed folders (it does it automatically)    ===> make dirs"
	@echo "clean the files after compilation ===> make clean"
