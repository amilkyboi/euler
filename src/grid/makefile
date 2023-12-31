# change terminal commands if the user is running windows
ifeq ($(OS),Windows_NT) 
RM = del /Q /F
CP = copy /Y
ifdef ComSpec
SHELL := $(ComSpec)
endif
ifdef COMSPEC
SHELL := $(COMSPEC)
endif
else
RM = rm -rf
CP = cp -f
endif

# compiler and linker flags used for debugging
DBUG_FLAGS = -Og -g3 -Wall -Wextra -pedantic -std=f2018 -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,underflow,overflow

# compiler and linker flags used for production
PROD_FLAGS = -O3 -std=f2018 -fimplicit-none -ffpe-trap=invalid,zero,underflow,overflow

# specify fortran compiler and executable filename
FC  = gfortran
EXE = main

# valid file suffixes
.SUFFIXES: .f90 .o

# source files, must be placed in order of variable use
SRC = mod_types.f90 input.f90 gridprop.f90 timing.f90 allocation.f90 functions.f90 algebraic.f90 \
	  elliptic.f90 main.f90

# compile .f90 files to .o files
OBJ = $(SRC:.f90=.o)

# options run with `make`
$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90 
	$(FC) $(FFLAGS) -c $<

# debug mode
dbug: FFLAGS  = $(DBUG_FLAGS)
dbug: LDFLAGS = $(DBUG_FLAGS)
dbug: clean $(EXE)

# production mode
prod: FFLAGS  = $(PROD_FLAGS)
prod: LDFLAGS = $(PROD_FLAGS)
prod: clean $(EXE)

# clean option
clean:
	$(RM) *.o
	$(RM) *.mod

# sets the default mode to debug, which allows `make dbug` to be run simply as `make`
.DEFAULT_GOAL := dbug
