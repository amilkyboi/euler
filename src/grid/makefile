# change terminal commands if OS=WIN
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

FC = gfortran # name of Fortran compiler
FFLAGS= -m64 -g -ffpe-trap=invalid,zero -fbacktrace -fdump-core
LDFLAGS = -m64 -g # -ffpe-trap=invalid,zero -fbacktrace -fdump-core
EXE = main # name of output .exe

.SUFFIXES: .f90 .o # valid file suffixes

SRC = mod_types.f90 functions.f90 gridprop.f90 algebraic.f90 elliptic.f90 main.f90 # all files $(wildcard *.f90)

OBJ = $(SRC:.f90=.o) # compile sources of object files

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90 
	$(FC) $(FFLAGS) -c $<

# delete output and mod files using make clean
clean:
	-$(RM) *.o
	-$(RM) *.mod
