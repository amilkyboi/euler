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
FFLAGS= -g# -fbacktrace -O -Wall -fcheck=all
LDFLAGS = -m64 -g # -ffpe-trap=invalid,zero -fbacktrace -fdump-core
EXE = main # name of output .exe

.SUFFIXES: .f90 .o # valid file suffixes

SRC = mod_types.f90 input.f90 timing.f90 gasprop.f90 freestream.f90 reference.f90 gridprop.f90 flowprop.f90 fluxes.f90 functions.f90 \
 	  startup.f90 allocation.f90 ghost.f90 area.f90 initializer.f90 bcwall.f90 bcinout.f90 flux.f90 residual.f90 dissipation.f90 \
	  timestep.f90 rk4.f90 plt.f90 main.f90 # all files $(wildcard *.f90)

OBJ = $(SRC:.f90=.o) # compile sources of object files

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90 
	$(FC) $(FFLAGS) -c $<

# delete output and mod files using make clean
clean:
	-$(RM) *.o
	-$(RM) *.mod
