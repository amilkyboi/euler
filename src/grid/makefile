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

# debugging flags:
# -Og -g3 -Wall -Wextra -pedantic -std=f2018 -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,underflow,overflow

# production flags:
# -O3 -std=f2018

FC      = gfortran
FFLAGS  = -Og -g3 -Wall -Wextra -pedantic -std=f2018 -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,underflow,overflow
LDFLAGS = -Og -g3 -Wall -Wextra -pedantic -std=f2018 -fimplicit-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,underflow,overflow
EXE     = main

.SUFFIXES: .f90 .o

# all files $(wildcard *.f90)
SRC = mod_types.f90 input.f90 gridprop.f90 allocate.f90 functions.f90 algebraic.f90 elliptic.f90 main.f90

OBJ = $(SRC:.f90=.o) # compile sources of object files

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LIBS) -o $(EXE)

%.o  : %.f90 
	$(FC) $(FFLAGS) -c $<

# delete output and mod files using make clean
clean:
	-$(RM) *.o
	-$(RM) *.mod
