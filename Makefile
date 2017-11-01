# The compiler
F_COMP = gfortran
F_DB_MOD = True
# F_COMP = ifort
# flags for debugging or for maximum performance, comment as necessary
FC_INTEL_DB = -g -check all -fpe0 -warn -traceback -debug extended
FC_INTEL_RL = -O2

FC_GNU_DB = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
FC_GNU_RL = -O2

ifeq ($(F_DB_MOD),True)
	FC_INTEL=$(FC_INTEL_DB)
	FC_GNU=$(FC_GNU_DB)
else
	FC_INTEL=$(FC_INTEL_RL)
	FC_GNU=$(FC_GNU_RL)
endif
ifeq ($(F_COMP),ifort)
	ifeq()
	FCFLAGS=$(FC_INTEL)
else
	FCFLAGS=$(FC_GNU)
endif


# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = main 

# Define python script to show spectr
PL_SPECTR = tests/plot_spectr.py

F_FILES := $(wildcard src/*.f90)
FT_FILES := $(wildcard tests/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F_FILES:.f90=.o)))

# "make" builds all
all: $(PROGRAMS)

main: $(OBJ_FILES)
	$(F_COMP) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.f90
	$(F_COMP) $(FCFLAGS) -c -o $@ $<

# tests/%.o: tests/%.f90
# 	$(F_COMP) $(FCFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_FILES)	$(PROGRAMS)

run:
	./main

spectr:
	python27 $(PL_SPECTR)