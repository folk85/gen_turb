# The compiler
# F_COMP = gfortran
F_DB_MOD = True
F_DB_MOD = False
F_COMP = ifort

# change these to proper directories where each file should be
SRCDIR   = src
BINDIR   = bin
OBJDIR   = obj
OBJDIR_DB   = $(BINDIR)/obj-db
OBJDIR_RL   = $(BINDIR)/obj-rl
TSTDIR   = tests
rm       = rm -f

# flags for debugging or for maximum performance, comment as necessary
USE_MKL = True
LD_MKL = -mkl
FC_INTEL_DB = -g -check all -fpe0 -warn -traceback -debug extended -module $(OBJDIR) $(LD_MKL)
FC_INTEL_RL = -O3 -sox -module $(OBJDIR) $(LD_MKL)

# FC_GNU_DB = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -ffpe-trap=zero,overflow,underflow -finit-real=nan -J$(OBJDIR)
FC_GNU_DB = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -J$(OBJDIR)
FC_GNU_RL = -O2 -J$(OBJDIR)

ifeq ($(F_DB_MOD),True)
	FC_INTEL=$(FC_INTEL_DB)
	FC_GNU=$(FC_GNU_DB)
	OBJDIR = $(OBJDIR_DB)
else
	FC_INTEL=$(FC_INTEL_RL)
	FC_GNU=$(FC_GNU_RL)
	OBJDIR = $(OBJDIR_RL)
endif
ifeq ($(F_COMP),ifort)
	FCFLAGS=$(FC_INTEL)
else
	FCFLAGS=$(FC_GNU)
endif


# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS =  $(LD_MKL) -lpthread -ldl -lm

# List of executables to be built within the package
TARGET = main 

# Define python script to show spectr
PL_SPECTR = $(TSTDIR)/plot_spectr.py

F_FILES := $(wildcard $(SRCDIR)/[^u][^s][^e]*.f90)
FT_FILES := $(wildcard $(TSTDIR)/*.f90)
OBJ_FILES := $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f90=.o)))
# OBJ_FILES += $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f=.o)))
MOD_FILES := $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f90=.mod)))

all: $(BINDIR)/$(TARGET)

build: $(OBJ_FILES)

$(BINDIR)/$(TARGET): $(OBJ_FILES)
	$(F_COMP) $(FCFLAGS) -o $@ $^ $(LDFLAGS) 
	@echo "Linking complete!"

$(OBJ_FILES): $(OBJDIR)/%.o : $(SRCDIR)/%.f90 ${OBJDIR}/.sentinel
	$(F_COMP) $(FCFLAGS) -c -o $@ $<

print:
	@echo "F_FILES    $(F_FILES)"
	@echo "OBJ_FILES  $(OBJ_FILES)"
	@echo "MOD_FILES  $(MOD_FILES)"

%/.sentinel:
	$(foreach d,$(subst /, ,$*),mkdir -p $d && cd $d && ):
	touch $@

.PHONY: clean
clean:
	@$(rm) $(OBJ_FILES) *.mod
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

run:
	$(BINDIR)/$(TARGET)

spectr:
	python $(PL_SPECTR)