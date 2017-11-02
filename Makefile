# The compiler
F_COMP = gfortran
F_DB_MOD = True
# F_COMP = ifort

# change these to proper directories where each file should be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin
TSTDIR   = tests
rm       = rm -f

# flags for debugging or for maximum performance, comment as necessary
FC_INTEL_DB = -g -check all -fpe0 -warn -traceback -debug extended -module $(OBJDIR)
FC_INTEL_RL = -O2 -module $(OBJDIR)

FC_GNU_DB = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -ffpe-trap=zero,overflow,underflow -finit-real=nan -J$(OBJDIR)
# FC_GNU_DB = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
FC_GNU_RL = -O2 -J$(OBJDIR)

ifeq ($(F_DB_MOD),True)
	FC_INTEL=$(FC_INTEL_DB)
	FC_GNU=$(FC_GNU_DB)
else
	FC_INTEL=$(FC_INTEL_RL)
	FC_GNU=$(FC_GNU_RL)
endif
ifeq ($(F_COMP),ifort)
	FCFLAGS=$(FC_INTEL)
else
	FCFLAGS=$(FC_GNU)
endif


# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS = 

# List of executables to be built within the package
TARGET = main 

# Define python script to show spectr
PL_SPECTR = $(TSTDIR)/plot_spectr.py

F_FILES := $(wildcard $(SRCDIR)/*.f90)
FT_FILES := $(wildcard $(TSTDIR)/*.f90)
OBJ_FILES := $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f90=.o)))
# OBJ_FILES += $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f=.o)))
MOD_FILES := $(addprefix $(OBJDIR)/,$(notdir $(F_FILES:.f90=.mod)))

# "make" builds all
# all: $(TARGET)

$(BINDIR)/$(TARGET): $(OBJ_FILES)
	# @mkdir -p $(BINDIR)
	$(F_COMP) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Linking complete!"

$(OBJ_FILES): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	# @mkdir -p $(OBJDIR)
	$(F_COMP) $(FCFLAGS) -c -o $@ $<

# $(OBJ_FILES): $(OBJDIR)/%.o : $(SRCDIR)/%.f
	# $(F_COMP) $(FCFLAGS) -c -o $@ $<

print:
	@echo "F_FILES    $(F_FILES)"
	@echo "OBJ_FILES  $(OBJ_FILES)"
	@echo "MOD_FILES  $(MOD_FILES)"


.PHONY: clean
clean:
	@$(rm) $(OBJ_FILES) *.mod
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"
# clean:
# 	@rm -rf $(OBJ_FILES)	$(PROGRAMS) *.mod *.MOD *~ 

run:
	$(BINDIR)/$(TARGET)

spectr:
	python27 $(PL_SPECTR)