# =========================
# Fortran Makefile (gfortran, debug)
# =========================

FC      := gfortran

# Executables
EXE_ONE    := trace_one_ray

# Build directories
BUILD   := build
OBJDIR  := $(BUILD)/obj
MODDIR  := $(BUILD)/mod

# Debug flags
FFLAGS  := -std=gnu -g -O0 -Wall -Wextra -fcheck=all -fbacktrace \
           -ffree-line-length-none \
           -J$(MODDIR) -I$(MODDIR)
#          -ffpe-trap=zero -ffpe-summary=all

# macOS SDK library path (fixes "library not found for -lSystem" error)
LDFLAGS := -L$(shell xcrun --show-sdk-path 2>/dev/null)/usr/lib

# =========================
# Source files (ORDER MATTERS)
# =========================

SRC_MOD := \
  mod1_internal_grids.f90 \
  mod2_constants.f90 \
  mod3_rootsfinding.f90 \
  mod4_ellfunction.f90

SRC_BLC :=

# Common library-ish sources (modules used by both drivers)
SRC_COMMON := \
  mod5_blcoordinate.f90 \
  mod6_raytrace_camera.f90

# Driver programs (each produces its own executable)
SRC_DRIVER_ONE    := trace_one_ray.f90

# All sources (for compilation)
SRCS := $(SRC_MOD) $(SRC_BLC) $(SRC_COMMON) $(SRC_DRIVER_ONE)

# Object files
OBJS_MOD_BLC_COMMON := $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC_MOD) $(SRC_BLC) $(SRC_COMMON))
OBJ_ONE             := $(OBJDIR)/$(SRC_DRIVER_ONE:.f90=.o)

# =========================
# Targets
# =========================

.PHONY: all
all: $(EXE_ONE)

.PHONY: one
one: $(EXE_ONE)

$(EXE_ONE): $(OBJS_MOD_BLC_COMMON) $(OBJ_ONE) | $(BUILD)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS_MOD_BLC_COMMON) $(OBJ_ONE)

# Compile rule
$(OBJDIR)/%.o: %.f90 | $(BUILD) $(OBJDIR) $(MODDIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

# Create build directories
$(BUILD) $(OBJDIR) $(MODDIR):
	@mkdir -p $@

# Clean
.PHONY: clean
clean:
	rm -rf $(BUILD) $(EXE_ONE) *.o *.mod

# Info
.PHONY: info
info:
	@echo "Compiler  : $(FC)"
	@echo "Flags     : $(FFLAGS)"
	@echo "Exe one   : $(EXE_ONE)"
