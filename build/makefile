# --------- toolchain & optimisation flags -----------------------------------
FC          = mpiifx

# Architecture & threading
ARCH_FLAG   = -xCORE-AVX2
THREAD_FLAG = -qopenmp

# Common optimization flags (used for both compile and link)
OPT_COMMON  = -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-streaming-stores=always \
              -qopt-mem-layout-trans=3 -align array64byte -heap-arrays 256 \
              -warn all -diag-error-limit=0

# Production build: fast & safe
FFLAGS_OPT  = -O3 -ipo $(ARCH_FLAG) $(THREAD_FLAG) $(OPT_COMMON) \
              -fp-model precise -fno-emit-llvm

# Debug/dev build: fast compile & debug symbols
FFLAGS_FAST = -O2 -g -traceback -debug minimal $(ARCH_FLAG) $(THREAD_FLAG) $(OPT_COMMON) \
              -fp-model source -heap-arrays 128 -fno-emit-llvm

# Consistent linker flags (link time optimization + threading + MKL)
LDFLAGS     = -qmkl=parallel $(THREAD_FLAG) $(ARCH_FLAG) -ipo

# --------- project layout ----------------------------------------------------
SRC_DIR     = ../src
TEST_DIR    = ../tests
BIN_DIR     = ../bin
OBJ_DIR     = .obj
DEP_DIR     = $(OBJ_DIR)

# --------- source lists (edit only these) ------------------------------------
COMMON_SRC  = mtdefs mt95 utils param_mod model solve_flow qn hmc markovchain
STEP1_SPEC  = update_bw
STEP2_SPEC  = generate_markov_chain
STEP3_SPEC  = virial
TEST1_SPEC  = $(TEST_DIR)/test
TEST2_SPEC  = $(TEST_DIR)/test2

# helper to prepend SRC_DIR and append .f90
define mk_src
$(addprefix $(SRC_DIR)/,$(addsuffix .f90,$(1)))
endef

SOURCES_COMMON := $(call mk_src,$(COMMON_SRC))
SOURCES_STEP1  := $(SOURCES_COMMON) $(call mk_src,$(STEP1_SPEC))
SOURCES_STEP2  := $(SOURCES_COMMON) $(call mk_src,$(STEP2_SPEC))
SOURCES_STEP3  := $(SOURCES_COMMON) $(call mk_src,$(STEP3_SPEC))
SOURCES_TEST1  := $(SOURCES_COMMON) $(TEST1_SPEC).f90
SOURCES_TEST2  := $(SOURCES_COMMON) $(TEST2_SPEC).f90

# map .f90 → .o inside $(OBJ_DIR)
OBJECTS_STEP1  := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SOURCES_STEP1)))
OBJECTS_STEP2  := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SOURCES_STEP2)))
OBJECTS_STEP3  := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SOURCES_STEP3)))
OBJECTS_TEST1  := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SOURCES_TEST1)))
OBJECTS_TEST2  := $(patsubst %.f90,$(OBJ_DIR)/%.o,$(notdir $(SOURCES_TEST2)))

# executables
BINARIES     = define_coefficients generate_markov_chain evaluate_expectations \
               test_program test_program2
EXECUTABLES  = $(addprefix $(BIN_DIR)/,$(BINARIES))

###############################################################################
#  build targets
###############################################################################
override FFLAGS := $(FFLAGS_OPT)   # <- now exists

.PHONY: all fast
all: $(EXECUTABLES)

fast: override FFLAGS := $(FFLAGS_FAST)
fast:  $(EXECUTABLES)

# compile .f90 -> .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -module $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -module $(OBJ_DIR)

# link rules
$(BIN_DIR)/define_coefficients: $(OBJECTS_STEP1) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/generate_markov_chain: $(OBJECTS_STEP2) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/evaluate_expectations: $(OBJECTS_STEP3) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/test_program: $(OBJECTS_TEST1) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

$(BIN_DIR)/test_program2: $(OBJECTS_TEST2) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $^

# create build directories if missing
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

###############################################################################
#  interactive menu
###############################################################################
define banner
	@printf "\n\033[1;36m%-60s\033[0m\n" "==============================="
	@printf "\033[1;36m%-60s\033[0m\n" "$(1)"
	@printf "\033[1;36m%-60s\033[0m\n\n" "==============================="
endef

mn: all
	$(call banner,HMC Simulation Steps)
	@echo "1. Define Coefficients"
	@echo "2. Generate Markov Chain"
	@echo "3. Evaluate Expectations"
	@echo "4. Run Test Program 1: Hamiltonian Conservation"
	@echo "5. Run Test Program 2: Action Derivatives"
	@echo "6. Exit"
	@echo "-----------------------------------------------"
	@read -p "Enter your choice: " choice; \
	case $$choice in \
	  1)  $(BIN_DIR)/define_coefficients ;; \
	  2)  $(BIN_DIR)/generate_markov_chain ;; \
	  3)  $(BIN_DIR)/evaluate_expectations ;; \
	  4)  $(BIN_DIR)/test_program ;; \
	  5)  $(BIN_DIR)/test_program2 ;; \
	  6)  echo "Exiting. Goodbye!" ;; \
	  *)  echo "Invalid option." ;; \
	esac

###############################################################################
#  housekeeping
###############################################################################
clean:
	$(RM) -r $(OBJ_DIR) *.mod

veryclean: clean
	$(RM) $(EXECUTABLES)