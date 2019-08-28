### Shell type ###
# REMEMBER: use hard tabs only in a makefile
UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
  CC     = g++
  SYSTEM = x86-64_linux
endif
ifeq ($(UNAME),Darwin)
  CC     = clang++
  SYSTEM = x86-64_osx
endif
RM = rm -f

### Build type ###
# Choose 'debug' or 'release'
# Can also be chosen through make "BUILD_CONFIG=XX" from command line 
# Or one can call make debug or make release directly
BUILD_CONFIG = release
BUILD_CONFIG = debug

### Variables user should set ###
PROJ_DIR=${PWD}
COIN_VERSION = 2.9
COIN_OR = $(PROJ_DIR)/lib/Cbc-$(COIN_VERSION)
ifeq ($(USER),otherperson)
  #COIN_OR = enter/dir/here

  # Optional (for testing branch and bound or enabling certain functions):
  #GUROBI_DIR = enter/dir/here
  #GUROBI_LINK="gurobi80"
  #CPLEX_DIR = enter/dir/here
endif

# For RF
ifeq ($(USER),"queen")
 CPLEX_DIR = /opt/ibm/ILOG/CPLEX_Studio128/cplex
endif

ifeq ($(USER),kazaalek)
  GUROBI_LINK = gurobi81
  #GUROBI_DIR = /opt/gurobi811/linux64
	#GUROBI_DIR = /home/gurobi/8.1.0/linux64
	GUROBI_DIR = ${HOME}/gurobi/linux64
  CPLEX_DIR = /home/ibm/cplex-studio/12.9.0.0/cplex
endif

ifeq ($(USER),akazachk)
  ifeq ($(UNAME),Linux)
  endif
  ifeq ($(UNAME),Darwin)
    GUROBI_LINK = gurobi81
    GUROBI_DIR = /Library/gurobi811/mac64
    CPLEX_DIR = /Applications/CPLEX_Studio129/cplex
  endif
endif

# Options for solvers
USE_COIN   = 1
USE_CLP    = 1
USE_CBC    = 1
USE_GUROBI = 1
USE_CPLEX  = 0
USE_CLP_SOLVER = 1
USE_CPLEX_SOLVER = 0

# Concerning executable
EXECUTABLE_STUB = vpc
SRC_DIR = src
SOURCES = main.cpp
DIR_LIST = $(SRC_DIR) $(SRC_DIR)/branch $(SRC_DIR)/cut $(SRC_DIR)/disjunction $(SRC_DIR)/utility

SOURCES += \
		branch/CbcBranchStrongDecision.cpp \
		branch/CbcCompareBFS.cpp \
		branch/OsiChooseStrongCustom.cpp \
    utility/nbspace.cpp \
    utility/OsiProblemData.cpp \
    utility/preprocess.cpp \
    utility/SolverHelper.cpp \
		utility/utility.cpp \
		utility/VPCSolverInterface.cpp \
		branch/VPCEventHandler.cpp \
		cut/CglVPC.cpp \
		cut/CutHelper.cpp \
    cut/PRLP.cpp \
    disjunction/Disjunction.cpp \
    disjunction/PartialBBDisjunction.cpp \
    disjunction/SplitDisjunction.cpp

# For running tests (need not include these or main if releasing code to others)
DIR_LIST += $(SRC_DIR)/test
SOURCES += test/analysis.cpp test/BBHelper.cpp test/DisjunctionHelper.cpp

### Set build values based on user variables ###
ifeq ($(BUILD_CONFIG),debug)
  # "Debug" build - no optimization, include debugging symbols, and keep inline functions
	SOURCES += utility/debug.cpp
  OUT_DIR = Debug
  DEBUG_FLAG = -g3
  OPT_FLAG = -O0
  DEFS = -DTRACE -DPRINT_LP_WITH_CUTS
  # message-length sets line wrapping for error messages; 0 = no line wrapping
  EXTRA_FLAGS = -fmessage-length=0
  ifeq ($(CC),g++)
    ifneq ($(USE_CPLEX),1)
      EXTRA_FLAGS += -fkeep-inline-functions 
    endif
  endif
endif
ifeq ($(BUILD_CONFIG),release)
  # "Release" build - maximum optimization, no debug symbols
  OUT_DIR = Release
  DEBUG_FLAG = 
  OPT_FLAG = -O3
  DEFS = 
  EXTRA_FLAGS = -fmessage-length=0 -ffast-math
endif
ifeq ($(USE_CLP),1)
  DEFS += -DUSE_CLP
endif
ifeq ($(USE_CLP_SOLVER),1)
  DEFS += -DUSE_CLP_SOLVER
endif
ifeq ($(USE_CBC),1)
  DEFS += -DUSE_CBC
endif
ifeq ($(USE_GUROBI),1)
  DEFS += -DUSE_GUROBI
  SOURCES += test/GurobiHelper.cpp
  GUROBI_INC="${GUROBI_DIR}/include"
  GUROBI_LIB="${GUROBI_DIR}/lib"
endif
ifeq ($(USE_CPLEX),1)
  DEFS += -DIL_STD -DUSE_CPLEX
endif
ifeq ($(USE_CPLEX_SOLVER),1)
  DEFS += -DUSE_CPLEX_SOLVER
endif
ifeq ($(COIN_VERSION),2.10)
  DEFS += -DCBC_VERSION_210PLUS
endif

EXECUTABLE = $(OUT_DIR)/$(EXECUTABLE_STUB)

# It is important that the only thing that changes about these directories
# is OUT_DIR depending on whether it is the debug or release build
# This is because later (building dependencies, archive file, cleaning, etc.)
# depends on this fact (when doing *_debug targets)
OBJ_DIR = $(OUT_DIR)/$(SRC_DIR)
OUT_DIR_LIST = $(addprefix $(OUT_DIR)/,$(DIR_LIST))
OBJECTS = $(SOURCES:.cpp=.o)
OUT_OBJECTS = $(addprefix $(OBJ_DIR)/,$(OBJECTS))

# Set includes
APPLINCLS = -Iinclude -Iinclude/test

APPLLIB = -lm -lz -lbz2 -lreadline

# Linker
CFLAGS = -Wall -MMD -MP
CFLAGS += -m64 $(DEBUG_FLAG) $(OPT_FLAG) $(EXTRA_FLAGS)
CXXFLAGS = $(CFLAGS) -std=c++14
#CXXFLAGS = $(CFLAGS) -std=c++14 -Wextra -Wpedantic
CXXLINKFLAGS += -std=c++14
ifeq ($(CC),clang++)
  CXXFLAGS += -Wno-gnu-zero-variadic-macro-arguments
  #CXXFLAGS += -stdlib=libc++ 
  #CXXLINKFLAGS += -stdlib=libc++ 
  APPLLIB += -framework Accelerate
endif
ifeq ($(CC),g++)
  ifneq (${ENV_BLAS_LIB},)
    APPLLIB += -L${ENV_BLAS_LIB} -lblas
  endif
  ifneq (${ENV_LAPACK_LIB},)
    APPLLIB += -L${ENV_LAPACK_LIB} -llapack
  endif
endif

# Set up COIN-OR stuff
ifeq ($(USE_COIN),1)
	# If not defined for the environment, define CBC / BCP here
	ifeq ($(BUILD_CONFIG),debug)
		CBC = $(COIN_OR)/buildg
	endif
	ifeq ($(BUILD_CONFIG),release)
		CBC = $(COIN_OR)/build
	endif
	CBClib = $(CBC)/lib
	CBCinc = $(CBC)/include/coin
	APPLINCLS += -isystem $(CBCinc)
	APPLLIB += -L$(CBClib)
  CXXLINKFLAGS += -Wl,-rpath $(CBClib)
	ifeq ($(USE_CBC),1)
		APPLLIB += -lCbcSolver -lCbc
	endif
  ifeq ($(USE_CLP),1)
    APPLLIB += -lOsiClp
  endif
  APPLLIB += -lCgl
  APPLLIB += -lOsi
  ifeq ($(USE_CLP),1)
    APPLLIB += -lClp
  endif
  APPLLIB += -lCoinUtils
endif
ifeq ($(USE_GUROBI),1)
  APPLINCLS += -isystem ${GUROBI_INC}
  APPLLIB   += -L${GUROBI_LIB} -lgurobi_c++ -l${GUROBI_LINK} -lm
endif
ifeq ($(USE_CPLEX),1)
  CPLEX_INC_DIR   =  $(CPLEX_DIR)/include
  CPLEX_LIB_DIR   =  $(CPLEX_DIR)/lib/$(SYSTEM)/static_pic
  APPLINCLS      += -isystem "$(CPLEX_INC_DIR)"
  APPLLIB        += -L${CPLEX_LIB_DIR} -lilocplex -lcplex -lm -lpthread -ldl
  CXXLINKFLAGS   += -Wl,-rpath $(CPLEX_LIB_DIR)
  #APPLLIB       += -lOsiCpx
endif

### Targets ###
all: | directories $(EXECUTABLE)
debug: FORCE
	@$(MAKE) "BUILD_CONFIG=debug"
release: FORCE
	@$(MAKE) "BUILD_CONFIG=release"

$(EXECUTABLE): $(OUT_OBJECTS)
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'linker'
		$(CC) $(DEFS) $(CXXLINKFLAGS) $(APPLINCLS) -o $@ $^ $(APPLLIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
		@echo ' '
		@echo 'Building target: $@'
		@echo 'Invoking' $(CC) 'compiler'
		$(CC) $(CXXFLAGS) $(DEFS) $(APPLINCLS) -c $< -o $@ 
		@echo 'Finished building: $@'

### Archive file ###
AR = ar
AR_FLAGS = -rcs
LIB_DIR=lib
archive_%:
	@$(MAKE) archive "BUILD_CONFIG=$*"

archive: $(LIB_DIR)/lib$(EXECUTABLE_STUB).a 

$(LIB_DIR)/lib$(EXECUTABLE_STUB).a: $(OUT_OBJECTS)
		@echo ' '
		@echo 'Making archive file: $@'
		@echo 'Invoking archiver'
		$(AR) $(AR_FLAGS) $@ $^
		@echo 'Finished making archive'

### Dependencies ###
# Dependencies (the -include says to ignore errors)
DEPENDENCIES = $(OUT_OBJECTS:.o=.d)
-include $(DEPENDENCIES)

### Phony ###
#.PHONY = \
#		all \
#		clean clean_debug distclean_debug clean_release distclean_release \
#		directories dir_debug dir_release dir_lib_debug dir_lib_release \
#		print print_dep \
#		archive_debug archive_release archive

.PHONY = all clean directories distclean docs print

### Docs ###
docs:
	@doxygen

### Cleaning ###
clean_%: FORCE
	@$(MAKE) clean "BUILD_CONFIG=$*"
distclean_%: FORCE
	@$(MAKE) distclean "BUILD_CONFIG=$*"

clean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE)

distclean: FORCE
	@$(RM) $(OUT_OBJECTS) $(EXECUTABLE) $(DEPENDENCIES) $(LIB_DIR)/lib$(EXECUTABLE_STUB).a

### Making directories that you need ###
MKDIR_P = mkdir -p

dir_%: FORCE
	@$(MAKE) directories "BUILD_CONFIG=$*"
directories: $(OUT_DIR_LIST)
$(OUT_DIR_LIST):
	$(MKDIR_P) $(OUT_DIR_LIST)

dir_lib_%: FORCE
	@$(MAKE) dir_lib "BUILD_CONFIG=$*"
dir_lib: FORCE
	$(MKDIR_P) $(LIB_DIR)
print: FORCE
	$(info UNAME: ${UNAME})
	$(info CC: ${CC})
	$(info COIN_OR: ${COIN_OR})
	$(info CPLEX_DIR: ${CPLEX_DIR})
	$(info GUROBI_DIR: ${GUROBI_DIR})
	$(info GUROBI_LINK: ${GUROBI_LINK})
	$(info USE_COIN: ${USE_COIN})
	$(info USE_CLP: ${USE_CLP})
	$(info USE_CBC: ${USE_CBC})
	$(info USE_GUROBI: ${USE_GUROBI})
	$(info USE_CPLEX: ${USE_CPLEX})
	$(info OUT_DIR: ${OUT_DIR})
	$(info DEBUG_FLAG: ${DEBUG_FLAG})
	$(info OPT_FLAG: ${OPT_FLAG})
	$(info DEFS: ${DEFS})
	$(info EXTRA_FLAGS: ${EXTRA_FLAGS})
	$(info LIB_DIR: ${LIB_DIR})
	$(info SOURCES: ${SOURCES})
	$(info OUT_OBJECTS: ${OUT_OBJECTS})

print_dep: FORCE
	$(info DEPENDENCIES: ${DEPENDENCIES})

FORCE: 
