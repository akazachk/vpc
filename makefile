### Shell type ###
# REMEMBER: use hard tabs only in a makefile
UNAME := $(shell uname)
ifeq ($(UNAME),Linux)
  CC = g++
endif
ifeq ($(UNAME),Darwin)
  CC = clang++
endif
RM = rm -f

### Build type ###
# Choose 'debug' or 'release'
# Can also be chosen through make "BUILD_CONFIG=XX" from command line 
# Or one can call make debug or make release directly
BUILD_CONFIG = release
BUILD_CONFIG = debug

USERNAME := ${USER}
### Variables user should set ###
ifeq ($(USERNAME),otherperson)
	#PROJ_DIR = enter/dir/here
	#COIN_OR = enter/dir/here
  #GUROBI_DIR = enter/dir/here
  #GUROBI_LINK="gurobi80"
else
	PROJ_DIR = ${REPOS_DIR}/vpc
  ifeq ($(UNAME),Linux)
	  COIN_OR = $(PROJ_DIR)/lib/Cbc-2.9
	  GUROBI_LINK = "gurobi80"
	else
  	GUROBI_LINK = "gurobi81"
	  #COIN_OR = $(PROJ_DIR)/coin-or/Cbc-2.9
    #GUROBI_DIR="/Library/gurobi810/mac64"
    #GUROBI_LINK="gurobi81"
	endif
endif

# Options for solvers
USE_COIN=1
USE_CLP=1
USE_CBC=1
USE_GUROBI=1
USE_CPLEX=0

# blas and lapack
#ENV_LAPACK_LIB = $(PROJ_DIR)/lib
#ENV_BLAS_LIB = $(PROJ_DIR)/lib

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
    utility/SolverHelper.cpp \
		utility/utility.cpp \
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
  DEFS = -DTRACE
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
  DEBUG_FLAG = -g3
  OPT_FLAG = -O3
  DEFS = 
  EXTRA_FLAGS = -fmessage-length=0 -ffast-math
endif

ifeq ($(USE_CLP),1)
  DEFS += -DUSE_CLP
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
INCL_SRC_DIRS = $(addprefix -I,$(DIR_LIST))
APPLINCLS = $(INCL_SRC_DIRS) -Iinclude -Iinclude/test

APPLLIB = -lm -lz -lbz2 -lreadline

# Linker
CFLAGS = -Wall -MMD -MP
CFLAGS += -m64 $(DEBUG_FLAG) $(OPT_FLAG) $(EXTRA_FLAGS)
CXXFLAGS = $(CFLAGS) -std=c++11
#CXXFLAGS = $(CFLAGS) -std=c++11 -Wextra -Wpedantic
CXXLINKFLAGS += -std=c++11
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
		APPLLIB += \
										-lCbcSolver \
										-lCbc
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
  APPLLIB += -L${GUROBI_LIB} -lgurobi_c++ -l${GUROBI_LINK} -lm
endif
ifeq ($(USE_CPLEX),1)
  APPLINCLS += -I${ENV_CPLEX_H}/.. -I${ENV_CONCERT_H}/..
  APPLINCLS += -I${ENV_CPLEX_H} -I${ENV_CONCERT_H}
  APPLLIB += -L${ENV_CPLEX_LIB} -L${ENV_CONCERT_LIB} -lilocplex -lconcert -lcplex -lm -lpthread
  CXXLINKFLAGS += -Wl,-rpath ${ENV_CPLEX_LIB}
  CXXLINKFLAGS += -Wl,-rpath ${ENV_CONCERT_LIB}
  #APPLLIB += -lOsiCpx
endif

### Targets ###
all: 	$(EXECUTABLE)
debug: FORCE
	@$(MAKE) "dir_debug"
	@$(MAKE) "BUILD_CONFIG=debug"
release: FORCE
	@$(MAKE) "dir_release"
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

.PHONY = all clean directories print

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
	@echo 'OUT_DIR: $(OUT_DIR)'
	@echo 'DEBUG_FLAG: $(DEBUG_FLAG)'
	@echo 'OPT_FLAG: $(OPT_FLAG)'
	@echo 'DEFS: $(DEFS)'
	@echo 'EXTRA_FLAGS: $(EXTRA_FLAGS)'
	@echo 'LIB_DIR: $(LIB_DIR)'
	@echo 'SOURCES: $(SOURCES)'
	@echo 'OUT_OBJECTS: $(OUT_OBJECTS)'

print_dep: FORCE
	@echo 'DEPENDENCIES: $(DEPENDENCIES)'

FORCE: 
