
include makefile.in
# executable name
BIN = exec 
TEST = test_bis
SYMBOLS = symbols
LIB = hedgehog
all:  ${LIB}  ${TEST} ${BIN}


# build profiler symbols: TODO figure out why this isn't seen in Instruments
# dsymutil -symbol-map=lib/libhedgehog.a build/test_bis
# build directory, mirriors source code directory structure
BUILD_DIR = build
LIB_DIR = lib
LIB_NAME = libhedgehog.a

# make a list of source files to compile
SRC_PREFIX= src/
DIRS = bie3d bdry3d common deps fmm3d 
SRC_DIRS = ${DIRS:%=${SRC_PREFIX}/%}

print-%:
	@echo '$*=$($*)'

# find all the source files
SRC = $(shell find ${SRC_DIRS} -not \( -type d -path src/fmm3d/kifmm -prune \) -name *.cpp )
OBJ = $(SRC:%.cpp=$(BUILD_DIR)/%.o)
DEP = $(OBJ:%.o=%.d)

# find all test source files
TEST_DIRS = tests
TEST_SRC = $(shell find ${TEST_DIRS}  -not \( -type d -path tests/exec -prune -o -type f -path tests/test_bis.cpp -prune \) -name *.cpp )
TEST_OBJ = $(TEST_SRC:%.cpp=$(BUILD_DIR)/%.o)
TEST_DEP = $(TEST_OBJ:%.o=%.d)

# find all executable source files
BIN_DIRS = tests/exec
BIN_SRC = $(shell find ${BIN_DIRS} -name *.cpp )
BIN_OBJ = $(BIN_SRC:%.cpp=$(BUILD_DIR)/%.o)
BIN_DEP = $(BIN_OBJ:%.o=%.d)
BIN_EXEC = $(BIN_OBJ:%.o=%)

##### WARNING SEGFAULT IN POINT MARKING WHEN NOT COMPILED WITH -g WITH INTEL COMPILER###########
# linking flags
LDFLAGS = -lstdc++  -O3  -mavx
LDLIBS = -L$(LIB_DIR) -l$(LIB) ${HEDGEHOG_LIB} 

# Compiling flags
INC = -I. ${HEDGEHOG_INC}
CXXFLAGS = -std=c++11  -Wall -O3 -mavx  -fPIC -fopenmp -DFFTW3  $(INC) 


$(BIN): $(BIN_EXEC) MOVE_EXEC
MOVE_EXEC: $(BIN_EXEC)
	mv $^ $(BUILD_DIR)
	
$(TEST):  $(OBJ) $(LIB) $(BUILD_DIR)/$(TEST) $(TEST_OBJ) 
#$(TEST): $(BUILD_DIR)/$(TEST) $(TEST_OBJ) $(LIB_DIR)/$(LIB_NAME)

$(LIB): $(OBJ)
	mkdir -p $(LIB_DIR)
	rm -f $(LIB_DIR)/$(LIB_NAME)
	$(AR) $(AR_FLAGS) $(LIB_DIR)/$(LIB_NAME) $?
	$(RANLIB) $(RANLIB_FLAGS) $(LIB_DIR)/$(LIB_NAME)

# link executable
# make build file structure copying src structure
$(BUILD_DIR)/$(BIN):  $(LIB)
	mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ ${LDLIBS} -o $@
# include dependecies
$(BUILD_DIR)/$(TEST): $(BUILD_DIR)/$(TEST_DIRS)/$(TEST).o  $(TEST_OBJ) $(OBJ)
	mkdir -p $(@D)
	$(CXX) $(LDFLAGS)  $(BUILD_DIR)/$(TEST_DIRS)/$(TEST).o -o $@ $(TEST_OBJ)  ${LDLIBS}
-include $(DEP)

# Build object files for every source .cpp file
# make build file structure copying src structure 
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@
docker-build: Dockerfile
	docker build -t hedgehog:base -f Dockerfile .
.PHONY: clean

clean: 
	-rm $(BUILD_DIR)/$(BIN) $(OBJ) $(DEP) $(TEST_OBJ) $(TEST_DEP) $(BUILD_DIR)/* $(BIN_OBJ) build/tests/test_bis.o
