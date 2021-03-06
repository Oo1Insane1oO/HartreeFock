INC = -I/usr/local/include -I/usr/local/include/eigen3
LIBS = -L/usr/local/lib

TESTS ?= OFF
ifeq ($(TESTS),ON)
	INC += -I/usr/local/include/UnitTest++
	LIBS += -L/usr/local/lib/libUnitTest++.a -lUnitTest++
endif

CC = mpic++
LDFLAGS	= -fopenmp -lpthread -lyaml-cpp -lboost_filesystem -lstdc++fs $(LIBS)
CXXFLAGS = -O3 -fopenmp -std=c++17 $(INC)

WRITE ?= OFF
ifeq ($(WRITE),ON)
	CXXFLAGS += -D$(WRITE)
endif

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -g -O0 -std=c++17 -Wall -Wextra $(INC)
endif

ifeq ($(TESTS),ON)
	CXXFLAGS += -DTESTS
endif

# WAVEFUNCTION = HYDROGENMOLECULE
# CXXFLAGS += -D$(WAVEFUNCTION)

EXEC = main
SRC = src

SUFFIXES += .d

INPUTS_DIR = inputs
DATA_DIR = data

ifneq ($(WAVEFUNCTION),)
	WAVE_DIR = $(SRC)/wavefunctions
	WAVE_LOWER = $(shell echo $(WAVEFUNCTION) | tr A-Z a-z)
	WAVE = $(WAVE_DIR)/$(WAVE_LOWER)
	WAVESRC = $(WAVE).cpp
	WAVEHDR = $(WAVE).h
endif
OBJ_DIR = bin
TEST_DIR = tests

BASIS = GAUSSHERMITE

ifeq ($(BASIS),GAUSSHERMITE)
	OTHER_BASIS_SOURCE := $(SRC)/basisfunctions/cartesian.cpp $(SRC)/hermite/hexpander.cpp $(SRC)/gaussianquadrature.cpp
	INTEGRALS := $(wildcard $(SRC)/integrals/gaussianintegrals.cpp)
	BASISSRC = $(SRC)/basis/gaussianbasis.cpp
endif

ifeq ($(BASIS),DOUBLEWELL)
	OTHER_BASIS_SOURCE := $(SRC)/basisfunctions/cartesian.cpp $(SRC)/hermite/hexpander.cpp $(SRC)/gaussianquadrature.cpp
	INTEGRALS := $(wildcard $(SRC)/integrals/gaussianintegrals.cpp) $(SRC)/integrals/doublewell.cpp
	BASISSRC = $(SRC)/basis/gaussianbasis.cpp
endif

ifeq ($(BASIS),STYPEGAUSSIAN)
	OTHER_BASIS_SOURCE := $(SRC)/basisfunctions/gaussianprimitivebasis.cpp $(SRC)/basisfunctions/gaussiancontractedbasis.cpp
	INTEGRALS := $(wildcard $(SRC)/integrals/gaussianstypeintegrals.cpp)
	BASISSRC = $(SRC)/basis/stypebasis.cpp
endif

BASISHDR = $(patsubst %.cpp,%.h,$(BASISSRC))
	
CXXFLAGS += -D$(BASIS)

MAINSRC = $(SRC)/main.cpp

ifeq ($(TESTS), ON)
	TEST_SRC = $(wildcard $(TEST_DIR)/*.cpp)
	TEST_OBJ += $(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC))
	TEST_DEPS = $(patsubst %.o,%.d,$(TEST_OBJ))
endif

TMPSRC = $(wildcard $(SRC)/*.cpp) $(WAVESRC) $(BASISSRC) $(OTHER_BASIS_SOURCE) $(HF_SOURCE) $(INTEGRALS)
SOURCES = $(filter-out $(MAINSRC), $(TMPSRC))
HEADERS = $(wildcard $(SRC)/*.h) $(WAVEHDR) $(BASISHDR) $($(OTHER_BASIS_SOURCE):.cpp=.h) $($(HF_SOURCE):.cpp=.h) $($(integrals):.cpp=.h)
SOURCES += $(MAINSRC)
OBJECTS = $(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

DEPS = $(patsubst %.o,$(OBJ_DIR)/%.d,$(OBJECTS))

build: $(EXEC)

clean:
	rm -rf $(OBJ_DIR)
	rm $(EXEC)

createDataDir:
	mkdir $(DATA_DIR)

cleanData:
	rm -rf $(DATA_DIR)

cleanInputs:
	rm $(wildcard $(INPUTS_DIR)/*)

.PHONY: clean build

NODEPS := clean tags svn
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
    -include $(DEPS) $(TEST_DEPS)
endif

$(EXEC): $(OBJECTS) | $(TEST_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

$(OBJ_DIR)/%.d: $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -MM -MT '$(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$<)' $< -MF $@

$(OBJ_DIR)/%.d: $(TEST_DIR)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -MM -MT '$(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/%.o,$<)' $< -MF $@

$(OBJ_DIR)/%.o: $(SRC)/%.cpp $(OBJ_DIR)/%.d $(HEADERS)
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp $(OBJ_DIR)/%.d
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ -c $<
