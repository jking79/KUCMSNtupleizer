# Compiler and flags

CC = g++
CXXFLAGS_COMMON = -O3 -march=native -ffast-math -funroll-loops -flto -fno-trapping-math -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
CXXFLAGS_DEBUG = -g -O0 -Wall -Wextra
CXXFLAGS_RELEASE = $(CXXFLAGS_COMMON)
CXXFLAGS ?= $(CXXFLAGS_RELEASE)

# Root and RestFrames

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)
RFCFLAGS = $(shell restframes-config --cxxflags)
RFGLIBS = $(shell restframes-config --libs)

# General libs

GLIBS = -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-f4591b847fcbe5753bfc5d2b02f57089/lib/ -lgmp -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-f76596f4b83666ac3468f34a5f342677/lib/ -lboost_thread

# BayesianClustering, Eigen, etc.

BCFLAGS = -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/ \
-I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/CGAL-5.6.1/include \
-I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/frugally-deep-master/include \
-I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/FunctionalPlus-master/include \
-I/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math

BCGLIBS = -L/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClustering $(GLIBS)

# Project include and object directories

INCDIR = ./include
OBJDIR = ./obj
$(shell mkdir -p $(OBJDIR))
MAINOBJDIR = ./mainobj
$(shell mkdir -p $(MAINOBJDIR))

IFLAGS = -I../interface

TCALIDIR = KUCMSTimeCaliFiles
TCALI_SRCDIR = $(TCALIDIR)/src
TCALI_INCDIR = $(TCALIDIR)/include

# Source and object files

SRCDIR = ./src

MAIN_CPP := $(SRCDIR)/runKUCMSAodSVSkimmer.cpp $(SRCDIR)/runLocalKUCMSAodSVSkimmer.cpp 
MAIN_CPP += $(TCALI_SRCDIR)/makeKUCMSTimeCalibration.cpp $(TCALI_SRCDIR)/makeKUCMSTimeResolution.cpp
SUPPORT_CPP := $(filter-out $(MAIN_CPP), $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*.cc) ) 
SUPPORT_CPP += $(filter-out $(MAIN_CPP), $(wildcard $(TCALI_SRCDIR)/*.cpp) $(wildcard $(TCALI_SRCDIR)/*.cc))
OBJ_FILES := $(patsubst %,$(OBJDIR)/%,$(notdir $(SUPPORT_CPP:.cpp=.o)))
OBJ_FILES := $(OBJ_FILES:.cc=.o) # handle .cc extension

# Default target

EXE_MAIN = runKUCMSAodSVSkimmer.obj
EXE_LOCAL = runLocalKUCMSAodSVSkimmer.obj
EXE_CALI = makeKUCMSTimeCalibration.obj
EXE_RES = makeKUCMSTimeResolution.obj

# Colored output

COLOR_YELLOW = \033[1;33m
COLOR_GREEN = \033[1;32m
COLOR_RED = \033[1;31m
COLOR_RESET = \033[0m

# Pattern rule for compiling .cpp and .cc

$(MAINOBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo -e "$(COLOR_YELLOW)Compiling main file $<$(COLOR_RESET)"
	$(CC) $(CXXFLAGS) $(CXXFLAGS_COMMON) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) -I$(INCDIR) -I$(TCALI_INCDIR) $(IFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@echo -e "$(COLOR_YELLOW)Compiling $<$(COLOR_RESET)"
	$(CC) $(CXXFLAGS) $(CXXFLAGS_COMMON) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) -I$(INCDIR) -I$(TCALI_INCDIR) $(IFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(TCALI_SRCDIR)/%.cpp
	@echo -e "$(COLOR_YELLOW)Compiling $<$(COLOR_RESET)"
	$(CC) $(CXXFLAGS) $(CXXFLAGS_COMMON) -I$(TCALI_INCDIR) -I$(INCDIR) $(RTCFLAGS) $(RTGLIBS) -c $< -o $@

$(MAINOBJDIR)/%.o: $(TCALI_SRCDIR)/%.cpp
	@echo -e "$(COLOR_YELLOW)Compiling $<$(COLOR_RESET)"
	$(CC) $(CXXFLAGS) $(CXXFLAGS_COMMON) -I$(TCALI_INCDIR) -I$(INCDIR) $(RTCFLAGS) $(RTGLIBS) -c $< -o $@


MAIN_OBJ := $(MAINOBJDIR)/runKUCMSAodSVSkimmer.o
LOCAL_OBJ := $(MAINOBJDIR)/runLocalKUCMSAodSVSkimmer.o
TCALI_OBJ := $(MAINOBJDIR)/makeKUCMSTimeCalibration.o
TRES_OBJ  := $(MAINOBJDIR)/makeKUCMSTimeResolution.o

# Executables

remote: $(EXE_MAIN)
local: $(EXE_LOCAL)
time: $(EXE_CALI) $(EXE_RES)
cali: $(EXE_CALI)
res: $(EXE_RES)

$(EXE_MAIN): $(OBJ_FILES) $(MAIN_OBJ)
	@echo -e "$(COLOR_GREEN)Linking $@$(COLOR_RESET)"
	$(CC) -o $@ $(OBJ_FILES) $(MAIN_OBJ) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCGLIBS) -pthread

$(EXE_LOCAL): $(OBJ_FILES) $(LOCAL_OBJ)
	@echo -e "$(COLOR_GREEN)Linking $@$(COLOR_RESET)"
	$(CC) -o $@ $(OBJ_FILES) $(LOCAL_OBJ) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCGLIBS) -pthread

$(EXE_CALI): $(TCALI_OBJ) $(OBJDIR)/KUCMS_TimeCalibration.o
	@echo -e "$(COLOR_GREEN)Linking $@$(COLOR_RESET)"
	$(CC) -o $@ $(TCALI_OBJ) $(OBJDIR)/KUCMS_TimeCalibration.o $(CXXFLAGS_COMMON) $(RTCFLAGS) $(RTGLIBS) -pthread

$(EXE_RES): $(TRES_OBJ) $(OBJDIR)/KUCMS_TimeCalibration.o
	@echo -e "$(COLOR_GREEN)Linking $@$(COLOR_RESET)"
	$(CC) -o $@ $(TRES_OBJ) $(OBJDIR)/KUCMS_TimeCalibration.o $(CXXFLAGS_COMMON) $(RTCFLAGS) $(RTGLIBS) -pthread


# Condor tarball

configtar: remote
	@echo -e "$(COLOR_GREEN)Creating condor tarball$(COLOR_RESET)"
	cp $(EXE_MAIN) config/
	cp -r ntuple_master_lists/ config/
	cp -r ecal_config/ config/
	tar -czf condor/config.tgz -C config .

all: remote configtar

# debug

print-main:
	@echo "MAINOBJDIR = '$(MAINOBJDIR)'"
	@echo "MAIN_CPP = '$(MAIN_CPP)'"
	@echo "MAIN_OBJ = '$(MAIN_OBJ)'"
	@echo "OBJ_FILES = '$(OBJ_FILES)'"

# Clean

clean:
	@echo -e "$(COLOR_RED)Cleaning objects and executables$(COLOR_RESET)"
	rm -rf $(OBJDIR)/*.o *.obj
	rm -rf condor/config.tgz

# Phony targets

.PHONY: all clean configtar

