##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

############################################################################################3
#
#  Need to place the following line into your bash_profile or run from command line
#  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/lib
#
###############################################################3

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

RFCFLAGS = $(shell restframes-config --cxxflags)
RFGLIBS = $(shell restframes-config --libs)

GLIBS8 += -L/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/boost/1.80.0-0eec8fa902b6ad913364ef1d6ac9c8d0/lib/ -lboost_thread
GLIBS8 += -L/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/gmp-static/6.2.1-b661550d99fe45fff5976e08ec5e3b5f/lib/ -lgmp
GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-f4591b847fcbe5753bfc5d2b02f57089/lib/ -lgmp
GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-f76596f4b83666ac3468f34a5f342677/lib/ -lboost_thread


BCFLAGS += -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/CGAL-5.6.1/include
BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/frugally-deep-master/include
BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/FunctionalPlus-master/include
BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math
BCGLIBS += -L/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClustering
#BCFLAGS += -I$(CMSSW_BASE)/src/BayesianClustering/include -frounding-math
#BCGLIBS += -L$(CMSSW_BASE)/src/BayesianClustering/lib -lBayesianClustering
BCGLIBS += $(GLIBS)

BCFLAGS8 += -I/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
#BCFLAGS8 += -I$(CMSSW_BASE)/src/CGAL-5.6.1/include
#BCFLAGS8 += -I$(CMSSW_BASE)/src/frugally-deep-master/include
#BCFLAGS8 += -I$(CMSSW_BASE)/src/FunctionalPlus-master/include
BCFLAGS8 += -I/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
BCFLAGS8 += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math
BCFLAGS8 += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/CGAL-5.6.1/include
BCFLAGS8 += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/frugally-deep-master/include
BCFLAGS8 += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/FunctionalPlus-master/include
BCGLIBS8 += -L/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClustering
BCGLIBS8 += $(GLIBS8)

CC=g++

OPTS = -O3 -march=native -ffast-math -funroll-loops -flto
OPTS += -fno-trapping-math -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno

p1=runKUCMSAodSVSkimmer
p2=fillKUCMSTimeCalibration
p3=runLocalKUCMSAodSVSkimmer
p4=runLocalCKUCMSAodSVSkimmer

TCOBJ = KUCMS_TimeCalibration.so
TCINC = -I/KUCMSTimeCaliFiles/include

# Automatically find all .cpp files under KUCMSTimeCaliFiles/src/
#SRC_FILES := $(wildcard KUCMSTimeCaliFiles/src/*.cpp)
#OBJ_FILES := $(patsubst KUCMSTimeCaliFiles/src/%.cpp, KUCMSTimeCaliFiles/obj/%.o, $(SRC_FILES))
#
## Ensure obj directory exists
#$(shell mkdir -p KUCMSTimeCaliFiles/obj)
#
## Compile rule for all .cpp in that directory
#KUCMSTimeCaliFiles/obj/%.o: KUCMSTimeCaliFiles/src/%.cpp
#	$(CC) $(RTCFLAGS) $(RTGLIBS) $(INC) -c $< -o $@

all: $(p1) 
TimeCali: KUCMS_TimeCalibration.so

KUCMS_TimeCalibration.so: KUCMSTimeCaliFiles/src/KUCMS_TimeCalibration.cpp KUCMSTimeCaliFiles/include/*.hh
	$(CC) $(RTCFLAGS) $(RTGLIBS) $(TCINC) -fPIC -c $< -o $@

$(p1): $(p1).cpp $(TCOBJ)
	$(CC) -o $(p1).obj $(p1).cpp $(TCOBJ) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS)  

$(p2): $(p2).cpp $(TCOBJ)
	$(CC) -o $(p2).obj $(p2).cpp $(TCOBJ) $(RTCFLAGS) $(RTGLIBS)

$(p3): $(p3).cpp $(TCOBJ)
	$(CC) -o $(p3).obj $(p3).cpp $(TCOBJ) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) $(OPTS)

$(p4): $(p4).cpp $(TCOBJ)
	$(CC) -o $(p4).obj $(p4).cpp $(TCOBJ) $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) $(OPTS) 

configtar:
	cp runKUCMSAodSVSkimmer.obj config/
	cp -r ntuple_master_lists/ config/
#	cp -r json/ config/
	cp -r ecal_config/ config/
	tar -czf condor/config.tgz config/

# Clean up
#clean:
#	rm -f *.o *.obj
