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

GLIBS += -L/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/boost/1.80.0-0eec8fa902b6ad913364ef1d6ac9c8d0/lib/ -lboost_thread
GLIBS += -L/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/gmp-static/6.2.1-b661550d99fe45fff5976e08ec5e3b5f/lib/ -lgmp
#GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-f4591b847fcbe5753bfc5d2b02f57089/lib/ -lgmp
#GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-f76596f4b83666ac3468f34a5f342677/lib/ -lboost_thread

#BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math
BCFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math
#BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/CGAL-5.6.1/include
BCFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/CGAL-5.6.1/include
#BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/frugally-deep-master/include
BCFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/frugally-deep-master/include
#BCFLAGS += -I/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
BCFLAGS += -I/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-f9c27fce684e89466e2ef07869cd264d/include/eigen3/
#BCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_0_13/src/FunctionalPlus-master/include
BCFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/FunctionalPlus-master/include
#BCGLIBS += -L/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClustering
BCGLIBS += -L/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClustering
BCGLIBS += $(GLIBS)

TCFLAGS += -I/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSTimeCaliFiles

CC = g++

OPTS = -O3 -march=native -ffast-math -funroll-loops -flto 
OPTS += -fno-trapping-math -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno

p1=runKUCMSAodSVSkimmer
p2=fillKUCMSTimeCalibration
p3=spaceholder

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) $(OPTS) 

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS) $(OPTS)

$(p3): $(p3).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)
