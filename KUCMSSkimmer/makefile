##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

export $LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/lib


RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

RFCFLAGS    = $(shell restframes-config --cxxflags)
RFGLIBS     = $(shell restframes-config --libs)

BCFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/include -frounding-math
BCGLIBS += -L/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/BayesianClustering/lib -lBayesianClusteringi

#CGFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/
CGFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/CGAL-5.6.1/include
FDFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/frugally-deep-master/include
EGFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/eigen
FPFLAGS += -I/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/FunctionalPlus-master/include

CC=g++

p1=runKUCMSAodSVSkimmer
p2=spaceholder

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCFLAGS) $(BCGLIBS) $(FDFLAGS) $(EGFLAGS) $(FPFLAGS) $(CGFLAGS)

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)
