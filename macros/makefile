##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

RFCFLAGS    = $(shell restframes-config --cxxflags)
RFGLIBS     = $(shell restframes-config --libs)

CC=g++

p1=runKUCMSAodSkimmer
p2=phoSelSkimHistMaker
p3=skimMultiHistMaker
p4=crStyleSkimHistMaker
p5=skimSVHistMaker
p6=wc_ku_InterCali_aveRecHit_isoot_kucc_fliter_red
p7=stackSkimHistMaker
p8=calibration_averhtime_isoot_filter_pho
p9=calibration_averhtime_isoot_filter_pho_ew
p10=skimRJRMultiHistMaker
p11=skimRJRMultiHistMaker2
p12=wc_ku_plot2dResolution_wGausErr_mf

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS)

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)

$(p3): $(p3).cpp
	$(CC) -o $(p3).obj $(p3).cpp $(RTCFLAGS) $(RTGLIBS)

$(p4): $(p4).cpp
	$(CC) -o $(p4).obj $(p4).cpp $(RTCFLAGS) $(RTGLIBS)

$(p5): $(p5).cpp
	$(CC) -o $(p5).obj $(p5).cpp $(RTCFLAGS) $(RTGLIBS)

$(p6): $(p6).cpp
	$(CC) -o $(p6).obj $(p6).cpp $(RTCFLAGS) $(RTGLIBS)

$(p7): $(p7).cpp
	$(CC) -o $(p7).obj $(p7).cpp $(RTCFLAGS) $(RTGLIBS)

$(p8): $(p8).cpp
	$(CC) -o $(p8).obj $(p8).cpp $(RTCFLAGS) $(RTGLIBS)

$(p9): $(p9).cpp
	$(CC) -o $(p9).obj $(p9).cpp $(RTCFLAGS) $(RTGLIBS)

$(p10): $(p10).cpp
	$(CC) -o $(p10).obj $(p10).cpp $(RTCFLAGS) $(RTGLIBS)

$(p11): $(p11).cpp
	$(CC) -o $(p11).obj $(p11).cpp $(RTCFLAGS) $(RTGLIBS)

$(p12): $(p12).cpp
	$(CC) -o $(p12).obj $(p12).cpp $(RTCFLAGS) $(RTGLIBS)

