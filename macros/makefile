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
p13=skimRJRMultiHistMaker3
p14=skimRJRMultiHistMaker3_batch
p15=skimRJRMultiHistMaker4_batch
p16=skimRJRMultiHistMaker5_batch
p17=runKUCMSAodSVSkimmer
p18=skimRJRMultiHistMakerZ5_batch
p19=skimRJRMultiHistMakerBatch6

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

$(p13): $(p13).cpp
	$(CC) -o $(p13).obj $(p13).cpp $(RTCFLAGS) $(RTGLIBS)

$(p14): $(p14).cpp
	$(CC) -o $(p14).obj $(p14).cpp $(RTCFLAGS) $(RTGLIBS)

$(p15): $(p15).cpp
	$(CC) -o $(p15).obj $(p15).cpp $(RTCFLAGS) $(RTGLIBS)

$(p16): $(p16).cpp
	$(CC) -o $(p16).obj $(p16).cpp $(RTCFLAGS) $(RTGLIBS)

$(p17): $(p17).cpp
	$(CC) -o $(p17).obj $(p17).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS)

$(p18): $(p18).cpp
	$(CC) -o $(p18).obj $(p18).cpp $(RTCFLAGS) $(RTGLIBS)

$(p19): $(p19).cpp
	$(CC) -o $(p19).obj $(p19).cpp $(RTCFLAGS) $(RTGLIBS)


