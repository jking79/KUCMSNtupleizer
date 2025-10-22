// -*- C++ -*-
//
// Package:    EventCounter
// Class:      EventCounter
//
/**\class EventCounter EventCounter.cc
 *
 *  A simple EDAnalyzer that counts how many events are processed by CMSSW.
 *  This count includes all events seen by the process, regardless of filtering
 *  or output selection.
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>
//#include "KUCMSBranchManager.hh"

#include <atomic>
#include <iostream>

typedef unsigned int uInt;

class KUCMSCounter : public edm::one::EDAnalyzer<> {

public:

  explicit KUCMSCounter(const edm::ParameterSet&) {}
  ~KUCMSCounter() override {}

  void beginJob(){ 

	edm::Service<TFileService> fs;
	ecTree = fs->make<TTree>("eventCountTree","KUCMSNtuple");	
	//ecBranches.makeBranch("nTotEvts","nTotEvts",UINT);
	//ecBranches.attachBranches(ecTree);
	
	eventCountBranch = ecTree->Branch( "nTotEvts", &eventCount_ );
	eventCount_ = 0;

  }

  void analyze(const edm::Event&, const edm::EventSetup&) override { ++eventCount_; }

  void endJob() override {

    //edm::LogVerbatim("EventCounter") << "Total events processed: " << eventCount_;
    //std::cout << "\n============================\n"
    //          << "Total events processed: " << eventCount_ << "\n"
    //          << "============================\n";

    //ecBranches.clearBranches();
	//ecBranches.fillBranch("nTotEvts",eventCount_);
	eventCountBranch->Fill();
	ecTree->Fill();

  }

private:

  //static std::atomic<unsigned long long> eventCount_;
  uInt eventCount_;
  TTree *ecTree;
  TBranch *eventCountBranch;
  //KUCMSBranchManager ecBranches;

};

//// Definition of static member
////std::atomic<unsigned long long> KUCMSEventCounter::eventCount_{0};

DEFINE_FWK_MODULE(KUCMSCounter);


