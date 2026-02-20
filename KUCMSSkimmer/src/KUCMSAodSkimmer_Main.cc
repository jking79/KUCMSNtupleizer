//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSHelperFunctions.hh"
#include <ctime>

#include "KUCMSAodSVSkimmer.hh"
/*
#include "KUCMSAodSkimmer_Helpers.hpp"
#include "KUCMSAodSkimmer_EvntMetVars.hpp"
#include "KUCMSAodSkimmer_Tracks.hpp"
#include "KUCMSAodSkimmer_Ecal.hpp"
#include "KUCMSAodSkimmer_Gen.hpp"
#include "KUCMSAodSkimmer_Photons.hpp"
#include "KUCMSAodSkimmer_RJR.hpp"
#include "KUCMSAodSkimmer_Electron.hpp"
#include "KUCMSAodSkimmer_Jets.hpp"
#include "KUCMSAodSkimmer_Muons.hpp"
#include "KUCMSAodSkimmer_SV.hpp"
#include "KUCMSAodSkimmer_BHC.hpp"
#include "KUCMSAodSkimmer_TimeSig.hpp"
*/

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
//// KUCMSAodSkimmer class ----------------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------------------------------------

//#define CLSTRMAPS true
#define CLSTRMAPS false

KUCMSAodSkimmer::KUCMSAodSkimmer(){


  ////    Base RJR ------------------------------------------

  LAB = new LabRecoFrame("LAB","LAB");
  S   = new DecayRecoFrame("S","#tilde{S}");

  X2a = new DecayRecoFrame("X2a","#tilde{#chi}_{2a}");
  X2b = new DecayRecoFrame("X2b","#tilde{#chi}_{2b}");	

  Ja = new DecayRecoFrame("Ja","jets_{a}");
  Jb = new DecayRecoFrame("Jb","jets_{b}");
  J1a = new VisibleRecoFrame("J1a","jets_{1a}");
  J2a = new VisibleRecoFrame("J2a","jets_{2a}");
  J1b = new VisibleRecoFrame("J1b","jets_{1b}");
  J2b = new VisibleRecoFrame("J2b","jets_{2b}");

  X1a = new InvisibleRecoFrame("X1a","#tilde{#chi}_{1a}");
  X1b = new InvisibleRecoFrame("X1b","#tilde{#chi}_{1b}");
	
  LAB->SetChildFrame(*S);
  S->AddChildFrame(*X2a);
  S->AddChildFrame(*X2b);
  X2a->AddChildFrame(*X1a);
  X2b->AddChildFrame(*X1b);
  X2a->AddChildFrame(*Ja);
  X2b->AddChildFrame(*Jb);
  Ja->AddChildFrame(*J1a);
  Ja->AddChildFrame(*J2a);
  Jb->AddChildFrame(*J1b);
  Jb->AddChildFrame(*J2b);

  if(!LAB->InitializeTree()){ std::cout << "Problem initializing tree" << std::endl; }
	
  INV = new InvisibleGroup("INV","Invisible System");
  INV->AddFrame(*X1a);
  INV->AddFrame(*X1b);
	
  InvM = new SetMassInvJigsaw("InvM", "Set inv. system mass");
  INV->AddJigsaw(*InvM);
	
  InvEta = new SetRapidityInvJigsaw("InvEta", "Set inv. system rapidity");
  INV->AddJigsaw(*InvEta);
  InvEta->AddVisibleFrames(S->GetListVisibleFrames());
	
  InvSplit = new MinMassesSqInvJigsaw("InvSplit", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}", 2);
  INV->AddJigsaw(*InvSplit);
  InvSplit->AddVisibleFrame(*J1a, 0);
  InvSplit->AddVisibleFrame(*J2a, 0);
  InvSplit->AddVisibleFrame(*J1b, 1);
  InvSplit->AddVisibleFrame(*J2b, 1);
  InvSplit->AddInvisibleFrame(*X1a, 0);
  InvSplit->AddInvisibleFrame(*X1b, 1);
	
  COMB_J =  new CombinatoricGroup("COMB_J", "Combinatoric System of Jets");
  CombSplit_J = new MinMassesSqCombJigsaw("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
  CombSplit_Ja = new MinMassesCombJigsaw("CombSplit_Ja", "Minimize M_{J1a} + M_{J2a} ");
  CombSplit_Jb = new MinMassesCombJigsaw("CombSplit_Jb", "Minimize M_{J1b} + M_{J2b} ");

  COMB_J->AddFrame(*J1a);
  COMB_J->AddFrame(*J2a);
  COMB_J->SetNElementsForFrame(*J1a, 1);
  COMB_J->SetNElementsForFrame(*J2a, 0);
  COMB_J->AddFrame(*J1b);
  COMB_J->AddFrame(*J2b);
  COMB_J->SetNElementsForFrame(*J1b, 1);
  COMB_J->SetNElementsForFrame(*J2b, 0);
	
  COMB_J->AddJigsaw(*CombSplit_J);
  CombSplit_J->AddCombFrames(Ja->GetListVisibleFrames(), 0);
  CombSplit_J->AddCombFrames(Jb->GetListVisibleFrames(), 1);
  //CombSplit_J->AddObjectFrame(*X2a, 0);// included vis & invis frames in split
  //CombSplit_J->AddObjectFrame(*X2b, 1);// check syntax from example
  CombSplit_J->AddObjectFrames(X2a->GetListVisibleFrames(), 0);
  CombSplit_J->AddObjectFrames(X2b->GetListVisibleFrames(), 1);// check syntax from example

  //add bloc do b
  COMB_J->AddJigsaw(*CombSplit_Ja);
  CombSplit_Ja->AddCombFrame(*J1a, 0);
  CombSplit_Ja->AddCombFrame(*J2a, 1);
  CombSplit_Ja->AddObjectFrame(*J1a, 0);
  CombSplit_Ja->AddObjectFrame(*J2a, 1);// check syntax from example
  COMB_J->AddJigsaw(*CombSplit_Jb);
  CombSplit_Jb->AddCombFrame(*J1b, 0);
  CombSplit_Jb->AddCombFrame(*J2b, 1);
  CombSplit_Jb->AddObjectFrame(*J1b, 0);
  CombSplit_Jb->AddObjectFrame(*J2b, 1);// check syntax from example

  if(!LAB->InitializeAnalysis()) std::cout << "Problem initializing analysis tree" << std::endl;

  // -----------  RJRISR -----------------------------------------

  LAB_c = new LabRecoFrame("LAB_c","LAB");
  CM_c  = new DecayRecoFrame("CM_c","CM");
  S_c  = new DecayRecoFrame("S_c","#tilde{S}");

  X2a_c = new DecayRecoFrame("X2a_c","#tilde{#chi}_{2a}");
  X2b_c = new DecayRecoFrame("X2b_c","#tilde{#chi}_{2b}");

  ISR_c = new VisibleRecoFrame("ISR_c","ISR");
  Ja_c = new DecayRecoFrame("Ja_c","jets_{a}");
  Jb_c = new DecayRecoFrame("Jb_c","jets_{b}");
  J1a_c = new VisibleRecoFrame("J1a_c","jets_{1a}");
  J2a_c = new VisibleRecoFrame("J2a_c","jets_{2a}");
  J1b_c = new VisibleRecoFrame("J1b_c","jets_{1b}");
  J2b_c = new VisibleRecoFrame("J2b_c","jets_{2b}");

  X1a_c = new InvisibleRecoFrame("X1a_c","#tilde{#chi}_{1a}");
  X1b_c = new InvisibleRecoFrame("X1b_c","#tilde{#chi}_{1b}");

  LAB_c->SetChildFrame(*CM_c);
  CM_c->AddChildFrame(*S_c);
  CM_c->AddChildFrame(*ISR_c);
  S_c->AddChildFrame(*X2a_c);
  S_c->AddChildFrame(*X2b_c);
  X2a_c->AddChildFrame(*X1a_c);
  X2b_c->AddChildFrame(*X1b_c);
  X2a_c->AddChildFrame(*Ja_c);
  X2b_c->AddChildFrame(*Jb_c);
  Ja_c->AddChildFrame(*J1a_c);
  Ja_c->AddChildFrame(*J2a_c);
  Jb_c->AddChildFrame(*J1b_c);
  Jb_c->AddChildFrame(*J2b_c);

  if(!LAB->InitializeTree()){ std::cout << "Problem initializing tree" << std::endl; }

  INV_c = new InvisibleGroup("INV_c","Invisible System");
  INV_c->AddFrame(*X1a_c);
  INV_c->AddFrame(*X1b_c);

  InvM_c = new SetMassInvJigsaw("InvM_c", "Set inv. system mass");
  INV_c->AddJigsaw(*InvM_c);

  InvEta_c = new SetRapidityInvJigsaw("InvEta_c", "Set inv. system rapidity");
  INV_c->AddJigsaw(*InvEta_c);
  InvEta_c->AddVisibleFrames(S_c->GetListVisibleFrames()); // CM or S ?

  InvSplit_c = new MinMassesSqInvJigsaw("InvSplit_c", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}", 2);
  INV_c->AddJigsaw(*InvSplit_c);
  InvSplit_c->AddVisibleFrame(*J1a_c, 0);
  InvSplit_c->AddVisibleFrame(*J2a_c, 0);
  InvSplit_c->AddVisibleFrame(*J1b_c, 1);
  InvSplit_c->AddVisibleFrame(*J2b_c, 1);
  InvSplit_c->AddInvisibleFrame(*X1a_c, 0);
  InvSplit_c->AddInvisibleFrame(*X1b_c, 1);

  COMB_J_c =  new CombinatoricGroup("COMB_J_c", "Combinatoric System of Jets");
  CombSplit_ISR_c = new MinMassesCombJigsaw("CombSplit_ISR_c", "Minimize M_{T}^{ISR} and M_{T}^{S}");
  CombSplit_J_c = new MinMassesSqCombJigsaw("CombSplit_J_c", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
  CombSplit_Ja_c = new MinMassesCombJigsaw("CombSplit_Ja_c", "Minimize M_{J1a} + M_{J2a} ");
  CombSplit_Jb_c = new MinMassesCombJigsaw("CombSplit_Jb_c", "Minimize M_{J1b} + M_{J2b} ");

  COMB_J_c->AddFrame(*ISR_c);
  COMB_J_c->SetNElementsForFrame(*ISR_c, 1);
  //COMB_J_c->AddFrame(*Ja_c);
  //COMB_J_c->SetNElementsForFrame(*Ja_c, 1);
  //COMB_J_c->AddFrame(*Jb_c);
  //COMB_J_c->SetNElementsForFrame(*Jb_c, 0);
  COMB_J_c->AddFrame(*J1a_c);
  COMB_J_c->AddFrame(*J2a_c);
  COMB_J_c->SetNElementsForFrame(*J1a_c, 0);
  COMB_J_c->SetNElementsForFrame(*J2a_c, 0);
  COMB_J_c->AddFrame(*J1b_c);
  COMB_J_c->AddFrame(*J2b_c);
  COMB_J_c->SetNElementsForFrame(*J1b_c, 0);
  COMB_J_c->SetNElementsForFrame(*J2b_c, 0);

  COMB_J_c->AddJigsaw(*CombSplit_ISR_c);
  CombSplit_ISR_c->SetTransverse();
  CombSplit_ISR_c->AddCombFrame(*ISR_c, 1);
  CombSplit_ISR_c->AddCombFrame(*J1a_c, 0);
  CombSplit_ISR_c->AddCombFrame(*J1b_c, 0);
  CombSplit_ISR_c->AddCombFrame(*J2a_c, 0);
  CombSplit_ISR_c->AddCombFrame(*J2b_c, 0);
  CombSplit_ISR_c->AddObjectFrame(*ISR_c, 1);
  CombSplit_ISR_c->AddObjectFrame(*S_c, 0);

  COMB_J_c->AddJigsaw(*CombSplit_J_c);
  CombSplit_J_c->AddCombFrames(Ja_c->GetListVisibleFrames(), 0);
  CombSplit_J_c->AddCombFrames(Jb_c->GetListVisibleFrames(), 1);
  //CombSplit_J_c->AddObjectFrame(*X2a, 0);// included vis & invis frames in split
  //CombSplit_J_c->AddObjectFrame(*X2b, 1);// check syntax from example
  CombSplit_J_c->AddObjectFrames(X2a_c->GetListVisibleFrames(), 0);
  CombSplit_J_c->AddObjectFrames(X2b_c->GetListVisibleFrames(), 1);// check syntax from example

  //add bloc do b
  COMB_J_c->AddJigsaw(*CombSplit_Ja_c);
  CombSplit_Ja_c->AddCombFrame(*J1a_c, 0);
  CombSplit_Ja_c->AddCombFrame(*J2a_c, 1);
  CombSplit_Ja_c->AddObjectFrame(*J1a_c, 0);
  CombSplit_Ja_c->AddObjectFrame(*J2a_c, 1);// check syntax from example
  COMB_J_c->AddJigsaw(*CombSplit_Jb_c);
  CombSplit_Jb_c->AddCombFrame(*J1b_c, 0);
  CombSplit_Jb_c->AddCombFrame(*J2b_c, 1);
  CombSplit_Jb_c->AddObjectFrame(*J1b_c, 0);
  CombSplit_Jb_c->AddObjectFrame(*J2b_c, 1);// check syntax from example

  if(!LAB_c->InitializeAnalysis()) std::cout << "Problem initializing analysis tree" << std::endl;

  //  ------------  Tree Ploter -----------------

  if( false ){
	
    TreePlot tree_plot("TreePlot","TreePlot");

    for(int t = 0; t < 2; t++){
      tree_plot.SetTree(*LAB);
      tree_plot.Draw("ANA_tree", "Reconstruction Tree");

      tree_plot.SetTree(*COMB_J);
      tree_plot.Draw("ANA_comb", "Combinatoric Jigsaws for jets/sig photons");

      //tree_plot.SetTree(*COMB_L);
      //tree_plot.Draw("ANA_comb_L", "Combinatoric Jigsaws for leps");

      tree_plot.SetTree(*INV);
      tree_plot.Draw("ANA_inv", "Invisible Jigsaws");
 
    }//<<>>for(int t = 0; t < 2; t++)

    tree_plot.WriteOutput("trees.root");

  }//<<>>treeplot on off

  // -------------  default root branch intilization -------------

  doGenInfoBase = false;
  doSVsBase = true;
  doNewSigBase = false;
  doHTLPathsBase = true;

  // -----------  Time Calibration -----------------

  // Cali Tags : Tags for calibrations to use 
  std::string r2EOY( "EG_EOY_MINI" ); 
  std::string r2Fall17AOD( "RunIIFall17DRPremix" ); 
  std::string r2Fall17MINIAOD( "RunIIFall17MiniAODv2" ); 
  std::string r2UL( "UL_R2_MINI" ); 

  timeCali = new KUCMS_TimeCalibration();
  timeCali->setTag(r2UL);
  std::map<UInt_t,kucms_DetIDStruct> detidmap = timeCali->getDetIDMapRef();
  for(auto it = detidmap.begin(); it != detidmap.end(); it++){
	  kucms_DetIDStruct detidstruct = it->second;
	  int ieta = detidstruct.i2;
	  int iphi = detidstruct.i1;
	  _detidmap[it->first] = std::make_pair(ieta, iphi);
   }

  // ------------------------------------------------------------------------------------------

  // setup lumi json map

  loadLumiJson("config/json/Cert_Collisions2024_378981_386951_Golden.json");
  loadLumiJson("config/json/Cert_Collisions2023_366442_370790_Golden.json");
  loadLumiJson("config/json/Cert_Collisions2022_355100_362760_Golden.json");
  loadLumiJson("config/json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.json");
  loadLumiJson("config/json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt");
  loadLumiJson("config/json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
  
  // condor event segmenting varibles : used to run over subset of events for condor jobs

  isLocal = false;

  _evti = -1;
  _evtj = -1;

  // event varibles

    dataSetKey = "single";
    xsctn = 1.f;
    gmass = -1.f;
    xmass = -1.f;
    mcwgt = 1.f;;
    mctype = 0;  // what type of input PD :  0 MC (AODSIM), 1 DATA (AOD), if we need fastSim ect, add new enrty here
    tctag = "none";

	// input tree names

	disphoTreeName = "tree/llpgtree";
	configTreeName = "tree/configtree";
	inFilePath = "none";
	inFileName = "none";
	listDirPath = "none";
	eosDirPath = "none";
	outFileName = "rjr_skimfile.root";

    // event varibles

	doSVs = true;
    doBHC = true;
	genSigPerfectFlag = false;
    noSVorPhoFlag = false;
    useEvtGenWgtFlag = true;
    hasGenInfoFlag = true;
    nEvents = 0, 
    nSelectedEvents = 0;
    sumEvtGenWgt = 0;
	setGenInfoBase(hasGenInfoFlag); 
	setDoSVsBase(doSVs);
	setNewSigBase(false);
	setHTLPathsBase(false);

	// BNC Intiation

    _ca.SetVerbosity(-1); //can turn on to see low-level warnings
    _ca.SetDetIDs(_detidmap);
    _ca.SetCNNModel("config/json/KU-CNN_detector_1000epochs_archsmall3_2017and2018_CMS.json");
    _ca.SetBarrelDNNModel("config/json/KU-DNN_photonID_JetHT18RunBC_EGamma18RunC_SMSGlGl_StandardNorm_shape_1000epochs_large8_barrelOnly.json");
    _ca.SetEndcapDNNModel("config/json/KU-DNN_photonID_JetHT18RunBC700chunks_EGamma18RunC_SMSGlGl_StandardNorm_isoShape_1000epochs_large8_endcapOnly.json");

}//<<>>KUCMSAodSkimmer::KUCMSAodSkimmer()

KUCMSAodSkimmer::~KUCMSAodSkimmer(){

    delete LAB;
    delete S;
    delete X2a;
    delete X2b;
    delete Ja;
    delete Jb;
    delete X1a;
    delete X1b;
    delete J1a;
    delete J2a;
    delete J1b;
    delete J2b;      

    delete INV;
    delete InvM;
    delete InvEta;
    delete InvSplit;
      
    delete COMB_J;
    delete CombSplit_J;
    delete CombSplit_Ja;
    delete CombSplit_Jb;

    delete LAB_c;
    delete CM_c;
    delete S_c;
    delete X2a_c;
    delete X2b_c;
	delete ISR_c;
    delete Ja_c;
    delete Jb_c;
    delete X1a_c;
    delete X1b_c;
    delete J1a_c;
    delete J2a_c;
    delete J1b_c;
    delete J2b_c;

    delete INV_c;
    delete InvM_c;
    delete InvEta_c;
    delete InvSplit_c;

    delete COMB_J_c;
    delete CombSplit_J_c;
    delete CombSplit_Ja_c;
    delete CombSplit_Jb_c;

    delete timeCali;
      
}//<<>>KUCMSAodSkimmer::~KUCMSAodSkimmer()

void KUCMSAodSkimmer::ProcessMainLoop( TChain* fInTree, TChain* fInConfigTree ){

    // --- setup for main loop
    TTree* fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
    TTree* fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");

	setGenInfoBase( hasGenInfoFlag );
	setDoSVsBase( doSVs );
	setNewSigBase( doNewSigBase );
	setHTLPathsBase( doHTLPathsBase );
	Init( fInTree );
	
    ////Init( fInTree, hasGenInfoFlag, doSVs );
    auto nEntries = fInTree->GetEntries();

    //loop over events
    if( _evti < 0 || _evti == _evtj){ _evti = 0; _evtj = nEntries; }
	else {
		if( _evtj > nEntries ){ _evtj = nEntries; } //cap at max number of entries
		if( _evti > nEntries ){ cout << "Starting event " << _evti << " above # of entries in tree " << nEntries << " returning." << endl; return; }
	}//<<>> if( _evti < 0 ) else
	//if( not doSVs ){ _evtj = 2500000; _evti = 1500000; }
	int nEventsProcessed = _evtj - _evti;
    initHists();
    setOutputBranches(fOutTree);

    startJobs(); // clear && init count varibles

    // ---   do config input

	ProcessConfigTree( fInConfigTree );
	if( not isLocal ) ProcessConfigFile();

    // --  main loop

    std::cout << "Setting up For Main Loop." << std::endl;
	int nForPrecent =  ( nEventsProcessed < 10 ) ? 1 : 10;
	if( not doSVs ){ std::cout << " ::::  Doing noSv for EGamma " << std::endl; } //nForPrecent = 100; }
    int loopCounter( nEventsProcessed / nForPrecent );
    std::cout << "Processing " << nEntries << " entries." << std::endl;
    nEvents = nEventsProcessed;//save # of events actually ran over so that when files are hadded, this should total nEntries in TChain
    cout << "Running over events " << _evti << " to " << _evtj << endl;
    cout << "loopCounter " << loopCounter << " nEventsProcessed " << nEventsProcessed << " nForPrecent " << nForPrecent << endl;
	int count_of_events = 0;
    for (Long64_t centry = _evti; centry < _evtj; centry++){
        if( centry%loopCounter == 0 ){
            std::time_t now = std::time(nullptr);
            std::string curtime = std::ctime(&now);
            curtime.pop_back();
            std::cout << "Proccessed " << centry << " of " << nEntries << " entries at " << curtime << std::endl;
        }//<<>>if( centry%loopCounter == 0 )
 	    auto entry = fInTree->LoadTree(centry);
        if(DEBUG) std::cout << " -- Getting Branches " << std::endl;
        getBranches( entry );
        ////getBranches( entry, hasGenInfoFlag, doSVs );
		//std::cout << " -- Check Valid Lumi with mctype " << mctype << " run " << Evt_run << " block " << Evt_luminosityBlock << std::endl;
        if( mctype==1 && not isValidLumisection( Evt_run, Evt_luminosityBlock ) ) continue;
		//std::cout << " --- Valid Lumi Processing Event " << std::endl;

        geCnts.clear();
        geVars.clear();
        hemBits.clear();
        hasHemObj = false;

        if( genSigPerfectFlag ) geVars.set( "genSigPerfect", 1 ); else geVars.set( "genSigPerfect", 0 );
        if( noSVorPhoFlag ) geVars.set( "noSVorPho", 1 ); else geVars.set( "noSVorPho", 0 );

        if(DEBUG) std::cout << " -- Event Loop " << std::endl;
		count_of_events++;
        auto saveToTree = eventLoop(entry);
        if( saveToTree ){ fOutTree->Fill(); }

    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
	if( count_of_events != nEvents ) std::cout << " ! nEvents " << nEvents << " not = count " << count_of_events << " !" << std::endl;

    // -- main end jobs

    fillConfigTree( fConfigTree );

    endJobs();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;
    std::cout << "Writing to " << outFileName  << endl;
    TFile* fOutFile = new TFile( outFileName.c_str(), "RECREATE" );
    fOutFile->cd();

    fOutTree->Write();
    fConfigTree->Write();

    // --- write & delete histograms

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

    if( CLSTRMAPS ){
        nMaps = 0;
        for( int it = 0; it < nEBEEMaps; it++ ){
            ebeeMapP[it]->Write(); delete ebeeMapP[it];
            ebeeMapT[it]->Write(); delete ebeeMapT[it];
            ebeeMapR[it]->Write(); delete ebeeMapR[it];
        }//<<>>for( int it = 0; it < nEBEEMaps; it++ )
    }//<<>>f( clsttrMaps )

    // ----- closeout open sys


    fOutFile->Close();

    delete fInTree;
    delete fOutTree;
    delete fConfigTree;
    delete fOutFile;

}//<<>>void KUCMSAodSkimmer::ProcessMainLoop()


int KUCMSAodSkimmer::ProcessFile(string infilename, TChain*& fInTree, TChain*& fInConfigTree){

    const std::string disphotreename = "tree/llpgtree";
    const std::string configtreename = "tree/configtree";
    if(gSystem->AccessPathName(infilename.c_str())){
      cout << "Error: File " << infilename << " not found" << endl;
      return -1;
    }
  
    if( DEBUG ) std:: cout << "InFile: " << infilename << std::endl;
    if( DEBUG ) std:: cout << "Key: " << dataSetKey << std::endl;
    if( DEBUG ) std:: cout << "XSec: " << xsctn << std::endl;
    if( DEBUG ) std:: cout << "GM: " << gmass << std::endl;
    if( DEBUG ) std:: cout << "XM: " << xmass << std::endl;
    if( DEBUG ) std:: cout << "MCw: " << mcwgt << std::endl;
    if( DEBUG ) std:: cout << "MCt: " << mctype << std::endl;
    if( DEBUG ) std:: cout << "tcTag: " << tctag << std::endl;

    std::cout << "Processing Events for : " << infilename << std::endl;
    fInTree = new TChain(disphotreename.c_str());
    fInConfigTree = new TChain(configtreename.c_str());
    fInTree->Add(infilename.c_str());
    fInConfigTree->Add(infilename.c_str());
    if(DEBUG) std::cout << "--  adding file: " << infilename << std::endl; //else std::cout << ".";
    return 0;


}



int KUCMSAodSkimmer::ProcessFilelist(string eosdir, string infilename, TChain*& fInTree, TChain*& fInConfigTree){

    const std::string disphotreename = "tree/llpgtree";
    const std::string configtreename = "tree/configtree";


    if(gSystem->AccessPathName(infilename.c_str())){
      cout << "Error: File " << infilename << " not found" << endl;
      return -1;
    }
  
    if( DEBUG ) std:: cout << "InFile: " << infilename << std::endl;
    if( DEBUG ) std:: cout << "Key: " << dataSetKey << std::endl;
    if( DEBUG ) std:: cout << "XSec: " << xsctn << std::endl;
    if( DEBUG ) std:: cout << "GM: " << gmass << std::endl;
    if( DEBUG ) std:: cout << "XM: " << xmass << std::endl;
    if( DEBUG ) std:: cout << "MCw: " << mcwgt << std::endl;
    if( DEBUG ) std:: cout << "MCt: " << mctype << std::endl;
    if( DEBUG ) std:: cout << "tcTag: " << tctag << std::endl;

    std::cout << "Processing Events for : " << infilename << std::endl;
    fInTree = new TChain(disphotreename.c_str());
    fInConfigTree = new TChain(configtreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilename << " >> " << fInTree << std::endl;
    std::string str;
    if( not DEBUG ) std::cout << "--  adding files";
    int nfiles = 0;	
    if( dataSetKey !=  "single" ){
      std::ifstream infile(infilename);
      while( std::getline( infile, str ) ){
	      if(str.find("#") != string::npos) continue;
	nfiles++;
	//if( skipCnt != 0 && ( nfiles%skipCnt != 0 ) ) continue;
    string sample_str = str.substr(0,str.find("/"));
	string prefix = str.substr(str.find("/")+1);
	string match = "kucmsntuple_";
	string proc = prefix.substr(match.size());
        proc = proc.substr(0,proc.find("_"));	
	string inpath = prefix.substr(0,prefix.find(sample_str,prefix.find(sample_str)+sample_str.size())-1);
	prefix = prefix.substr(0,prefix.rfind(proc));
	prefix.replace(prefix.rfind("_"),1,"/");
	inpath = prefix+str;
	auto tfilename = eosdir + inpath;
	fInTree->Add(tfilename.c_str());
	fInConfigTree->Add(tfilename.c_str());
	if(DEBUG) std::cout << "--  adding file: " << tfilename << std::endl; //else std::cout << ".";
	//if(DEBUG) break;
      }//<<>>while (std::getline(infile,str))
    } else { // single infile and not a list of infile
    	string sample_str = str.substr(0,str.find("/"));
	string prefix = str.substr(str.find("/")+1);
	string inpath = prefix.substr(0,prefix.find(sample_str,prefix.find(sample_str)+sample_str.size())-1)+"/";
      auto tfilename = eosdir + inpath + infilename;
      fInTree->Add(tfilename.c_str());
      fInConfigTree->Add(tfilename.c_str());
      nfiles++;			
    }//<<>>if( key !=  "test" )
    if( not DEBUG ) std::cout << std::endl;
    if( nfiles == 0 ){ std::cout << " !!!!! no input files !!!!! " << std::endl; return -1; }
    return 0;    
}//<<>>int KUCMSAodSkimmer::ProcessFilelist(string eosdir, string infilename, TChain*& fInTree, TChain*& fInConfigTree)

int KUCMSAodSkimmer::ProcessFilelistOfLists(string eosdir, vector<string> processed_strings, TChain*& fInTree, TChain*& fInConfigTree){

    const std::string disphotreename = "tree/llpgtree";
    const std::string configtreename = "tree/configtree";

    auto inpath = processed_strings[0];
    cout << "inpath " << inpath << endl;
    auto infile_str = processed_strings[1];
    dataSetKey = processed_strings[2];
    xsctn = std::stof( processed_strings[3] );
    gmass = std::stof( processed_strings[4] );
    xmass = std::stof( processed_strings[5] );
    mcwgt = std::stof( processed_strings[6] );
    mctype = std::stoi( processed_strings[7] );
    tctag = processed_strings[8];
    if( DEBUG ) std:: cout << "InFile: " << infile_str << std::endl;
    if( DEBUG ) std:: cout << "XSec: " << xsctn << std::endl;
    if( DEBUG ) std:: cout << "GM: " << gmass << std::endl;
    if( DEBUG ) std:: cout << "XM: " << xmass << std::endl;
    if( DEBUG ) std:: cout << "MCw: " << mcwgt << std::endl;
    if( DEBUG ) std:: cout << "MCt: " << mctype << std::endl;
    if( DEBUG ) std:: cout << "tcTag: " << tctag << std::endl;
		
    if(gSystem->AccessPathName(infile_str.c_str())){
      cout << "Error: File " << infile_str << " not found" << endl;
      return -1;
    }
  
    std::cout << "Processing Events for : " << infile_str << std::endl;
    fInTree = new TChain(disphotreename.c_str());
    fInConfigTree = new TChain(configtreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infile_str << " >> " << fInTree << std::endl;
    std::string str;
    if( not DEBUG ) std::cout << "--  adding files";
    int nfiles = 0;	
    if( dataSetKey !=  "single" ){
      std::ifstream infile(infile_str);
      while( std::getline( infile, str ) ){
	      if(str.find("#") != string::npos) continue;
	nfiles++;
	//if( skipCnt != 0 && ( nfiles%skipCnt != 0 ) ) continue;
	auto tfilename = eosdir + inpath + str;
	fInTree->Add(tfilename.c_str());
	fInConfigTree->Add(tfilename.c_str());
	if(DEBUG) std::cout << "--  adding file: " << tfilename << std::endl; //else std::cout << ".";
	//if(DEBUG) break;
      }//<<>>while (std::getline(infile,str))
    } else { // single infile and not a list of infile
      auto tfilename = eosdir + inpath + infile_str;
      fInTree->Add(tfilename.c_str());
      fInConfigTree->Add(tfilename.c_str());
      nfiles++;			
    }//<<>>if( key !=  "test" )
    if( not DEBUG ) std::cout << std::endl;
    if( nfiles == 0 ){ std::cout << " !!!!! no input files !!!!! " << std::endl; return -1; }
    return 0;    
}//<<>>int KUCMSAodSkimmer::ProcessFilelistOfLists(string eosdir, vector<string> processed_strings, TChain*& fInTree, TChain*& fInConfigTree)

void KUCMSAodSkimmer::ProcessConfigTree( TChain* fInConfigTree ){

    int nTotEvts;
    int nFltrdEvts;
    float sumEvtWgt;
    float sumFltrdEvtWgt;
    int nMetFltrdEvts;
    int nPhoFltrdEvts;
   
	TBranch *b_nTotEvts;
    TBranch *b_nFltrdEvts;
    TBranch *b_sumEvtWgt;
    TBranch *b_sumFltrdEvtWgt;
    TBranch *b_nMetFltrdEvts;
    TBranch *b_nPhoFltrdEvts;

    fInConfigTree->SetBranchAddress("nTotEvts", &nTotEvts, &b_nTotEvts);
    fInConfigTree->SetBranchAddress("nFltrdEvts", &nFltrdEvts, &b_nFltrdEvts);
    fInConfigTree->SetBranchAddress("sumEvtWgt", &sumEvtWgt, &b_sumEvtWgt);
    fInConfigTree->SetBranchAddress("sumFltrdEvtWgt", &sumFltrdEvtWgt, &b_sumFltrdEvtWgt);
    fInConfigTree->SetBranchAddress("nMetFltrdEvts", &nMetFltrdEvts, &b_nMetFltrdEvts);
    fInConfigTree->SetBranchAddress("nPhoFltrdEvts", &nPhoFltrdEvts, &b_nPhoFltrdEvts);

    auto nConfigEntries = fInConfigTree->GetEntries();
    std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
    configCnts.clear();
    configWgts.clear();

	if( _evti == 0 ){

      for (Long64_t centry = 0; centry < nConfigEntries; centry++){

        auto entry = fInConfigTree->LoadTree(centry);

        b_nTotEvts->GetEntry(entry);   //!
        b_nFltrdEvts->GetEntry(entry);   //!
        b_sumEvtWgt->GetEntry(entry);   //!  
        b_sumFltrdEvtWgt->GetEntry(entry);   //!
        b_nMetFltrdEvts->GetEntry(entry);   //!
        b_nPhoFltrdEvts->GetEntry(entry);   //!

        cutflow["nTotEvts"] += nTotEvts;
        cutflow["nFltrdEvts"] += nFltrdEvts;

        configCnts["nTotEvts"] += nTotEvts;
        configCnts["nFltrdEvts"] += nFltrdEvts;
        configWgts["sumEvtWgt"] += useEvtGenWgtFlag ? sumEvtWgt : nTotEvts;
        configWgts["sumFltrdEvtWgt"] += useEvtGenWgtFlag ? sumFltrdEvtWgt : nFltrdEvts;
        configCnts["nMetFltrdEvts"] += nMetFltrdEvts;
        configCnts["nPhoFltrdEvts"] += nPhoFltrdEvts;


      }//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)

	} else {

      cutflow["nTotEvts"] = 0;
      cutflow["nFltrdEvts"] = 0;

      configCnts["nTotEvts"] = 0;
      configCnts["nFltrdEvts"] = 0;
      configWgts["sumEvtWgt"] = 0;
      configWgts["sumFltrdEvtWgt"] = 0;
      configCnts["nMetFltrdEvts"] = 0;
      configCnts["nPhoFltrdEvts"] = 0;
	
	}//<<>if( _evti == 0 )
	

    std::cout << "configCnts ( ";
    for( auto item : configCnts ){
      std::cout << item.first <<  " " << item.second << " ";
    }//<<>>for( auto item : configInfo )
    std::cout << ")" << std::endl;

    std::cout << "configWgts ( ";
    for( auto item : configWgts ){
      std::cout << item.first <<  " " << item.second << " ";
    }//<<>>for( auto item : configInfo )
    std::cout << ")" << std::endl;

}//<<>>void KUCMSAodSkimmer::ProcessConfigTree(TChain* fInConfigTree)

void KUCMSAodSkimmer::ProcessConfigFile(){

	std::cout << " -- reading config/EventCount.txt " << std::endl;
    std::ifstream infile("config/EventCount.txt");
    if( !infile ){
        std::cerr << "  --  Could not open config/EventCount.txt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
        return;
    }//<<>>if( !infile )

	configData.clear();

    std::string key;
    float total, sum;
    while( infile >> key >> total >> sum ){ configData[key] = {total, sum}; }

    // Print everything
    for( const auto& kv : configData ){
        std::cout << " -- ConfigFile : " << kv.first << " : total=" << kv.second.first << ", sum=" << kv.second.second << std::endl;
    }//<<>>for (const auto& kv : data)//<<>>if( DEBUG )

}//<<>>void KUCMSAodSkimmer::ProcessConfigTreeAndFile(TChain* fInConfigTree)

void KUCMSAodSkimmer::kucmsAodSkimmer( std::string infile, std::string outfilename){

    std::cout << "Processing Input Lists for : " << infile << std::endl;

    //add path to input file
    //TTree* fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
    //TTree* fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");
    TChain* fInTree = nullptr;
    TChain* fInConfigTree = nullptr;
    int ret = ProcessFile( infile, fInTree, fInConfigTree);
    if(ret < 0) return;

	SetOutFileName(outfilename);
    ProcessMainLoop(fInTree, fInConfigTree);
    std::cout << "Finished processing events for : " << infile << std::endl;
  	std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

}//<<>>void KUCMSAodSkimmer::kucmsAodSkimmer( std::string eosdir, std::string infilelist, std::string outfilename)

void KUCMSAodSkimmer::kucmsAodSkimmer_Filelist( std::string eosdir, std::string infilelist, std::string outfilename){

    std::cout << "Processing Input Lists for : " << infilelist << std::endl;

    //add path to input file
    //TTree* fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
    //TTree* fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");
    TChain* fInTree = nullptr;
    TChain* fInConfigTree = nullptr;
    int ret = ProcessFilelist(eosdir, infilelist, fInTree, fInConfigTree);
    if(ret < 0) return;

	SetOutFileName(outfilename);
	cout << "is in tree null " << (fInTree == nullptr) << endl; 
    ProcessMainLoop(fInTree, fInConfigTree);
    std::cout << "Finished processing events for : " << infilelist << std::endl;
  	std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

}//<<>>void KUCMSAodSkimmer::kucmsAodSkimmer( std::string eosdir, std::string infilelist, std::string outfilename)


void KUCMSAodSkimmer::kucmsAodSkimmer_listsOfLists( std::string eosdir, std::string infilelist, std::string outfilename){
	//std::string inpath, infiles, key, 
	std::string masterstr; 
	//int mct;
	//float crossSection, gmsblam, gmsbct, mcw;
	std::cout << "Processing Input Lists for : " << infilelist << std::endl;
	std::ifstream masterInfile(infilelist);
	
	string listdir = infilelist.substr(0,infilelist.find("/")+1); 
	while( std::getline( masterInfile, masterstr ) ){
	
	  if( DEBUG ) std:: cout << masterstr << std::endl;
	  if( masterstr[0] == '#' ) continue;
	  if( masterstr == " " ) continue;
	  auto instrs = splitString( masterstr, " " );
	  if( DEBUG ) std:: cout << instrs.size() << std::endl;
	  if( instrs.size() < 9 ) continue;
	  //add path to input file
	  instrs[1] = listdir+instrs[1];
	  TTree* fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
	  TTree* fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");
	  TChain* fInTree = nullptr;
	  TChain* fInConfigTree = nullptr;
	  int ret = ProcessFilelistOfLists(eosdir, instrs, fInTree, fInConfigTree);
	  if(ret < 0) return;
	  SetOutFileName(outfilename);
	  ProcessMainLoop(fInTree, fInConfigTree);	
    	std::cout << "Finished processing events for : " << instrs[1] << std::endl;
	}//<<>>while (std::getline(infile,str))

  std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;


}//<<>>void KUCMSAodSkimmer::kucmsAodSkimmer_listsOfLists( std::string eosdir, std::string infilelist, std::string outfilename)

void KUCMSAodSkimmer::kucmsAodSkimmer_local( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename ){

  // set "meta" parameters

  eosDirPath = eosdir;
  listDirPath = listdir;

  //isLocal = true;

  // ----- parse input batch file ---------------------------------------------------------------------

  std::string masterstr;
  std::cout << "Processing Input Lists for : " << listdir+infilelist << std::endl;
  std::ifstream masterInfile(listdir+infilelist);
  //while( masterInfile >> inFilePath >> inFileName >> key >> crossSection >> gmsblam >> gmsbct >> mcwgt >> mctype ){
  while( std::getline( masterInfile, masterstr ) ){

    // ---- parse input params

    //if( DEBUG ) std:: cout << masterstr << std::endl;
    if( masterstr[0] == '#' ) continue;
    if( masterstr == "" ) continue;
    auto instrs = splitString( masterstr, " " );
    if( DEBUG ) std::cout << instrs.size() << std::endl;
    if( instrs.size() != 9 ) continue;

    inFilePath = instrs[0];
    inFileName = instrs[1];

    auto key = instrs[2];
    auto crossSection = std::stof( instrs[3] );
    auto gmsbgm = std::stof( instrs[4] );
    auto gmsbxm = std::stof( instrs[5] );
    auto mcw = std::stof( instrs[6] );
    auto mct = std::stoi( instrs[7] );
    auto tct = instrs[8];

    if( DEBUG ) std::cout << "InPath: " << inFilePath << std::endl;
    if( DEBUG ) std::cout << "InFile: " << inFileName << std::endl;
    if( DEBUG ) std::cout << "Key: " << key << std::endl;
    if( DEBUG ) std::cout << "XSec: " << crossSection << std::endl;
    if( DEBUG ) std::cout << "GM: " << gmsbgm << std::endl;
    if( DEBUG ) std::cout << "XM: " << gmsbxm << std::endl;
    if( DEBUG ) std::cout << "MCw: " << mcw << std::endl;
    if( DEBUG ) std::cout << "MCt: " << mct << std::endl;
    if( DEBUG ) std::cout << "tcTag: " << tctag << std::endl;

    dataSetKey = key;
    xsctn = crossSection;
    gmass = gmsbgm; // = 0 if not gmsb
    xmass = gmsbxm; // = 0 if not gmsb
    mcwgt = mcw; // default 1
    mctype = mct; // 0 for MC, 1 for data
    tctag = tct; // r2_ul18 Data, r2_ul18_mc MC

    // ---- load files

    std::cout << "Processing Events for : " << inFileName << std::endl;
    TChain* fInTree = new TChain(disphoTreeName.c_str());
    TChain* fInConfigTree = new TChain(configTreeName.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << listDirPath+inFileName << " >> " << fInTree << std::endl;
    std::string str;
    if( not DEBUG ) std::cout << "--  adding files";
    int nfiles = 0;
	int skipCnt = 20;
    if( skipCnt != 0 ) std::cout << "-- !! Skipping every " << skipCnt << std::endl;
    if( dataSetKey !=  "Single" ){
        std::ifstream infile(listDirPath+inFileName);
        while( std::getline( infile, str ) ){
            if(str.find("#") != string::npos) continue;
            nfiles++;
            if( skipCnt != 0 && ( nfiles%skipCnt != 0 ) ) continue;
            auto tfilename = eosDirPath + inFilePath + str;
            fInTree->Add(tfilename.c_str());
            fInConfigTree->Add(tfilename.c_str());
            if(DEBUG) std::cout << "--  adding file: " << tfilename << std::endl; else std::cout << ".";
        }//<<>>while (std::getline(infile,str))
    } else { // single infile and not a list of inFileName
        auto tfilename = eosDirPath + inFilePath + inFileName;
        fInTree->Add(tfilename.c_str());
        fInConfigTree->Add(tfilename.c_str());
        nfiles++;
    }//<<>>if( key !=  "test" )
    if( not DEBUG ) std::cout << std::endl;
    if( nfiles == 0 ){ std::cout << " !!!!! no input files !!!!! " << std::endl; return; }

    // ------ Do Main Loop
	_evti = -1; _evti = -1;
    auto ext = splitString( inFileName, "." );
    std::string extOutFileName( ext[0] + outfilename );
	SetOutFileName( extOutFileName );
    ProcessMainLoop( fInTree, fInConfigTree );

  // ------- while loop -----------------------------------------------------
  }//<<>>while (std::getline(infile,str))

  std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

}//<<>>void KUCMSAodSkimmer::kucmsAodSkimmer_local(std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename)

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// event loop and startup jobs 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree ){

	setEvtVarMetBranches( fOutTree );
    setTrackBranches( fOutTree );
    setEcalBranches( fOutTree );
    setGenBranches( fOutTree );
    setPhotonBranches( fOutTree );
    setRJRBranches( fOutTree );
	setRJRISRBranches( fOutTree );
    setElectronBranches( fOutTree );
    setJetsBranches( fOutTree );
    setMuonsBranches( fOutTree );
    setSVBranches( fOutTree );
    //setBCBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setOutputBranches(fOutTree)

void KUCMSAodSkimmer::startJobs(){

    std::cout << "Runing StartJobs." << std::endl;

    nEvents = 0;
    nSelectedEvents = 0;

    //sumEvtGenWgt = 0.0;

    configCnts.clear();
    configWgts.clear();
    configCnts["nTotEvts"] = 0;
    configCnts["nFltrdEvts"] = 0;
    configWgts["sumEvtWgt"] = 0;
    configWgts["sumFltrdEvtWgt"] = 0;
    configCnts["nMetFltrdEvts"] = 0;
    configCnts["nPhoFltrdEvts"] = 0;

    cutflow.clear();
    cutflow["nTotEvts"] = 0;
    cutflow["nFltrdEvts"] = 0;
    cutflow["met150"] = 0;
    cutflow["m_gt2jets"] = 0;
    cutflow["mj_gt1phos"] = 0;
    cutflow["mjp_leadPhoPt30"] = 0;
    cutflow["sel_m"] = 0;
    cutflow["sel_j"] = 0;
    cutflow["sel_p"] = 0;
    cutflow["sel_ppt"] = 0;

	std::cout << "Finished StartJobs." << std::endl;

};//<<>>void KUCMSAodSkimmer::startJobs()

bool KUCMSAodSkimmer::eventLoop( Long64_t entry ){

  // counts events and saves event varibles
  // --------------------------------------

  bool debugl1 = false;
  
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processRechits" << std::endl;
  processRechits();// must be done before rechitID to Iter map used
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processMet" << std::endl;
  processMet();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processPhotons" << std::endl;
  processPhotons();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processElectrons" << std::endl;
  processElectrons();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processMuons" << std::endl;
  if( doSVs ) processMuons();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processJets" << std::endl;
  processJets();
  if( debugl1 && doSVs ) std::cout << "KUCMSAodSkimmer : processSV" << std::endl;
  if( doSVs ) processSV();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processTracks" << std::endl;
  if( doSVs ) processTracks();
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : processEvntVars" << std::endl;
  processEvntVars();// process last to catch Hem issue
  if( debugl1 && hasGenInfoFlag ) std::cout << "KUCMSAodSkimmer : processGenParticles" << std::endl;
  if( hasGenInfoFlag ){ processGenParticles(); }
  
  // select events to process and store
  //--------------------------------------
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : eventSelection" << std::endl;
  auto saveToTree = eventSelection();
  if( saveToTree ){ 

    //if( doBHC ){ //processBHCPhotons();
    if( debugl1 && doBHC ) std::cout << "KUCMSAodSkimmer : processBHCJets" << std::endl;
	if( doBHC ){ processBHCJets(); }
    if( debugl1 ) std::cout << "KUCMSAodSkimmer : processRJR 0" << std::endl;
	processRJR(0,true); 
    //if( debugl1 ) std::cout << "KUCMSAodSkimmer : processRJR 1" << std::endl;
	//processRJR(1,false); 
    if( debugl1 ) std::cout << "KUCMSAodSkimmer : processRJR ISR" << std::endl;
	processRJRISR();

  }//<<>>if( saveToTree )
  if( debugl1 ) std::cout << "KUCMSAodSkimmer : Return saveToTree" << std::endl;
  return saveToTree;
  
}//<<>>void KUCMSAodSkimmer::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
// decide which events to save to tree
//------------------------------------------------------------------------------------------------------------

bool KUCMSAodSkimmer::eventSelection(){
  // select which events to save and fill output branches

  if( DEBUG ) std::cout << "Event selection." << std::endl;
  // determine if we want to save event
  bool dobase = ( geVars("noSVorPho") == 1 ) ? true : false;

  float evtMet = geVars("cmet");

  int nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
  int nQJets = geCnts("nQJets");
  int vetoJets = geCnts("jetEventVeto");

  float nSelPhotons = geCnts("nSelPhotons"); //selPhotons.getUIBranchValue("nSelPhotons");
  auto lSigPhoIndx = geVars("lSigPhoIndx");
  auto slSigPhoIndx = geVars("slSigPhoIndx");
  float leadPhoPt = ( lSigPhoIndx > -1 ) ? geVects("selPhoPt").at(lSigPhoIndx) : 0;
  float subLeadPhoPt = ( slSigPhoIndx > -1 ) ? geVects("selPhoPt").at(slSigPhoIndx) : 0;

  bool hasLepSV = doSVs ? geVars("nSVLep") > 0 : false;
  bool hasHadSV = doSVs ? geVars("nSVHad") > 0 : false;
  bool hasSV = hasLepSV || hasHadSV;

  bool met100 = evtMet >= 100;
  bool met150 = evtMet >= 150;
  bool gt1phos = lSigPhoIndx > -1;
  bool gt2jets = nSelJets >= 2;
  bool gt2qjets = nQJets >= 2;
  bool allgjets = vetoJets < 1;

  bool gt2phos = lSigPhoIndx > -1 and slSigPhoIndx > -1;
  bool leadPhoPt70 = leadPhoPt >= 70;
  bool leadPhoPt30 = leadPhoPt >= 30;
  bool subLeadPhoPt40 = subLeadPhoPt >= 40; 

  bool basesel = met100 && gt2jets && allgjets && gt2qjets;
  //bool svsel = basesel && hasSV;	
  bool phosel = basesel && ( ( gt1phos && leadPhoPt30 ) || hasSV );
  //auto evtSelected = leadPhoPt70 && subLeadPhoPt40 && gt2jets && gt2phos;

  bool evtSelected = dobase ? basesel : phosel;
  //if( hasHemObj ) evtSelected = false;
  if( not doSVs ) evtSelected = gt2qjets && gt2jets; 

  if( met150 ) cutflow["met150"]++;
  if( met150 && gt2jets ) cutflow["m_gt2jets"]++;
  if( met150 && gt2jets && gt1phos ){ cutflow["mj_gt1phos"]++; cutflow["sel_ppt"]++; }
  if( met150 && gt2jets && gt1phos && leadPhoPt30 ) cutflow["mjp_leadPhoPt30"]++;

  if( met150 && gt2jets && leadPhoPt30 ) cutflow["sel_p"]++;
  if( met150 && gt1phos && leadPhoPt30 ) cutflow["sel_j"]++;
  if( gt2jets && gt1phos && leadPhoPt30 ) cutflow["sel_m"]++;

  if( DEBUG ){ 
    if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl; 
    else std::cout << " --------------------- Event Failed." << std::endl;}
  //if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl;
  if( evtSelected ){ nSelectedEvents++; return true; } 
  else return false;

}//<<>>void KUCMSAodSkimmer::eventSelection()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::endJobs(){ 

  std::cout << "Running EndJobs." << std::endl;
  hist1d[11]->Divide(hist1d[10]);
  hist1d[12]->Divide(hist1d[10]);
  hist1d[13]->Divide(hist1d[10]);
  hist1d[14]->Divide(hist1d[10]);

}//void KUCMSAodSkimmer::endJobs()

void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree ){

  std::cout << "Filling ConfigTree." << std::endl;
  TBranch *nEventBranch = fConfigTree->Branch( "nEvents", &nEvents );
  nEventBranch->Fill();
  TBranch *nSelectedEventsBranch = fConfigTree->Branch( "nSelectedEvents", &nSelectedEvents );
  nSelectedEventsBranch->Fill();
  TBranch *sKeyBranch = fConfigTree->Branch( "sKey", &dataSetKey ); 
  sKeyBranch->Fill();
  TBranch *sCrossSectionBranch = fConfigTree->Branch( "sCrossSection", &xsctn );
  sCrossSectionBranch->Fill();
  TBranch *sLambdaBranch = fConfigTree->Branch( "sGMSBGravMass", &gmass );
  sLambdaBranch->Fill();
  TBranch *sCTauBranch = fConfigTree->Branch( "sGMSBChi1Mass", &xmass );
  sCTauBranch->Fill();
  TBranch *sMCWgtBranch = fConfigTree->Branch( "sMCWgt", &mcwgt );
  sMCWgtBranch->Fill();
  TBranch *sMCTypeBranch = fConfigTree->Branch( "sMCType", &mctype );
  sMCTypeBranch->Fill();
  //TBranch *sumEvtGenWgtBranch = fConfigTree->Branch( "sumEvtGenWgt", &sumEvtGenWgt );
  //sumEvtGenWgtBranch->Fill();

  for( auto item : cutflow ){
    std::string bname = "cf_" + item.first;
    TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
    cfBranch->Fill();
  }//<<>>for( auto item : configInfo )

  for( auto item : configCnts ){
	//std::cout << " -- configCnts : " << item.first << " local " << isLocal << std::endl;
    if( item.first == "nTotEvts" && not isLocal ){
		int nTotEvtsKey = configData[dataSetKey].first;
		//std::cout << " -- Filling nTotEvts : " << dataSetKey << " " << nTotEvtsKey << " for " << item.first << std::endl;
		TBranch *cfBranch = fConfigTree->Branch( "nTotEvts", &nTotEvtsKey );
		cfBranch->Fill();
    } else {
    	std::string bname = item.first;
    	TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
    	cfBranch->Fill();
	}//<<>>
  }//<<>>for( auto item : configInfo )

  for( auto item : configWgts ){
	//std::cout << " -- configWgts : " << item.first << " local " << isLocal << std::endl;
    if( item.first == "sumEvtWgt" && not isLocal ){
		float sumEvtWgtKey = useEvtGenWgtFlag ? configData[dataSetKey].second : configData[dataSetKey].first;
        //std::cout << " -- Filling sumEvtWgt : " << dataSetKey << " " << sumEvtWgtKey << " for " << item.first << std::endl;
		TBranch *cfBranch = fConfigTree->Branch( "sumEvtWgt", &sumEvtWgtKey );
		cfBranch->Fill();
    } else {
    	std::string bname = item.first;
    	TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
    	cfBranch->Fill();
    }//<<>>
  }//<<>>for( auto item : configInfo )

  fConfigTree->Fill();

}//<<>>void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree, std::string key )

//------------------------------------------------------------------------------------------------------------
// Init histograms ------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = new TH1D("ecalrhenergy0", "RecHit Energy;rechit E [GeV]",1000,0,1000);
    hist1d[1] = new TH1D("ecalrhenergy1", "RecHit Energy;rechit E [GeV]",200,0,20);
    hist1d[2] = new TH1D("ecalrhenergy2", "RecHit Energy;rechit E [GeV]",200,0,2);
    hist1d[3] = new TH1D("sctype","SC !Orig, Orig, OOT, Excl",4,0,4);
    hist1d[4] = new TH1D("scorigtype","Orig+OOT, Orig+!OOT, Orig+Exc, !Orig+OOT, !Orig+Exc, !Orig+!OOT",6,0,6);
    hist1d[5] = new TH1D("ecalrhtime", "RecHit Time;rechit t [ns]",200,-10,10);
    hist1d[6] = new TH1D("scenergy","SC rawEnergy;energy [GeV]",200,0,1000);
    hist1d[7] = new TH1D("trckpt","Track pt;Track pt [GeV]",500,0,1000);
    hist1d[8] = new TH1D("NRecHitsInSCs","NRecHitsInSCs;NRecHitsInSCs",100,0,100);
    hist1d[9] = new TH1D("sceta","SC Eta;Eta",200,-10,10);

	hist1d[10] = new TH1D("elept","Electron Pt [GeV];Electron Pt [GeV];a.u.",100,0,1000);
    hist1d[11] = new TH1D("elept_lv","LooseVeto;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[12] = new TH1D("elept_v","Veto;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[13] = new TH1D("elept_l","Loose;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[14] = new TH1D("elept_m","Medium;Electron Pt [GeV];Eff",100,0,1000);

    hist1d[20] = new TH1D("erhnamps","Pho SC Nomilized RH Ampsig;Norm Ampsig",100,0,1);

    ////hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    ////hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[1] = new TH2D("scorigcross", "SC Types;Mian-Orig;Std-OOT-Excluded", 2, 0, 2, 3, 0, 3); 
    hist2d[2] = new TH2D("rjrPtS_v_rjrdPhiSI_t0", "rjrDPhiSI_v_rjrPtS Type 0;rjrDphiSI;rjrPtS", 32, 0, 3.2, 60, 0, 3000);
    hist2d[3] = new TH2D("rjrPtS_v_rjrdPhiSI_t1", "rjrDPhiSI_v_rjrPtS Type 1;rjrDphiSI;rjrPtS", 32, 0, 3.2, 60, 0, 3000);

	//------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

    //hist3d[0] = new TH3D("phoNRH_ClR9_phoID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",200,0,200,100,0,1,8,0,4);

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------

    if( CLSTRMAPS ){
		nMaps = 0;
		for(int it=0; it<nEBEEMaps; it++){
			fMap[it] = false;
			std::string label(";iEta;iPhi");
    	    std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
    	    ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
    	    std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
    	    ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
			std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
    	    ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
		}//<<>>for(int it=0; it<nEBEEMaps; it++)
	}//<<>>if( CLSTRMAPS )

}//<<>>void KUCMSAodSkimmer::initHists()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
