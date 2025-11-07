//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSVSkimmer.hh"
#include "KUCMSHelperFunctions.hh"

//#define DEBUG true
#define DEBUG false

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//  do any processing and calulations for objects and save values to output varibles 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// void KUCMSAodSkimmer::processTemplate(){
// 	
// 		Clear out branch vector varibles &/or initilize other out branch vars
//	------------------------------------------------
// 		Do any calculations / cuts
//  ------------------------------------------------
//		Fill out branch varibles
//
//}//<<>>void KUCMSAodSkimmer::processTemplate()

enum FinalStateType {kNHad1, kNHadGe2, kNLep1, kNLepGe2, kMixed, kNoSV};

void KUCMSAodSkimmer::processSV(){


  bool doEVSVs = true;
  if( not doEVSVs ) std::cout << " ------ !!!!!!!!!! Fill SV Branches turned off !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  
  selSV.clearBranches(); // <<<<<<<   must do
  
  // SVs
  if( DEBUG ) std::cout << "Finding SV Vars" << std::endl;
  
  int nSVs = 0;
  if( doEVSVs )  nSVs = Vertex_mass->size();
  int nLsv = 0;
  int nHsv = 0;
  int nEle = 0;
  int nMuon = 0;
  //Z-Window Count
  int nHadPassZWindow(0);
  int nElePassZWindow(0);
  int nMuPassZWindow(0);
  //Final State Counters
  int nLepPassCR(0), nLepPassSR(0), nLepPassCRTight(0), nLepPassSRTight(0);
  int nHadPassCR(0), nHadPassSR(0), nHadPassCRTight(0), nHadPassSRTight(0);
  
  if( doEVSVs ){
    for( int svit = 0; svit < nSVs; svit++ ){
      
      const float mass((*Vertex_mass)[svit]);
      const float dxy((*Vertex_dxy)[svit]);
      const float dxyError((*Vertex_dxyError)[svit]);
      const float dxySig(dxy/dxyError);
      const float x((*Vertex_x)[svit]),  y((*Vertex_y)[svit]), z((*Vertex_z)[svit]);
      const float cosTheta((*Vertex_cosTheta)[svit]);
      const uInt ntracks((*Vertex_nTracks)[svit]);
      const bool peleid((*Vertex_passLooseElectronID)[svit]);
      const bool pmuonid((*Vertex_passLooseMuonID)[svit]);

      ROOT::Math::XYZVectorF vec3D(x, y, z);
      const float etaOrigin(vec3D.eta());

      // Apply baseline selection
      if(dxy < 2 || (dxy < 15 && fabs(z) > 15) || (dxy > 15 && fabs(etaOrigin) > 1.) || cosTheta > 0.995)
	continue;
      
      MaterialVeto vetoTool("fullMask60cm_Data.root");
      const bool tightOnZSelection((peleid && mass > 60) || (pmuonid && mass > 70) || (ntracks >= 5 && mass > 20 && mass/ntracks > 1));

      string svType;
      
      if(ntracks == 2) {
	svType = "Leptonic";
	if( peleid || pmuonid ) {
	  nLsv++;

	  if(dxySig < 1000) {
	    nLepPassCR++;
	    if(tightOnZSelection)
	      nLepPassCRTight++;
	  }
	  else {
	    nLepPassSR++;
	    if(tightOnZSelection)
              nLepPassSRTight++;
	  }
	}
	if(peleid) {
	  nEle++;
	  selSV.fillBranch("LeptonicSV_electronIndex", int(svit));
	  if(tightOnZSelection)
	    nElePassZWindow++;
	}
	if(pmuonid) {
	  nMuon++;
	  selSV.fillBranch("LeptonicSV_muonIndex", int(svit));
	  if(tightOnZSelection)
            nMuPassZWindow++;
	}

	selSV.fillBranch("LeptonicSV_passMuonID", pmuonid);
	selSV.fillBranch("LeptonicSV_passElectronID", peleid);
	if(hasGenInfoFlag) 
	  selSV.fillBranch( "LeptonicSV_isGold", (*Vertex_isGold)[svit]);
      }
      else if( ( ntracks >= 5 ) && ( mass/ntracks > 1 ) ) {
	svType = "Hadronic"; 
	nHsv++;

	if(dxySig < 1000) {
	  nHadPassCR++;
	  if(tightOnZSelection)
	    nHadPassCRTight++;
	}
	else {
	  nHadPassSR++;
	  if(tightOnZSelection)
	    nHadPassSRTight++;
	}

	selSV.fillBranch( "HadronicSV_massOverNtracks", mass/ntracks);
	if(hasGenInfoFlag)
	  selSV.fillBranch( "HadronicSV_matchRatio", (*Vertex_matchRatio)[svit]);
	if(tightOnZSelection)
	  nHadPassZWindow++;

      }
      else continue;

      const float p((*Vertex_p)[svit]), decayAngle((*Vertex_decayAngle)[svit]);
      
      selSV.fillBranch( svType+"SV_nTracks", int(ntracks));
      selSV.fillBranch( svType+"SV_mass", mass);
      selSV.fillBranch( svType+"SV_pOverE", float(p/sqrt(p*p + mass*mass)));
      selSV.fillBranch( svType+"SV_decayAngle", decayAngle);
      selSV.fillBranch( svType+"SV_cosTheta", cosTheta);
      selSV.fillBranch( svType+"SV_x", x);
      selSV.fillBranch( svType+"SV_y", y);
      selSV.fillBranch( svType+"SV_z", z);
      selSV.fillBranch( svType+"SV_dxy", dxy);
      selSV.fillBranch( svType+"SV_dxySig", dxy/dxyError);
      selSV.fillBranch( svType+"SV_passMaterialVeto",  vetoTool.PassVeto(x, y));
      selSV.fillBranch( svType+"SV_tightZWindowSelection", tightOnZSelection);
    }//<<>>for( svit = 0; svit < nSVs; scit++ )
  }//<<>>if( doEVSVs )
  
  selSV.fillBranch( "SV_nLeptonic", nLsv );
  selSV.fillBranch( "SV_nHadronic", nHsv );
  selSV.fillBranch( "SV_nElectron", nEle );
  selSV.fillBranch( "SV_nMuon", nMuon );

  selSV.fillBranch( "SV_nLeptonPassZWindow", unsigned(nMuPassZWindow + nElePassZWindow) );
  selSV.fillBranch( "SV_nHadronPassZWindow", unsigned(nHadPassZWindow) );
  selSV.fillBranch( "SV_nElectronPassZWindow", unsigned(nElePassZWindow) );
  selSV.fillBranch( "SV_nMuonPassZWindow", unsigned(nMuPassZWindow) );

  const bool anyInSR(nHadPassSR>0 || nLepPassSR>0);
  //Flags final state CR Loose
  selSV.fillBranch( "passNHad1SelectionCRLoose", bool(nHadPassCR==1 && nLepPassCR==0 && !anyInSR));
  selSV.fillBranch( "passNHadGe2SelectionCRLoose", bool(nHadPassCR>1 && nLepPassCR==0 && !anyInSR));
  selSV.fillBranch( "passNLep1SelectionCRLoose", bool(nHadPassCR<1 && nLepPassCR==1 && !anyInSR));
  selSV.fillBranch( "passNLepGe2SelectionCRLoose", bool(nHadPassCR<1 && nLepPassCR>1 && !anyInSR));
  selSV.fillBranch( "passGe1HadAndLepSelectionCRLoose", bool(nHadPassCR>0 && nLepPassCR>0 && !anyInSR));
  //Flags final	state CR Tight
  selSV.fillBranch( "passNHad1SelectionCRTight", bool(nHadPassCRTight==1 && nLepPassCRTight<1 && !anyInSR));
  selSV.fillBranch( "passNHadGe2SelectionCRTight", bool(nHadPassCRTight>1 && nLepPassCR<1 && !anyInSR));
  selSV.fillBranch( "passNLep1SelectionCRTight", bool(nHadPassCRTight<1 && nLepPassCRTight==1 && !anyInSR));
  selSV.fillBranch( "passNLepGe2SelectionCRTight", bool(nHadPassCRTight<1 && nLepPassCRTight>1 && !anyInSR));
  selSV.fillBranch( "passGe1HadAndLepSelectionCRTight", bool(nHadPassCRTight>0 && nLepPassCRTight>0 && !anyInSR));
  //Flags final	state SR Loose
  selSV.fillBranch( "passNHad1SelectionSRLoose", bool(nHadPassSR==1 && nLepPassSR<1));
  selSV.fillBranch( "passNHadGe2SelectionSRLoose", bool(nHadPassSR>1 && nLepPassSR<1));
  selSV.fillBranch( "passNLep1SelectionSRLoose", bool(nHadPassSR<1 && nLepPassSR==1));
  selSV.fillBranch( "passNLepGe2SelectionSRLoose", bool(nHadPassSR<1 && nLepPassSR>1));
  selSV.fillBranch( "passGe1HadAndLepSelectionSRLoose", bool(nHadPassSR>0 && nLepPassSR>0));
  //Flags final state SR Tight
  selSV.fillBranch( "passNHad1SelectionSRTight", bool(nHadPassSRTight==1 && nLepPassSRTight<1));
  selSV.fillBranch( "passNHadGe2SelectionSRTight", bool(nHadPassSRTight>1 && nLepPassSR<1));
  selSV.fillBranch( "passNLep1SelectionSRTight", bool(nHadPassSRTight<1 && nLepPassSRTight==1));
  selSV.fillBranch( "passNLepGe2SelectionSRTight", bool(nHadPassSRTight<1 && nLepPassSRTight>1));
  selSV.fillBranch( "passGe1HadAndLepSelectionSRTight", bool(nHadPassSRTight>0 && nLepPassSRTight>0));
  
  geVars.set("nSVLep", nLsv );
  geVars.set("nSVHad", nHsv );
  
}//<<>>void KUCMSAodSkimmer::processSV()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setSVBranches( TTree* fOutTree ){

  for(const string svType : {"Hadronic", "Leptonic"}) {
    selSV.makeBranch( svType+"SV_nTracks", VINT);
    selSV.makeBranch( svType+"SV_pOverE", VFLOAT);
    selSV.makeBranch( svType+"SV_decayAngle", VFLOAT);
    selSV.makeBranch( svType+"SV_cosTheta", VFLOAT);
    selSV.makeBranch( svType+"SV_mass", VFLOAT);
    selSV.makeBranch( svType+"SV_x", VFLOAT);
    selSV.makeBranch( svType+"SV_y", VFLOAT);
    selSV.makeBranch( svType+"SV_z", VFLOAT);
    selSV.makeBranch( svType+"SV_dxy", VFLOAT);
    selSV.makeBranch( svType+"SV_dxySig", VFLOAT);
    selSV.makeBranch( svType+"SV_passMaterialVeto", VBOOL );
    selSV.makeBranch( svType+"SV_tightZWindowSelection", VBOOL );
  }

  selSV.makeBranch( "LeptonicSV_nTotal", UINT);
  selSV.makeBranch( "LeptonicSV_electronIndex", VINT);
  selSV.makeBranch( "LeptonicSV_muonIndex", VINT);
  selSV.makeBranch( "LeptonicSV_passMuonID", VBOOL);
  selSV.makeBranch( "LeptonicSV_passElectronID", VBOOL);
    
  selSV.makeBranch( "HadronicSV_nTotal", UINT);
  selSV.makeBranch( "HadronicSV_massOverNtracks", VFLOAT);
  selSV.makeBranch( "SV_nLeptonic", INT );
  selSV.makeBranch( "SV_nHadronic", INT );
  selSV.makeBranch( "SV_nElectron", INT );
  selSV.makeBranch( "SV_nMuon", INT );

  selSV.makeBranch( "SV_nLeptonPassZWindow", UINT );
  selSV.makeBranch( "SV_nHadronPassZWindow", UINT );
  selSV.makeBranch( "SV_nElectronPassZWindow", UINT );
  selSV.makeBranch( "SV_nMuonPassZWindow", UINT );

  for(const string stringency : {"Loose", "Tight"}) {
    for(const string binType : {"SR", "CR"}) {
      selSV.makeBranch( "passNHad1Selection"+binType+stringency, BOOL );
      selSV.makeBranch( "passNHadGe2Selection"+binType+stringency, BOOL );
      selSV.makeBranch( "passNLep1Selection"+binType+stringency, BOOL );
      selSV.makeBranch( "passNLepGe2Selection"+binType+stringency, BOOL );
      selSV.makeBranch( "passGe1HadAndLepSelection"+binType+stringency, BOOL );
    }
  }
    
  if( hasGenInfoFlag ){
    selSV.makeBranch( "HadronicSV_matchRatio", VFLOAT);
    selSV.makeBranch( "LeptonicSV_isGold", VBOOL);
  }
  selSV.attachBranches( fOutTree );
    
}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

