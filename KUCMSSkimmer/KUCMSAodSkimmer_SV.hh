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
  
  if( doEVSVs ){
    for( int svit = 0; svit < nSVs; svit++ ){
      
      const float mass((*Vertex_mass)[svit]);
      const float dxy((*Vertex_dxy)[svit]);
      const float x((*Vertex_x)[svit]),  y((*Vertex_y)[svit]), z((*Vertex_z)[svit]);
      const float cosTheta((*Vertex_cosTheta)[svit]);
      const uInt ntracks((*Vertex_nTracks)[svit]);
      const bool peleid((*Vertex_passLooseElectronID)[svit]);
      const bool pmuonid((*Vertex_passLooseMuonID)[svit]);

      ROOT::Math::XYZVectorF vec3D(x, y, z);
      const float etaOrigin(vec3D.eta());

      if(dxy < 2 || (dxy < 15 && fabs(z) > 15) || (dxy > 15 && fabs(etaOrigin) > 1.) || cosTheta > 0.995)
	continue;
      
      MaterialVeto vetoTool("fullMask60cm_Data.root");
      const bool tightOnZSelection((peleid && mass > 60) || (pmuonid && mass > 70) || (ntracks >= 5 && mass > 20 && mass/ntracks > 1));

      string svType;
      
      if(ntracks == 2) {
	svType = "Leptonic";
	if( peleid || pmuonid ) nLsv++;
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

	selSV.fillBranch( "HadronicSV_massOverNtracks", mass/ntracks);
	if(hasGenInfoFlag)
	  selSV.fillBranch( "HadronicSV_matchRatio", (*Vertex_matchRatio)[svit]);
	if(tightOnZSelection)
	  nHadPassZWindow++;

      }
      else continue;

      const float p((*Vertex_p)[svit]), dxyError((*Vertex_dxyError)[svit]), decayAngle((*Vertex_decayAngle)[svit]);
      
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
  
  if( hasGenInfoFlag ){
    selSV.makeBranch( "HadronicSV_matchRatio", VFLOAT);
    selSV.makeBranch( "LeptonicSV_isGold", VBOOL);
  }
  selSV.attachBranches( fOutTree );
    
}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

