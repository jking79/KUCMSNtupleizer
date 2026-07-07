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

  if( DEBUG ) std::cout << "Finding SV Vars" << std::endl;

  int nSVs = 0;
  if( doEVSVs ) nSVs = Vertex_mass->size();
  int nLsv = 0;
  int nHsv = 0;
  //Analysis channel flag accumulators (from AnalysisConfig YAMLs)
  bool anyHadMassGeq15{false}, anyHadMassGeq15DxySigGeq800{false}, anyHadDxySigGeq800{false};
  bool anyLepMassGeq15{false}, anyLepMassGeq10{false}, anyLepMassGeq10DxySigGeq500{false}, anyLepDxySigGeq500{false};
  //Validation region accumulators
  bool anyHadMassGeq15DxySigLt100{false}, anyHadMassGeq15DxySigGeq100Lt800{false}, anyHadDxySigLt200{false};
  bool anyLepMassGeq10DxySigLt100{false}, anyLepMassGeq10DxySigGeq100Lt500{false}, anyLepDxySigLt50{false};

  if( doEVSVs ){
    for( int svit = 0; svit < nSVs; svit++ ){

      const float mass((*Vertex_mass)[svit]);
      const float dxy((*Vertex_dxy)[svit]);
      const float dxyError((*Vertex_dxyError)[svit]);
      if (!std::isfinite(dxyError) || dxyError <= 0.0f) continue;
      const float dxySig(dxy/dxyError);
      const float x((*Vertex_x)[svit]),  y((*Vertex_y)[svit]), z((*Vertex_z)[svit]);
      const float cosTheta((*Vertex_cosTheta)[svit]);
      const uInt ntracks((*Vertex_nTracks)[svit]);
      const bool peleid((*Vertex_passLooseElectronID)[svit]);
      const bool pmuonid((*Vertex_passLooseMuonID)[svit]);
      const bool hadronicId(mass/ntracks > 1 && ntracks>3);

      ROOT::Math::XYZVectorF vec3D(x, y, z);
      const float etaOrigin(vec3D.eta());

      // Apply baseline selection
      if(dxy < 1 || (dxy < 15 && fabs(z) > 15) || (dxy > 15 && fabs(etaOrigin) > 1.) || cosTheta > 0.995 || cosTheta <= -0.8)
	continue;

      MaterialVeto vetoTool("fullMask60cm_Data.root");
      const bool tightOnZSelection((peleid && mass > 60) || (pmuonid && mass > 70) || (ntracks >= 5 && mass > 20 && hadronicId));

      string svType;

      if(ntracks == 2) {
	svType = "Leptonic";
	if( peleid || pmuonid ) {
	  nLsv++;
	  if(mass >= 15)                  anyLepMassGeq15 = true;
	  if(mass >= 10)                  anyLepMassGeq10 = true;
	  if(mass >= 10 && dxySig >= 500) anyLepMassGeq10DxySigGeq500 = true;
	  if(dxySig >= 500)               anyLepDxySigGeq500 = true;
	  if(mass >= 10 && dxySig < 100)  anyLepMassGeq10DxySigLt100 = true;
	  if(mass >= 10 && dxySig >= 100 && dxySig < 500) anyLepMassGeq10DxySigGeq100Lt500 = true;
	  if(dxySig < 50)                 anyLepDxySigLt50 = true;
	}
	if(peleid)
	  selSV.fillBranch("LeptonicSV_electronIndex", int(svit));
	if(pmuonid)
	  selSV.fillBranch("LeptonicSV_muonIndex", int(svit));

	selSV.fillBranch("LeptonicSV_passMuonID", pmuonid);
	selSV.fillBranch("LeptonicSV_passElectronID", peleid);
	if(hasGenInfoFlag)
	  selSV.fillBranch("LeptonicSV_isGold", (*Vertex_isGold)[svit]);
      }
      else if(hadronicId) {
	svType = "Hadronic";
	nHsv++;
	if(mass >= 15)                  anyHadMassGeq15 = true;
	if(mass >= 15 && dxySig >= 800) anyHadMassGeq15DxySigGeq800 = true;
	if(dxySig >= 800)               anyHadDxySigGeq800 = true;
	if(mass >= 15 && dxySig < 100)  anyHadMassGeq15DxySigLt100 = true;
	if(mass >= 15 && dxySig >= 100 && dxySig < 800) anyHadMassGeq15DxySigGeq100Lt800 = true;
	if(dxySig < 200)                anyHadDxySigLt200 = true;

	selSV.fillBranch("HadronicSV_massOverNtracks", mass/ntracks);
	if(hasGenInfoFlag)
	  selSV.fillBranch("HadronicSV_matchRatio", (*Vertex_matchRatio)[svit]);
      }
      else continue;

      const float p((*Vertex_p)[svit]), decayAngle((*Vertex_decayAngle)[svit]);

      selSV.fillBranch(svType+"SV_nTracks",              int(ntracks));
      selSV.fillBranch(svType+"SV_mass",                 mass);
      selSV.fillBranch(svType+"SV_pOverE",               float(p/sqrt(p*p + mass*mass)));
      selSV.fillBranch(svType+"SV_decayAngle",           decayAngle);
      selSV.fillBranch(svType+"SV_cosTheta",             cosTheta);
      selSV.fillBranch(svType+"SV_x",                    x);
      selSV.fillBranch(svType+"SV_y",                    y);
      selSV.fillBranch(svType+"SV_z",                    z);
      selSV.fillBranch(svType+"SV_dxy",                  dxy);
      selSV.fillBranch(svType+"SV_dxySig",               dxy/dxyError);
      selSV.fillBranch(svType+"SV_passMaterialVeto",     vetoTool.PassVeto(x, y));
      selSV.fillBranch(svType+"SV_tightZWindowSelection",tightOnZSelection);
    }//<<>>for( svit = 0; svit < nSVs; svit++ )
  }//<<>>if( doEVSVs )

  selSV.fillBranch("SV_nLeptonic", nLsv);
  selSV.fillBranch("SV_nHadronic", nHsv);

  //Analysis channel flags — LittleGuy compressed regions
  selSV.fillBranch("passNSVGe1Selection",
      bool((nHsv >= 1 && anyHadMassGeq15) || (nLsv >= 1 && anyLepMassGeq15)));
  selSV.fillBranch("passNSVGe1SelectionLowDxySigCR",
      bool(((nHsv >= 1 && anyHadMassGeq15) || (nLsv >= 1 && anyLepMassGeq10))
           && !anyHadMassGeq15DxySigGeq800 && !anyLepMassGeq10DxySigGeq500));
  selSV.fillBranch("passNSVGe1SelectionHighDxySigSR",
      bool(anyHadMassGeq15DxySigGeq800 || anyLepMassGeq10DxySigGeq500));
  selSV.fillBranch("passNHadEq1SelectionLowDxySigCR",
      bool(nHsv == 1 && anyHadMassGeq15 && !anyHadMassGeq15DxySigGeq800));
  selSV.fillBranch("passNSVEq0Selection",
      bool(nHsv == 0 && nLsv == 0));
  //Analysis channel flags — BigGuy non-compressed regions
  selSV.fillBranch("passNLepGe1SelectionHighDxySigSR",
      bool(nLsv >= 1 && anyLepDxySigGeq500));
  selSV.fillBranch("passNLepGe1SelectionLowDxySigCR",
      bool(nLsv >= 1 && !anyLepDxySigGeq500));
  selSV.fillBranch("passNHadGe1SelectionHighDxySigSR",
      bool(nLsv == 0 && nHsv >= 1 && anyHadDxySigGeq800));
  selSV.fillBranch("passNHadGe1SelectionLowDxySigCR",
      bool(nLsv == 0 && nHsv >= 1 && !anyHadDxySigGeq800));
  //Analysis channel flags — LittleGuy compressed validation regions
  selSV.fillBranch("passNSVGe1SelectionLowDxySigValCR",
      bool((anyHadMassGeq15DxySigLt100 || anyLepMassGeq10DxySigLt100)
           && !anyHadMassGeq15DxySigGeq100Lt800 && !anyLepMassGeq10DxySigGeq100Lt500
           && !anyHadMassGeq15DxySigGeq800 && !anyLepMassGeq10DxySigGeq500));
  selSV.fillBranch("passNSVGe1SelectionMidDxySigValSR",
      bool((anyHadMassGeq15DxySigGeq100Lt800 || anyLepMassGeq10DxySigGeq100Lt500)
           && !anyHadMassGeq15DxySigGeq800 && !anyLepMassGeq10DxySigGeq500));
  //Analysis channel flags — BigGuy non-compressed validation regions
  selSV.fillBranch("passNLepGe1SelectionLowDxySigValCR",
      bool(nLsv >= 1 && anyLepDxySigLt50));
  selSV.fillBranch("passNHadGe1SelectionLowDxySigValCR",
      bool(nLsv == 0 && nHsv >= 1 && anyHadDxySigLt200));

  geVars.set("nSVLep", nLsv);
  geVars.set("nSVHad", nHsv);

}//<<>>void KUCMSAodSkimmer::processSV()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setSVBranches( TTree* fOutTree ){

  std::cout << " - Making Branches for SV." << std::endl;
  for(const string svType : {"Hadronic", "Leptonic"}) {
    selSV.makeBranch(svType+"SV_nTracks",               VINT);
    selSV.makeBranch(svType+"SV_pOverE",                VFLOAT);
    selSV.makeBranch(svType+"SV_decayAngle",            VFLOAT);
    selSV.makeBranch(svType+"SV_cosTheta",              VFLOAT);
    selSV.makeBranch(svType+"SV_mass",                  VFLOAT);
    selSV.makeBranch(svType+"SV_x",                     VFLOAT);
    selSV.makeBranch(svType+"SV_y",                     VFLOAT);
    selSV.makeBranch(svType+"SV_z",                     VFLOAT);
    selSV.makeBranch(svType+"SV_dxy",                   VFLOAT);
    selSV.makeBranch(svType+"SV_dxySig",                VFLOAT);
    selSV.makeBranch(svType+"SV_passMaterialVeto",      VBOOL);
    selSV.makeBranch(svType+"SV_tightZWindowSelection", VBOOL);
  }

  selSV.makeBranch("LeptonicSV_electronIndex",  VINT);
  selSV.makeBranch("LeptonicSV_muonIndex",      VINT);
  selSV.makeBranch("LeptonicSV_passMuonID",     VBOOL);
  selSV.makeBranch("LeptonicSV_passElectronID", VBOOL);

  selSV.makeBranch("HadronicSV_massOverNtracks", VFLOAT);

  selSV.makeBranch("SV_nLeptonic", INT);
  selSV.makeBranch("SV_nHadronic", INT);

  //Analysis channel flags — LittleGuy compressed regions
  selSV.makeBranch("passNSVGe1Selection",              BOOL);
  selSV.makeBranch("passNSVGe1SelectionLowDxySigCR",   BOOL);
  selSV.makeBranch("passNSVGe1SelectionHighDxySigSR",  BOOL);
  selSV.makeBranch("passNHadEq1SelectionLowDxySigCR",  BOOL);
  selSV.makeBranch("passNSVEq0Selection",              BOOL);
  //Analysis channel flags — BigGuy non-compressed regions
  selSV.makeBranch("passNLepGe1SelectionHighDxySigSR", BOOL);
  selSV.makeBranch("passNLepGe1SelectionLowDxySigCR",  BOOL);
  selSV.makeBranch("passNHadGe1SelectionHighDxySigSR", BOOL);
  selSV.makeBranch("passNHadGe1SelectionLowDxySigCR",  BOOL);
  //Analysis channel flags — validation regions
  selSV.makeBranch("passNSVGe1SelectionLowDxySigValCR",  BOOL);
  selSV.makeBranch("passNSVGe1SelectionMidDxySigValSR",  BOOL);
  selSV.makeBranch("passNLepGe1SelectionLowDxySigValCR", BOOL);
  selSV.makeBranch("passNHadGe1SelectionLowDxySigValCR", BOOL);

  if( hasGenInfoFlag ){
    selSV.makeBranch("HadronicSV_matchRatio", VFLOAT);
    selSV.makeBranch("LeptonicSV_isGold",     VBOOL);
  }
  selSV.attachBranches(fOutTree);

}//<<>>void KUCMSAodSkimmer::setSVBranches( TTree* fOutTree )
