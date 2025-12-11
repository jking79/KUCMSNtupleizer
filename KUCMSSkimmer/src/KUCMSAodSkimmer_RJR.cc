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

//#define RJRDEBUG true
#define RJRDEBUG false

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

//------------------------------------------------------------------------------------------------------------
// process RJR for event
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::processRJR( int type, bool newEvent ){

  //bool verbose = true;
  bool verbose = false;

  if( RJRDEBUG || verbose ) std::cout << "Processing RJR event varibles for type " << type << std::endl;
  // ------------------  Loading &  Processing RJR event varibles  ---------------------------------------------------------

  // intilize out branches
  if( newEvent ) selRjrVars.clearBranches(); // <<<<<<<   must do

  // process event
  LAB->ClearEvent();

  auto sPhoPt = geVects( "selPhoPt");
  auto sPhoEta = geVects( "selPhoEta");
  auto sPhoPhi = geVects( "selPhoPhi");
  auto sPhoMx = geVects( "selPhoEMx");
  auto sPhoMy = geVects( "selPhoEMy");

  auto nSelPhotons = geCnts("nSelPhotons");

  auto phoRMetCPx = geVars("metPx"); //selMet.getFLBranchValue("metPx");
  auto phoRMetCPy = geVars("metPy"); //selMet.getFLBranchValue("metPy");
  float unCorMet = hypo(phoRMetCPx,phoRMetCPy);

  auto phoBMetPx = geVars("metBPx"); //selMet.getFLBranchValue("metPx");
  auto phoBMetPy = geVars("metBPy"); //selMet.getFLBranchValue("metPy");
  float unCorBMet = hypo(phoBMetPx,phoBMetPy);

  if( RJRDEBUG ) std::cout << " - Loading Photons." << std::endl;

  int nRJRPhos = 0;
  std::vector<RFKey> jetID;
  std::vector<TLorentzVector> pho4vec;
  std::vector<TLorentzVector> jet4vec;
  bool zsig = ( geVars("noSVorPho") == 1 ) ? true : false;	
  if( zsig ) std::cout << " ------------ noSVorPho is True !!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
  if( true ){ // old : if( not zsig ){ ?? never do this ??
	nRJRPhos = ( nSelPhotons > 1 ) ? 2 : ( nSelPhotons > 0 ) ? 1 : 0;
  	for( int spidx = 0; spidx < nSelPhotons; spidx++ ){
    
    	if( RJRDEBUG ) std::cout << " - Loading Pho " << spidx << std::endl;
        if( type == 0 ){
            phoBMetPx += sPhoMx[spidx]; 
            phoBMetPy += sPhoMy[spidx];
        }//<<>>if( type == 0 )
        TLorentzVector phojet;
        phojet.SetPtEtaPhiM( sPhoPt[spidx], sPhoEta[spidx], sPhoPhi[spidx], 0 );
        jetID.push_back( COMB_J->AddLabFrameFourVector(phojet) );
        if( spidx < nRJRPhos ) pho4vec.push_back(phojet);
		else jet4vec.push_back(phojet);		

  	}//<<>>for( spidx = 0; spidx < nSelPhotons; spidx++ )
/*
	if( nSelPhotons != 0 ){
    	nRJRPhos = 1;
      	if( RJRDEBUG ) std::cout << " - Loading Lead/SubLead Pho" << std::endl;
      	if( type == 0 ){
			phoBMetPx += sPhoMx[0]; 
			phoBMetPy += sPhoMy[0];
      	}//<<>>if( type == 0 )
		TLorentzVector phojet;
		phojet.SetPtEtaPhiM( sPhoPt[0], sPhoEta[0], sPhoPhi[0], 0 );
		jetID.push_back( COMB_J->AddLabFrameFourVector(phojet) );
		pho4vec.push_back(phojet);
    }//<<>>if( nSelPhotons > 0 )
    if( nSelPhotons > 1 ){ 
    	nRJRPhos = 2;
      	if( type == 0 ){
			phoBMetPx += sPhoMx[1];
			phoBMetPy += sPhoMy[1];
      	}//<<>>if( type == 0 )
		TLorentzVector sphojet;
		sphojet.SetPtEtaPhiM( sPhoPt[1], sPhoEta[1], sPhoPhi[1], 0 );
		jetID.push_back( COMB_J->AddLabFrameFourVector(sphojet) );
		pho4vec.push_back(sphojet);
	}//<<>>if( nSelPhotons > 1 )
*/

  }//<<>>if( not zsig )
  else std::cout << " --- !!!  Photons not loaded into RJR !!!!!!!!!!!!!!!!!!" << std::endl;

  if( RJRDEBUG ) std::cout << " - Loading MET." << std::endl;
  if( verbose ){
    std::cout << " - Loading MET lPt: " << sPhoPt[0] << " lPhi: " << sPhoPhi[0] << std::endl;
    std::cout << " - Loading MET slPt: " << sPhoPt[1] << " slPhi: " << sPhoPhi[1] << std::endl; 
    std::cout << " - Loading MET x: " << geVars("metPx") << " -> " << phoRMetCPx;
    std::cout << " y: " << geVars("metPy") << " -> " << phoRMetCPy << std::endl;
  }//<<>>if( verbose )

  TVector3 ETBMiss(phoBMetPx,phoBMetPy,0);
  TVector3 ETMiss(phoRMetCPx,phoRMetCPy,0);
  if( type == 0 ) INV->SetLabFrameThreeVector(ETBMiss);
  else if( type == 1 ) INV->SetLabFrameThreeVector(ETMiss);
  else std::cout << " !!!!!!!! Valid RJR Photon Processing Option Not Specified !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

  auto nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
  auto selJetPt = geVects( "selJetPt");
  auto selJetEta = geVects( "selJetEta");
  auto selJetPhi = geVects( "selJetPhi");
  auto selJetMass = geVects( "selJetMass");
  std::vector<RFKey> leadJetKey;
  if( RJRDEBUG ) std::cout << " - Loading Jets." << std::endl;
  for( uInt it = 0; it < nSelJets; it++ ){
    auto sjetPt = selJetPt[it]; //selJets.getFLBranchValue( "selJetPt", it );
    auto sjetEta = selJetEta[it]; //selJets.getFLBranchValue( "selJetEta", it );
    auto sjetPhi = selJetPhi[it]; //selJets.getFLBranchValue( "selJetPhi", it ); 
    auto sjetMass = ( selJetMass[it] > 0 ) ? selJetMass[it] : 0; //selJets.getFLBranchValue( "selJetMass", it );   
    TLorentzVector jet;
    jet.SetPtEtaPhiM( sjetPt, sjetEta, sjetPhi, sjetMass );
    if( verbose ) std::cout << " - Loading Jet Pt: " << sjetPt << " Eta: " << sjetEta;
    if( verbose ) std::cout << " Phi: " << sjetPhi << " M: " << sjetMass << std::endl;
    if( it == 0 ){ leadJetKey.push_back( COMB_J->AddLabFrameFourVector(jet) ); jetID.push_back(leadJetKey[0]); } 
    else jetID.push_back(COMB_J->AddLabFrameFourVector(jet));
    jet4vec.push_back(jet); 
  }//<<>>for( int i = 0; i < nSelJets; i++ )

  if( !LAB->AnalyzeEvent() ) std::cout << "Something went wrong with tree event analysis" << std::endl;
	
  if( RJRDEBUG ) std::cout << " - Getting RJR Object Location Information." << std::endl;
  // ---------  Finished Processing RJR varibles --------------------------------------------------------------------

  uInt firstJet = 0;
  if( ( type < 2 ) && ( nRJRPhos > 0 ) ) firstJet = ( nRJRPhos == 2 ) ? 2 : 1 ;
  int nJetsJa = 0;
  int nJetsJb = 0;
  bool isALeadPhoSide = true;	
  int subPhoLocation = 0;
  int leadPhoLocation = ( type < 2 ) ? 1 : 0;
  int xtraPhos = ( nSelPhotons > nRJRPhos ) ? nSelPhotons - nRJRPhos : 0;
  int nVisObjects = jetID.size();
  if( nVisObjects != ( (type < 2 ) ? nRJRPhos + nSelJets + xtraPhos : nSelJets + xtraPhos ) ){ 
    std::cout << " !!!!!!    nVisObjects != ( nRJRPhos + nSelJets )  !!!!!!!!! " << std::endl;
    std::cout << " !!!!!!    " << nVisObjects << " != ( " << nRJRPhos + xtraPhos << " + " << nSelJets << " )  !!!!!!!!! " << std::endl;
  }//<<>>if( nVisObjects != ( nRJRPhos + nSelJets ) )
  if( COMB_J->GetFrame(jetID[0]) == *J1b || COMB_J->GetFrame(jetID[0]) == *J2b ){ isALeadPhoSide = false; if( type < 2 ) leadPhoLocation = 2; }
  //  -- redifine A & B side for jets to be : is jet on lead pho side -> A ; is not on lead pho side -> B // done
  for( uInt it = firstJet; it < nVisObjects; it++ ){	

    if( COMB_J->GetFrame(jetID[it]) == *J1a || COMB_J->GetFrame(jetID[it]) == *J2a  ){ isALeadPhoSide ? nJetsJa++ : nJetsJb++; } // one for each frame 
    else { isALeadPhoSide ? nJetsJb++ : nJetsJa++; }

  }//<<>>for( int i = 0; i < nSelJets; i++ )
  float abDiffSide = isALeadPhoSide ? 1 : -1;
  if( ( type < 2 ) && ( nRJRPhos == 2 ) ) subPhoLocation = ( COMB_J->GetFrame(jetID[1]) == *J1a || COMB_J->GetFrame(jetID[1]) == *J2a ) ? 1 : 2;

  selRjrVars.fillBranch( "rjrType", type );
  selRjrVars.fillBranch( "rjrABSide", isALeadPhoSide );
  selRjrVars.fillBranch( "rjrNJetsJa", nJetsJa );
  selRjrVars.fillBranch( "rjrNJetsJb", nJetsJb );
  selRjrVars.fillBranch( "rjrNJets", int(geCnts("nSelJets")) );
  selRjrVars.fillBranch( "rjrNPhotons", nRJRPhos );
  selRjrVars.fillBranch( "rjrNVisObjects", nVisObjects );
  selRjrVars.fillBranch( "rjrLeadPhoLocation", leadPhoLocation );
  selRjrVars.fillBranch( "rjrSubPhoLocation", subPhoLocation );

  std::vector<int> phoside;
  std::vector<int> jetside;
  for( uInt it = 0; it < nVisObjects; it++ ){

    bool onAside = true;
    int whichside = -1;
    ////if( COMB_J->GetFrame(jetID[it]) == *J1b || COMB_J->GetFrame(jetID[it]) == *J2b  ) onAside = false;
    if( COMB_J->GetFrame(jetID[it]) == *J1a ) whichside = 0;
    if( COMB_J->GetFrame(jetID[it]) == *J2a ) whichside = 1;
    if( COMB_J->GetFrame(jetID[it]) == *J1b ) whichside = 2;
    if( COMB_J->GetFrame(jetID[it]) == *J2b ) whichside = 3;
    if( whichside < 0 ) std::cout << " !!!!!!!  Visible object is not on a side !!!!!!! " << std::endl;
    if( whichside > 1 ) onAside = false;
    if( it < firstJet ){ 
      selRjrVars.fillBranch( "rjrVisPhoSide", onAside );
      phoside.push_back(whichside);
    } else { 
      selRjrVars.fillBranch( "rjrVisJetSide", onAside );
      jetside.push_back(whichside);
    }//<<>>if( it < nSelPhotons )

  }//<<>>for( uInt it = 0; it < nVisObjects; it++ )

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles." << std::endl;

  float m_MS = S->GetMass();
  float m_cosS  = S->GetCosDecayAngle();
  float m_dphiS = S->GetDeltaPhiDecayAngle();
  float m_dphiSI  = S->GetDeltaPhiBoostVisible();
  
  float m_MX2a = X2a->GetMass();
  float m_cosX2a = X2a->GetCosDecayAngle();
  float m_MX2b = X2b->GetMass();
  float m_cosX2b = X2b->GetCosDecayAngle();

  float m_PTS = S->GetFourVector().Pt();
  float m_PZS = S->GetFourVector().Pz();
  float m_PS = S->GetMomentum(*LAB);

  selRjrVars.fillBranch( "rjrX2aMass", m_MX2a );
  selRjrVars.fillBranch( "rjrX2aCosA", m_cosX2a );
  selRjrVars.fillBranch( "rjrX2bMass", m_MX2b );
  selRjrVars.fillBranch( "rjrX2bCosA", m_cosX2b );

  selRjrVars.fillBranch( "rjrSMass", m_MS );
  selRjrVars.fillBranch( "rjrSCosA", m_cosS );

  selRjrVars.fillBranch( "rjrSdphiDA", m_dphiS  );
  selRjrVars.fillBranch( "rjrSdphiBV", m_dphiSI ); // !!!!!!!!!!!!!  dphiCMI
  selRjrVars.fillBranch( "rjrDPhiSI", m_dphiSI ); // !!!!!!!!!!!!!  dphiCMI
  selRjrVars.fillBranch( "rjrPTS", m_PTS ); // !!!!!!!!!!!!!!!!!!  P_{t}^{S}
  selRjrVars.fillBranch( "rjrPZS", m_PZS );
  selRjrVars.fillBranch( "rjrPS", m_PS );

 
  if( type == 0 ) hist2d[2]->Fill(m_dphiSI,m_PTS);
  if( type == 1 ) hist2d[3]->Fill(m_dphiSI,m_PTS); 

  bool rjrCleaningVeto0 =  ( m_PTS > 150.0 ) ? true : false; 
  selRjrVars.fillBranch( "rjrCleaningVeto0", rjrCleaningVeto0 );

  bool rjrCleaningVeto1 = false;
  if( ( m_dphiSI <= 0.6 ) && ( m_PTS > (250.0 * m_dphiSI) ) ) rjrCleaningVeto1 = true;
  if( ( m_dphiSI > 0.6 ) && ( m_dphiSI < 2.4 ) && ( m_PTS > 150.0 ) ) rjrCleaningVeto1 = true; 
  if( ( m_dphiSI >= 2.4 ) && ( m_PTS > (150.0 - (187.5 * (m_dphiSI - 2.4))) ) ) rjrCleaningVeto1 = true;
  selRjrVars.fillBranch( "rjrCleaningVeto1", rjrCleaningVeto1 );  

  bool rjrCleaningVeto2 = false; 
  if( ( m_dphiSI <= 0.3 ) && ( m_PTS > (500.0 * m_dphiSI) ) ) rjrCleaningVeto2 = true; 
  if( ( m_dphiSI > 0.3 ) && ( m_dphiSI < 2.8 ) && ( m_PTS > 150.0 ) ) rjrCleaningVeto2 = true;  
  if( ( m_dphiSI >= 2.8 ) && ( m_PTS > (150.0 - (375.0 * (m_dphiSI - 2.8))) ) ) rjrCleaningVeto2 = true; 
  selRjrVars.fillBranch( "rjrCleaningVeto2", rjrCleaningVeto2 );   

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 2." << std::endl;

  std::vector< TLorentzVector > hp4;
  std::vector< TLorentzVector > hxp4;
	
	int nPhov = pho4vec.size();
	if( nPhov != phoside.size() ) std::cout << " !!!!!!!  pho4vec - phoside count mismatch !!!!!!!" << std::endl; 
	int nJetv = jet4vec.size();
    if( nJetv != jetside.size() ) std::cout << " !!!!!!!  jet4vec - jetside count mismatch !!!!!!!" << std::endl;
	//std::cout << " nPhov " << nPhov << " phoside " << phoside.size() << " nJetv " << nJetv << " jetside " << jetside.size() << std::endl; 

	TLorentzVector nulltv(0,0,0,0);
	std::vector< TLorentzVector > pJetSideSum4Vec{nulltv,nulltv,nulltv,nulltv};
	std::vector< TLorentzVector > phoSideSum4Vec{nulltv,nulltv,nulltv,nulltv};
	for ( int pit = 0; pit < nPhov; pit++ ){ phoSideSum4Vec[phoside[pit]] += pho4vec[pit]; pJetSideSum4Vec[phoside[pit]] += pho4vec[pit];}
    std::vector< TLorentzVector > jetSideSum4Vec{nulltv,nulltv,nulltv,nulltv};
    for ( int jit = 0; jit < nJetv; jit++ ){ jetSideSum4Vec[jetside[jit]] += jet4vec[jit]; pJetSideSum4Vec[jetside[jit]] += jet4vec[jit];}

    hp4.push_back(S->GetFourVector(pJetSideSum4Vec[0]));//0 J1a
    hp4.push_back(S->GetFourVector(pJetSideSum4Vec[1]));//1 J2a
    hp4.push_back(S->GetFourVector(pJetSideSum4Vec[2]));//2 J1b
    hp4.push_back(S->GetFourVector(pJetSideSum4Vec[3]));//3 J2b

    selRjrVars.fillBranch( "rjr_sJ1a_px", float(hp4[0].Px()) );
    selRjrVars.fillBranch( "rjr_sJ1a_py", float(hp4[0].Py()) );
    selRjrVars.fillBranch( "rjr_sJ1a_pz", float(hp4[0].Pz()) );
    selRjrVars.fillBranch( "rjr_sJ1a_e", float(hp4[0].E()) );

    selRjrVars.fillBranch( "rjr_sJ2a_px", float(hp4[1].Px()) );
    selRjrVars.fillBranch( "rjr_sJ2a_py", float(hp4[1].Py()) );
    selRjrVars.fillBranch( "rjr_sJ2a_pz", float(hp4[1].Pz()) );
    selRjrVars.fillBranch( "rjr_sJ2a_e", float(hp4[1].E()) );

    selRjrVars.fillBranch( "rjr_sJ1b_px", float(hp4[2].Px()) );
    selRjrVars.fillBranch( "rjr_sJ1b_py", float(hp4[2].Py()) );
    selRjrVars.fillBranch( "rjr_sJ1b_pz", float(hp4[2].Pz()) );
    selRjrVars.fillBranch( "rjr_sJ1b_e", float(hp4[2].E()) );

    selRjrVars.fillBranch( "rjr_sJ2b_px", float(hp4[3].Px()) );
    selRjrVars.fillBranch( "rjr_sJ2b_py", float(hp4[3].Py()) );
    selRjrVars.fillBranch( "rjr_sJ2b_pz", float(hp4[3].Pz()) );
    selRjrVars.fillBranch( "rjr_sJ2b_e", float(hp4[3].E()) );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 3." << std::endl;

	float x1sp = ( X1a->GetFourVector(*S) + X1b->GetFourVector(*S) ).P();
	float x1asp = ( X1a->GetFourVector(*S) ).P();
    float x1bsp = ( X1b->GetFourVector(*S) ).P();

    hp4.push_back( hp4[0] + hp4[1] );//4 Ja
    hp4.push_back( hp4[2] + hp4[3] );//5 Jb

    float phs20 = hp4[4].P() + hp4[5].P();
    float phs40 = hp4[0].P() + hp4[1].P() + hp4[2].P() + hp4[3].P();
    float phs21 = phs20 + x1sp;
    float phs41 = phs40 + x1sp;

    float phs21a = hp4[0].P() + hp4[1].P() + x1asp;
    float phs21b = hp4[2].P() + hp4[3].P() + x1bsp;
    float phs11a = hp4[4].P() + x1asp;
    float phs11b = hp4[5].P() + x1bsp;

    selRjrVars.fillBranch( "rjr_pHs20", phs20 );
    selRjrVars.fillBranch( "rjr_pHs40", phs40 );
    selRjrVars.fillBranch( "rjr_pHs21", phs21 );
    selRjrVars.fillBranch( "rjr_pHs41", phs41 );
    selRjrVars.fillBranch( "rjr_pHs21a", phs21a );
    selRjrVars.fillBranch( "rjr_pHs21b", phs21b );
    selRjrVars.fillBranch( "rjr_pHs11a", phs11a );
    selRjrVars.fillBranch( "rjr_pHs11b", phs11b );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 4." << std::endl;

    hxp4.push_back(X2a->GetFourVector(pJetSideSum4Vec[0]));//0
    hxp4.push_back(X2a->GetFourVector(pJetSideSum4Vec[1]));//1
    hxp4.push_back(X1a->GetFourVector(*X2a));//2
    hxp4.push_back(X2b->GetFourVector(pJetSideSum4Vec[2]));//3
    hxp4.push_back(X2b->GetFourVector(pJetSideSum4Vec[3]));//4
    hxp4.push_back(X1b->GetFourVector(*X2b));//5

    hxp4.push_back( hxp4[0] + hxp4[1] );//6
    hxp4.push_back( hxp4[3] + hxp4[4] );//7

    float phxa20 = hxp4[0].P() + hxp4[1].P();
    float phxb20 = hxp4[3].P() + hxp4[4].P();
    float phxa21 = phxa20 + hxp4[2].P();
    float phxb21 = phxb20 + hxp4[5].P();
    float phxa10 = hxp4[6].P();
    float phxb10 = hxp4[7].P();
    float phxa11 = hxp4[6].P() + hxp4[2].P();
    float phxb11 = hxp4[7].P() + hxp4[5].P();

    selRjrVars.fillBranch( "rjr_pHxa20", phxa20 );
    selRjrVars.fillBranch( "rjr_pHxb20", phxb20 );
    selRjrVars.fillBranch( "rjr_pHxa21", phxa21 );
    selRjrVars.fillBranch( "rjr_pHxb21", phxb21 );
    selRjrVars.fillBranch( "rjr_pHxa10", phxa10 );
    selRjrVars.fillBranch( "rjr_pHxb10", phxb10 );
    selRjrVars.fillBranch( "rjr_pHxa11", phxa11 );
    selRjrVars.fillBranch( "rjr_pHxb11", phxb11 );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 5." << std::endl;

    float pj1apt = S->GetTransverseMomentum(pJetSideSum4Vec[0]);
    float pj2apt = S->GetTransverseMomentum(pJetSideSum4Vec[1]);
    float pj1bpt = S->GetTransverseMomentum(pJetSideSum4Vec[2]);
    float pj2bpt = S->GetTransverseMomentum(pJetSideSum4Vec[3]);
    float pjapt = S->GetTransverseMomentum(pJetSideSum4Vec[0]+pJetSideSum4Vec[1]);
    float pjbpt = S->GetTransverseMomentum(pJetSideSum4Vec[2]+pJetSideSum4Vec[3]);

	float x1spt = S->GetTransverseMomentum( X1a->GetFourVector(*LAB) + X1b->GetFourVector(*LAB) );
    float x1aspt = S->GetTransverseMomentum( X1a->GetFourVector(*LAB) );
    float x1bspt = S->GetTransverseMomentum( X1b->GetFourVector(*LAB) );

    float phts20 = pjapt + pjbpt;
    float phts40 = pj1apt + pj2apt + pj1bpt + pj2bpt;
    float phts21 = phts20 + x1spt;
    float phts41 = phts40 + x1spt;

    float phts21a = pj1apt + pj2apt + x1aspt;
    float phts21b = pj1bpt + pj2bpt + x1bspt;
    float phts11a = pjapt + x1aspt;
    float phts11b = pjbpt + x1bspt;

    selRjrVars.fillBranch( "rjr_pHts20", pjapt + pjbpt );
    selRjrVars.fillBranch( "rjr_pHts40", pj1apt + pj2apt + pj1bpt + pj2bpt );
    selRjrVars.fillBranch( "rjr_pHts21", phts21 );
    selRjrVars.fillBranch( "rjr_pHts41", phts41 );
    selRjrVars.fillBranch( "rjr_pHts21a", phts21a );
    selRjrVars.fillBranch( "rjr_pHts21b", phts21b );
    selRjrVars.fillBranch( "rjr_pHts11a", phts11a );
    selRjrVars.fillBranch( "rjr_pHts11b", phts11b );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 6." << std::endl;

    float ptj1apt = X2a->GetTransverseMomentum(pJetSideSum4Vec[0]);
    float ptj2apt = X2a->GetTransverseMomentum(pJetSideSum4Vec[1]);
    float ptx1apt = X2a->GetTransverseMomentum(X1a->GetFourVector(*LAB));
    float ptj1bpt = X2b->GetTransverseMomentum(pJetSideSum4Vec[2]);
    float ptj2bpt = X2b->GetTransverseMomentum(pJetSideSum4Vec[3]);
    float ptx1bpt = X2b->GetTransverseMomentum(X1b->GetFourVector(*LAB));
    float ptj12apt = X2a->GetTransverseMomentum(pJetSideSum4Vec[0]+pJetSideSum4Vec[1]);
    float ptj12bpt = X2b->GetTransverseMomentum(pJetSideSum4Vec[2]+pJetSideSum4Vec[3]);

    float phtxa20 = ptj1apt + ptj2apt;
    float phtxb20 = ptj1bpt + ptj2bpt;
    float phtxa21 = phtxa20 + ptx1apt;
    float phtxb21 = phtxb20 + ptx1bpt;
    float phtxa10 = ptj12apt;
    float phtxb10 = ptj12bpt;
    float phtxa11 = ptj12apt + ptx1apt;
    float phtxb11 = ptj12bpt + ptx1bpt;

    selRjrVars.fillBranch( "rjr_pHtxa20", phtxa20 );
    selRjrVars.fillBranch( "rjr_pHtxb20", phtxb20 );
    selRjrVars.fillBranch( "rjr_pHtxa21", phtxa21 );
    selRjrVars.fillBranch( "rjr_pHtxb21", phtxb21 );
    selRjrVars.fillBranch( "rjr_pHtxa10", phtxa10 );
    selRjrVars.fillBranch( "rjr_pHtxb10", phtxb10 );
    selRjrVars.fillBranch( "rjr_pHtxa11", phtxa11 );
    selRjrVars.fillBranch( "rjr_pHtxb11", phtxb11 );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 7." << std::endl;

	float phopt1 = 0;
    float phopt2 = 0;

	if( nPhov > 0 ) phopt1 = S->GetTransverseMomentum( pho4vec[0] );
    if( nPhov > 1 ) phopt2 = S->GetTransverseMomentum( pho4vec[1] );
    //if( nPhov > 0 ) phopt1 = S->GetFourVector( pho4vec[0] ).P();
    //if( nPhov > 1 ) phopt2 = S->GetFourVector( pho4vec[1] ).P();

    selRjrVars.fillBranch( "rjr_p1Pts", phopt1 );
    selRjrVars.fillBranch( "rjr_p2Pts", phopt2 );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 8." << std::endl;

    float phopt1xa = 0;
    float phopt1xb = 0;
    float phopt2xa = 0;
    float phopt2xb = 0;

    if( nPhov > 0 ){ 
		if( phoside[0] < 2 ) phopt1xa = X2a->GetTransverseMomentum( pho4vec[0] );
		else phopt1xb = X2b->GetTransverseMomentum( pho4vec[0] );
	}//<<>>if( nPhov > 0 )
    if( nPhov > 1 ){ 
		if( phoside[1] < 2 ) phopt2xa = X2a->GetTransverseMomentum( pho4vec[1] );
        else phopt2xb = X2b->GetTransverseMomentum( pho4vec[1] );
	}//<<>>if( nPhov > 1 )

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 9." << std::endl;

    float phop1xa = 0;
    float phop1xb = 0;
    float phop2xa = 0;
    float phop2xb = 0;

    if( nPhov > 0 ){
        if( phoside[0] < 2 ) phop1xa = X2a->GetFourVector( pho4vec[0] ).P();
        else phop1xb = X2b->GetFourVector( pho4vec[0] ).P();
    }//<<>>if( nPhov > 0 )
    if( nPhov > 1 ){
        if( phoside[1] < 2 ) phop2xa = X2a->GetFourVector( pho4vec[1] ).P();
        else phop2xb = X2b->GetFourVector( pho4vec[1] ).P();
    }//<<>>if( nPhov > 1 )


    selRjrVars.fillBranch( "rjr_p1Ptxa", phopt1xa );
    selRjrVars.fillBranch( "rjr_p1Ptxb", phopt1xb );
    selRjrVars.fillBranch( "rjr_p2Ptxa", phopt2xa );
    selRjrVars.fillBranch( "rjr_p2Ptxb", phopt2xb );

    selRjrVars.fillBranch( "rjr_p1Pxa", phop1xa );
    selRjrVars.fillBranch( "rjr_p1Pxb", phop1xb );
    selRjrVars.fillBranch( "rjr_p2Pxa", phop2xa );
    selRjrVars.fillBranch( "rjr_p2Pxb", phop2xb );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 10." << std::endl;

	std::vector< TLorentzVector > p4;
	p4.push_back(Ja->GetFourVector(*S));
    p4.push_back(Jb->GetFourVector(*S));
    p4.push_back(X1a->GetFourVector(*S));
    p4.push_back(X1b->GetFourVector(*S));

	float pHs22 = p4[0].P() + p4[1].P() + p4[2].P() + p4[3].P();
    float pHs11 =  ( p4[0] + p4[1] ).P() + ( p4[2] + p4[3] ).P();

    float pHJas = Ja->GetTransverseMomentum(*S);
    float pHJbs = Jb->GetTransverseMomentum(*S);
    float pHX1as = X1a->GetTransverseMomentum(*S);
    float pHX1bs = X1b->GetTransverseMomentum(*S);
	float pHJs = S->GetTransverseMomentum( Ja->GetFourVector(*LAB) + Jb->GetFourVector(*LAB) );
    float pHX1s = S->GetTransverseMomentum( X1a->GetFourVector(*LAB) + X1b->GetFourVector(*LAB) );

	float pHts22 = pHJas + pHJbs + pHX1as + pHX1bs;
	float pHts11 = pHJs + pHX1s;

	float Rs = (phxa11+phxb11)/pHs22;
	float Rxa = phxa11/phxa21;
    float Rxb = phxb11/phxb21;
	float Rx = hypo(phxa11/phxa21,phxb11/phxb21);

    selRjrVars.fillBranch( "rjr_pHs22", pHs22 );
    selRjrVars.fillBranch( "rjr_pHs11", pHs11 );
    selRjrVars.fillBranch( "rjr_pHts22", pHts22 );
    selRjrVars.fillBranch( "rjr_pHts11", pHts11 );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 11." << std::endl;

	for( int i = 0; i < 4; i++ ){

		float px = p4[i].Px();
		float py = p4[i].Py(); 
		float pz = p4[i].Pz(); 
		float e = sqrt( px*px + py*py + pz*pz );
		p4[i].SetPxPyPzE(px,py,pz,e);	

	}//<<>>for( int i = 0; i < 4; i++ )

	float a_MX2a = ( p4[0] + p4[2] ).M();
	float a_MX2b = ( p4[1] + p4[3] ).M();
	float a_MS = ( p4[0] + p4[2] + p4[1] + p4[3] ).M();// -- !! Mr	

    float AX2QSum = std::sqrt((sq2(a_MX2a)+sq2(a_MX2b))/2);
    float AX2NQSum = 2*AX2QSum/a_MS;//*2 -- !! R
    float AX2Diff = abDiffSide*(a_MX2a-a_MX2b)/(a_MX2a+a_MX2b);

    float pf_pX2a = X2a->GetMomentum(*S);
    float pf_pX2b = X2b->GetMomentum(*S);

    selRjrVars.fillBranch( "rjr_Mr", a_MS );
    selRjrVars.fillBranch( "rjr_R", AX2NQSum );
    selRjrVars.fillBranch( "rjrAX2Diff", AX2Diff );

    selRjrVars.fillBranch( "rjr_Ms", phs41 );
    selRjrVars.fillBranch( "rjr_Rs", Rs );
    selRjrVars.fillBranch( "rjr_Rx", Rx );
    selRjrVars.fillBranch( "rjr_Rxa", Rxa );
    selRjrVars.fillBranch( "rjr_Rxb", Rxb );

    selRjrVars.fillBranch( "rjrX2aPs", pf_pX2a );
    selRjrVars.fillBranch( "rjrX2bPs", pf_pX2b );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 12." << std::endl;

  	float m_EVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).E();
  	float m_EVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).E();
  	float m_PVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).P();
  	float m_PVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).P();

  	float m_MX1a = X1a->GetMass();
  	float m_cosX1a = X1a->GetCosDecayAngle();
  	float m_MX1b = X1b->GetMass();
  	float m_cosX1b = X1b->GetCosDecayAngle();  
    float x1a_pt = X1a->GetFourVector().Pt();
    float x1a_eta = X1a->GetFourVector().Eta();
    float x1a_phi = X1a->GetFourVector().Phi();
    float x1b_pt = X1b->GetFourVector().Pt();
    float x1b_eta = X1b->GetFourVector().Eta();
    float x1b_phi = X1b->GetFourVector().Phi();

    float m_MVa = Ja->GetMass();
    float m_MVb = Jb->GetMass();

    float m_MVSum = std::sqrt((sq2(m_MVa)+sq2(m_MVb))/2);
	float m_MVNSum = 2*m_MVSum/a_MS;// -- !! Rv
    float m_MVab = std::sqrt(sq2(m_MVa)+sq2(m_MVb))/phs40;

    selRjrVars.fillBranch( "rjr_Rv", m_MVNSum );
    selRjrVars.fillBranch( "rjr_Mva", m_MVa );
    selRjrVars.fillBranch( "rjr_Mvb", m_MVb );
    selRjrVars.fillBranch( "rjr_Rm", m_MVab );

  	float m_PV_lab = S->GetListVisibleFrames().GetFourVector().P();
  	float m_dphiMET_V = 0.f; 
	if( type == 0 ) m_dphiMET_V = S->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETBMiss);
    else if( type == 1 ) m_dphiMET_V = S->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

    selRjrVars.fillBranch( "rjrX1aMass", m_MX1a );
    selRjrVars.fillBranch( "rjrX1aCosA", m_cosX1a );
    selRjrVars.fillBranch( "rjrX1bMass", m_MX1b );
    selRjrVars.fillBranch( "rjrX1bCosA", m_cosX1b );
    selRjrVars.fillBranch( "rjrX1a_Pt", x1a_pt );
    selRjrVars.fillBranch( "rjrX1a_Eta", x1a_eta );
    selRjrVars.fillBranch( "rjrX1a_Phi", x1a_phi );
    selRjrVars.fillBranch( "rjrX1b_Pt", x1b_pt );
    selRjrVars.fillBranch( "rjrX1b_Eta", x1b_eta );
    selRjrVars.fillBranch( "rjrX1b_Phi", x1b_phi );

    selRjrVars.fillBranch( "rjrEVa", m_EVa );
    selRjrVars.fillBranch( "rjrEVb", m_EVb );
    selRjrVars.fillBranch( "rjrPVa", m_PVa );
    selRjrVars.fillBranch( "rjrPVb", m_PVb );

    selRjrVars.fillBranch( "rjrPVlab", m_PV_lab );
    selRjrVars.fillBranch( "rjrDphiMETV", m_dphiMET_V );

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles 13." << std::endl;

	float rjrMet = hypo(phoRMetCPx,phoRMetCPy);

    selRjrVars.fillBranch( "rjrMET", rjrMet ); 

	// attempt to use RJR to guess p vector of N2
    TLorentzVector x1alab = X1a->GetFourVector(*LAB);
    TLorentzVector x1blab = X1b->GetFourVector(*LAB);
	float n2apx = x1alab.Px();
    float n2apy = x1alab.Py();
    float n2apz = x1alab.Pz();
	float n2bpx = x1blab.Px();
    float n2bpy = x1blab.Py();
    float n2bpz = x1blab.Pz();
	float n2ap = hypo( n2apx, n2apy, n2apz );
    float n2bp = hypo( n2bpx, n2bpy, n2bpz );
	bool aorb = true;
	//bool aorb = n2ap > n2bp;
	//bool aorb = isALeadPhoSide;
	float n2px = ( aorb ) ? n2apx : n2bpx;
    float n2py = ( aorb ) ? n2apy : n2bpy;
    float n2pz = ( aorb ) ? n2apz : n2bpz;
    selRjrVars.fillBranch( "rjrN2Px", n2px );
    selRjrVars.fillBranch( "rjrN2Py", n2py );
    selRjrVars.fillBranch( "rjrN2Pz", n2pz );

}//<<>>void KUCMSAodSkimmer::processRJR( int type, bool newEvent )

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setRJRBranches( TTree* fOutTree ){

	std::cout << " - Making Branches for RJR." << std::endl;
    selRjrVars.makeBranch( "rjrType", VINT );

    selRjrVars.makeBranch( "rjrX1aMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX1aCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX1bMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX1bCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Pt", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Eta", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Phi", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Pt", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Eta", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Phi", VFLOAT );

    selRjrVars.makeBranch( "rjrX2aMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX2aCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX2aPs", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bPs", VFLOAT );

    selRjrVars.makeBranch( "rjrSMass", VFLOAT );
    selRjrVars.makeBranch( "rjrSCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrSdphiDA", VFLOAT );
    selRjrVars.makeBranch( "rjrSdphiBV", VFLOAT );
    selRjrVars.makeBranch( "rjrPTS", VFLOAT );
    selRjrVars.makeBranch( "rjrPZS", VFLOAT );
    selRjrVars.makeBranch( "rjrPS", VFLOAT );

    selRjrVars.makeBranch( "rjrEVa", VFLOAT );
    selRjrVars.makeBranch( "rjrEVb", VFLOAT );
    selRjrVars.makeBranch( "rjrPVa", VFLOAT );
    selRjrVars.makeBranch( "rjrPVb", VFLOAT );

    selRjrVars.makeBranch( "rjrNPhotons", VINT );
    selRjrVars.makeBranch( "rjrMET", VFLOAT );
    selRjrVars.makeBranch( "rjrNJets", VINT );
    selRjrVars.makeBranch( "rjrLeadPhoLocation", VINT);
    selRjrVars.makeBranch( "rjrSubPhoLocation", VINT);
    selRjrVars.makeBranch( "rjrNVisObjects", VINT );
    selRjrVars.makeBranch( "rjrABSide", VBOOL );

    selRjrVars.makeBranch( "rjrVisPhoSide", VBOOL );
    selRjrVars.makeBranch( "rjrVisJetSide", VBOOL );

    selRjrVars.makeBranch( "rjrPVlab", VFLOAT );
    selRjrVars.makeBranch( "rjrDphiMETV", VFLOAT );
    selRjrVars.makeBranch( "rjrDPhiSI", VFLOAT );
    selRjrVars.makeBranch( "rjrCleaningVeto0", VBOOL );
    selRjrVars.makeBranch( "rjrCleaningVeto1", VBOOL );
    selRjrVars.makeBranch( "rjrCleaningVeto2", VBOOL );

    selRjrVars.makeBranch( "rjrNJetsJb", VINT );
    selRjrVars.makeBranch( "rjrNJetsJa", VINT );

    selRjrVars.makeBranch( "rjr_pHts20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts40", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts41", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts21a", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts21b", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts11a", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts11b", VFLOAT );

    selRjrVars.makeBranch( "rjr_pHs20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs40", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs41", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs21a", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs21b", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs11a", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs11b", VFLOAT );

    selRjrVars.makeBranch( "rjr_pHxa20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxb20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxa21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxb21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxa10", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxb10", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxa11", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHxb11", VFLOAT );

    selRjrVars.makeBranch( "rjr_pHtxa20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxb20", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxa21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxb21", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxa10", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxb10", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxa11", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHtxb11", VFLOAT );

    selRjrVars.makeBranch( "rjr_pHs22", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHs11", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts22", VFLOAT );
    selRjrVars.makeBranch( "rjr_pHts11", VFLOAT );

    selRjrVars.makeBranch( "rjr_p1Pts", VFLOAT );
    selRjrVars.makeBranch( "rjr_p2Pts", VFLOAT );

    selRjrVars.makeBranch( "rjr_p1Ptxa", VFLOAT );
    selRjrVars.makeBranch( "rjr_p1Ptxb", VFLOAT );
    selRjrVars.makeBranch( "rjr_p2Ptxa", VFLOAT );
    selRjrVars.makeBranch( "rjr_p2Ptxb", VFLOAT );

    selRjrVars.makeBranch( "rjr_p1Pxa", VFLOAT );
    selRjrVars.makeBranch( "rjr_p1Pxb", VFLOAT );
    selRjrVars.makeBranch( "rjr_p2Pxa", VFLOAT );
    selRjrVars.makeBranch( "rjr_p2Pxb", VFLOAT );

    selRjrVars.makeBranch( "rjrAX2Diff","rjr_Rdiff", VFLOAT );

    selRjrVars.makeBranch( "rjr_Mr", VFLOAT );
    selRjrVars.makeBranch( "rjr_R", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rv", VFLOAT );
    selRjrVars.makeBranch( "rjr_Ms", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rs", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rx", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rm", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rxa", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rxb", VFLOAT );

    selRjrVars.makeBranch( "rjr_sJ1a_px", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1a_py", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1a_pz", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1a_e", VFLOAT );

    selRjrVars.makeBranch( "rjr_sJ2a_px", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2a_py", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2a_pz", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2a_e", VFLOAT );

    selRjrVars.makeBranch( "rjr_sJ1b_px", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1b_py", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1b_pz", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ1b_e", VFLOAT );

    selRjrVars.makeBranch( "rjr_sJ2b_px", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2b_py", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2b_pz", VFLOAT );
    selRjrVars.makeBranch( "rjr_sJ2b_e", VFLOAT );

    selRjrVars.makeBranch( "rjr_Mva", VFLOAT );
    selRjrVars.makeBranch( "rjr_Mvb", VFLOAT );

    selRjrVars.makeBranch( "rjrN2Px", VFLOAT );
    selRjrVars.makeBranch( "rjrN2Py", VFLOAT );
    selRjrVars.makeBranch( "rjrN2Pz", VFLOAT );

    selRjrVars.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

