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

void KUCMSAodSkimmer::processRJRISR(){

  //bool verbose = true;
  bool verbose = false;

  if( RJRDEBUG || verbose ) std::cout << "Processing RJR ISR event varibles " << std::endl;
  // ------------------  Loading &  Processing RJR event varibles  ---------------------------------------------------------

  // intilize out branches
  selRjrIsrVars.clearBranches(); // <<<<<<<   must do

  // process event
  LAB_c->ClearEvent();

  auto sPhoPt = geVects( "selPhoPt");
  auto sPhoEta = geVects( "selPhoEta");
  auto sPhoPhi = geVects( "selPhoPhi");
  auto sPhoMx = geVects( "selPhoEMx");
  auto sPhoMy = geVects( "selPhoEMy");

  auto nSelPhotons = geCnts("nSelPhotons");

  auto phoRMetCPx = geVars("metPx"); //selMet.getFLBranchValue("metPx");
  auto phoRMetCPy = geVars("metPy"); //selMet.getFLBranchValue("metPy");
  float unCorMet = hypo(phoRMetCPx,phoRMetCPy);

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
        TLorentzVector phojet;
        phojet.SetPtEtaPhiM( sPhoPt[spidx], sPhoEta[spidx], sPhoPhi[spidx], 0 );
        jetID.push_back( COMB_J_c->AddLabFrameFourVector(phojet) );
        if( spidx < nRJRPhos ) pho4vec.push_back(phojet);
		else jet4vec.push_back(phojet);		

  	}//<<>>for( spidx = 0; spidx < nSelPhotons; spidx++ )
  }//<<>>if( not zsig )
  else std::cout << " --- !!!  Photons not loaded into RJR !!!!!!!!!!!!!!!!!!" << std::endl;

  if( RJRDEBUG ) std::cout << " - Loading MET." << std::endl;
  if( verbose && nSelPhotons > 0 ){
    std::cout << " - Loading MET lPt: " << sPhoPt[0] << " lPhi: " << sPhoPhi[0] << std::endl;
    std::cout << " - Loading MET slPt: " << sPhoPt[1] << " slPhi: " << sPhoPhi[1] << std::endl; 
    std::cout << " - Loading MET x: " << geVars("metPx") << " -> " << phoRMetCPx;
    std::cout << " y: " << geVars("metPy") << " -> " << phoRMetCPy << std::endl;
  }//<<>>if( verbose )

  TVector3 ETMiss(phoRMetCPx,phoRMetCPy,0);
  INV_c->SetLabFrameThreeVector(ETMiss);

  auto nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
  auto selJetPt = geVects( "selJetPt");
  auto selJetEta = geVects( "selJetEta");
  auto selJetPhi = geVects( "selJetPhi");
  auto selJetMass = geVects( "selJetMass");

  if( RJRDEBUG ) std::cout << " - Loading Jets." << std::endl;

  std::vector< TLorentzVector > rjr_jets;
  for( uInt it = 0; it < nSelJets; it++ ){
	float sjetPt = selJetPt[it]; //selJets.getFLBranchValue( "selJetPt", it );
    float sjetEta = selJetEta[it]; //selJets.getFLBranchValue( "selJetEta", it );
    float sjetPhi = selJetPhi[it]; //selJets.getFLBranchValue( "selJetPhi", it ); 
    float sjetMass = ( selJetMass[it] > 0 ) ? selJetMass[it] : 0; //selJets.getFLBranchValue( "selJetMass", it );   
    TLorentzVector jet;
    jet.SetPtEtaPhiM( sjetPt, sjetEta, sjetPhi, sjetMass );
	rjr_jets.push_back(jet);
    if( verbose ) std::cout << " - Loading Jet Pt: " << sjetPt << " Eta: " << sjetEta;
    if( verbose ) std::cout << " Phi: " << sjetPhi << " M: " << sjetMass << std::endl;
    ////if( it == 0 ){ leadJetKey.push_back( COMB_J_c->AddLabFrameFourVector(jet) ); jetID.push_back(leadJetKey[0]); } 
    ////else jetID.push_back(COMB_J_c->AddLabFrameFourVector(jet));
    ////jet4vec.push_back(jet); 
  }//<<>>for( int i = 0; i < nSelJets; i++ )

  if( RJRDEBUG ) std::cout << " -- nRJRJets Pre : " << nSelJets << std::endl;
  std::vector<RFKey> leadJetKey;
  //std::vector< TLorentzVector > merged_rjr_jets = BinaryMerge( rjr_jets, 13 -  nRJRPhos );
  BinaryMergeInPlace( rjr_jets, 16 -  nRJRPhos );
  nSelJets = rjr_jets.size();
  if( RJRDEBUG ) std::cout << " -- nRJRJets Post : " << nSelJets << std::endl;
  for( uInt it = 0; it < nSelJets; it++ ){
    ////auto sjetPt = selJetPt[it]; //selJets.getFLBranchValue( "selJetPt", it );
    ////auto sjetEta = selJetEta[it]; //selJets.getFLBranchValue( "selJetEta", it );
    ////auto sjetPhi = selJetPhi[it]; //selJets.getFLBranchValue( "selJetPhi", it ); 
    ////auto sjetMass = ( selJetMass[it] > 0 ) ? selJetMass[it] : 0; //selJets.getFLBranchValue( "selJetMass", it );   
    TLorentzVector jet = rjr_jets[it];
    ////jet.SetPtEtaPhiM( sjetPt, sjetEta, sjetPhi, sjetMass );
    if( verbose ) std::cout << " - Loading Jet Pt: " << jet.Pt() << " Eta: " << jet.Eta();
    if( verbose ) std::cout << " Phi: " << jet.Phi() << " M: " << jet.M() << std::endl;
    if( it == 0 ){ leadJetKey.push_back( COMB_J_c->AddLabFrameFourVector(jet) ); jetID.push_back(leadJetKey[0]); }
    else jetID.push_back(COMB_J_c->AddLabFrameFourVector(jet));
    jet4vec.push_back(jet);
  }//<<>>for( int i = 0; i < nSelJets; i++ )

  if( RJRDEBUG ) std::cout << " -- Processing RJR Tree ------ " << std::endl;

  if( !LAB_c->AnalyzeEvent() ) std::cout << "Something went wrong with tree event analysis" << std::endl;
	
  if( RJRDEBUG ) std::cout << " - Getting RJR Object Location Information." << std::endl;
  // ---------  Finished Processing RJR varibles --------------------------------------------------------------------

  uInt firstJet = 0;
  if( nRJRPhos > 0 ) firstJet = ( nRJRPhos == 2 ) ? 2 : 1 ;
  int nJetsJa = 0;
  int nJetsJb = 0;
  int nJetsISR = 0;
  int nPhosISR = 0;
  bool isALeadPhoSide = true;	
  int subPhoSide = ( nRJRPhos > 1 ) ? 0 : -1;
  int leadPhoSide = ( nRJRPhos > 0 ) ? 0 : -1;
  int leadPhoIndx = ( nRJRPhos > 0 ) ? 0 : -1;
  int xtraPhos = ( nSelPhotons > nRJRPhos ) ? nSelPhotons - nRJRPhos : 0;
  int nVisObjects = jetID.size();
  if( nVisObjects != nRJRPhos + nSelJets + xtraPhos ){ 
    std::cout << " !!!!!!    nVisObjects != ( nRJRPhos + nSelJets )  !!!!!!!!! " << std::endl;
    std::cout << " !!!!!!    " << nVisObjects << " != ( " << nRJRPhos + xtraPhos << " + " << nSelJets << " )  !!!!!!!!! " << std::endl;
  }//<<>>if( nVisObjects != ( nRJRPhos + nSelJets ) )

  if( nRJRPhos > 0 ){
	if( COMB_J_c->GetFrame(jetID[0]) == *ISR_c ){ nPhosISR++; isALeadPhoSide = false; leadPhoSide = 2; }	
	else if( ( COMB_J_c->GetFrame(jetID[0]) == *J1b_c || COMB_J_c->GetFrame(jetID[0]) == *J2b_c ) ){ 
	isALeadPhoSide = false; 
	leadPhoSide = 1;
	}//<<>>if( ( COMB_J_c->GetFrame(jetID[0]) == *J1b_c || COMB_J_c->GetFrame(jetID[0]) == *J2b_c ) ) 
  }//<<>>if( leadPhoIndx > -1 &&  COMB_J_c->GetFrame(jetID[0]) == *J1b_c || COMB_J_c->GetFrame(jetID[0]) == *J2b_c )
  //  -- redifine A & B side for jets to be : is jet on lead pho side -> A ; is not on lead pho side -> B // done
  for( uInt it = firstJet; it < nVisObjects; it++ ){	

	if( COMB_J_c->GetFrame(jetID[it]) == *ISR_c ){ nJetsISR++; continue; }
    if( COMB_J_c->GetFrame(jetID[it]) == *J1a_c || COMB_J_c->GetFrame(jetID[it]) == *J2a_c  ){ isALeadPhoSide ? nJetsJa++ : nJetsJb++; } // one per frame 
    else { isALeadPhoSide ? nJetsJb++ : nJetsJa++; }

  }//<<>>for( int i = 0; i < nSelJets; i++ )
  float abDiffSide = isALeadPhoSide ? 1 : -1;
  if( ( nRJRPhos > 1 ) ){ 
	if( COMB_J_c->GetFrame(jetID[1]) == *ISR_c ){ subPhoSide = 2; nPhosISR++; }
	else subPhoSide = ( COMB_J_c->GetFrame(jetID[1]) == *J1a_c || COMB_J_c->GetFrame(jetID[1]) == *J2a_c ) ? 0 : 1;
  }//<<>>if( ( nRJRPhos > 1 ) )

  int nIsrVisObjects = nPhosISR + nJetsISR;
  int nSVisObjects = nVisObjects - nIsrVisObjects;

  selRjrIsrVars.fillBranch( "rjrIsrNISRPhos", nPhosISR );
  selRjrIsrVars.fillBranch( "rjrIsrABSide", isALeadPhoSide );
  selRjrIsrVars.fillBranch( "rjrIsrNJetsJa", nJetsJa );
  selRjrIsrVars.fillBranch( "rjrIsrNJetsJb", nJetsJb );
  selRjrIsrVars.fillBranch( "rjrIsrNJetsISR", nJetsISR );
  selRjrIsrVars.fillBranch( "rjrIsrNJets", int(geCnts("nSelJets")) );
  selRjrIsrVars.fillBranch( "rjrIsrNPhotons", nRJRPhos );
  selRjrIsrVars.fillBranch( "rjrIsrNVisObjects", nVisObjects );
  selRjrIsrVars.fillBranch( "rjrIsrLeadPhoLocation", leadPhoSide );
  selRjrIsrVars.fillBranch( "rjrIsrSubPhoLocation", subPhoSide );
  selRjrIsrVars.fillBranch( "rjrIsr_nSVisObjects", nSVisObjects );
  selRjrIsrVars.fillBranch( "rjrIsr_nIsrVisObjects", nIsrVisObjects );

  std::vector<int> phoside;
  std::vector<int> jetside;
  for( uInt it = 0; it < nVisObjects; it++ ){

    int onAside = -1;
    int whichside = -1;
    ////if( COMB_J_c->GetFrame(jetID[it]) == *J1b || COMB_J_c->GetFrame(jetID[it]) == *J2b  ) onAside = false;
    if( COMB_J_c->GetFrame(jetID[it]) == *J1a_c ) whichside = 0;
    if( COMB_J_c->GetFrame(jetID[it]) == *J2a_c ) whichside = 1;
    if( COMB_J_c->GetFrame(jetID[it]) == *J1b_c ) whichside = 2;
    if( COMB_J_c->GetFrame(jetID[it]) == *J2b_c ) whichside = 3;
	if( COMB_J_c->GetFrame(jetID[it]) == *ISR_c ) whichside = 4;
    if( whichside < 0 ) std::cout << " !!!!!!!  Visible object is not on a side !!!!!!! " << std::endl;
    if( whichside != 4 ) onAside = ( whichside > 1 ) ? 0 : 1;
    if( it < firstJet ){ 
      	selRjrIsrVars.fillBranch( "rjrIsrVisPhoSide", onAside );
      	phoside.push_back(whichside);
        //std::cout << " -- pho Side : " << whichside << " for # " << nRJRPhos << std::endl;
    } else { 
      	selRjrIsrVars.fillBranch( "rjrIsrVisJetSide", onAside );
      	jetside.push_back(whichside);
		//std::cout << " -- jet Side : " << whichside << " for # " << nVisObjects - nRJRPhos << std::endl;
    }//<<>>if( it < nSelPhotons )
  }//<<>>for( uInt it = 0; it < nVisObjects; it++ )

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles." << std::endl;

  float m_MS = CM_c->GetMass();
  float m_cosS  = CM_c->GetCosDecayAngle();
  float m_dphiS = CM_c->GetDeltaPhiDecayAngle();
  float m_dphiSI  = CM_c->GetDeltaPhiBoostVisible();
  
  float m_MX2a = X2a_c->GetMass();
  float m_cosX2a = X2a_c->GetCosDecayAngle();
  float m_MX2b = X2b_c->GetMass();
  float m_cosX2b = X2b_c->GetCosDecayAngle();

  float m_PTS = CM_c->GetFourVector().Pt();
  float m_PZS = CM_c->GetFourVector().Pz();
  float m_PS = CM_c->GetMomentum(*LAB_c);

  selRjrIsrVars.fillBranch( "rjrIsrX2aMass", m_MX2a );
  selRjrIsrVars.fillBranch( "rjrIsrX2aCosA", m_cosX2a );
  selRjrIsrVars.fillBranch( "rjrIsrX2bMass", m_MX2b );
  selRjrIsrVars.fillBranch( "rjrIsrX2bCosA", m_cosX2b );

  selRjrIsrVars.fillBranch( "rjrIsrSMass", m_MS );
  selRjrIsrVars.fillBranch( "rjrIsrSCosA", m_cosS );

  selRjrIsrVars.fillBranch( "rjrIsrSdphiDA", m_dphiS  );
  selRjrIsrVars.fillBranch( "rjrIsrSdphiBV", m_dphiSI ); // !!!!!!!!!!!!!  dphiCMI
  selRjrIsrVars.fillBranch( "rjrIsrDPhiSI", m_dphiSI ); // !!!!!!!!!!!!!  dphiCMI
  selRjrIsrVars.fillBranch( "rjrIsrPTS", m_PTS ); // !!!!!!!!!!!!!!!!!!  P_{t}^{S}
  selRjrIsrVars.fillBranch( "rjrIsrPZS", m_PZS );
  selRjrIsrVars.fillBranch( "rjrIsrPS", m_PS );

  //if( type == 0 ) hist2d[2]->Fill(m_dphiSI,m_PTS);
  //if( type == 1 ) hist2d[3]->Fill(m_dphiSI,m_PTS); 

  bool rjrCleaningVeto0 =  ( m_PTS > 150.0 ) ? true : false; 
  selRjrIsrVars.fillBranch( "rjrIsrCleaningVeto0", rjrCleaningVeto0 );

  bool rjrCleaningVeto1 = false;
  if( ( m_dphiSI <= 0.6 ) && ( m_PTS > (250.0 * m_dphiSI) ) ) rjrCleaningVeto1 = true;
  if( ( m_dphiSI > 0.6 ) && ( m_dphiSI < 2.4 ) && ( m_PTS > 150.0 ) ) rjrCleaningVeto1 = true; 
  if( ( m_dphiSI >= 2.4 ) && ( m_PTS > (150.0 - (187.5 * (m_dphiSI - 2.4))) ) ) rjrCleaningVeto1 = true;
  selRjrIsrVars.fillBranch( "rjrIsrCleaningVeto1", rjrCleaningVeto1 );  

  bool rjrCleaningVeto2 = false; 
  if( ( m_dphiSI <= 0.3 ) && ( m_PTS > (500.0 * m_dphiSI) ) ) rjrCleaningVeto2 = true; 
  if( ( m_dphiSI > 0.3 ) && ( m_dphiSI < 2.8 ) && ( m_PTS > 150.0 ) ) rjrCleaningVeto2 = true;  
  if( ( m_dphiSI >= 2.8 ) && ( m_PTS > (150.0 - (375.0 * (m_dphiSI - 2.8))) ) ) rjrCleaningVeto2 = true; 
  selRjrIsrVars.fillBranch( "rjrIsrCleaningVeto2", rjrCleaningVeto2 );   


  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 2." << std::endl;

  std::vector< TLorentzVector > hp4;
  std::vector< TLorentzVector > hxp4;
	
	int nPhov = pho4vec.size();
	if( nPhov != phoside.size() ) std::cout << " !!!!!!!  pho4vec - phoside count mismatch !!!!!!!" << std::endl; 
	int nJetv = jet4vec.size();
    if( nJetv != jetside.size() ) std::cout << " !!!!!!!  jet4vec - jetside count mismatch !!!!!!!" << std::endl;
	//std::cout << " nPhov " << nPhov << " phoside " << phoside.size() << " nJetv " << nJetv << " jetside " << jetside.size() << std::endl; 

	TLorentzVector nulltv(0,0,0,0);
	std::vector< TLorentzVector > pJetSideSum4Vec{nulltv,nulltv,nulltv,nulltv,nulltv};
	std::vector< TLorentzVector > phoSideSum4Vec{nulltv,nulltv,nulltv,nulltv,nulltv};
    std::vector< TLorentzVector > jetSideSum4Vec{nulltv,nulltv,nulltv,nulltv,nulltv};
	for ( int pit = 0; pit < nPhov; pit++ ){ phoSideSum4Vec[phoside[pit]] += pho4vec[pit]; pJetSideSum4Vec[phoside[pit]] += pho4vec[pit];}
    for ( int jit = 0; jit < nJetv; jit++ ){ jetSideSum4Vec[jetside[jit]] += jet4vec[jit]; pJetSideSum4Vec[jetside[jit]] += jet4vec[jit];}

    hp4.push_back(S_c->GetFourVector(pJetSideSum4Vec[0]));//0 J1a
    hp4.push_back(S_c->GetFourVector(pJetSideSum4Vec[1]));//1 J2a
    hp4.push_back(S_c->GetFourVector(pJetSideSum4Vec[2]));//2 J1b
    hp4.push_back(S_c->GetFourVector(pJetSideSum4Vec[3]));//3 J2b

    selRjrIsrVars.fillBranch( "rjrIsr_sJ1a_px", float(hp4[0].Px()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1a_py", float(hp4[0].Py()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1a_pz", float(hp4[0].Pz()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1a_e", float(hp4[0].E()) );

    selRjrIsrVars.fillBranch( "rjrIsr_sJ2a_px", float(hp4[1].Px()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2a_py", float(hp4[1].Py()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2a_pz", float(hp4[1].Pz()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2a_e", float(hp4[1].E()) );

    selRjrIsrVars.fillBranch( "rjrIsr_sJ1b_px", float(hp4[2].Px()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1b_py", float(hp4[2].Py()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1b_pz", float(hp4[2].Pz()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ1b_e", float(hp4[2].E()) );

    selRjrIsrVars.fillBranch( "rjrIsr_sJ2b_px", float(hp4[3].Px()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2b_py", float(hp4[3].Py()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2b_pz", float(hp4[3].Pz()) );
    selRjrIsrVars.fillBranch( "rjrIsr_sJ2b_e", float(hp4[3].E()) );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 3." << std::endl;

	float x1sp = ( X1a_c->GetFourVector(*S_c) + X1b_c->GetFourVector(*S_c) ).P();
	float x1asp = ( X1a_c->GetFourVector(*S_c) ).P();
    float x1bsp = ( X1b_c->GetFourVector(*S_c) ).P();

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

//
//  // ISR related variables
//  if(m_Njet_ISR > 0){
//    TVector3 vPISR = S->GetFourVector(*CM).Vect();
//    TVector3 vPINV = (X1a->GetFourVector(*CM)+X1b->GetFourVector(*CM)).Vect();
//
//    m_PISR = vPISR.Mag();
//    m_MISR = ISR->GetMass();
//    m_RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
//
//    // transverse
//    TVector3 vPTISR = S->GetTransverseFourVector(*CM).Vect();
//    TVector3 vPTINV = (X1a->GetTransverseFourVector(*CM)+X1b->GetTransverseFourVector(*CM)).Vect();
//
//    m_PTISR  = vPTISR.Mag();
//    m_RISRT  = fabs(vPTINV.Dot(vPTISR.Unit())) / vPTISR.Mag();
//    
//    if((std::isnan(m_RISR) || std::isnan(m_Mperp))){
//      cout << "VAR NAN " << vPTISR.Mag() << " " << vPTINV.Mag() << " NjetS=" << m_Njet_S << " Njeta=" << m_Njet_a << " Njetb=" << m_Njet_b << " Nlep=" << m_Nlep << " Nlepa=" << m_Nlep_a << " Nlepb=" << m_Nlep_b << " " << m_Mperp << endl;
//      cout << m_RISR << endl;
//    }
//
//    TVector3 isr_t = ISR->GetTransverseFourVector(*S).Vect();
//    TVector3 isr   = ISR->GetFourVector(*S).Vect();
//


	float msperp = -1;
	float m_PISR = -1;
    float m_PTISR = -1;
	float m_RISR = -1;
	if( nJetsJa > 0 || nJetsJb > 0 ){ // must have jets in ja or jb 

        TVector3 vPISR = S_c->GetFourVector(*CM_c).Vect();
        TVector3 vPINV = (X1a_c->GetFourVector(*CM_c)+X1b_c->GetFourVector(*CM_c)).Vect();
        TVector3 vPTISR = S_c->GetTransverseFourVector(*CM_c).Vect();
        m_PISR = vPISR.Mag();

		m_PTISR = vPTISR.Mag();
		m_RISR = std::fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();

		TLorentzVector j1as =  S_c->GetFourVector(pJetSideSum4Vec[0]);
		TLorentzVector j2as =  S_c->GetFourVector(pJetSideSum4Vec[1]);
        TLorentzVector j1bs =  S_c->GetFourVector(pJetSideSum4Vec[2]);
        TLorentzVector j2bs =  S_c->GetFourVector(pJetSideSum4Vec[3]);
		TLorentzVector x1s =  X1a_c->GetFourVector(*S_c) + X1b_c->GetFourVector(*S_c);

  		TVector3 boostVis = (j1as+j1bs+j2as+j2bs).BoostVector();
  		TVector3 boostInv = (x1s).BoostVector();
  		TVector3 daBoost = vPISR.Unit();

  		boostVis = (boostVis.Dot(daBoost))*daBoost;
  		boostInv = (boostInv.Dot(daBoost))*daBoost;

  		if((!std::isnan(boostVis.Mag())) && (boostVis.Mag() < 1)){
        	j1as.Boost(-boostVis);
        	j1bs.Boost(-boostVis);
        	j2as.Boost(-boostVis);
        	j2bs.Boost(-boostVis);
  		} else {
        	j1as.SetVectM(TVector3(0.,0.,0.),std::max(0.,j1as.M()));
        	j1bs.SetVectM(TVector3(0.,0.,0.),std::max(0.,j1bs.M()));
        	j2as.SetVectM(TVector3(0.,0.,0.),std::max(0.,j2as.M()));
        	j2bs.SetVectM(TVector3(0.,0.,0.),std::max(0.,j2bs.M()));
  		}//<<>>if((!std::isnan(boostVis.Mag())) && (boostVis.Mag() < 1))

  		if( (!std::isnan(boostInv.Mag())) && (boostInv.Mag() < 1)){ x1s.Boost(-boostInv); } 
		else { x1s.SetVectM(TVector3(0.,0.,0.),std::max(0.,x1s.M())); }

        msperp = j1as.P() + j1bs.P() + j2as.P() + j2bs.P() + x1s.P();

	}//<<>>if( true )

    selRjrIsrVars.fillBranch( "rjrIsr_MsPerp", msperp );
    selRjrIsrVars.fillBranch( "rjrIsr_PIsr", m_PISR );
    selRjrIsrVars.fillBranch( "rjrIsr_PtIsr", m_PTISR );
    selRjrIsrVars.fillBranch( "rjrIsr_RIsr", m_RISR );

    selRjrIsrVars.fillBranch( "rjrIsr_pHs20", phs20 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs40", phs40 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs21", phs21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs41", phs41 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs21a", phs21a );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs21b", phs21b );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs11a", phs11a );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs11b", phs11b );

  	//if( RJRDEBUG ) std::cout << " - Getting RJR varibles 4." << std::endl;


    hxp4.push_back(X2a_c->GetFourVector(pJetSideSum4Vec[0]));//0
    hxp4.push_back(X2a_c->GetFourVector(pJetSideSum4Vec[1]));//1
    hxp4.push_back(X1a_c->GetFourVector(*X2a_c));//2
    hxp4.push_back(X2b_c->GetFourVector(pJetSideSum4Vec[2]));//3
    hxp4.push_back(X2b_c->GetFourVector(pJetSideSum4Vec[3]));//4
    hxp4.push_back(X1b_c->GetFourVector(*X2b_c));//5

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

    selRjrIsrVars.fillBranch( "rjrIsr_pHxa20", phxa20 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxb20", phxb20 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxa21", phxa21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxb21", phxb21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxa10", phxa10 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxb10", phxb10 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxa11", phxa11 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHxb11", phxb11 );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 5." << std::endl;

    float pj1apt = S_c->GetTransverseMomentum(pJetSideSum4Vec[0]);
    float pj2apt = S_c->GetTransverseMomentum(pJetSideSum4Vec[1]);
    float pj1bpt = S_c->GetTransverseMomentum(pJetSideSum4Vec[2]);
    float pj2bpt = S_c->GetTransverseMomentum(pJetSideSum4Vec[3]);
    float pjapt = S_c->GetTransverseMomentum(pJetSideSum4Vec[0]+pJetSideSum4Vec[1]);
    float pjbpt = S_c->GetTransverseMomentum(pJetSideSum4Vec[2]+pJetSideSum4Vec[3]);

	float x1spt = S_c->GetTransverseMomentum( X1a_c->GetFourVector(*LAB_c) + X1b_c->GetFourVector(*LAB_c) );
    float x1aspt = S_c->GetTransverseMomentum( X1a_c->GetFourVector(*LAB_c) );
    float x1bspt = S_c->GetTransverseMomentum( X1b_c->GetFourVector(*LAB_c) );

    float phts20 = pjapt + pjbpt;
    float phts40 = pj1apt + pj2apt + pj1bpt + pj2bpt;
    float phts21 = phts20 + x1spt;
    float phts41 = phts40 + x1spt;

    float phts21a = pj1apt + pj2apt + x1aspt;
    float phts21b = pj1bpt + pj2bpt + x1bspt;
    float phts11a = pjapt + x1aspt;
    float phts11b = pjbpt + x1bspt;

    selRjrIsrVars.fillBranch( "rjrIsr_pHts20", pjapt + pjbpt );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts40", pj1apt + pj2apt + pj1bpt + pj2bpt );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts21", phts21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts41", phts41 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts21a", phts21a );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts21b", phts21b );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts11a", phts11a );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts11b", phts11b );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 6." << std::endl;

    float ptj1apt = X2a_c->GetTransverseMomentum(pJetSideSum4Vec[0]);
    float ptj2apt = X2a_c->GetTransverseMomentum(pJetSideSum4Vec[1]);
    float ptx1apt = X2a_c->GetTransverseMomentum(X1a_c->GetFourVector(*LAB_c));
    float ptj1bpt = X2b_c->GetTransverseMomentum(pJetSideSum4Vec[2]);
    float ptj2bpt = X2b_c->GetTransverseMomentum(pJetSideSum4Vec[3]);
    float ptx1bpt = X2b_c->GetTransverseMomentum(X1b_c->GetFourVector(*LAB_c));
    float ptj12apt = X2a_c->GetTransverseMomentum(pJetSideSum4Vec[0]+pJetSideSum4Vec[1]);
    float ptj12bpt = X2b_c->GetTransverseMomentum(pJetSideSum4Vec[2]+pJetSideSum4Vec[3]);

    float phtxa20 = ptj1apt + ptj2apt;
    float phtxb20 = ptj1bpt + ptj2bpt;
    float phtxa21 = phtxa20 + ptx1apt;
    float phtxb21 = phtxb20 + ptx1bpt;
    float phtxa10 = ptj12apt;
    float phtxb10 = ptj12bpt;
    float phtxa11 = ptj12apt + ptx1apt;
    float phtxb11 = ptj12bpt + ptx1bpt;

    selRjrIsrVars.fillBranch( "rjrIsr_pHtxa20", phtxa20 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxb20", phtxb20 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxa21", phtxa21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxb21", phtxb21 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxa10", phtxa10 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxb10", phtxb10 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxa11", phtxa11 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHtxb11", phtxb11 );


  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 7." << std::endl;

	float phopt1 = 0;
    float phopt2 = 0;

	if( nPhov > 0 ) phopt1 = S_c->GetTransverseMomentum( pho4vec[0] );
    if( nPhov > 1 ) phopt2 = S_c->GetTransverseMomentum( pho4vec[1] );
    //if( nPhov > 0 ) phopt1 = S_c->GetFourVector( pho4vec[0] ).P();
    //if( nPhov > 1 ) phopt2 = S_c->GetFourVector( pho4vec[1] ).P();

    selRjrIsrVars.fillBranch( "rjrIsr_p1Pts", phopt1 );
    selRjrIsrVars.fillBranch( "rjrIsr_p2Pts", phopt2 );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 8." << std::endl;

    float phopt1xa = 0;
    float phopt1xb = 0;
    float phopt2xa = 0;
    float phopt2xb = 0;

    if( nPhov > 0 ){ 
		if( phoside[0] < 2 ) phopt1xa = X2a_c->GetTransverseMomentum( pho4vec[0] );
		else phopt1xb = X2b_c->GetTransverseMomentum( pho4vec[0] );
	}//<<>>if( nPhov > 0 )
    if( nPhov > 1 ){ 
		if( phoside[1] < 2 ) phopt2xa = X2a_c->GetTransverseMomentum( pho4vec[1] );
        else phopt2xb = X2b_c->GetTransverseMomentum( pho4vec[1] );
	}//<<>>if( nPhov > 1 )

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 9." << std::endl;

    float phop1xa = 0;
    float phop1xb = 0;
    float phop2xa = 0;
    float phop2xb = 0;

    if( nPhov > 0 ){
        if( phoside[0] < 2 ) phop1xa = X2a_c->GetFourVector( pho4vec[0] ).P();
        else phop1xb = X2b_c->GetFourVector( pho4vec[0] ).P();
    }//<<>>if( nPhov > 0 )
    if( nPhov > 1 ){
        if( phoside[1] < 2 ) phop2xa = X2a_c->GetFourVector( pho4vec[1] ).P();
        else phop2xb = X2b_c->GetFourVector( pho4vec[1] ).P();
    }//<<>>if( nPhov > 1 )


    selRjrIsrVars.fillBranch( "rjrIsr_p1Ptxa", phopt1xa );
    selRjrIsrVars.fillBranch( "rjrIsr_p1Ptxb", phopt1xb );
    selRjrIsrVars.fillBranch( "rjrIsr_p2Ptxa", phopt2xa );
    selRjrIsrVars.fillBranch( "rjrIsr_p2Ptxb", phopt2xb );

    selRjrIsrVars.fillBranch( "rjrIsr_p1Pxa", phop1xa );
    selRjrIsrVars.fillBranch( "rjrIsr_p1Pxb", phop1xb );
    selRjrIsrVars.fillBranch( "rjrIsr_p2Pxa", phop2xa );
    selRjrIsrVars.fillBranch( "rjrIsr_p2Pxb", phop2xb );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 10." << std::endl;

	std::vector< TLorentzVector > p4;
	p4.push_back(Ja_c->GetFourVector(*S_c));
    p4.push_back(Jb_c->GetFourVector(*S_c));
    p4.push_back(X1a_c->GetFourVector(*S_c));
    p4.push_back(X1b_c->GetFourVector(*S_c));

	float pHs22 = p4[0].P() + p4[1].P() + p4[2].P() + p4[3].P();
    float pHs11 =  ( p4[0] + p4[1] ).P() + ( p4[2] + p4[3] ).P();

    float pHJas = Ja_c->GetTransverseMomentum(*S_c);
    float pHJbs = Jb_c->GetTransverseMomentum(*S_c);
    float pHX1as = X1a_c->GetTransverseMomentum(*S_c);
    float pHX1bs = X1b_c->GetTransverseMomentum(*S_c);
	float pHJs = S_c->GetTransverseMomentum( Ja_c->GetFourVector(*LAB_c) + Jb_c->GetFourVector(*LAB_c) );
    float pHX1s = S_c->GetTransverseMomentum( X1a_c->GetFourVector(*LAB_c) + X1b_c->GetFourVector(*LAB_c) );

	float pHts22 = pHJas + pHJbs + pHX1as + pHX1bs;
	float pHts11 = pHJs + pHX1s;

	float Rs = (phxa11+phxb11)/pHs22;
	float Rxa = phxa11/phxa21;
    float Rxb = phxb11/phxb21;
	float Rx = hypo(phxa11/phxa21,phxb11/phxb21);

    selRjrIsrVars.fillBranch( "rjrIsr_pHs22", pHs22 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHs11", pHs11 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts22", pHts22 );
    selRjrIsrVars.fillBranch( "rjrIsr_pHts11", pHts11 );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 11." << std::endl;

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

    float pf_pX2a = X2a_c->GetMomentum(*S_c);
    float pf_pX2b = X2b_c->GetMomentum(*S_c);

    selRjrIsrVars.fillBranch( "rjrIsr_Mr", a_MS );
    selRjrIsrVars.fillBranch( "rjrIsr_R", AX2NQSum );
    selRjrIsrVars.fillBranch( "rjrIsrAX2Diff", AX2Diff );

    selRjrIsrVars.fillBranch( "rjrIsr_Ms", phs41 );
    selRjrIsrVars.fillBranch( "rjrIsr_Rs", Rs );
    selRjrIsrVars.fillBranch( "rjrIsr_Rx", Rx );
    selRjrIsrVars.fillBranch( "rjrIsr_Rxa", Rxa );
    selRjrIsrVars.fillBranch( "rjrIsr_Rxb", Rxb );

    selRjrIsrVars.fillBranch( "rjrIsrX2aPs", pf_pX2a );
    selRjrIsrVars.fillBranch( "rjrIsrX2bPs", pf_pX2b );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 12." << std::endl;

  	float m_EVa = X2a_c->GetListVisibleFrames().GetFourVector(*X2a_c).E();
  	float m_EVb = X2b_c->GetListVisibleFrames().GetFourVector(*X2b_c).E();
  	float m_PVa = X2a_c->GetListVisibleFrames().GetFourVector(*X2a_c).P();
  	float m_PVb = X2b_c->GetListVisibleFrames().GetFourVector(*X2b_c).P();

  	float m_MX1a = X1a_c->GetMass();
  	float m_cosX1a = X1a_c->GetCosDecayAngle();
  	float m_MX1b = X1b_c->GetMass();
  	float m_cosX1b = X1b_c->GetCosDecayAngle();  
    float x1a_pt = X1a_c->GetFourVector().Pt();
    float x1a_eta = X1a_c->GetFourVector().Eta();
    float x1a_phi = X1a_c->GetFourVector().Phi();
    float x1b_pt = X1b_c->GetFourVector().Pt();
    float x1b_eta = X1b_c->GetFourVector().Eta();
    float x1b_phi = X1b_c->GetFourVector().Phi();

    float m_MVa = Ja_c->GetMass();
    float m_MVb = Jb_c->GetMass();

    float m_MVSum = std::sqrt((sq2(m_MVa)+sq2(m_MVb))/2);
	float m_MVNSum = 2*m_MVSum/a_MS;// -- !! Rv
    float m_MVab = std::sqrt(sq2(m_MVa)+sq2(m_MVb))/phs40;

    selRjrIsrVars.fillBranch( "rjrIsr_Rv", m_MVNSum );
    selRjrIsrVars.fillBranch( "rjrIsr_Mva", m_MVa );
    selRjrIsrVars.fillBranch( "rjrIsr_Mvb", m_MVb );
    selRjrIsrVars.fillBranch( "rjrIsr_Rm", m_MVab );

  	float m_PV_lab = S_c->GetListVisibleFrames().GetFourVector().P();
  	float m_dphiMET_V = S_c->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

    selRjrIsrVars.fillBranch( "rjrIsrX1aMass", m_MX1a );
    selRjrIsrVars.fillBranch( "rjrIsrX1aCosA", m_cosX1a );
    selRjrIsrVars.fillBranch( "rjrIsrX1bMass", m_MX1b );
    selRjrIsrVars.fillBranch( "rjrIsrX1bCosA", m_cosX1b );
    selRjrIsrVars.fillBranch( "rjrIsrX1a_Pt", x1a_pt );
    selRjrIsrVars.fillBranch( "rjrIsrX1a_Eta", x1a_eta );
    selRjrIsrVars.fillBranch( "rjrIsrX1a_Phi", x1a_phi );
    selRjrIsrVars.fillBranch( "rjrIsrX1b_Pt", x1b_pt );
    selRjrIsrVars.fillBranch( "rjrIsrX1b_Eta", x1b_eta );
    selRjrIsrVars.fillBranch( "rjrIsrX1b_Phi", x1b_phi );

    selRjrIsrVars.fillBranch( "rjrIsrEVa", m_EVa );
    selRjrIsrVars.fillBranch( "rjrIsrEVb", m_EVb );
    selRjrIsrVars.fillBranch( "rjrIsrPVa", m_PVa );
    selRjrIsrVars.fillBranch( "rjrIsrPVb", m_PVb );

    selRjrIsrVars.fillBranch( "rjrIsrPVlab", m_PV_lab );
    selRjrIsrVars.fillBranch( "rjrIsrDphiMETV", m_dphiMET_V );

  //if( RJRDEBUG ) std::cout << " - Getting RJR varibles 13." << std::endl;

	float rjrMet = hypo(phoRMetCPx,phoRMetCPy);

    selRjrIsrVars.fillBranch( "rjrIsrMET", rjrMet ); 

/*
	// attempt to use RJR to guess p vector of N2
    TLorentzVector x1alab = X1a_c->GetFourVector(*LAB_c);
    TLorentzVector x1blab = X1b_c->GetFourVector(*LAB_c);
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
    selRjrIsrVars.fillBranch( "rjrIsrN2Px", n2px );
    selRjrIsrVars.fillBranch( "rjrIsrN2Py", n2py );
    selRjrIsrVars.fillBranch( "rjrIsrN2Pz", n2pz );
*/

  if( RJRDEBUG ) std::cout << " - Getting RJR varibles finished." << std::endl;

}//<<>>void KUCMSAodSkimmer::processRJR( int type, bool newEvent )

std::vector< TLorentzVector > KUCMSAodSkimmer::BinaryMerge( const std::vector< TLorentzVector > & rjr_jets, int Nmax ){

	if( rjr_jets.size() <= Nmax ) return rjr_jets;

  	const int N = (int)rjr_jets.size();
  	double Rmin = 2.;
  	int imin = -1;
  	int jmin = -1;
  	for( int i = 0; i < N-1; i++ ){
    	for( int j = i+1; j < N; j++ ){
      		const double M2 = ( rjr_jets[i] + rjr_jets[j] ).M2();
			const double sumM = rjr_jets[i].M() + rjr_jets[j].M();
      		const double R = ( M2 > 0. ) ? ( M2 - sumM*sumM ) / M2 : 0.;
      		if( R < Rmin ){ Rmin = R; imin = i; jmin = j; }//<<>>if(R < Rmin)
    	}//<<>>for(int j = i+1; j < N; j++)
  	}//<<>>for(int i = 0; i < N-1; i++)

  	std::vector< TLorentzVector > merged_rjr_jets;
    merged_rjr_jets.reserve(N - 1);
	merged_rjr_jets.push_back( rjr_jets[imin] + rjr_jets[jmin] );
  	for(int i = 0; i < N; i++){ if( i != imin && i != jmin ) merged_rjr_jets.push_back( rjr_jets[i] ); }

  return BinaryMerge( merged_rjr_jets, Nmax );

}//<<>>void KUCMSAodSkimmer::BinaryMerge( std::vector< TLorentzVector > rjr_jets )

void KUCMSAodSkimmer::BinaryMergeInPlace( std::vector<TLorentzVector>& jets, int Nmax ){

    if( jets.size() < 2 ) return;

    auto ptDesc = []( const TLorentzVector& a, const TLorentzVector& b ){ return a.Pt() > b.Pt(); };

    // Ensure initial ordering: lead jet first
    std::sort(jets.begin(), jets.end(), ptDesc);

    
    while( (int)jets.size() > Nmax ){

        const int N = (int)jets.size();
        if( N < Nmax ) break;

		//std::cout << " - Jet Merge for RJR ISR " << std::endl;
  		//for( uInt it = 0; it < N; it++ ){
    	//	TLorentzVector jet = jets[it];
    	//	std::cout << " - Jet " << it << " Pt: " << jet.Pt() << " Eta: " << jet.Eta();
    	//	std::cout << " Phi: " << jet.Phi() << " M: " << jet.M() << std::endl;
  		//}//<<>>for( int i = 0; i < nSelJets; i++ )

        // Default: merge two smallest-pT jets (last two entries)
        int imin = N - 2;
        int jmin = N - 1;

        double Rmin = 100000000.0;
        // Search for best pair with R < 2
        for( int i = 0; i < N - 1; ++i ){
            for( int j = i + 1; j < N; ++j ){

                const double M2 = (jets[i] + jets[j]).M2();
                const double sumM = jets[i].M() + jets[j].M();
                const double R = ( M2 > 0.0 ) ? ( M2 - sumM * sumM ) : 0.0;
                if( R < Rmin ){ Rmin = R; imin = i; jmin = j; }

            }//<<>>for( int j = i + 1; j < N; ++j )
        }//<<>>for( int i = 0; i < N - 1; ++i )

        // Merge jmin into imin
        jets[imin] += jets[jmin];
        // Remove the merged-away jet. Use erase to keep list clean.
        jets.erase(jets.begin() + jmin);
        // Restore pT ordering so lead jet is first for next iteration
        std::sort(jets.begin(), jets.end(), ptDesc);

    }//while( (int)jets.size() > Nmax )

}//<<>>BinaryMergeInPlace(std::vector<TLorentzVector>& jets, int Nmax)

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setRJRISRBranches( TTree* fOutTree ){

	std::cout << " - Making Branches for RJRISR." << std::endl;

    selRjrIsrVars.makeBranch( "rjrIsrX1aMass", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1aCosA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1bMass", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1bCosA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1a_Pt", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1a_Eta", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1a_Phi", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1b_Pt", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1b_Eta", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX1b_Phi", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsrX2aMass", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX2aCosA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX2aPs", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX2bMass", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX2bCosA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrX2bPs", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsrSMass", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrSCosA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrSdphiDA", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrSdphiBV", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrPTS", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrPZS", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrPS", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsrEVa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrEVb", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrPVa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrPVb", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_nSVisObjects", INT );
    selRjrIsrVars.makeBranch( "rjrIsr_nIsrVisObjects", INT );
    selRjrIsrVars.makeBranch( "rjrIsrNPhotons", INT );
    selRjrIsrVars.makeBranch( "rjrIsrMET", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrNJets", INT );
    selRjrIsrVars.makeBranch( "rjrIsrLeadPhoLocation", INT);
    selRjrIsrVars.makeBranch( "rjrIsrSubPhoLocation", INT);
    selRjrIsrVars.makeBranch( "rjrIsrNVisObjects", INT );
    selRjrIsrVars.makeBranch( "rjrIsrABSide", BOOL );

    selRjrIsrVars.makeBranch( "rjrIsrVisPhoSide", INT );
    selRjrIsrVars.makeBranch( "rjrIsrVisJetSide", INT );

    selRjrIsrVars.makeBranch( "rjrIsrPVlab", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrDphiMETV", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrDPhiSI", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsrCleaningVeto0", BOOL );
    selRjrIsrVars.makeBranch( "rjrIsrCleaningVeto1", BOOL );
    selRjrIsrVars.makeBranch( "rjrIsrCleaningVeto2", BOOL );

    selRjrIsrVars.makeBranch( "rjrIsrNISRPhos", INT );
    selRjrIsrVars.makeBranch( "rjrIsrNJetsJb", INT );
    selRjrIsrVars.makeBranch( "rjrIsrNJetsJa", INT );
    selRjrIsrVars.makeBranch( "rjrIsrNJetsISR",INT );

    selRjrIsrVars.makeBranch( "rjrIsr_pHts20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts40", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts41", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts21a", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts21b", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts11a", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts11b", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_pHs20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs40", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs41", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs21a", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs21b", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs11a", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs11b", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_pHxa20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxb20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxa21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxb21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxa10", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxb10", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxa11", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHxb11", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_pHtxa20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxb20", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxa21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxb21", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxa10", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxb10", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxa11", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHtxb11", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_pHs22", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHs11", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts22", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_pHts11", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_p1Pts", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p2Pts", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_p1Ptxa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p1Ptxb", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p2Ptxa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p2Ptxb", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_p1Pxa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p1Pxb", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p2Pxa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_p2Pxb", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsrAX2Diff","rjr_Rdiff", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_Mr", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_R", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rv", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Ms", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rs", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rx", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rm", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rxa", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Rxb", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_sJ1a_px", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1a_py", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1a_pz", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1a_e", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_sJ2a_px", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2a_py", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2a_pz", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2a_e", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_sJ1b_px", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1b_py", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1b_pz", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ1b_e", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_sJ2b_px", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2b_py", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2b_pz", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_sJ2b_e", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_Mva", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_Mvb", FLOAT );

    selRjrIsrVars.makeBranch( "rjrIsr_MsPerp", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_PIsr", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_PtIsr", FLOAT );
    selRjrIsrVars.makeBranch( "rjrIsr_RIsr", FLOAT );

    //selRjrIsrVars.makeBranch( "rjrIsrN2Px", FLOAT );
    //selRjrIsrVars.makeBranch( "rjrIsrN2Py", FLOAT );
    //selRjrIsrVars.makeBranch( "rjrIsrN2Pz", FLOAT );

    selRjrIsrVars.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

