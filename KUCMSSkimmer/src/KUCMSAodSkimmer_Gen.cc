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

void KUCMSAodSkimmer::processGenParticles(){

  // initilize
  selGenPart.clearBranches(); // <<<<<<<   must do

  // calc
  if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
  //------------  genparts ------------------------

  int nLSPfXfSg = 0;
  int nQfSg = 0;
  int nQfSqk = 0;
  int nPHOfX = 0;
  int nZfX = 0;
  int nSGlue = 0;
  int nSQuark = 0;
  int nX234 = 0;
  int nLZX = 0;
  int nQfZ = 0;
  int nN0fsqk = 0;
  int nN0fsg = 0;

  //std::cout << "New Event ----------------------------------" << std::endl;
  int nGenParts = Gen_pdgId->size();
  for( int it = 0; it < nGenParts; it++ ){

    float displacment = (*Gen_momDisplacment)[it];
    float energy = (*Gen_energy)[it];
    float eta = (*Gen_eta)[it];
    float phi = (*Gen_phi)[it];
    float pt = (*Gen_pt)[it];
    uInt  pdgId = (*Gen_pdgId)[it];
    int momIndx = (*Gen_motherIdx)[it];
    int   susId = (*Gen_susId)[it];
    float charge = (*Gen_charge)[it];
    float mass = (*Gen_mass)[it];
    int status = (*Gen_status)[it];
    float vx = (*Gen_vx)[it];
    float vy = (*Gen_vy)[it];
    float vz = (*Gen_vz)[it];
    float px = (*Gen_px)[it];
    float py = (*Gen_py)[it];
    float pz = (*Gen_pz)[it];

    //if( status != 1 ) continue;
    if( pdgId > 40 && pdgId < 1000000 ) continue;

    float mommass = ( momIndx > -1 ) ? (*Gen_mass)[momIndx] : -1;
    uInt mompdg = ( momIndx > -1 ) ? (*Gen_pdgId)[momIndx] : 0;
    float mompx = ( momIndx > -1 ) ? (*Gen_px)[momIndx] : -1;
    float mompy = ( momIndx > -1 ) ? (*Gen_py)[momIndx] : -1;
    float mompz = ( momIndx > -1 ) ? (*Gen_pz)[momIndx] : -1;
    float genmomp = hypo( mompx, mompy, mompz );
    //float beta = ( mommass > 0 ) ? genmomp/mommass : -1;
    //float gama = ( beta >= 0 ) ? 1/std::sqrt( 1 - beta*beta ) : -1; 
    float gbeta = ( mommass > 0 ) ? genmomp/mommass : -1;
    //float gbeta = ( gama >= 0 && beta >= 0 ) ? gama*beta : -1;
    float ctau = ( ( gbeta >= 0 ) && ( displacment >= 0 ) ) ? displacment/gbeta : -1;
    if( mompdg == 1000023 ){ selGenPart.fillBranch( "genXMomCTau", ctau ); }
    //if( mompdg != 0 ){
    if( false ){
      std::cout << " ctau for : " << pdgId << " mother: " << mompdg << " with mommass " << mommass;
      std::cout << " genp: " << genmomp << " gbeta: " << gbeta << " dis: " << displacment;
      std::cout << " ctau: " << ctau << std::endl; 
    }//<<>>if( mompdg != 0 )

    if( ( pdgId > 1000000 ) && ( pdgId < 1000007 ) ){ nSQuark++; selGenPart.fillBranch( "genSQMass", mass ); }
    if( pdgId == 1000021 ){ nSGlue++; selGenPart.fillBranch( "genSGMass", mass ); }
    if( pdgId == 1000023 ){ selGenPart.fillBranch( "genLLPMass", mass ); }
    if( pdgId == 1000022 ){ selGenPart.fillBranch( "genLSPMass", mass ); }
    if( pdgId == 1000039 ){ selGenPart.fillBranch( "genGrvtinoMass", mass ); }

    //if( pdgId < 7 ) continue;
    bool hasMom( ( momIndx < nGenParts ) && ( momIndx > -1 ) );
    int gMomIndx = hasMom ? Gen_motherIdx->at(momIndx) : -1;
    bool hasGrandMom( gMomIndx > -1 );
    bool lsp( ( pdgId == 1000039 ) || ( ( pdgId == 1000022 ) && ( status == 1 ) ) );
    bool hasX234( ( pdgId > 1000022 ) && ( pdgId < 1000038 ) );
    bool fromX( hasMom && ( ( Gen_pdgId->at(momIndx) == 1000022 ) || ( Gen_pdgId->at(momIndx) == 1000023 ) ) );
    bool fromSg( hasMom && ( Gen_pdgId->at(momIndx) == 1000021 ) );
    bool fromSqkL( hasMom && ( Gen_pdgId->at(momIndx) > 1000000 ) && ( Gen_pdgId->at(momIndx) < 1000009 ) );
    bool fromSqkR( hasMom && ( Gen_pdgId->at(momIndx) > 2000000 ) && ( Gen_pdgId->at(momIndx) < 2000009 ) );
    bool fromSqk( fromSqkL || fromSqkR );
    bool fromZ( hasMom && ( Gen_pdgId->at(momIndx) == 23 ) );
    bool momFromX( hasGrandMom && ( ( Gen_pdgId->at(gMomIndx) == 1000022 ) || ( Gen_pdgId->at(gMomIndx) == 1000023 ) ) );
    bool momFromSg( hasGrandMom && ( Gen_pdgId->at(gMomIndx) == 1000021 ) );
    bool momFromSqkL( hasGrandMom && ( Gen_pdgId->at(gMomIndx) > 1000000 ) && ( Gen_pdgId->at(gMomIndx) < 1000009 ) );
    bool momFromSqkR( hasGrandMom && ( Gen_pdgId->at(gMomIndx) > 2000000 ) && ( Gen_pdgId->at(gMomIndx) < 2000009 ) );
    bool momFromSqk( momFromSqkL || momFromSqkR );
    bool quark( pdgId < 9 );
    bool photon( pdgId == 22 ); 
    bool zee( pdgId == 23 );
    bool lept( ( pdgId > 10 ) && ( pdgId < 19 ) );
    bool N0( pdgId == 1000022 );

    //if( lsp && fromX && momFromSg ) nLSPfXfSg++;
    //if( lsp && fromX ) nLSPfXfSg++;
    if( lsp ) nLSPfXfSg++;
    if( quark && fromSg ) nQfSg++;
    if( quark && fromSqk ) nQfSqk++;
    if( photon && fromX ) nPHOfX++;
    if( zee && fromX ) nZfX++;
    if( hasX234 ) nX234++;
    if( lept && fromZ && momFromX ) nLZX++;
    if( quark && fromSqk ) nQfSqk++;
    if( quark && fromZ && momFromX ) nQfZ++;
    if( N0 && fromSqk ) nN0fsqk++;
    if( N0 && fromSg ) nN0fsg++;

    selGenPart.fillBranch( "genPartEnergy", energy );
    selGenPart.fillBranch( "genPartEta", eta );
    selGenPart.fillBranch( "genPartPhi", phi );
    selGenPart.fillBranch( "genPartPt", pt );
    selGenPart.fillBranch( "genPartPdgId", pdgId );
    selGenPart.fillBranch( "genPartSusId", susId );
    selGenPart.fillBranch( "genPartMomIdx", momIndx );
    selGenPart.fillBranch( "genMomCTau", ctau );
    selGenPart.fillBranch( "genCharge", charge );
    selGenPart.fillBranch( "genMass", mass );
    selGenPart.fillBranch( "genStatus", status );
    selGenPart.fillBranch( "genVx", vx );
    selGenPart.fillBranch( "genVy", vy );
    selGenPart.fillBranch( "genVz", vz );
    selGenPart.fillBranch( "genPx", px );
    selGenPart.fillBranch( "genPy", py );
    selGenPart.fillBranch( "genPz", pz );

  }//<<>>for( int it = 0; it < nGenParts; it++ )

  //bool hasLSP( nLSPfXfSg > 0 );
  bool hasLSP( true );
  bool noX234( nX234 == 0 );

  bool has2Nfsqk( nN0fsqk == 2 );
  bool has2Nfsg( nN0fsg == 2 );
  bool has4QfSg( nQfSg > 3 );
  bool has2QfSqk( nQfSqk > 1 );

  bool has2PfN( nPHOfX == 2 );
  bool has2ZfN( nZfX == 2 );

  bool hasLLfz( nLZX == 4 );
  bool hasLfz( nLZX == 2 );
  bool hasNoLfz( nLZX == 0 );

  //if( nZfX > 1 && nLZX > 0 && nX234 == 0 ){
  //if( nX234 == 0 ){
  //if( nQfSg || nQfSqk ){ 
  if( false ){
    std::cout << "GenEventType::";
    std::cout << " nPHOfX: " << nPHOfX << " nZfX: " << nZfX;
    std::cout << " nX234: " << nX234 << " nN0fsqk: " << nN0fsqk << " nN0fsg: " << nN0fsg;
    std::cout << " nLZX: " << nLZX  << " nQfZ: " << nQfZ; 
    std::cout << " nQfSg: " << nQfSg << " nQfSqk: " << nQfSqk;
    std::cout << std::endl;
  }//<<>>if(

  bool isSTqp( has2Nfsqk && has2QfSqk && has2PfN && noX234 );
  bool isSTqqp( has2Nfsg && has4QfSg && has2PfN && noX234 );
  bool isSTqqzll( has2Nfsg && has4QfSg && has2ZfN && hasLLfz && noX234 ); 
  bool isSTqqzl( has2Nfsg && has4QfSg && has2ZfN && hasLfz && noX234 );
  bool isSTqqz( has2Nfsg && has4QfSg && has2ZfN && hasNoLfz && noX234 );	

  //if( nZfX > 1 && nLZX > 0 && nX234 == 0 ){
  //if( nX234 == 0 ){ 
  //if( isSTqp || isSTqqp || isSTqqzll || isSTqqzl || isSTqqz ){
  //  std::cout << "GenSTFlag:";
  //  std::cout << " isSTqqp: " << isSTqqp << " isSTqqzll: " << isSTqqzll << " isSTqqzl: " << isSTqqzl;
  //  std::cout << " isSTqqz: " << isSTqqz << " isSTqp: " << isSTqp;
  //  std::cout << std::endl;
  //}//<<>>if(

	
  selGenPart.fillBranch( "genSTFlagQQP", isSTqqp );
  selGenPart.fillBranch( "genSTFlagQQZLL", isSTqqzll );
  selGenPart.fillBranch( "genSTFlagQQZL", isSTqqzl );
  selGenPart.fillBranch( "genSTFlagQQZ", isSTqqz );
  selGenPart.fillBranch( "genSTFlagQP", isSTqp );
  selGenPart.fillBranch( "genSigType", Gen_susEvtType->at(0) );

}//<<>>void KUCMSAodSkimmer::processGenParticles()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setGenBranches( TTree* fOutTree ){

  selGenPart.makeBranch( "genPartEnergy", VFLOAT );
  selGenPart.makeBranch( "genPartEta", VFLOAT );
  selGenPart.makeBranch( "genPartPhi", VFLOAT );
  selGenPart.makeBranch( "genPartPt", VFLOAT );
  selGenPart.makeBranch( "genPartPdgId", VUINT );
  selGenPart.makeBranch( "genPartSusId", VINT );
  selGenPart.makeBranch( "genPartMomIdx", VINT );
  selGenPart.makeBranch( "genXMomCTau", VFLOAT );
  selGenPart.makeBranch( "genMomCTau", VFLOAT );
  selGenPart.makeBranch( "genCharge", VINT );   //!
  selGenPart.makeBranch( "genMass", VFLOAT );   //!   
  selGenPart.makeBranch( "genStatus", VINT );   //!
  selGenPart.makeBranch( "genVx", VFLOAT );   //!
  selGenPart.makeBranch( "genVy", VFLOAT );   //!
  selGenPart.makeBranch( "genVz", VFLOAT );   //!
  selGenPart.makeBranch( "genPx", VFLOAT );   //!
  selGenPart.makeBranch( "genPy", VFLOAT );   //!
  selGenPart.makeBranch( "genPz", VFLOAT );   //!
  selGenPart.makeBranch( "genSigType", INT );   //!
  selGenPart.makeBranch( "genSTFlagQQP", BOOL );   //!
  selGenPart.makeBranch( "genSTFlagQQZLL", BOOL );   //!
  selGenPart.makeBranch( "genSTFlagQQZL", BOOL );   //!
  selGenPart.makeBranch( "genSTFlagQQZ", BOOL );   //!
  selGenPart.makeBranch( "genSTFlagQP", BOOL );   //!
  selGenPart.makeBranch( "genSQMass", VFLOAT );   //! 
  selGenPart.makeBranch( "genSGMass", VFLOAT );   //! 
  selGenPart.makeBranch( "genLLPMass", VFLOAT );   //! 
  selGenPart.makeBranch( "genLSPMass", VFLOAT );   //! 
  selGenPart.makeBranch( "genGrvtinoMass", VFLOAT );   //! 

  selGenPart.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

