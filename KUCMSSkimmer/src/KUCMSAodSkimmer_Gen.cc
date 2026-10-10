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

#include <algorithm>
#include <cmath>
#include <vector>

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

  struct X2LifetimeInfo {
    int genIndex;
    float energy;
    float eta;
    float phi;
    float mass;
    float displacement;
    float momentum;
    uInt pdgId;
    float pt;
    float vx;
    float vy;
    float vz;
    float beta;
    float ctau;
  };

  //std::cout << "New Event ----------------------------------" << std::endl;
  int nGenParts = Gen_pdgId->size();
  std::vector<X2LifetimeInfo> validX2Lifetimes;
  std::vector<bool> hasValidLifetimeForMother(nGenParts, false);
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

    const bool hasMom = momIndx >= 0 && momIndx < nGenParts;

    //if( status != 1 ) continue;
    if( pdgId > 40 && pdgId < 1000000 ) continue;

    float mommass = hasMom ? (*Gen_mass)[momIndx] : -1;
    uInt mompdg = hasMom ? (*Gen_pdgId)[momIndx] : 0;
    float mompx = hasMom ? (*Gen_px)[momIndx] : -1;
    float mompy = hasMom ? (*Gen_py)[momIndx] : -1;
    float mompz = hasMom ? (*Gen_pz)[momIndx] : -1;
    float genmomp = hypo( mompx, mompy, mompz );
    //float beta = ( mommass > 0 ) ? genmomp/mommass : -1;
    //float gama = ( beta >= 0 ) ? 1/std::sqrt( 1 - beta*beta ) : -1; 
    float gbeta = ( mommass > 0 ) ? genmomp/mommass : -1;
    //float gbeta = ( gama >= 0 && beta >= 0 ) ? gama*beta : -1;
    float ctau = ( ( gbeta > 0 ) && ( displacment >= 0 ) ) ? displacment/gbeta : -1;
    if( mompdg == 1000023 ){ selGenPart.fillBranch( "genXMomCTau", ctau ); }

    // Reconstruct each physical X2 once from a retained decay child.  The
    // ntuplizer's normalized mother index skips same-PDG copies, so a photon,
    // Z, or LSP with an X2 mother identifies the terminal/decaying X2 copy.
    const bool isX2DecayChild = pdgId == 22 || pdgId == 23 || pdgId == 1000022;
    if( hasMom && mompdg == 1000023 && isX2DecayChild &&
        !hasValidLifetimeForMother[momIndx] ){
      const float momenergy = (*Gen_energy)[momIndx];
      const float mometa = (*Gen_eta)[momIndx];
      const float momphi = (*Gen_phi)[momIndx];
      const float mompt = (*Gen_pt)[momIndx];
      const float momvx = (*Gen_vx)[momIndx];
      const float momvy = (*Gen_vy)[momIndx];
      const float momvz = (*Gen_vz)[momIndx];
      const float mombeta = momenergy > 0.f ? genmomp/momenergy : -1.f;
      const float momctau = ( genmomp > 0.f && mommass > 0.f ) ?
        displacment*mommass/genmomp : -1.f;

      const bool validLifetime = std::isfinite(displacment) && displacment > 0.f &&
        std::isfinite(genmomp) && genmomp > 0.f &&
        std::isfinite(mommass) && mommass > 0.f &&
        std::isfinite(momenergy) && momenergy > 0.f &&
        std::isfinite(mometa) && std::isfinite(momphi) && std::isfinite(mompt) &&
        std::isfinite(momvx) && std::isfinite(momvy) && std::isfinite(momvz) &&
        std::isfinite(mombeta) && mombeta > 0.f &&
        std::isfinite(momctau) && momctau > 0.f;

      if( validLifetime ){
        validX2Lifetimes.push_back({
          momIndx, momenergy, mometa, momphi, mommass, displacment, genmomp,
          mompdg, mompt, momvx, momvy, momvz, mombeta, momctau
        });
        hasValidLifetimeForMother[momIndx] = true;
      }
    }
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
    int gMomIndx = hasMom ? Gen_motherIdx->at(momIndx) : -1;
    bool hasGrandMom( gMomIndx >= 0 && gMomIndx < nGenParts );
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
    //if( zee && fromX ) nZfX++;
    if( zee ) nZfX++;
    if( hasX234 ) nX234++;
    if( lept && fromZ && momFromX ) nLZX++;
    if( quark && fromZ && momFromX ) nQfZ++;
    if( N0 && fromSqk ) nN0fsqk++;
    if( N0 && fromSg ) nN0fsg++;

    selGenPart.fillBranch( "genPartEnergy", energy );
    selGenPart.fillBranch( "genPartEta", eta );
    selGenPart.fillBranch( "genPartPhi", phi );
    selGenPart.fillBranch( "genPartPt", pt );
    selGenPart.fillBranch( "genPartPdgId", pdgId );
    selGenPart.fillBranch( "genPartSusId", susId );
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

  std::sort(validX2Lifetimes.begin(), validX2Lifetimes.end(),
    []( const X2LifetimeInfo& lhs, const X2LifetimeInfo& rhs ){
      return lhs.genIndex < rhs.genIndex;
    });
  const int nValidXs = static_cast<int>(validX2Lifetimes.size());

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

  bool evtIsZZ = ( nZfX == 2 ) ? true : false;
  bool evtIsZG = ( nZfX == 1 ) ? true : false;
  bool evtIsGG = ( nZfX == 0 ) ? true : false;
	
  selGenPart.fillBranch( "genSTFlagQQP", isSTqqp );
  selGenPart.fillBranch( "genSTFlagQQZLL", isSTqqzll );
  selGenPart.fillBranch( "genSTFlagQQZL", isSTqqzl );
  selGenPart.fillBranch( "genSTFlagQQZ", isSTqqz );
  selGenPart.fillBranch( "genSTFlagQP", isSTqp );
  selGenPart.fillBranch( "genSigType", Gen_susEvtType->at(0) );

  if( doNewSigBase ){

    const auto fillX2Lifetime = [&]( const X2LifetimeInfo& x2, const std::string& label ){
      selGenPart.fillBranch( label + "_energy", x2.energy );
      selGenPart.fillBranch( label + "_phi", x2.phi );
      selGenPart.fillBranch( label + "_mass", x2.mass );
      selGenPart.fillBranch( label + "_Displacment", x2.displacement );
      selGenPart.fillBranch( label + "_p", x2.momentum );
      selGenPart.fillBranch( label + "_pdgId", x2.pdgId );
      selGenPart.fillBranch( label + "_pt", x2.pt );
      selGenPart.fillBranch( label + "_vx", x2.vx );
      selGenPart.fillBranch( label + "_vy", x2.vy );
      selGenPart.fillBranch( label + "_vz", x2.vz );
      selGenPart.fillBranch( label + "_beta", x2.beta );
      selGenPart.fillBranch( label + "_ctau", x2.ctau );
      geVects.set( label == "Xa" ? "xa5vec" : "xb5vec",
        { x2.vx, x2.vy, x2.vz, x2.beta, x2.displacement, x2.eta, x2.phi } );
    };

    geVects.set( "xa5vec", std::vector<float>(7, -1.f) );
    geVects.set( "xb5vec", std::vector<float>(7, -1.f) );
    if( nValidXs > 0 ) fillX2Lifetime(validX2Lifetimes[0], "Xa");
    if( nValidXs > 1 ) fillX2Lifetime(validX2Lifetimes[1], "Xb");

    //selGenPart.fillBranch( "Evt_isGG", evtIsZZ );
    //selGenPart.fillBranch( "Evt_isGZ", evtIsZG );
    //selGenPart.fillBranch( "Evt_isZZ", evtIsGG );
    selGenPart.fillBranch( "Evt_isGG", Evt_isGG );
    selGenPart.fillBranch( "Evt_isGZ", Evt_isGZ );
    selGenPart.fillBranch( "Evt_isZZ", Evt_isZZ );
    // Preserve the ntuplizer's raw X2-record count as the reconstruction flag.
    selGenPart.fillBranch( "Evt_nXs", Evt_nXs );
    selGenPart.fillBranch( "Evt_nValidXs", nValidXs );

  }//<<>>if( doNewSigBase )

}//<<>>void KUCMSAodSkimmer::processGenParticles()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setGenBranches( TTree* fOutTree ){

  std::cout << " - Making Branches for Gen." << std::endl;
  selGenPart.makeBranch( "genPartEnergy", VFLOAT );
  selGenPart.makeBranch( "genPartEta", VFLOAT );
  selGenPart.makeBranch( "genPartPhi", VFLOAT );
  selGenPart.makeBranch( "genPartPt", VFLOAT );
  selGenPart.makeBranch( "genPartPdgId", VUINT );
  selGenPart.makeBranch( "genPartSusId", VINT );
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

  selGenPart.makeBranch( "Xa_energy", FLOAT );
  selGenPart.makeBranch( "Xa_phi", FLOAT ); // opps -> missed eta : actually did phi twice :(
  selGenPart.makeBranch( "Xa_mass", FLOAT );
  selGenPart.makeBranch( "Xa_Displacment", FLOAT );
  selGenPart.makeBranch( "Xa_p", FLOAT );
  selGenPart.makeBranch( "Xa_pdgId", UINT );
  selGenPart.makeBranch( "Xa_pt", FLOAT );
  selGenPart.makeBranch( "Xa_vx", FLOAT );
  selGenPart.makeBranch( "Xa_vy", FLOAT );
  selGenPart.makeBranch( "Xa_vz", FLOAT );
  selGenPart.makeBranch( "Xa_beta", FLOAT );
  selGenPart.makeBranch( "Xa_ctau", FLOAT );

  selGenPart.makeBranch( "Xb_energy", FLOAT );
  selGenPart.makeBranch( "Xb_phi", FLOAT );
  selGenPart.makeBranch( "Xb_mass", FLOAT );
  selGenPart.makeBranch( "Xb_Displacment", FLOAT );
  selGenPart.makeBranch( "Xb_p", FLOAT );
  selGenPart.makeBranch( "Xb_pdgId", UINT );
  selGenPart.makeBranch( "Xb_pt", FLOAT );
  selGenPart.makeBranch( "Xb_vx", FLOAT );
  selGenPart.makeBranch( "Xb_vy", FLOAT );
  selGenPart.makeBranch( "Xb_vz", FLOAT );
  selGenPart.makeBranch( "Xb_beta", FLOAT );
  selGenPart.makeBranch( "Xb_ctau", FLOAT );

  selGenPart.makeBranch( "Evt_isGG", BOOL );
  selGenPart.makeBranch( "Evt_isGZ", BOOL );
  selGenPart.makeBranch( "Evt_isZZ", BOOL );
  selGenPart.makeBranch( "Evt_nXs", INT );
  selGenPart.makeBranch( "Evt_nValidXs", INT,
    "unique terminal X2 mothers with a positive finite reconstructed lifetime" );

  selGenPart.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )
