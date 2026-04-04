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

void KUCMSAodSkimmer::processJets(){

  // intilize
  selJets.clearBranches(); // <<<<<<<   must do

  // calc
  //if( DEBUG ) std::cout << "Finding jets" << std::endl;
  //--------- jets --------------------------
  isSelJet.clear(); 

  int nSelJetsWithWPixTime = 0;
  int nJetsWithWPixTime = 0;
  int nSelJetsWithWTime = 0;
  int nJetsWithWTime = 0;
  int jetEventVeto = 0;
  uInt nSelJets = 0;
  int nQJets = 0;
  uInt nJets = Jet_energy->size();
  std::vector<float> seljetpt;
  std::vector<float> seljeteta;
  std::vector<float> seljetphi;
  std::vector<float> seljetmass;
  std::vector<float> alljetwtimevar;
  std::vector<float> alljetwtenergy;
  std::vector<float> alljetpixwtimevar;
  std::vector<float> alljetpixwtenergy;

  uInt nSCBase = SuperCluster_excluded->size();
  std::vector<bool> usedsc( nSCBase, false );
  //std::vector<bool> scskip( nSCBase, false );
  std::vector<bool> sckeep(nSCBase, false);
  std::vector<uInt> kept;
  kept.reserve(nSCBase);
  std::vector<uInt> order(nSCBase);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](uInt a, uInt b){ float ea = (*SuperCluster_energy)[a]; float eb = (*SuperCluster_energy)[b];
        	if( ea != eb ) return ea > eb; return a < b; } );

  //for( uInt pit = 0; pit < nSCBase; pit++ ){ scskip.push_back(false); usedsc.push_back(false);}
  for( uInt idx : order ){
	bool overlaps_kept = false;
    for( uInt jdx : kept ){
        if( getOverLapCnt((*SuperCluster_rhIds)[idx], (*SuperCluster_rhIds)[jdx]) > 0 ){
            overlaps_kept = true;
            break;
        }//<<>>if( getOverLapCnt((*SuperCluster_rhIds)[idx], (*SuperCluster_rhIds)[jdx]) > 0 )
    }//<<>>for( uInt jdx : kept )
    if( !overlaps_kept ){ sckeep[idx] = true; kept.push_back(idx);}
  }//<<>>for( uInt idx : order )

/*
  for( uInt pit = 0; pit < nSCBase; pit++ ){
	if( scskip[pit] ) continue;
    float peng = (*SuperCluster_energy)[pit];
	for( uInt olit = pit + 1; olit < nSCBase; ++olit ){
		if( scskip[olit] ) continue;
        float oeng = (*SuperCluster_energy)[olit];
        int nol = getOverLapCnt( (*SuperCluster_rhIds)[pit], (*SuperCluster_rhIds)[olit] );
		//int oid = (*SuperCluster_ObjectPdgId)[olit];
        //if( nol > 0 ){ std::cout << " -- ! SC Overlap w/ " << olit << " of " << nol << " e " << oeng << " t " << oid << std::endl; }
        if( nol > 0 ){ if(  peng > oeng ) scskip[olit] = true; else scskip[pit] = true; }
    }//<<>>for( uInt pitt = 0; pitt < nPhotons; pitt++ )
  }//<<>>for( uInt pit = 0; pit < nSCBase; pit++ )
*/

  alljetwtime.clear();
  if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets" << std::endl;
  for( uInt it = 0; it < nJets; it++ ){

    // pull values ---------------------------------------------------
    float energy = (*Jet_energy)[it];
    float mass = (*Jet_mass)[it];
    //if( DEBUG ) std::cout << " - Finding Jet Quality" << std::endl;
    int quality = getJetQuality(it);
    float pt = (*Jet_pt)[it];
    float eta = (*Jet_eta)[it];
    float phi = (*Jet_phi)[it];
    auto rhids = (*Jet_drRhIds)[it];
    if( DEBUG ) std::cout << " -- jettDrId size: " << rhids.size() << std::endl;		

    float area = (*Jet_area)[it];
    float chEmEF = (*Jet_chEmEF)[it];
    float chHEF = (*Jet_chHEF)[it];
    float chHM = (*Jet_chHM)[it];
    float muEF = (*Jet_muEF)[it];
    float neEmEF = (*Jet_neEmEF)[it];
    float neHEF = (*Jet_neHEF)[it];
    float neHM = (*Jet_neHM)[it];

  	bool nhfup( neHEF  <= 0.90 );   // delpho
  	bool nemfup2( neEmEF <= 0.99 );
  	bool shm1( ( neHM + chHM ) >= 1 );
  	bool mufllp( muEF > 0.80 );
 	bool chf0( chHEF  >= 0 );       // delpho
  	bool chm0( chHM >= 0 );       // delpho
  	bool cemfup1( chEmEF <= 0.99 );
  	bool isLLPJet = nhfup && nemfup2 && shm1 && chf0 && chm0 && ( cemfup1 || mufllp );
	//if( quality == 4 && isLLPJet == false ) std::cout << "!!!!!!!!!!!!!!!!!!! bad jet llp qualty match !!!!!!!!!!!!!!!!!!!!!" << std::endl;

    int qrkllpId = 0; //*Jet_genLlpId)[it];
    int gjllpId = 0; //*Jet_genLlpId)[it];
    float jgdpt = 0; //*Jet_genDptMatch)[it];
    float jgdr = 0; //*Jet_genDrMatch)[it];
    float jge = 0; //*Jet_genEnergy)[it];
    float jgeta = 0; //*Jet_genEta)[it];
    float jgimpang = 0; //*Jet_genImpactAngle)[it];
    float llpdp = 0; //*Jet_genLlpDp)[it];
    float llpdr = 0; //*Jet_genLlpDr)[it];
    float jgpt = 0; //*Jet_genPt)[it];
    float jgtof = 0; //*Jet_genTOF)[it];
    float jgt = 0; //(*Jet_genTime)[it];
    float jgllpt = 0; //*Jet_genTimeLLP)[it];
    if( hasGenInfoFlag ){
      //qrkllpId = int((*Jet_genQrkLlpId)[it]);
      gjllpId = int((*Jet_genJetLlpId)[it]);
      jgdpt = (*Jet_genDptMatch)[it];
      jgdr = (*Jet_genDrMatch)[it];
      jge = (*Jet_genEnergy)[it];
      jgeta = (*Jet_genEta)[it];
      jgimpang = (*Jet_genImpactAngle)[it];
      //llpdp = (*Jet_genQrkLlpDp)[it];
      //llpdr = (*Jet_genQrkLlpDr)[it];
      jgpt = (*Jet_genPt)[it];
      jgtof = (*Jet_genTOF)[it];
      jgt = (*Jet_genTime)[it];
      jgllpt = (*Jet_genTimeLLP)[it];
    }//if( doGenInfo )

    // get cacluated values-------------------------------------------
    //if( DEBUG ) 
	//std::cout << " - Finding Jet Time." << std::endl;


    auto rhenergies = getRhGrpEnergies( rhids );
    auto rhtimes = getRhGrpTimes( rhids );
    auto timedist = getDistStats( rhtimes, rhenergies );
    auto time = timedist[6];

    float sumwte = 0;
    float sumpixwte = 0;
	std::vector<int> scIndexs;
    std::vector<int> scPixIndexs;
	std::vector<float> eledelay;
    //uInt nSCBase = Photon_excluded->size();
    //uInt nSCBase = SuperCluster_excluded->size();
	//std::vector<bool> scskip;
	//for( uInt pit = 0; pit < nSCBase; pit++ ){ scskip.push_back(false); }
	int nPhoInJet = ( nSCBase > 0 ) ? 0 : -1;
    int nPhoInPixJet = ( nSCBase > 0 ) ? 0 : -1;
    for( uInt pit = 0; pit < nSCBase; pit++ ){ 

		//if( scskip[pit] ) continue;
        if( not sckeep[pit] ) continue;
        if( usedsc[pit] ) continue;

        if( (*SuperCluster_isOot)[pit] ) continue;
        if( (*SuperCluster_excluded)[pit] ) continue;
		int phoindx = (*SuperCluster_PhotonIndx)[pit];
		if( phoindx > -1 && isBaseLinePho[phoindx] ) continue;

        float peta = (*SuperCluster_eta)[pit];
        float pphi = (*SuperCluster_phi)[pit];
        float peng = (*SuperCluster_energy)[pit];
        int pid = (*SuperCluster_ObjectPdgId)[pit];
        auto trkindxs = (*SuperCluster_TrackIndx)[pit];
        int trkindx = ( trkindxs.size() > 0 ) ? trkindxs[0] : -1;

		// skip "soft" leptons with tracks
        if( pid == 0 && trkindx > -1 && peng < 20*std::cosh(peta) ) continue;
		//std::cout << " -- Scanning : e " << peng << " eta " << peta << " phi " << pphi << " : " << pit << std::endl; 
		float dr = dR1(peta, pphi, eta, phi);
      	if( dr < 0.4 ){ 

			usedsc[pit] = true;
			//std::cout << " -- In Jet !!! t " << pid << std::endl;
            scPixIndexs.push_back( pit );
            nPhoInPixJet++; 
            sumpixwte += peng; 
            if( pid > -1 ){
				if( pid == 11 ){ // correction for curved path w/ pt estimation
					float dx = (*SuperCluster_clcx)[pit] - PV_x;
					float dy = (*SuperCluster_clcy)[pit] - PV_y;
					float dz = (*SuperCluster_clcz)[pit] - PV_z;
					int eleindx = (*SuperCluster_ElectronIndx)[pit]; 
					float rad = (*Electron_pt)[eleindx]/1.14;
					float dtan = hypo( dx, dy );
					float arg = dtan/(2.0*rad);
				    arg = std::max(-1.0f, std::min(1.0f, arg));
					float dP = 2.0*std::asin(arg);
					float length = hypo( rad*dP, dz );
					float beta = (*Electron_pt)[eleindx]/(*Electron_energy)[eleindx];
					float delay = ( hypo( dx, dy, dz ) - length )/(beta*SOL);
					delay = std::min( 0.0f, delay );
					eledelay.push_back( delay );
				} else { eledelay.push_back( 0 ); }//<<>>if( pid == 11 ) 
				//scIndexs.push_back( (*Photon_scIndex)[pit] ); 
				scIndexs.push_back( pit );
				nPhoInJet++; 
				sumwte += peng; 
			}//<<>>if( (*Photon_pixelSeed)[pit] )
	
		}//<<>>if( minDr )

	}//<<>>for(  uInt pit = 0; pit < nSCBase; pit++  )

	int nSCIndexs = scIndexs.size();
	float jetTime = -40;
	float jetTimeRes = 10;
	float jetTimeSig = ( nSCIndexs > 0 ) ? getTimeSig( scIndexs, jetTime, jetTimeRes, eledelay ) : -40; 
	alljetwtime.push_back( jetTime );
    alljetwtimevar.push_back( jetTimeRes*jetTimeRes );
    alljetwtenergy.push_back( sumwte );

    int nPixSCIndexs = scPixIndexs.size();
    float jetPixTime = -40;
    float jetPixTimeRes = 10;
    float jetPixTimeSig = ( nPixSCIndexs > 0 ) ? getTimeSig( scPixIndexs, jetPixTime, jetPixTimeRes ) : -40;
    alljetpixwtime.push_back( jetPixTime );
    alljetpixwtimevar.push_back( jetPixTimeRes*jetPixTimeRes );
    alljetpixwtenergy.push_back( sumpixwte );

    if( DEBUG ) std::cout << " - Jet Obj selection." << std::endl;
    // jet object selection ------------------------------------------
    bool overMinPt( pt > 20 );
    bool overMinQaulityPt( pt > 30 ); 
    ////bool underMaxEta( std::abs(eta) <= 3.2 );
    bool underMaxEta( std::abs(eta) <= 2.4 );
    bool isMinQuality( quality > 1 ); // 2 = "tight" 3 = "tighter"
    //bool isMinQuality( quality == 4 );
    bool isNotPhoJet = not phoJetVeto[it];
    //if( isMinQuality != isLLPJet ) std::cout << "!!!!!!!!!!!!!!!!!!! bad jet llp qualty match !!!!!!!!!!!!!!!!!!!!!" << std::endl;

    bool hemEligible1( pt > 10.0 );
    bool hemEligible2( pt > 20.0 );
	bool isInHemRegion = inHEMRegion( eta, phi );
    hemBits.set( "jet1hvl", isInHemRegion && hemEligible1 );
    hemBits.set( "jet2hvm", isInHemRegion && hemEligible2 );
    //if( hemEligible && inHEMRegion( eta, phi ) ) hasHemObj = true;

    std::vector<uInt> allJetSCIndexs;
	for( auto index : scIndexs ){ uInt idx = index;  allJetSCIndexs.push_back( idx ); }

	if( isMinQuality ){

    	selJets.fillBranch( "allJetPt", pt);
    	selJets.fillBranch( "allJetMass", mass);
    	selJets.fillBranch( "allJetEnergy", energy);
    	selJets.fillBranch( "allJetEta", eta);
    	selJets.fillBranch( "allJetPhi", phi);
        selJets.fillBranch( "allJetNPho", nPhoInJet );
        selJets.fillBranch( "allJetWTime", jetTime );
        selJets.fillBranch( "allJetWTimeSig", jetTimeSig );
        selJets.fillBranch( "allJetWTimeRes", jetTimeRes );
        selJets.fillBranch( "allJetPixNPho", nPhoInPixJet );
        selJets.fillBranch( "allJetSCIndexs", allJetSCIndexs );
        selJets.fillBranch( "allJetPixWTime", jetPixTime );
        selJets.fillBranch( "allJetPixWTimeSig", jetPixTimeSig );
        selJets.fillBranch( "allJetPixWTimeRes", jetPixTimeRes );
		if( jetTimeRes > 0 ) nJetsWithWTime++;
		if( jetPixTimeRes > 0 ) nJetsWithWPixTime++;

	}//<<>>if( isMinQuality ){

	bool inAcceptance = underMaxEta && overMinPt && isNotPhoJet;
	if( inAcceptance && not isLLPJet ) jetEventVeto++;

    bool jetQualityGood = underMaxEta && isMinQuality && overMinQaulityPt && isNotPhoJet;
	bool jetSelected = underMaxEta && isLLPJet && overMinPt && isNotPhoJet;

	if( jetQualityGood ){ nQJets++; isSelJet.push_back(true); }
	else isSelJet.push_back(false);

	if( not jetSelected ) continue;
    nSelJets++;

    seljetpt.push_back(pt);
    seljeteta.push_back(eta);
    seljetphi.push_back(phi);
    seljetmass.push_back(mass);

    if( jetTimeRes > 0 ) nSelJetsWithWTime++;
    if( jetPixTimeRes > 0 ) nSelJetsWithWPixTime++;

    // fill vectors
    selJets.fillBranch( "selJetGenLlpId", gjllpId );
    selJets.fillBranch( "selJetQrkLlpId", qrkllpId );
    selJets.fillBranch( "selJetQuality", quality );
    selJets.fillBranch( "selJetPt", pt);
    selJets.fillBranch( "selJetMass", mass);
    selJets.fillBranch( "selJetEnergy", energy);
    selJets.fillBranch( "selJetEta", eta);
    selJets.fillBranch( "selJetPhi", phi);
    selJets.fillBranch( "selJetTime", time);
    selJets.fillBranch( "selJetWTime", jetTime );
    selJets.fillBranch( "selJetWTimeSig", jetTimeSig );
    selJets.fillBranch( "selJetWTimeRes", jetTimeRes );
    selJets.fillBranch( "selJetPixWTime", jetPixTime );
    selJets.fillBranch( "selJetPixWTimeSig", jetPixTimeSig );
    selJets.fillBranch( "selJetPixWTimeRes", jetPixTimeRes );

    selJets.fillBranch( "selJetArea", area ); //*Jet_area)[it]; 
    selJets.fillBranch( "selJetChEmEF", chEmEF ); //*Jet_chEmEF)[it]; 
    selJets.fillBranch( "selJetchHEF", chHEF ); //*Jet_chHEF)[it]; 
    selJets.fillBranch( "selJetChHM", chHM ); //*Jet_chHM)[it]; 
    selJets.fillBranch( "selJetMuEF", muEF ); //*Jet_muEF)[it]; 
    selJets.fillBranch( "selJetNeEmEF", neEmEF ); //*Jet_neEmEF)[it]; 
    selJets.fillBranch( "selJetNeHEF", neHEF ); //*Jet_neHEF)[it]; 
    selJets.fillBranch( "selJetNeHM", neHM ); //*Jet_neHM)[it]; 

    selJets.fillBranch( "selGenJetDpt", jgdpt ); //*Jet_genDptMatch)[it]; 
    selJets.fillBranch( "selGenJetdr", jgdr ); //*Jet_genDrMatch)[it]; 
    selJets.fillBranch( "selGenJetEnergy", jge ); //*Jet_genEnergy)[it]; 
    selJets.fillBranch( "selGenJeteta", jgeta ); //*Jet_genEta)[it]; 
    selJets.fillBranch( "selGenJetImpAng", jgimpang ); //*Jet_genImpactAngle)[it]; 
    selJets.fillBranch( "selJetLlpDp", llpdp ); //*Jet_genLlpDp)[it]; 
    selJets.fillBranch( "selJetLlpDr", llpdr ); //*Jet_genLlpDr)[it]; 
    selJets.fillBranch( "selGenJetPt", jgpt ); //*Jet_genPt)[it]; 
    selJets.fillBranch( "selGenJetTof", jgtof ); //*Jet_genTOF)[it]; 
    selJets.fillBranch( "selGenJetTime", jgt ); // (*Jet_genTime)[it]; 
    selJets.fillBranch( "selGenJetLlpTime", jgllpt ); //*Jet_genTimeLLP)[it]; 

  }//<<>>for( int it = 0; it < nJets; it++ )
  if( DEBUG ) std::cout << " - Finished Jet loop." << std::endl;
  if( nJets == 0 ){ hemBits.set( "jet1hvl", false ); hemBits.set( "jet2hvm", false ); }

  geVects.set( "selJetPt", seljetpt );
  geVects.set( "selJetEta", seljeteta );
  geVects.set( "selJetPhi", seljetphi );
  geVects.set( "selJetMass", seljetmass );

  geCnts.set( "jetEventVeto", jetEventVeto ); 
  geCnts.set( "nQJets", nQJets );
  geCnts.set( "nSelJets", nSelJets );

  selJets.fillBranch( "nJets", nJets);
  selJets.fillBranch( "nSelJets", nSelJets );
  selJets.fillBranch( "allJetsWithWTime", nJetsWithWTime );
  selJets.fillBranch( "allJetsWithWPixTime", nJetsWithWPixTime );
  selJets.fillBranch( "nSelJetsWithWTime", nSelJetsWithWTime );
  selJets.fillBranch( "nSelJetsWithWPixTime", nSelJetsWithWPixTime );

  float sumt = 0;
  float sumv = 0;
  float nsumed = 0;
  float sume = 0;
  int nAllJets = alljetwtime.size();
  for( int ijet = 0; ijet < nAllJets; ijet++ ){

    if( std::abs(alljetwtime[ijet]) > 20.0 ) continue;
    float invres = 1/alljetwtimevar[ijet];
    sumt += invres*alljetwtime[ijet];
    sume += alljetwtenergy[ijet];
    sumv += invres;
    nsumed++;

  }//<<>>for( ijet = 0; ijet < nAllJets; ijet++ )

  int hasPVTime = 1;
  float pvTRes = 0;
  if( nsumed == 0 ){
	sumv = 1; 
	sumt = -30.0; 
	sume = -1; 	
	hasPVTime = 0;
  }//<<>>if( nsumed == 0 )
  else pvTRes = std::sqrt(1/sumv);
  float pvTime = sumt/sumv;
  float pvPEnergy = sume;
  float pvHTVar = std::sqrt(nsumed/sumv);

  selJets.fillBranch( "pv_wtime", pvTime );
  selJets.fillBranch( "pv_wtimeres", pvTRes );
  selJets.fillBranch( "pv_wpenergy", pvPEnergy );
  selJets.fillBranch( "pv_wphartres", pvHTVar );
  selJets.fillBranch( "pv_haspvtime", hasPVTime );

  GetDiJetsSel();

  float gammatime = ( gammaJetIndex[0] > -1 ) ? allphowtime[gammaJetIndex[0]] : -50; 
  bool gammabaseline = ( gammaJetIndex[0] > -1 ) ? allphoBaseline[gammaJetIndex[0]] : false;
  float gammaminjetdr = ( gammaJetIndex[0] > -1 ) ? allphominjetdr[gammaJetIndex[0]] : -1;
  if( gammaJetIndex[0] > -1 && gammaminjetdr < 0.4 ) gammatime = -40;
  float dijet0time = ( diJetIndex[0] > -1 ) ? alljetwtime[diJetIndex[0]] : -50;
  float dijet1time = ( diJetIndex[1] > -1 ) ? alljetwtime[diJetIndex[1]] : -50;

  selJets.fillBranch( "pv_gjGammaTime", gammatime );
  selJets.fillBranch( "pv_gjJetTime", float(gammaJetIndex[1]) );
  float digjettime = ( gammatime > -20 && pvTime > -20 ) ? gammatime - pvTime : -50;
  selJets.fillBranch( "pv_dGJTime", digjettime );
  selJets.fillBranch( "pv_dGJPEnergy", pvPEnergy );
  selJets.fillBranch( "pv_dGJMinDr", gammaminjetdr );
  selJets.fillBranch( "pv_gjPhoIndex", gammaJetIndex[0] );

  float totenergy = ( diJetIndex[0] > -1 && diJetIndex[1] > -1 ) ? alljetwtenergy[diJetIndex[0]] + alljetwtenergy[diJetIndex[1]] : -1;
  selJets.fillBranch( "pv_diJet1Time", dijet0time );
  selJets.fillBranch( "pv_diJet2Time", dijet1time );
  float dijettime = ( dijet0time > -20 && dijet1time > -20 ) ? dijet0time - dijet1time : -50;
  selJets.fillBranch( "pv_dDiJetTime", dijettime );
  selJets.fillBranch( "pv_dDiJetPEnergy", totenergy );

}//<<>>void KUCMSAodSkimmer::processJets()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setJetsBranches( TTree* fOutTree ){

	std::cout << " - Making Branches for Jets." << std::endl;
    //selJets.makeBranch( "JetHt", &JetHt );
    selJets.makeBranch( "nJets", UINT );
    selJets.makeBranch( "nSelJets", UINT );
    selJets.makeBranch( "selJetGenLlpId", VINT ); 
    selJets.makeBranch( "selJetQrkLlpId", VINT );
    selJets.makeBranch( "selJetQuality", VINT ); 
    selJets.makeBranch( "selJetPt", VFLOAT ); 
    selJets.makeBranch( "selJetEnergy", VFLOAT );
    selJets.makeBranch( "selJetEta", VFLOAT ); 
    selJets.makeBranch( "selJetPhi", VFLOAT );
    selJets.makeBranch( "selJetTime", VFLOAT ); 
    selJets.makeBranch( "selJetMass", VFLOAT );
    selJets.makeBranch( "selJetWTime", VFLOAT );
    selJets.makeBranch( "selJetWTimeSig", VFLOAT );
    selJets.makeBranch( "selJetWTimeRes", VFLOAT );
    selJets.makeBranch( "selJetPixWTime", VFLOAT );
    selJets.makeBranch( "selJetPixWTimeSig", VFLOAT );
    selJets.makeBranch( "selJetPixWTimeRes", VFLOAT );
    selJets.makeBranch( "nSelJetsWithWTime", VINT );
    selJets.makeBranch( "nSelJetsWithWPixTime", VINT );

    selJets.makeBranch( "selJetArea", VFLOAT ); //*Jet_area)[it]; 
    selJets.makeBranch( "selJetChEmEF", VFLOAT ); //*Jet_chEmEF)[it]; 
    selJets.makeBranch( "selJetchHEF", VFLOAT ); //*Jet_chHEF)[it]; 
    selJets.makeBranch( "selJetChHM", VFLOAT ); //*Jet_chHM)[it]; 
    selJets.makeBranch( "selJetMuEF", VFLOAT ); //*Jet_muEF)[it]; 
    selJets.makeBranch( "selJetNeEmEF", VFLOAT ); //*Jet_neEmEF)[it]; 
    selJets.makeBranch( "selJetNeHEF", VFLOAT ); //*Jet_neHEF)[it]; 
    selJets.makeBranch( "selJetNeHM", VFLOAT ); //*Jet_neHM)[it]; 

    selJets.makeBranch( "selGenJetDpt", VFLOAT ); //*Jet_genDptMatch)[it]; 
    selJets.makeBranch( "selGenJetdr", VFLOAT ); //*Jet_genDrMatch)[it]; 
    selJets.makeBranch( "selGenJetEnergy", VFLOAT ); //*Jet_genEnergy)[it]; 
    selJets.makeBranch( "selGenJeteta", VFLOAT ); //*Jet_genEta)[it]; 
    selJets.makeBranch( "selGenJetImpAng", VFLOAT ); //*Jet_genImpactAngle)[it]; 
    selJets.makeBranch( "selJetLlpDp", VFLOAT ); //*Jet_genLlpDp)[it]; 
    selJets.makeBranch( "selJetLlpDr", VFLOAT ); //*Jet_genLlpDr)[it]; 
    selJets.makeBranch( "selGenJetPt", VFLOAT ); //*Jet_genPt)[it]; 
    selJets.makeBranch( "selGenJetTof", VFLOAT ); //*Jet_genTOF)[it]; 
    selJets.makeBranch( "selGenJetTime", VFLOAT ); // (*Jet_genTime)[it]; 
    selJets.makeBranch( "selGenJetLlpTime", VFLOAT ); //*Jet_genTimeLLP)[it]; 

    selJets.makeBranch( "allJetPt", VFLOAT);
    selJets.makeBranch( "allJetMass", VFLOAT);
    selJets.makeBranch( "allJetEnergy", VFLOAT);
    selJets.makeBranch( "allJetEta", VFLOAT);
    selJets.makeBranch( "allJetPhi", VFLOAT);
    selJets.makeBranch( "allJetNPho", VINT );
    selJets.makeBranch( "allJetSCIndexs", VVUINT );
    selJets.makeBranch( "allJetWTime", VFLOAT );
    selJets.makeBranch( "allJetWTimeSig", VFLOAT );
    selJets.makeBranch( "allJetWTimeRes", VFLOAT );
    selJets.makeBranch( "allJetsWithWTime", VINT );
    selJets.makeBranch( "allJetsWithWPixTime", VINT );
    selJets.makeBranch( "allJetPixNPho", VINT );
    selJets.makeBranch( "allJetPixWTime", VFLOAT );
    selJets.makeBranch( "allJetPixWTimeSig", VFLOAT );
    selJets.makeBranch( "allJetPixWTimeRes", VFLOAT );

  	selJets.makeBranch( "pv_gjGammaTime", VFLOAT );
  	selJets.makeBranch( "pv_gjJetTime", VFLOAT );
 	selJets.makeBranch( "pv_dGJTime", VFLOAT );
    selJets.makeBranch( "pv_gjPhoIndex", VINT );
    selJets.makeBranch( "pv_dGJMinDr", VFLOAT );
  	selJets.makeBranch( "pv_diJet1Time", VFLOAT );
  	selJets.makeBranch( "pv_diJet2Time", VFLOAT );
  	selJets.makeBranch( "pv_dDiJetTime", VFLOAT );
    selJets.makeBranch( "pv_wtime", VFLOAT );
    selJets.makeBranch( "pv_wtimeres", VFLOAT );
    selJets.makeBranch( "pv_wphartres", VFLOAT );
    selJets.makeBranch( "pv_wpenergy", VFLOAT );
    selJets.makeBranch( "pv_dGJPEnergy", VFLOAT );
    selJets.makeBranch( "pv_dDiJetPEnergy", VFLOAT );
    selJets.makeBranch( "pv_haspvtime", VINT );

    // add new jet branches above
    selJets.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getJetQuality( int it ){

  const auto eta  = std::abs((*Jet_eta)[it]);
  const auto NHF  = (*Jet_neHEF)[it];
  const auto NEMF = (*Jet_neEmEF)[it];
  const auto CHF  = (*Jet_chHEF)[it];
  const auto CEMF = (*Jet_chEmEF)[it];
  const auto NHM  = (*Jet_neHM)[it];
  const auto CHM  = (*Jet_chHM)[it];
  const auto SHM  = NHM + CHM;
  const auto MUF  = (*Jet_muEF)[it];

  int tighter = 3;
  int tightLepVeto = 0;
  int tight = 2;
  int loose = 1;

  bool nhfup  = NHF  <= 0.90;   // delpho
  bool nhfup9  = NHF  <= 0.99;
  bool nhflw  = NHF  >= 0.2;

  bool nemfup1 = NEMF <= 0.90; // delpho
  bool nemfup2 = NEMF <= 0.99;
  bool nemf80 = NEMF <= 0.80;
  bool nemflw = NEMF >= 0.01;
  bool nemf10 = NEMF >= 0.1;

  bool shm1  = SHM  >= 1;
  bool muf8  = MUF  <= 0.80;
  bool mufllp  = MUF > 0.80;
  bool chf0  = CHF  >= 0;       // delpho
  bool chf10  = CHF  >= 0.10;
  bool chm0  = CHM  >= 0;       // delpho
  bool cemf8 = CEMF >= 0.80;
  bool cemfup1 = CEMF <= 0.99;
  bool nhm2  = NHM  >= 1;
  bool nhm10 = NHM  >= 10;

  bool eta0 = eta <= 2.4;
  bool eta1 = eta <= 2.6;
  bool eta2 = eta <= 2.7;
  bool eta3 = eta <= 3.0;

  if( nhfup && nemfup2 && shm1 && chf0 && chm0 && ( cemfup1 || mufllp ) ) return 4;
  if (eta1){
    if      (nhfup && nemfup1 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
    else if (nhfup && nemf80 && shm1 && chf10 && chm0) return tighter;
    else if (nhfup && nemfup1 && shm1 && chf0 && chm0) return tight;
    else    return loose;
  } else if (!eta1 && eta2 ){ //<<>>if (eta1)
    if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
    else if (nhfup && nemf80 && chm0) return tighter;
    else if (nhfup && nemfup2 && chm0) return tight;
    else    return loose;
  } else if (!eta2 && eta3){ //<<>>if (eta1) : else if
    if      (nemf10 && nemf80 && nhm2) return tighter;
    else if (nemflw && nemfup2 && nhm2) return tight;
    else    return loose;
  } else { //<<>>if (eta1) : else if : else if
    if      (nhflw && nemfup1 && nhm10) return tight;
    else    return loose;
  }//<<>>if (eta1) : else if : else if : else

  return -1; // should not happen

}//<<>>int KUCMSAodSkimmer::getJetQuality( int iter )


