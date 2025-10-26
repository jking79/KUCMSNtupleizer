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

void KUCMSAodSkimmer::processPhotons(){

  //bool verbose = true;
  bool verbose = false;

  // intilize
  selPhotons.clearBranches(); // <<<<<<<   must do
  isSelPho.clear();

  // calc
  if( DEBUG ) std::cout << "Finding photons" << std::endl;
  //----------------- photons ------------------

  std::vector<int> phoOrderIndx;
  uInt nSelPhotons = 0;
  //uint nEleVeto = 0;
  uInt nPhotons = Photon_excluded->size();	
  if( DEBUG || verbose ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
  for( uInt it = 0; it < nPhotons; it++ ){

    //---------------------------------------------------
    ///////////  pho selection ////////////////////////////////////////////////////////////////////
    if( DEBUG ) std::cout << " -- looping photons : getting pho isEB, has min pt, not excluded, electron veto " << std::endl;
    auto isExcluded = (*Photon_excluded)[it];
    auto scIndx = (*Photon_scIndex)[it];
    auto isEB = (*SuperCluster_seedIsEB)[scIndx];
    auto hasEleVeto = (*Photon_electronVeto)[it];
    bool hasPixSeed = (*Photon_pixelSeed)[it];

    //if( geVars("genSigPerfect") == 1 ) std::cout << " -- pho sel susid " << (*Gen_susId)[(*Photon_genIdx)[it]] << std::endl;
    bool isGenSig = hasGenInfoFlag ? ( (*Gen_susId)[(*Photon_genIdx)[it]] == 22 )  : 0;

    auto pt = (*Photon_pt)[it];
    bool underMinPt = pt < 30;
    auto eta = (*Photon_eta)[it];
    auto overMaxEta = std::abs(eta) > 1.479;
    auto phi = (*Photon_phi)[it];

	//if( ( eta > -3.2 && eta < -1.77 ) && ( phi > -1.77 && phi < -0.67 ) ) hasHemObj = true;

    if( DEBUG ) std::cout << " -- looping photons : getting pho iso " << std::endl;
    auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
    bool passHcalSum = true;
    auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
    bool passTrkSum = tspscdr4 < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
    auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
    bool passsEcalRhSum = erhsecdr4 < 10.0;
    auto htoem = (*Photon_hadTowOverEM)[it];
    bool passHOE = htoem < 0.02;
    bool isoPho = passHOE && passsEcalRhSum && passTrkSum && passHcalSum;
    bool failPhoIso = not isoPho;

    /*      // moved to after photon selection
	    if( DEBUG ) std::cout << " -- looping photons : getting phojet iso " << std::endl;
	    bool isJetPhoton = false;
	    phoJetVeto.clear();
	    for( int jit = 0; jit < Jet_energy->size(); jit++ ){ phoJetVeto.push_back(false); }
	    for( int jit = 0; jit < Jet_energy->size(); jit++ ){

            bool underMinJPt = (*Jet_pt)[jit] < 30.0;
            bool underMinJQual = getJetQuality(it)  < 2;
	    auto jeta = (*Jet_eta)[jit];
            auto jphi = (*Jet_phi)[jit];
	    auto overMaxJEta = std::abs(jeta) > 2.4;
            if( underMinJPt || underMinJQual || overMaxJEta ) continue;
			
            float dpjeta = jeta - eta;
            float dpjphi = dPhi( jphi, phi );
            float dr = hypo( dpjeta, dpjphi );
	    bool minDr = dr < 0.4;
            if( minDr ) isJetPhoton = true; 
	    if( minDr && isoPho ) phoJetVeto[jit] = true;

	    } // for( int jit = 0; jit < nSelJets; jit++ )
    */

    //  change to skip jets and keep all photons regardless of photon iso with jet
    bool phoskip = isExcluded || hasPixSeed || overMaxEta || underMinPt || failPhoIso;
	bool hemEligible = not underMinPt && not isExcluded; 
    //if( geVars("genSigPerfect") == 1 && isGenSig ) std::cout << " -- pho sel: phoskip " << phoskip <<  " isGenSig " << isGenSig << std::endl;  
    //if( geVars("genSigPerfect") == 1 &&  phoskip && isGenSig  ){ 
    //		std::cout << "   -- xsepji: " << isExcluded  << hasPixSeed << overMaxEta << underMinPt << isJetPhoton << failPhoIso << std::endl; }

    if( hemEligible && inHEMRegion( eta, phi ) ) hasHemObj = true;
	
    if( ( geVars("genSigPerfect") != 1 ) && phoskip ){ isSelPho.push_back(false); continue; }		
    if( ( geVars("genSigPerfect") == 1 ) &&  ( not isGenSig ) ){ isSelPho.push_back(false);  continue; }
	isSelPho.push_back(true);

    ///////////  pho selection ////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << " -- pho pull SC info" << std::endl;
    auto scSize = SuperCluster_seedIsEB->size();
    auto rhids = (*SuperCluster_rhIds)[scIndx];
    uInt nrh = rhids.size();

    if( DEBUG ) std::cout << " -- pho pull info" << std::endl;
    auto isOOT = (*Photon_isOot)[it];
    auto time = (*Photon_seedTOFTime)[it];
    auto smaj = (*SuperCluster_smaj)[scIndx];
    auto smin = (*SuperCluster_smin)[scIndx];
    auto r9 = (*Photon_r9)[it];
    auto sieie = 0.f;//(*Photon_sieie)[it];
    auto energy = (*Photon_energy)[it];
    auto cee = (*SuperCluster_covEtaEta)[scIndx];
    auto cep = (*SuperCluster_covEtaPhi)[scIndx];
    auto cpp = (*SuperCluster_covPhiPhi)[scIndx];
    auto s4 = (*Photon_s4)[it];
    auto salp = (*SuperCluster_salp)[scIndx];
    auto etawide = (*SuperCluster_etaWidth)[scIndx];
    auto phiwide = (*SuperCluster_phiWidth)[scIndx];
    auto sigmaieie = (*Photon_SigmaIEtaIEta)[it];   //!
    auto scx = (*SuperCluster_x_calo)[scIndx];
    auto scy = (*SuperCluster_y_calo)[scIndx];
    auto scz = (*SuperCluster_z_calo)[scIndx];

    int genIdx = -1;
    int momIdx = -1;
    float susId = -1;
    float genpt = -10;
    float momEnergy = -10;   //!
    float momEta = -4;   //!
    float momMass = -10;   //!
    float momPhi = -4;   //!
    float momPt = -10;   //!
    float momPx = -10;   //!
    float momPy = -10;   //!
    float momPz = -10;   //!
    float momVx = -1;   //!
    float momVy = -1;   //!
    float momVz = -100;   //!

    float ctau = -10;

    if( hasGenInfoFlag ){

      genIdx = (*Photon_genIdx)[it];
      momIdx = (*Gen_motherIdx)[genIdx];
      //momIdx = (*Photon_genSigXMomId)[it];
      susId = (*Gen_susId)[genIdx];

      if( momIdx > -1.0 ){

		momEnergy = (*Gen_energy)[momIdx];   //!
		momEta = (*Gen_eta)[momIdx];   //!
		momMass = (*Gen_mass)[momIdx];   //!
		momPhi = (*Gen_phi)[momIdx];   //!
		momPt = (*Gen_pt)[momIdx];   //!
		momPx = (*Gen_px)[momIdx];   //!
		momPy = (*Gen_py)[momIdx];   //!
		momPz = (*Gen_pz)[momIdx];   //!
		momVx = (*Gen_vx)[momIdx];   //!
		momVy = (*Gen_vy)[momIdx];   //!
		momVz = (*Gen_vz)[momIdx];   //!

      }//<<>>if( momIdx > -1.0 )

    }//if( doGenInfo )

    float rtime = time;
    float rt2 = sq2(time);
    float re = energy;
    float s1 = 0.99736-0.0035*rtime-0.00319*rt2;
    float s2 = 1.205-0.079*rtime+0.00093*rt2;//<100, <3 // 
    float s3 = std::exp( -1*sq2(rtime-0.2037)/144.942338);
    //float core = ( rtime < -0.549 ) ? 1 : ( rtime < 3.37 ) ? s1 : s2;
    float core = ( rtime < -0.549 ) ? 1 : s1;
    if( ( re < 325 ) && ( rtime > 7.0 ) ) core = s2;
    if( ( re < 200 ) && ( rtime > 6.0 ) ) core = s2;
    if( ( re < 150 ) && ( rtime > 5.0 ) ) core = s2;
    if( ( re < 100 ) && ( rtime > 4.0 ) ) core = s2;
    if( rtime > 16.0 ) core = s3;
    ////float core = ( rtime < 0.204 ) ? 1 : std::exp(-1*sq2(rtime-0.204)/144.942338);
    float ce = energy/core;
    float cpt = pt/core;		

    float cms000 = hypo(scx,scy,scz);
    float calcor = cms000/SOL;
    float tofPVtoRH = hypo(scx-PV_x,scy-PV_y,scz-PV_z);
    float pvtof = tofPVtoRH/SOL;
	float tofcor = calcor - pvtof;
    float crtime = rtime - pvtof;
    float tmeasured = rtime*SOL;
    float m1 = std::sqrt(sq2(tofPVtoRH)+8*sq2(tmeasured));
    float m2 = (tofPVtoRH+m1)*(tofPVtoRH/(4*sq2(tmeasured)));
    float m3 = tofPVtoRH/tmeasured;
    float m2_phys = ( m2 > 1 ) ? 1 : m2;
    float m3_phys = ( m3 > 1 ) ? 1 : m3;
    float MBetaEqual = 2*ce*std::sqrt(1-sq2(m2_phys));
    //float MBetaEqual = 2*ce*m2_phys;
    float MBetaPrompt = 2*ce*std::sqrt(1-sq2(m3_phys));
    //float MBetaPrompt = 2*ce*m3_phys;

    //int nSCRhids = scrhids.size();
    float sumtw = 0;
	float sumtrw = 0;
    float sumw = 0;
    float leadE = 0;
    float leadTime = -99;
    float leadAres = 0;
    float leadTres = 2000000;
	float leadSX = -1;
	bool leadWried = false;
    uInt leadRHID = 0;
    float seedE = 0;
    float seedTime = -99;
    float seedAres = 0;
    float seedTres = 2000000;
    float seedSX = -1;
    bool seedWried = false;
	bool seedGS = false;
	uInt seedRHID = 0;
	std::vector<float> erhamps;
	float sumerha = 0;
    //auto nRecHits = ECALRecHit_ID->size();
    for( int sciter = 0; sciter < nrh; sciter++  ){
        uInt pscrhid = rhids[sciter];
		int erhiter = ( rhIDtoIterMap.find(pscrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[pscrhid] : -1;
        if( erhiter != -1 ){
			float swcrss = (*ECALRecHit_swCross)[erhiter];
			bool isWeird = (*ECALRecHit_isWrd)[erhiter] || (*ECALRecHit_isDiWrd)[erhiter];
			float gainwt = 1;
            float erhe = (*ECALRecHit_energy)[erhiter];
			float erampres = (*ECALRecHit_ampres)[erhiter];
            float erhct = erh_corTime[erhiter];
			float erx = (*ECALRecHit_rhx)[erhiter];
            float ery = (*ECALRecHit_rhy)[erhiter];
            float erz = (*ECALRecHit_rhz)[erhiter];
			float cor_cms000 = hypo(erx,ery,erz)/SOL;
			float cor_tofPVtoRH = hypo(erx-PV_x,ery-PV_y,erz-PV_z)/SOL;
			float ertoftime = erhct - cor_cms000 + cor_tofPVtoRH;
			float ertres = erh_timeRes[erhiter];
			bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];
			if( hasGainSwitch ) gainwt = 0;
			float invertres = 1/ertres;
			erhamps.push_back(invertres);
			sumerha += invertres;
			float erhar = invertres*gainwt;
            sumtw += erhar*ertoftime;
			sumtrw += erhar*ertres;
            sumw += erhar;
            if( erhe*gainwt > leadE ){ 
				leadE = erhe; leadTime = ertoftime; leadAres = erampres; leadRHID = pscrhid; 
				leadTres = ertres; leadSX = swcrss; leadWried = isWeird;
			}//<<>>if( erhe*gainwt > leadE )
			if( erhe > seedE ){ 
				seedE = erhe; seedTime = ertoftime; seedAres = erampres; seedRHID = pscrhid; 
				seedTres = ertres; seedSX = swcrss; seedWried = isWeird; seedGS = hasGainSwitch;
			}//<<>>if( erhe > seedE ) 
        }//<<>>if( scrhid == rhid )
    }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
	if( sumw == 0 ){ sumw = 1; sumtw = -100; sumtrw = -1000; }
    float phoWTime = sumtw/sumw;
	float phoWRes = sumtrw/sumw;
    //float ltimeres = timeCali->getTimeResoltuion( leadAres, leadRHID, Evt_run, tctag, mctype );
	//if( ltimeres != leadTres ) std::cout << " !!!!!!!   lead res mis : " << ltimeres << " v " << leadTres << std::endl;
    //float stimeres = timeCali->getTimeResoltuion( seedAres, seedRHID, Evt_run, tctag, mctype );
    //if( stimeres != seedTres ) std::cout << " !!!!!!!   seed res mis : " << stimeres << " v " << seedTres << std::endl;
	leadTres = std::sqrt(leadTres/2);
	seedTres = std::sqrt(seedTres/2);
	phoWRes = std::sqrt(phoWRes/2);
	float leadtimesig = leadTime/leadTres;
    float seedtimesig = seedTime/seedTres;
	float wttimesig = phoWTime/phoWRes;

	for( auto as : erhamps ){ hist1d[20]->Fill(as/sumerha); }

    selPhotons.fillBranch( "selPhoLTimeSig", leadtimesig );
    selPhotons.fillBranch( "selPhoSTimeSig", seedtimesig );
    selPhotons.fillBranch( "selPhoWTimeSig", wttimesig );
    selPhotons.fillBranch( "selPhoWTime", phoWTime );
    selPhotons.fillBranch( "selPhoLTime", leadTime );
    selPhotons.fillBranch( "selPhoSTime", seedTime );
    selPhotons.fillBranch( "selPhoWTRes", phoWRes );
    selPhotons.fillBranch( "selPhoLTRes", leadTres );
    selPhotons.fillBranch( "selPhoSTRes", seedTres );
    selPhotons.fillBranch( "selPhoLSCross", leadSX );
    selPhotons.fillBranch( "selPhoSSCross", seedSX );
    selPhotons.fillBranch( "selPhoLWeird", leadWried );
    selPhotons.fillBranch( "selPhoSWeird", seedWried );
    selPhotons.fillBranch( "selPhoShasGS", seedGS );

    auto htsebcdr4 = (*Photon_hcalTowerSumEtBcConeDR04)[it];
    auto tsphcdr3 = (*Photon_trkSumPtHollowConeDR03)[it];
    auto tsphcdr4 = (*Photon_trkSumPtHollowConeDR04)[it];
    //auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
    //auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
    //auto htoem = (*Photon_hadTowOverEM)[it];
    auto hoem = (*Photon_hadOverEM)[it];
    auto sieip = 0.f;//(*Photon_sieip)[it];
    auto sipip = 0.f;//(*Photon_sipip)[it];
    auto phoOOT = (*Photon_isOot)[it];
    //auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
    auto nthcdr4 = (*Photon_nTrkHollowConeDR04)[it];   //!
    auto ntscdr4 = (*Photon_nTrkSolidConeDR04)[it];   //!
    auto epfciso = (*Photon_ecalPFClusterIso)[it];   //!
    auto hct = (*hasConversionTracks)[it];   //!
    auto hpfciso = (*Photon_hcalPFClusterIso)[it];   //!
    auto hoepuc = (*Photon_Hoe_PUcorr)[it];   //!
    auto pfciso = (*Photon_pfChargedIso)[it];   //!
    auto pscisopfpv = (*Photon_pfChargedIsoPFPV)[it];   //!
    auto pfcisowv = (*Photon_pfChargedIsoWorstVtx)[it];   //!
    auto pfphoisso = (*Photon_pfPhoIso03)[it];   //!
    auto pfriso3aq = (*pfRelIso03_all_quadratic)[it];   //!
    auto pfrtso3cq = (*pfRelIso03_chg_quadratic)[it];   //!

    //--------------------------------------------------------------
    if( DEBUG ) std::cout << " -- pho get calclated values" << std::endl;
    int phoQuality = getPhoQuality(it);
    auto phoClstrR9 = clstrR9( rhids );
    //auto phoEigens2D = getRhGrpEigenFromAngles( rhids );
    //auto evaluegeo = phoEigens2D[2];
    //auto geosmaj = phoEigens2D[3];
    //auto geosmin = phoEigens2D[4];

    float phoPhoIsoDr = 10.0;
    for( uInt it2 = 0; it2 < nPhotons; it2++ ){
      if( it == it2 ) continue;
      float eta2 = (*Photon_eta)[it2];
      float phi2 = (*Photon_phi)[it2];
      auto dphi12 = dPhi( phi, phi2 );
      auto pho2dr = hypo( eta-eta2, dphi12 ); 
      //std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << pho2dr << std::endl;
      if( pho2dr < phoPhoIsoDr ) phoPhoIsoDr = pho2dr;
    }//for( uInt it2 = it+1; it2 < nPhotons; it2++ )

    /*
    // pho object selection ------------------------------------------
    if( DEBUG ) std::cout << " -- pho obj selection" << std::endl;
    auto isMinMedQuality = phoQuality > 1;
    auto underMaxSMaj = smaj == 1.3;
    auto underMaxSMin = smin <= 0.4;
    auto overMinR9 = r9 >= 0.9;
    auto underMaxSieie = sieie <= 0.014;
    auto overMinRhCnt = nrh >= 20;

    //auto phoSelected = isMinMedQuality;
    auto phoSelected = true;
    if( not phoSelected ) continue;
    */
    //float ordpt = pt;
    //while( phoOrderIndx.count(ordpt) > 0 ){ ordpt -= 0.000001; }
    //if( DEBUG ) std::cout << " -- setting pho index : " << it << " for pt : " << ordpt << std::endl;
    phoOrderIndx.push_back(it);
    //if( DEBUG ) std::cout << " -- pho index set : " << phoOrderIndx[ordpt] << std::endl;

    // fill ( vectors )
    if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
    selPhotons.fillBranch( "selPhoOOT", phoOOT );
    selPhotons.fillBranch( "selPhoEleVeto", hasEleVeto );
    selPhotons.fillBranch( "selPhoQuality", phoQuality );
    selPhotons.fillBranch( "selPhoTime", time );
    selPhotons.fillBranch( "selPhoSCx", scx );
    selPhotons.fillBranch( "selPhoSCy", scy );
    selPhotons.fillBranch( "selPhoSCz", scz );

    selPhotons.fillBranch( "selPhoCorEnergy", ce );
    selPhotons.fillBranch( "selPhoCorPt", cpt );
    selPhotons.fillBranch( "selPhoMBetaEql", MBetaEqual );
    selPhotons.fillBranch( "selPhoMBetaPmt", MBetaPrompt );
    selPhotons.fillBranch( "selPhoPVtoRH", tofPVtoRH );

    //selPhotons.fillBranch( "selPhoGeoEgnVal", evaluegeo );
    selPhotons.fillBranch( "selPhoEta", eta );
    selPhotons.fillBranch( "selPhoPhi", phi );
    selPhotons.fillBranch( "selPhoPt", pt );
    selPhotons.fillBranch( "selPhoGenPt", genpt );
    selPhotons.fillBranch( "selPhoSMaj", smaj );
    selPhotons.fillBranch( "selPhoSMin", smin );
    selPhotons.fillBranch( "selPhoClstrRn", phoClstrR9 );
    selPhotons.fillBranch( "selPhoR9", r9 );
    selPhotons.fillBranch( "selPhoSieie", sieie );
    selPhotons.fillBranch( "selPhoSieip", sieip );
    selPhotons.fillBranch( "selPhoSipip", sipip );
    selPhotons.fillBranch( "selPhoCovEtaEta", cee );
    selPhotons.fillBranch( "selPhoCovEtaPhi", cep );
    selPhotons.fillBranch( "selPhoCovPhiPhi", cpp );
    selPhotons.fillBranch( "selPhoS4", s4 );
    selPhotons.fillBranch( "selPhoSAlp", salp );
    selPhotons.fillBranch( "selPhoEtaWidth", etawide );
    selPhotons.fillBranch( "selPhoPhiWidth", phiwide );
    selPhotons.fillBranch( "selPhoNrh", nrh );
    selPhotons.fillBranch( "selPhoEnergy" , energy );
    selPhotons.fillBranch( "selPhoSusyId" , susId );
    selPhotons.fillBranch( "selPhoSigmaIEtaIEta" , sigmaieie );
    //selPhotons.fillBranch( "selPhoGeoSMaj", geosmaj );
    //selPhotons.fillBranch( "selPhoGeoSMin", geosmin );
    selPhotons.fillBranch( "selPhoHcalTowerSumEtBcConeDR04", htsebcdr4 );
    selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR03", tsphcdr3 );
    selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR04", tsphcdr4 );
    selPhotons.fillBranch( "selPhoTrkSumPtSolidConeDR04", tspscdr4 );
    selPhotons.fillBranch( "selPhoPixelSeed", hasPixSeed );
    selPhotons.fillBranch( "selPhoEcalRHSumEtConeDR04", erhsecdr4 );
    selPhotons.fillBranch( "selPhoHadTowOverEM", htoem );
    selPhotons.fillBranch( "selPhoHadOverEM", hoem );
    selPhotons.fillBranch( "selPhoHcalTowerSumEtConeDR04", htsecdr4  );   //!
    selPhotons.fillBranch( "selPhoNTrkHollowConeDR04", nthcdr4  );   //!
    selPhotons.fillBranch( "selPhoNTrkSolidConeDR04", ntscdr4  );   //!
    selPhotons.fillBranch( "selPhoEcalPFClusterIso", epfciso  );   //!
    selPhotons.fillBranch( "selPhoHasConversionTracks", hct  );   //!
    selPhotons.fillBranch( "selPhoHcalPFClusterIso", hpfciso  );   //!
    selPhotons.fillBranch( "selPhoHoe_PUcorr", hoepuc  );   //!
    selPhotons.fillBranch( "selPhoPfChargedIso", pfciso  );   //!
    selPhotons.fillBranch( "selPhoPfChargedIsoPFPV", pscisopfpv  );   //!
    selPhotons.fillBranch( "selPhoPfChargedIsoWorstVtx", pfcisowv  );   //!
    selPhotons.fillBranch( "selPhoPfPhoIso03", pfphoisso  );   //!
    selPhotons.fillBranch( "selPhoPfRelIso03_all_quadratic", pfriso3aq  );   //!
    selPhotons.fillBranch( "selPhoPfRelIso03_chg_quadratic", pfrtso3cq  );   //!

    selPhotons.fillBranch( "selPhoPhoIsoDr", phoPhoIsoDr );
    selPhotons.fillBranch( "selPhoGenIdx", genIdx );
    selPhotons.fillBranch( "selPhoGenSigMomEnergy", momEnergy  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomEta", momEta  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomMass", momMass  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomPhi", momPhi  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomPt", momPt  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomPx", momPx  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomPy", momPy  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomPz", momPz  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomVx", momVx  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomVy", momVy  );   //!
    selPhotons.fillBranch( "selPhoGenSigMomVz", momVz  );   //!
    //if( verbose ) std::cout << " -- selPho Pt: " << pt << " phi: " << phi << " geo: " << evaluegeo << " clrn: " << phoClstrR9;
    if( verbose ) std::cout << " nrh: " << nrh << " quality: " << phoQuality << std::endl;

    nSelPhotons++;
  }//<<>>for( int it = 0; it < nPhotons; it++ )
  if( DEBUG ) std::cout << " -- pho loop finished" << std::endl;
  // fill ( other )

  /*
    std::vector<unsigned int> phoPtOrder;
    if( phoOrderIndx.size() > 0 ){ 
    for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ ){ 
    phoPtOrder.push_back(phoptit->second); 
    if( DEBUG || verbose ) std::cout << " phoptit pushback : " << phoptit->second << std::endl;
    }//for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ )
    }//if( phoOrderIndx.size() > 0 )
  */

  int loopEnd = ( nSelPhotons > 0 ) ? ( nSelPhotons < 2  ) ? 1 : 2 : 0;
  phoJetVeto.clear();
  for( int jit = 0; jit < Jet_energy->size(); jit++ ){ phoJetVeto.push_back(false); } 
  for( int pit = 0; pit < loopEnd; pit++ ){ 

    float eta = (*Photon_eta)[phoOrderIndx[pit]];
    float phi = (*Photon_phi)[phoOrderIndx[pit]];

    if( DEBUG ) std::cout << " -- looping photons : getting phojet iso " << std::endl;
    //bool isJetPhoton = false;
    for( int jit = 0; jit < Jet_energy->size(); jit++ ){

      bool underMinJPt = (*Jet_pt)[jit] < 30.0;
      bool underMinJQual = getJetQuality(jit)  < 2;
      auto jeta = (*Jet_eta)[jit];
      auto jphi = (*Jet_phi)[jit];
      auto overMaxJEta = std::abs(jeta) > 2.4;
      if( underMinJPt || underMinJQual || overMaxJEta ) continue;

      float dpjeta = jeta - eta;
      float dpjphi = dPhi( jphi, phi );
      float dr = hypo( dpjeta, dpjphi );
      bool minDr = dr < 0.4;
      //if( minDr ) isJetPhoton = true;
      if( minDr ) phoJetVeto[jit] = true;

    } // for( int jit = 0; jit < nSelJets; jit++ )

  }//<<>>for( int pit = 0; pit < loopEnd; pit++ )

  if( DEBUG || verbose ) std::cout << " - Selected " << nSelPhotons << " photons" << std::endl;
  geCnts.set( "nSelPhotons", nSelPhotons );	
  selPhotons.fillBranch( "nSelPhotons",  nSelPhotons );
  selPhotons.fillBranch( "nPhotons", nPhotons );	
  //selPhotons.fillBranch( "nEleVeto", nEleVeto );
  //selPhotons.fillBranch( "selPhoPtOrder", phoPtOrder );

  //if( DEBUG || verbose ) std::cout << " - Setting Leading/Subleading Photon with " << phoPtOrder.size() << std::endl;
  uInt leadPhoIdx = ( nSelPhotons >= 1 ) ? phoOrderIndx[0]  : 9999;
  uInt subLeadPhoIdx = ( nSelPhotons >= 2 ) ? phoOrderIndx[1] : 9999;
  geCnts.set("leadPho",leadPhoIdx);
  geCnts.set("subLeadPho",subLeadPhoIdx);
  //selPhotons.fillBranch( "leadSelPho", leadPho );
  //selPhotons.fillBranch( "subLeadSelPho", subLeadPho );
  //if( verbose ) std::cout << " -- pho lead & sublead selected : " << leadPho << " - " << subLeadPho << std::endl;
  if( DEBUG ) std::cout << " -- pho lead & sublead idx selected : " << leadPhoIdx << " - " << subLeadPhoIdx << std::endl;	
  float lPhoPt = ( nSelPhotons > 0 ) ? (*Photon_pt)[leadPhoIdx] : 0.f;
  geVars.set( "leadPhoPt", lPhoPt );
  float slPhoPt = ( nSelPhotons > 1 ) ? (*Photon_pt)[subLeadPhoIdx] : 0.f;
  geVars.set( "subLeadPhoPt", slPhoPt);
  float lPhoPhi = ( nSelPhotons > 0 ) ? (*Photon_phi)[leadPhoIdx] : 0.f;
  geVars.set( "leadPhoPhi", lPhoPhi );
  float slPhoPhi = ( nSelPhotons > 1 ) ? (*Photon_phi)[subLeadPhoIdx] : 0.f;
  geVars.set( "subLeadPhoPhi", slPhoPhi );
  float lPhoEta = ( nSelPhotons > 0 ) ? (*Photon_eta)[leadPhoIdx] : 0.f;
  geVars.set( "leadPhoEta", lPhoEta );
  float slPhoEta = ( nSelPhotons > 1 ) ? (*Photon_eta)[subLeadPhoIdx] : 0.f;
  geVars.set( "subLeadPhoEta", slPhoEta );

}//<<>>void KUCMSAodSkimmer::processPhoton(){

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setPhotonBranches( TTree* fOutTree ){

  selPhotons.makeBranch( "nPhotons", UINT );
  selPhotons.makeBranch( "nSelPhotons", UINT ); 
  selPhotons.makeBranch( "selPhoEleVeto", UINT );
  selPhotons.makeBranch( "selPhoOOT", VBOOL );
  selPhotons.makeBranch( "selPhoSusyId", VFLOAT );
  selPhotons.makeBranch( "selPhoQuality", VINT ); 
  selPhotons.makeBranch( "selPhoTime", VFLOAT ); 
 
  selPhotons.makeBranch( "selPhoLTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoSTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoWTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoWTime", VFLOAT );
  selPhotons.makeBranch( "selPhoLTime", VFLOAT );
  selPhotons.makeBranch( "selPhoSTime", VFLOAT );
  selPhotons.makeBranch( "selPhoWTRes", VFLOAT ); 
  selPhotons.makeBranch( "selPhoLTRes", VFLOAT );
  selPhotons.makeBranch( "selPhoSTRes", VFLOAT );

  selPhotons.makeBranch( "selPhoShasGS", VBOOL );
  selPhotons.makeBranch( "selPhoLSCross", VFLOAT );
  selPhotons.makeBranch( "selPhoSSCross", VFLOAT );
  selPhotons.makeBranch( "selPhoLWeird", VBOOL );
  selPhotons.makeBranch( "selPhoSWeird", VBOOL );

  selPhotons.makeBranch( "selPhoEnergy", VFLOAT );
  selPhotons.makeBranch( "selPhoEta", VFLOAT ); 
  selPhotons.makeBranch( "selPhoPhi", VFLOAT );     
  selPhotons.makeBranch( "selPhoPt", VFLOAT ); 
  selPhotons.makeBranch( "selPhoSCx", VFLOAT );
  selPhotons.makeBranch( "selPhoSCy", VFLOAT );
  selPhotons.makeBranch( "selPhoSCz", VFLOAT );
  selPhotons.makeBranch( "selPhoCorEnergy", VFLOAT );
  selPhotons.makeBranch( "selPhoCorPt", VFLOAT );
  selPhotons.makeBranch( "selPhoMBetaEql", VFLOAT );
  selPhotons.makeBranch( "selPhoMBetaPmt", VFLOAT );
  selPhotons.makeBranch( "selPhoPVtoRH", VFLOAT );

  selPhotons.makeBranch( "selPhoHcalTowerSumEtBcConeDR04", VFLOAT );
  selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR03", VFLOAT );
  selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR04", VFLOAT );
  selPhotons.makeBranch( "selPhoTrkSumPtSolidConeDR04", VFLOAT );
  selPhotons.makeBranch( "selPhoPixelSeed", VBOOL );
  selPhotons.makeBranch( "selPhoEcalRHSumEtConeDR04", VFLOAT );
  selPhotons.makeBranch( "selPhoHadTowOverEM", VFLOAT );
  selPhotons.makeBranch( "selPhoHadOverEM", VFLOAT );

  selPhotons.makeBranch( "selPhoHcalTowerSumEtConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoNTrkHollowConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoNTrkSolidConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoEcalPFClusterIso", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoHasConversionTracks", VBOOL );   //!
  selPhotons.makeBranch( "selPhoHcalPFClusterIso", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoHoe_PUcorr", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfChargedIso", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfChargedIsoPFPV", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfChargedIsoWorstVtx", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfPhoIso03", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfRelIso03_all_quadratic", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoPfRelIso03_chg_quadratic", VFLOAT );   //!

  selPhotons.makeBranch( "selPhoSMaj", VFLOAT );
  selPhotons.makeBranch( "selPhoSMin", VFLOAT );
  selPhotons.makeBranch( "selPhoSAlp", VFLOAT );
  selPhotons.makeBranch( "selPhoClstrRn", VFLOAT );
  selPhotons.makeBranch( "selPhoR9", VFLOAT );
  selPhotons.makeBranch( "selPhoNrh", VUINT );
  selPhotons.makeBranch( "selPhoCovEtaEta", VFLOAT );
  selPhotons.makeBranch( "selPhoCovEtaPhi", VFLOAT );
  selPhotons.makeBranch( "selPhoCovPhiPhi", VFLOAT );
  selPhotons.makeBranch( "selPhoSieie", VFLOAT );
  selPhotons.makeBranch( "selPhoSieip", VFLOAT );
  selPhotons.makeBranch( "selPhoSipip", VFLOAT );
  selPhotons.makeBranch( "selPhoEtaWidth", VFLOAT );
  selPhotons.makeBranch( "selPhoPhiWidth", VFLOAT );
  selPhotons.makeBranch( "selPhoS4", VFLOAT );
  selPhotons.makeBranch( "selPhoSigmaIEtaIEta", VFLOAT );   //!

  selPhotons.makeBranch( "selPhoGenPt", VFLOAT );
  selPhotons.makeBranch( "selPhoPhoIsoDr", VFLOAT );
  selPhotons.makeBranch( "selPhoGenIdx", VINT );

  selPhotons.makeBranch( "selPhoGenSigMomEnergy", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomEta", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomMass", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomPhi", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomPt", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomPx", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomPy", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomPz", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomVx", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomVy", VFLOAT );   //!
  selPhotons.makeBranch( "selPhoGenSigMomVz", VFLOAT );   //!

  // add new photon branches above
  selPhotons.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getPhoQuality( int it ){

  // determine pog quality class
  // -----------------------------------------------------
  if( DEBUG ) std::cout << " -- pho id" << std::endl;
  //bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= ( 0.006*(*Photon_pt)[it] + 4.2 );
  bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= 10.0;
  //bool hcalTowIso = (*Photon_hcalTowerSumEtBcConeDR04)[it] <= ( 0.0025*(*Photon_pt)[it] + 2.2 );
  //if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
  //bool hcTrkIsoL = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 3.5 ); //hallow cone track iso
  //bool hcTrkIsoT = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 2 ); //hallow cone track iso
  bool hcTrkIso = (*Photon_trkSumPtSolidConeDR04)[it] <= 6.0;
  //if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
  //bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.05;
  bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.02;
  //auto sieie = (*Photon_sieie)[it];
  //bool sigmaIeieEE = sieie <= 0.03; // tight only EE
  //bool sigmaIeieEB = sieie <= 0.013; // tight only EB

  //if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
  //auto baseCut = rhIso && hcalTowIso && hadOverE;
  //auto looseCut = baseCut && hcTrkIsoL;
  //auto tightCut = baseCut && hcTrkIsoT;
  //auto tightEB = tightCut && sigmaIeieEB;
  //auto tightEE = tightCut && sigmaIeieEE;
  bool isocut = rhIso && hcTrkIso && hadOverE;

  auto phoClass = isocut ? 1 : 0;
  //auto phoClass = tightCut ? 3 : looseCut ? 2 : baseCut ? 1 : 0;

  return phoClass;

}//<<>>int KUCMSAodSkimmer::getPhoQuality( int iter )



