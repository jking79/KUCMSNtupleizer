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

  bool isfastsim = ( mctype == 2 ) ? true : false;

  // intilize
  selPhotons.clearBranches(); // <<<<<<<   must do
  isBaseLinePho.clear();

  // calc
  if( DEBUG ) std::cout << "Finding photons" << std::endl;
  //----------------- photons ------------------
  // prep fpr main photon loop
  // counts 
  std::vector<int> phoOrderIndx;
  std::vector<bool> isRJRPho;
  std::vector<int> phoExcIndx;
  uInt nBaseLinePhotons = 0; 
  uInt nPhotons = Photon_excluded->size();

  diJetIndex[0] = -1;
  diJetIndex[1] = -1;
  gammaJetIndex[0] = -1;
  gammaJetIndex[1] = -1;

  //tracking (first 5) selected photons 
  // -- photon region def varibles
  //Discriminant Score Cut Values
  float bhCutVal = 0.917252;
  float notBHCutVal = 0.18523645;//1-0.81476355;
  float EEtightIsoCutVal = 0.5220398;
  float EBtightIsoCutVal = 0.5209812;
  float EELnTIsoCutVal = 0.4;
  float EBLnTIsoCutVal = 0.4;
  float EELnTIso1CutVal = 0.45;
  float EBLnTIso1CutVal = 0.45;
  float EEVeryLooseIsoCutVal = 0.21944034;//loosest isolation selection on baseline/selected photons
  float EBVeryLooseIsoCutVal = 0.2919292;
  float tagged_lead_timesig = -999; //time sig of lead ID'd photon (either as bh or not bh)
  //counting number of baseline photons passing a certain ID
  uInt nBHphotons = 0; 
  uInt nNotBHphotons = 0;
  uInt nTightIsophotons = 0;
  uInt nLnTightIsophotons = 0;
  uInt nLnTightIso1photons = 0;
  uInt nLnTightIso2photons = 0;
  uInt nPromptphotons = 0;
  uInt nNonPromptphotons = 0;
  allphoBaseline.clear();
  allphowtime.clear();
  //Time sigma window cut values
  float earlyTimeCut = -2.5;
  float lateTimeCut = 2.5;
  if( DEBUG || verbose ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;

  //-------------------------------------------------------------------------------------------
  ///////////  Loop over raw photons  /////////////////////////////////////////////////////////
  //-------------------------------------------------------------------------------------------

  // count of # of photons that pass very loose base selection
  int nBaseRJRPhos = 0; 
  //---  main photon loop 
  for( uInt it = 0; it < nPhotons; it++ ){

    //---------------------------------------------------
    ///////////  pho selection base information //////////////////////////////////////////////
    //---------------------------------------------------

    if( DEBUG ) std::cout << " -- looping photons : getting pho isEB, has min pt, not excluded, electron veto " << std::endl;
    auto isExcluded = (*Photon_excluded)[it];
    auto scIndx = (*Photon_scIndex)[it];
    auto isEB = (*SuperCluster_seedIsEB)[scIndx];
    auto hasEleVeto = (*Photon_electronVeto)[it];
    bool hasPixSeed = (*Photon_pixelSeed)[it];

    //if( geVars("genSigPerfect") == 1 ) std::cout << " -- pho sel susid " << (*Gen_susId)[(*Photon_genIdx)[it]] << std::endl;
    bool isGenSig = hasGenInfoFlag ? ( (*Gen_susId)[(*Photon_genIdx)[it]] == 22 )  : 0;

    auto energy = (*Photon_energy)[it];
    auto pt = (*Photon_pt)[it];
    bool underMinPtEB = pt < 30;
    bool underMinPtEE = pt < 30;
    bool overMinPt = pt >= 30;
    auto eta = (*Photon_eta)[it];
    auto overMaxEta = std::abs(eta) > 1.479;
    auto phi = (*Photon_phi)[it];

    //---------------------------------------------------
    /////////// skip all excluded photons ///////////////////////////////////
    //---------------------------------------------------

    if( isExcluded ) continue;

    //---------------------------------------------------
    ///////////  processing SC rechit info for  time/res/sig ////////////////////////////////////////////////////////////
    //---------------------------------------------------

    if( DEBUG ) std::cout << " -- pho pull SC info" << std::endl;
    auto scSize = SuperCluster_seedIsEB->size();
    auto rhids = (*SuperCluster_rhIds)[scIndx];
    uInt nrh = rhids.size();

    //int nSCRhids = scrhids.size();
    float sumtw = 0;
    float sumtrw = 0;
    float sumw = 0;
    float sumtw1 = 0;
    float sumtrw1 = 0;
    float sumw1 = 0;
    float leadE = 0;
    float leadTime = -99;
    float leadAres = 0;
    float leadTvar = 2000000;
    float leadSX = -1;
    bool leadWried = false;
    uInt leadRHID = 0;
    float seedE = 0;
    float seedTime = -99;
    float rawTime = -99;
    float seedAres = 0;
    float seedTvar = 2000000;
    float seedSX = -1;
    bool seedWried = false;
    bool seedGS = false;
    uInt seedRHID = 0;
    std::vector<float> erhamps;
    float sumerha = 0;
    //auto nRecHits = ECALRecHit_ID->size();
    // -------   check isTimeValid on rechit?  -------------------------------------
    for( int sciter = 0; sciter < nrh; sciter++  ){
        uInt pscrhid = rhids[sciter];
        int erhiter = ( rhIDtoIterMap.find(pscrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[pscrhid] : -1;
        if( erhiter != -1 ){

            float gainwt = 1;
            float gainwt1 = 1;

            bool isValid = (*ECALRecHit_isTimeValid)[erhiter];
            float swcrss = (*ECALRecHit_swCross)[erhiter];
            bool isWeird = (*ECALRecHit_isWrd)[erhiter] || (*ECALRecHit_isDiWrd)[erhiter];
            float erhe = (*ECALRecHit_energy)[erhiter];
            float erampres = (*ECALRecHit_ampres)[erhiter];
            float erhct = erh_corTime[erhiter];
			float rawerhct = (*ECALRecHit_time)[erhiter];
            float erx = (*ECALRecHit_rhx)[erhiter];
            float ery = (*ECALRecHit_rhy)[erhiter];
            float erz = (*ECALRecHit_rhz)[erhiter];

            float cor_cms000 = hypo(erx,ery,erz)/SOL;
            float cor_tofPVtoRH = hypo(erx-PV_x,ery-PV_y,erz-PV_z)/SOL;
            float ertoftime = erhct + cor_cms000 - cor_tofPVtoRH;
            float ertres = erh_timeRes[erhiter];

            bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];
            if( hasGainSwitch ) gainwt = 0;
			if( (*ECALRecHit_hasGS1)[erhiter] ) gainwt1 = 0;
            if( not isValid ){ gainwt = 0; gainwt1 = 0; }

            //float terror = (*ECALRecHit_timeError)[erhiter];
            //bool isEE = fabs((*ECALRecHit_eta)[erhiter]) > 1.479;
            //bool isValidT = true;
            //if( terror > 900 ) isValidT = false;
            //if( hasGainSwitch ){
            //    if( isEE ){ if( terror < 1.0 ) isValidT = false; }
            //    else { if( terror < 0.6 ) isValidT = false; }
            //} else {
            //    if( isEE ){ if( terror < 1.05 ) isValidT = false; }
            //    else { if( terror < 0.65 ) isValidT = false; }
            //}//<<>>if( hasGainSwitch )
            //if( not isValidT ){ gainwt = 0; gainwt1 = 0; }

            float invertres = 1/ertres;
            erhamps.push_back(invertres);
            sumerha += invertres;
            float erhar = invertres*gainwt;
            sumtw += erhar*ertoftime;
            sumtrw += erhar*ertres;
            sumw += erhar;
            float erhar1 = invertres*gainwt1;
            sumtw1 += erhar1*ertoftime;
            sumtrw1 += erhar1*ertres;
            sumw1 += erhar;
            if( erhe*gainwt1 > leadE ){
                leadE = erhe; leadTime = ertoftime; leadAres = erampres; leadRHID = pscrhid;
                leadTvar = ertres; leadSX = swcrss; leadWried = isWeird;
            }//<<>>if( erhe*gainwt > leadE )
            if( erhe > seedE ){
                seedE = erhe; seedTime = ertoftime; seedAres = erampres; seedRHID = pscrhid; rawTime = rawerhct;
                seedTvar = ertres; seedSX = swcrss; seedWried = isWeird; seedGS = hasGainSwitch;
            }//<<>>if( erhe > seedE ) 

			bool hasBadTime = hasGainSwitch || !(*ECALRecHit_isTimeValid)[erhiter];
			//cout << "adding rh with e " << erhe << endl;
			_ca.AddRecHit(erx, ery, erz, erhe, erhct, pscrhid, hasBadTime);

        }//<<>>if( scrhid == rhid )
    }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

    //---------------------------------------------------
    // calc wieghted times and time sigs
    //===================================================
    if( sumw == 0 ){ sumw = 1; sumtw = -100; sumtrw = -1000; }
    float phoWTimeOld = sumtw/sumw;
    float phoWHVarOld = sumtrw/sumw;
    float phoWVarOld = 1/sumw;
    float phoWTime1 = sumtw1/sumw1;
    float phoWHVar1 = sumtrw1/sumw1;
    float phoWVar1 = 1/sumw1;
    float phoWResOld = std::sqrt(phoWVarOld);
	float phoWRes1 = std::sqrt(phoWVar1);
	float leadTres = std::sqrt(leadTvar);
	float seedTres = std::sqrt(seedTvar);
    //float ltimeres = timeCali->getTimeResoltuion( leadAres, leadRHID, Evt_run, tctag, mctype );
    //if( ltimeres != leadTres ) std::cout << " !!!!!!!   lead res mis : " << ltimeres << " v " << leadTres << std::endl;
    //float stimeres = timeCali->getTimeResoltuion( seedAres, seedRHID, Evt_run, tctag, mctype );
    //if( stimeres != seedTres ) std::cout << " !!!!!!!   seed res mis : " << stimeres << " v " << seedTres << std::endl;
    float leadtimesig = leadTime/leadTres;// leadTres == variance
    float seedtimesig = seedTime/seedTres;
    float wttimesig = phoWTimeOld/phoWResOld;
    float wttimesig1 = phoWTime1/phoWRes1;

    float phoWTime, phoWRes;
	float phoWTimeSig = getTimeSig( scIndx, phoWTime, phoWRes );
    allphowtime.push_back( phoWTime );

    float scx = (*SuperCluster_clcx)[scIndx];
    float scy = (*SuperCluster_clcy)[scIndx];
    float scz = (*SuperCluster_clcz)[scIndx];

    //----------------------------------------------------------------------------------
    //  getting pho GEN info
    //---------------------------------------------------------------------------------

    int genIdx = -1;
    int momIdx = -1;
    float susId = -1;

    float genEnergy = -10;   //!
    float genpt = -10;
    float genPt = -10;   //!
    float genPx = -10;   //!
    float genPy = -10;   //!
    float genPz = -10;   //!
    float genVx = -100;   //!
    float genVy = -100;   //!
    float genVz = -100;   //!

    float momEnergy = -10;   //!
    float momEta = -4;   //!
    float momMass = -10;   //!
    float momPhi = -4;   //!
    float momPt = -10;   //!
    float momPx = -10;   //!
    float momPy = -10;   //!
    float momPz = -10;   //!
    float momVx = -100;   //!
    float momVy = -100;   //!
    float momVz = -100;   //!

	float distPho = 0;
	float distMom = 0;
	float betamom = 1;

	float cor_gtofPVtoSC = 0;
    float labtime = -100;
    float gentime = -100;
    float labtimesig = -100;
	float gentimesig = -100;

    if( hasGenInfoFlag ){

        genIdx = (*Photon_genIdx)[it];
        momIdx = (*Gen_motherIdx)[genIdx];
        susId = (*Gen_susId)[genIdx];

        genEnergy = (*Gen_energy)[genIdx];   //!
        genPt = (*Gen_pt)[genIdx];   //!
        genPx = (*Gen_px)[genIdx];   //!
        genPy = (*Gen_py)[genIdx];   //!
        genPz = (*Gen_pz)[genIdx];   //!
        genVx = (*Gen_vx)[genIdx];   //!
        genVy = (*Gen_vy)[genIdx];   //!
        genVz = (*Gen_vz)[genIdx];   //!

		distPho = hypo( scx - genVx, scy - genVy, scz - genVz );
        cor_gtofPVtoSC = hypo(scx-PV_x,scy-PV_y,scz-PV_z);

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

			distMom = hypo( genVx - momVx, genVy - momVy, genVz - momVz );
			betamom = hypo( momPx, momPy, momPz )/momEnergy;

        }//<<>>if( momIdx > -1.0 )

		float cor_gtofPVtoSCSOL = cor_gtofPVtoSC/SOL;
		labtime = ( distPho + distMom/betamom )/SOL;
		gentime = timeCali->getSmearedTime( labtime, phoWRes/3.0);
		labtime -= cor_gtofPVtoSCSOL;
		gentime -= cor_gtofPVtoSCSOL;
		labtimesig = labtime/phoWRes;
		gentimesig = gentime/phoWRes;		

		if( susId == 22 && isfastsim ){

			phoWTime = gentime;
			phoWTimeSig = gentimesig;

		}//<<>>if( susId == 22 )

    }//if( doGenInfo )

    //---------------------------------------------------
    ///////////  saving hem region info  ////////////////////////////////////////////////////////////////////
    //---------------------------------------------------

    bool hemEligible1 = pt > 20 && not isExcluded && not hasPixSeed;
    bool hemEligible2 = pt > 30 && not isExcluded && not hasPixSeed;

    bool isInHemRegion = inHEMRegion( eta, phi );
    hemBits.set( "pho1hvl", isInHemRegion && hemEligible1 );
    hemBits.set( "pho2hvm", isInHemRegion && hemEligible2 );

    //---------------------------------------------------
    ///////////  pho disriminate ids ////////////////////////////////////////////////////////////////////
    //---------------------------------------------------

	float isobkg_score = 0;
	float nonisobkg_score = 0;
	float physbkg_score = 0;
	float bh_score = 0;
	bool bhtag = false;
	bool pbtag = false;
	bool tightisotag = false;
	bool lntisotag = false;
	bool lntiso1tag = false;
	bool lntiso2tag = false;
	bool gjetscr = false;
	bool dijetscr = false;
	bool spikecr = false;
	bool bhcr = false;
    	ClusterObj phoobj;
    	_ca.NoClusterRhs(phoobj, true);
    	vector<float> photonIDscores;
     	map<string, double> isomap;
        MakePhotonIsoMap(it, isomap);
	if(overMaxEta){ //endcap
		phoobj.CalculateEndcapPhotonIDScores(pt, isomap);
	} else { //barrel
		phoobj.CalculateBarrelPhotonIDScores(pt, isomap);
	}//<<>>if(overMaxEta)
	phoobj.GetPhotonIDScores(photonIDscores);
    isobkg_score = photonIDscores[0];
    nonisobkg_score = photonIDscores[1];
	//do det bkg scores
	vector<vector<float>> detbkgScores;
	phoobj.CalculateDetBkgScores(true);
	phoobj.GetDetBkgScores(detbkgScores);
	//binary classifier
	physbkg_score = detbkgScores[0][0];
	bh_score = detbkgScores[0][1];

	if(bh_score > bhCutVal){
		bhtag = true;
	}
	if(bh_score <= notBHCutVal){
		pbtag = true;
	}
	if(overMaxEta){
		if(isobkg_score > EELnTIsoCutVal && isobkg_score <= EEtightIsoCutVal ){
			lntisotag = true;
			if(isobkg_score > EELnTIso1CutVal)
				lntiso2tag = true;
			else
				lntiso1tag = true;
		}
		if(isobkg_score > EEtightIsoCutVal )
			tightisotag = true;
	}
	else{
		if(isobkg_score > EBLnTIsoCutVal && isobkg_score <= EBtightIsoCutVal ){
			lntisotag = true;
			if(isobkg_score > EBLnTIso1CutVal)
				lntiso2tag = true;
			else
				lntiso1tag = true;
		}
		if(isobkg_score > EBtightIsoCutVal )
			tightisotag = true;
	}

	gjetscr = GetDiJetsCR(it);
	dijetscr = GetGJetsCR(it);
	int detbkgcr = doSVs ? GetDetBkgCR(it) : -1;
	if(detbkgcr == 2) bhcr = true;
	if(detbkgcr == 3) spikecr = true;

    if(DEBUG) cout << "photon " << it << " iso score " << isobkg_score << " noniso score " << nonisobkg_score << endl;
    if( DEBUG ) std::cout << " -- looping photons : getting pho iso " << std::endl;
    auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
    auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
    auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
    auto htoem = (*Photon_hadTowOverEM)[it];
    bool isoPho;
    if(overMaxEta){ //endcap
    	if(isobkg_score >= EEVeryLooseIsoCutVal) isoPho = true;
	    else isoPho = false;
    } else { //barrel - turn off any iso presel (for now) to study the discriminator in v48 skims
    	if(isobkg_score >= EBVeryLooseIsoCutVal) isoPho = true;
	    else isoPho = false;
    	//bool passHcalSum = true;
    	//bool passTrkSum = tspscdr4/pt < 0.12; // using rel = 6.0/50  ( abs = 6.0 fir pt of 50 GeV )
    	//bool passsEcalRhSum = erhsecdr4/pt < 0.2; // using rel = 10/50 ( abs cut = 10, for pt of 50 GeV )
    	//bool passHOE = htoem < 0.02; // using abs value 0.02
    	//isoPho = passHOE && passsEcalRhSum && passTrkSum && passHcalSum;
    }

    //---------------------------------------------------
    //adding criteria for 'selected' or 'very loose base selection' photons - need at least 2 rechits that are within [-20, 20] ns
    bool badNRechits  = (_ca.GetNRecHits() < 2);

  	_ca.ClearRecHitList();

    //---------------------------------------------------------------------------------
    ///////////  Determining photon Base selection   ////////////////////////////////////////////////////////////////////
    //------------------------------------------------------------------------------

    // Only the 2 leading photons that pass very loose criteria should be in Base selction 
    //  change to skip jets and keep all photons regardless of photon iso with jet

    bool pass_very_loose_id = isoPho && overMinPt; // ?  this is place holder - was not in place in v37  and prior to 2/17/26 

	bool standard_selction = not badNRechits and not hasPixSeed and not isExcluded;
    bool underNMaxBasePhos = nBaseRJRPhos < 2;// checks how many RJR photons already found - need to have found less then 2 already
    bool in_base_selection = standard_selction and pass_very_loose_id and underNMaxBasePhos;

    //---------------------------------------------------
    ///////////  saving info on EB/EE photon information ///////////////////////////////////////////////////////
    //---------------------------------------------------
    bool isEESig = not isExcluded && not hasPixSeed && overMaxEta;
    if( isEESig ){

        selPhotons.fillBranch( "endcap_photon_baseline", in_base_selection );   //!
        selPhotons.fillBranch( "endcap_photon_energy", (*Photon_energy)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_eta", (*Photon_eta)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_phi", (*Photon_phi)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_tspscdr4", (*Photon_trkSumPtSolidConeDR04)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_erhsecdr4", (*Photon_ecalRHSumEtConeDR04)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_htoem", (*Photon_hadTowOverEM)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_pt", (*Photon_pt)[it] );   //!
        selPhotons.fillBranch( "endcap_photon_isoANNScore", isobkg_score );   //!

    }//if( doGenInfo )

    bool isEBSig = not isExcluded && not hasPixSeed && not overMaxEta;
    if( isEBSig ){

        selPhotons.fillBranch( "barrel_photon_baseline", in_base_selection );   //!
        selPhotons.fillBranch( "barrel_photon_energy", (*Photon_energy)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_eta", (*Photon_eta)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_phi", (*Photon_phi)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_tspscdr4", (*Photon_trkSumPtSolidConeDR04)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_erhsecdr4", (*Photon_ecalRHSumEtConeDR04)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_htoem", (*Photon_hadTowOverEM)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_pt", (*Photon_pt)[it] );   //!
        selPhotons.fillBranch( "barrel_photon_isoANNScore", isobkg_score );   //!

    }//if( doGenInfo )

    //---------------------------------------------------
    ///////////  saving info for all non excluded photons  ////////////////////////////////////////////////////////////////////
    //---------------------------------------------------

	allphoBaseline.push_back( in_base_selection );
    selPhotons.fillBranch( "photon_baseline", in_base_selection );   //!
    selPhotons.fillBranch( "photon_WTimeSig", phoWTimeSig );
    selPhotons.fillBranch( "photon_WTime", phoWTime );
    selPhotons.fillBranch( "photon_LeadTime", leadTime );
    selPhotons.fillBranch( "photon_LeadTimeSig", leadtimesig );
    selPhotons.fillBranch( "photon_Eta", eta );
    selPhotons.fillBranch( "photon_Phi", phi );
    selPhotons.fillBranch( "photon_Pt", pt );
    selPhotons.fillBranch( "photon_E", energy );
    selPhotons.fillBranch( "photon_PixSeed", hasPixSeed );
    selPhotons.fillBranch( "photon_isoANNScore", isobkg_score );
    selPhotons.fillBranch( "photon_beamHaloCNNScore", bh_score );

    ///////////  Very Loose Base Photon selection ////////////////////////////////////////////////////////////////////
    //---------------------------------------------------
    //  phoskip all photons that do not pass selection criteria
    //  moved to here so that EE and all photon info could be saved with WTSig info
    //---------------------------------------------------
    ///////////  pho selection ////////////////////////////////////////////////////////////////////

    if( ( geVars("genSigPerfect") != 1 ) && not in_base_selection ){ isBaseLinePho.push_back(false); continue; }
    if( ( geVars("genSigPerfect") == 1 ) &&  ( not isGenSig ) ){ isBaseLinePho.push_back(false);  continue; }
    isBaseLinePho.push_back(true);

    ///////////////////////////////  depreciated but still used in process
    // setting flag for photons that go into rjr - should be all that get to here now.
    // we want to reject (endcap) photons that are at the most noniso edge of the noniso sideband. 
    bool isRjrSigPho = true; // overMaxEta ? nonisobkg_score < EEVeryLnTIsoCutVal : isoPho;
    if( isRjrSigPho ){ isRJRPho.push_back(true); nBaseRJRPhos++; }
    else isRJRPho.push_back(false);

	//----------------------------------------------------------------------------------
	//  getting pho SC info
	//---------------------------------------------------------------------------------

    if( DEBUG ) std::cout << " -- pho pull info" << std::endl;
    auto isOOT = (*Photon_isOot)[it];
    auto time = (*Photon_seedTOFTime)[it];
    auto smaj = (*SuperCluster_smaj)[scIndx];
    auto smin = (*SuperCluster_smin)[scIndx];
    auto r9 = (*Photon_r9)[it];
    auto sieie = 0.f;//(*Photon_sieie)[it];
    auto cee = (*SuperCluster_covEtaEta)[scIndx];
    auto cep = (*SuperCluster_covEtaPhi)[scIndx];
    auto cpp = (*SuperCluster_covPhiPhi)[scIndx];
    auto s4 = (*Photon_s4)[it];
    auto salp = (*SuperCluster_salp)[scIndx];
    auto etawide = (*SuperCluster_etaWidth)[scIndx];
    auto phiwide = (*SuperCluster_phiWidth)[scIndx];
    auto sigmaieie = (*Photon_SigmaIEtaIEta)[it];   //!
	// moved to gen time info secction
    //auto scx = (*SuperCluster_x_calo)[scIndx];
    //auto scy = (*SuperCluster_y_calo)[scIndx];
    //auto scz = (*SuperCluster_z_calo)[scIndx];

	//---------------------------------------------------------------------------
	//  Calulate time based energy and pt bais corrections
	//---------------------------------------------------------------------------

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

    //---------------------------------------------------------------------------
    //  Calulate beta values photon mother without gen info - R&D project
    //---------------------------------------------------------------------------

    float cms000 = hypo(scx,scy,scz);
    float calcor = cms000/SOL;
    float tofPVtoRH = hypo(scx-PV_x,scy-PV_y,scz-PV_z);
    float pvtof = tofPVtoRH/SOL;
	float tofcor = calcor - pvtof;
    float crtime = rtime - tofcor;
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

    //---------------------------------------------------------------------------
    //  Fill Photon branch time and sc information
    //---------------------------------------------------------------------------

	for( auto as : erhamps ){ hist1d[20]->Fill(as/sumerha); }

    selPhotons.fillBranch( "baseLinePhoton_RawTime", rawTime );
    selPhotons.fillBranch( "baseLinePhoton_WTimeSig", phoWTimeSig );
    selPhotons.fillBranch( "baseLinePhoton_LTimeSig", leadtimesig );
    selPhotons.fillBranch( "baseLinePhoton_STimeSig", seedtimesig );
    selPhotons.fillBranch( "baseLinePhoton_WTimeSigOld", wttimesig );
    selPhotons.fillBranch( "baseLinePhoton_WTime", phoWTime );
    selPhotons.fillBranch( "baseLinePhoton_WTimeOld", phoWTimeOld );
    selPhotons.fillBranch( "baseLinePhoton_LTime", leadTime );
    selPhotons.fillBranch( "baseLinePhoton_STime", seedTime );
    selPhotons.fillBranch( "baseLinePhoton_WTRes", phoWRes );
    selPhotons.fillBranch( "baseLinePhoton_WTResOld", phoWResOld );
    selPhotons.fillBranch( "baseLinePhoton_LTRes", leadTres );
    selPhotons.fillBranch( "baseLinePhoton_STRes", seedTres );
    selPhotons.fillBranch( "baseLinePhoton_LSCross", leadSX );
    selPhotons.fillBranch( "baseLinePhoton_SSCross", seedSX );
    selPhotons.fillBranch( "baseLinePhoton_LWeird", leadWried );
    selPhotons.fillBranch( "baseLinePhoton_SWeird", seedWried );
    selPhotons.fillBranch( "baseLinePhoton_ShasGS", seedGS );

    selPhotons.fillBranch("baseLinePhoton_isoANNScore",isobkg_score);
    selPhotons.fillBranch("baseLinePhoton_nonIsoANNScore",nonisobkg_score);
    selPhotons.fillBranch("baseLinePhoton_physBkgCNNScore",physbkg_score);
    selPhotons.fillBranch("baseLinePhoton_beamHaloCNNScore",bh_score);
    selPhotons.fillBranch("baseLinePhoton_beamHaloID",bhtag);
    selPhotons.fillBranch("baseLinePhoton_notBeamHaloID",pbtag);
    selPhotons.fillBranch("baseLinePhoton_tightIsoID",tightisotag);
    selPhotons.fillBranch("baseLinePhoton_looseNotTightIsoID",lntisotag);
    selPhotons.fillBranch("baseLinePhoton_looseNotTightIso1ID",lntiso1tag);
    selPhotons.fillBranch("baseLinePhoton_looseNotTightIso2ID",lntiso2tag);
    selPhotons.fillBranch("baseLinePhoton_beamHaloCR",bhcr);
    selPhotons.fillBranch("baseLinePhoton_spikeCR",spikecr);
    selPhotons.fillBranch("baseLinePhoton_GJetsCR",gjetscr);
    selPhotons.fillBranch("baseLinePhoton_DiJetsCR",dijetscr);

    //---------------------------------------------------------------------------
    //  get Photon isolation information
    //---------------------------------------------------------------------------

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

    //---------------------------------------------------------------------------
    //  check photon v photon overlap for met correction
    //---------------------------------------------------------------------------

    float phoPhoIsoDr = 10.0;
	int phoIsoIndx = -1;
	if( (*Photon_isOot)[it] ){
    for( uInt it2 = 0; it2 < nPhotons; it2++ ){
    	if( it == it2 ) continue;
		if( not (*Photon_excluded)[it2] ) continue; 
      	float eta2 = (*Photon_eta)[it2];
      	float phi2 = (*Photon_phi)[it2];
      	auto dphi12 = dPhi( phi, phi2 );
      	auto pho2dr = hypo( eta-eta2, dphi12 ); 
      	//std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << pho2dr << std::endl;
      	if( pho2dr < phoPhoIsoDr ){ phoPhoIsoDr = pho2dr; phoIsoIndx = it2; }
    }//for( uInt it2 = it+1; it2 < nPhotons; it2++ )
	}//<<>>if( (*Photon_isOot)[it] )
	//bool matchPhoExc = ( phoIsoIndx > -1 ) ? (*Photon_excluded)[it] : false;
	if( phoPhoIsoDr < 0.3 ) phoExcIndx.push_back( phoIsoIndx ); else phoExcIndx.push_back( -1 );
	
    //---------------------------------------------------------------------------
    //  fill photon region defintion information 
    //---------------------------------------------------------------------------


    //if( DEBUG ) std::cout << " -- setting pho index : " << it << " for pt : " << ordpt << std::endl;
    phoOrderIndx.push_back(it);
    //phoOrderId.push_back(phoObjId);
    //if( DEBUG ) std::cout << " -- pho index set : " << phoOrderIndx[ordpt] << std::endl;


    //---------------------------------------------------------------------------
    //  fill photon branches 
    //---------------------------------------------------------------------------

    // fill ( vectors )
    if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
    selPhotons.fillBranch( "baseLinePhoton_OOT", phoOOT );
    selPhotons.fillBranch( "baseLinePhoton_EleVeto", hasEleVeto );
    selPhotons.fillBranch( "baseLinePhoton_Quality", phoQuality );
    selPhotons.fillBranch( "baseLinePhoton_Time", time );
    selPhotons.fillBranch( "baseLinePhoton_SCx", scx );
    selPhotons.fillBranch( "baseLinePhoton_SCy", scy );
    selPhotons.fillBranch( "baseLinePhoton_SCz", scz );

    selPhotons.fillBranch( "baseLinePhoton_CorEnergy", ce );
    selPhotons.fillBranch( "baseLinePhoton_CorPt", cpt );
    selPhotons.fillBranch( "baseLinePhoton_MBetaEql", MBetaEqual );
    selPhotons.fillBranch( "baseLinePhoton_MBetaPmt", MBetaPrompt );
    selPhotons.fillBranch( "baseLinePhoton_PVtoRH", tofPVtoRH );

    //selPhotons.fillBranch( "baseLinePhoton_GeoEgnVal", evaluegeo );
    selPhotons.fillBranch( "baseLinePhoton_Eta", eta );
    selPhotons.fillBranch( "baseLinePhoton_Phi", phi );
    selPhotons.fillBranch( "baseLinePhoton_Pt", pt );
    selPhotons.fillBranch( "baseLinePhoton_GenPt", genpt );
    selPhotons.fillBranch( "baseLinePhoton_SMaj", smaj );
    selPhotons.fillBranch( "baseLinePhoton_SMin", smin );
    selPhotons.fillBranch( "baseLinePhoton_ClstrRn", phoClstrR9 );
    selPhotons.fillBranch( "baseLinePhoton_R9", r9 );
    selPhotons.fillBranch( "baseLinePhoton_Sieie", sieie );
    selPhotons.fillBranch( "baseLinePhoton_Sieip", sieip );
    selPhotons.fillBranch( "baseLinePhoton_Sipip", sipip );
    selPhotons.fillBranch( "baseLinePhoton_CovEtaEta", cee );
    selPhotons.fillBranch( "baseLinePhoton_CovEtaPhi", cep );
    selPhotons.fillBranch( "baseLinePhoton_CovPhiPhi", cpp );
    selPhotons.fillBranch( "baseLinePhoton_S4", s4 );
    selPhotons.fillBranch( "baseLinePhoton_SAlp", salp );
    selPhotons.fillBranch( "baseLinePhoton_EtaWidth", etawide );
    selPhotons.fillBranch( "baseLinePhoton_PhiWidth", phiwide );
    selPhotons.fillBranch( "baseLinePhoton_Nrh", nrh );
    selPhotons.fillBranch( "baseLinePhoton_Energy" , energy );
    selPhotons.fillBranch( "baseLinePhoton_SusyId" , susId );
    selPhotons.fillBranch( "baseLinePhoton_SigmaIEtaIEta" , sigmaieie );
    //selPhotons.fillBranch( "baseLinePhoton_GeoSMaj", geosmaj );
    //selPhotons.fillBranch( "baseLinePhoton_GeoSMin", geosmin );
    selPhotons.fillBranch( "baseLinePhoton_HcalTowerSumEtBcConeDR04", htsebcdr4 );
    selPhotons.fillBranch( "baseLinePhoton_TrkSumPtHollowConeDR03", tsphcdr3 );
    selPhotons.fillBranch( "baseLinePhoton_TrkSumPtHollowConeDR04", tsphcdr4 );
    selPhotons.fillBranch( "baseLinePhoton_TrkSumPtSolidConeDR04", tspscdr4 );
    selPhotons.fillBranch( "baseLinePhoton_PixelSeed", hasPixSeed );
    selPhotons.fillBranch( "baseLinePhoton_EcalRHSumEtConeDR04", erhsecdr4 );
    selPhotons.fillBranch( "baseLinePhoton_HadTowOverEM", htoem );
    selPhotons.fillBranch( "baseLinePhoton_HadOverEM", hoem );
    selPhotons.fillBranch( "baseLinePhoton_HcalTowerSumEtConeDR04", htsecdr4  );   //!
    selPhotons.fillBranch( "baseLinePhoton_NTrkHollowConeDR04", nthcdr4  );   //!
    selPhotons.fillBranch( "baseLinePhoton_NTrkSolidConeDR04", ntscdr4  );   //!
    selPhotons.fillBranch( "baseLinePhoton_EcalPFClusterIso", epfciso  );   //!
    selPhotons.fillBranch( "baseLinePhoton_HasConversionTracks", hct  );   //!
    selPhotons.fillBranch( "baseLinePhoton_HcalPFClusterIso", hpfciso  );   //!
    selPhotons.fillBranch( "baseLinePhoton_Hoe_PUcorr", hoepuc  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfChargedIso", pfciso  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfChargedIsoPFPV", pscisopfpv  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfChargedIsoWorstVtx", pfcisowv  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfPhoIso03", pfphoisso  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfRelIso03_all_quadratic", pfriso3aq  );   //!
    selPhotons.fillBranch( "baseLinePhoton_PfRelIso03_chg_quadratic", pfrtso3cq  );   //!

    selPhotons.fillBranch( "baseLinePhoton_GenLabTime", labtime );
    selPhotons.fillBranch( "baseLinePhoton_GenLabTimeSig", labtimesig );
    selPhotons.fillBranch( "baseLinePhoton_GenTime", gentime );
    selPhotons.fillBranch( "baseLinePhoton_GenTimeSig", gentimesig );

    selPhotons.fillBranch( "baseLinePhoton_PhoIsoDr", phoPhoIsoDr );
    selPhotons.fillBranch( "baseLinePhoton_GenIdx", genIdx );
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomEnergy", momEnergy  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomEta", momEta  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomMass", momMass  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomPhi", momPhi  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomPt", momPt  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomPx", momPx  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomPy", momPy  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomPz", momPz  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomVx", momVx  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomVy", momVy  );   //!
    selPhotons.fillBranch( "baseLinePhoton_GenSigMomVz", momVz  );   //!
    //if( verbose ) std::cout << " -- selPho Pt: " << pt << " phi: " << phi << " geo: " << evaluegeo << " clrn: " << phoClstrR9;
    if( verbose ) std::cout << " nrh: " << nrh << " quality: " << phoQuality << std::endl;

	if(nBHphotons == 0 && bhtag)
		tagged_lead_timesig = phoWTimeSig;
	if(nBHphotons == 0 && nNotBHphotons == 0 && pbtag)
		tagged_lead_timesig = phoWTimeSig;
	//count photons passing ID
    if(bhtag)
	    nBHphotons++;
    if(pbtag)
	    nNotBHphotons++;
    if(tightisotag)
	    nTightIsophotons++;
    if(lntisotag){
	    nLnTightIsophotons++;
    	if(lntiso1tag)
		nLnTightIso1photons++;
	if(lntiso2tag)
		nLnTightIso2photons++;
    }
    if(phoWTimeSig < earlyTimeCut || phoWTimeSig > lateTimeCut)
	    nNonPromptphotons++;
    else
	    nPromptphotons++;

    nBaseLinePhotons++;

  }//<<>>for( int it = 0; it < nPhotons; it++ )
  // --------------    End of Photon Loop ----------------------------------------------------------------------


  if( nPhotons == 0 ){ hemBits.set( "pho1hvl", false ); hemBits.set( "pho2hvm", false ); }
  if( DEBUG ) std::cout << " -- pho loop finished" << std::endl;

  //---------------------------------------------------------------------------
  //  setting jet veto for for overlapping photons
  //---------------------------------------------------------------------------
  
  int loopEnd = ( nBaseLinePhotons > 0 ) ? ( nBaseLinePhotons < 2  ) ? 1 : 2 : 0;
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

  //---------------------------------------------------------------------------
  //  photon ID counting
  //---------------------------------------------------------------------------

  if( DEBUG || verbose ) std::cout << " - Selected " << nBaseLinePhotons << " photons" << std::endl;
 
  geCnts.set( "nSelPhotons", nBaseLinePhotons );	
  selPhotons.fillBranch( "nBaseLinePhotons",  nBaseLinePhotons );
  selPhotons.fillBranch( "nPhotons", nPhotons );

  //---------------------------------------------------------------------------
  //  Fill region index
  //---------------------------------------------------------------------------
  selPhotons.fillBranch( "nblBHPhotons", nBHphotons );
  selPhotons.fillBranch( "nblNotBHPhotons", nNotBHphotons );
  selPhotons.fillBranch( "nblTightIsoPhotons", nTightIsophotons );
  selPhotons.fillBranch( "nblLooseNotTightIsoPhotons", nLnTightIsophotons );
  selPhotons.fillBranch( "nblNonpromptPhotons", nNonPromptphotons );
  selPhotons.fillBranch( "nblPromptPhotons", nPromptphotons );


  //---------------------------------------------------------------------------
  //  end of region defintions fill - saving region defintion branches
  //---------------------------------------------------------------------------

  //delayed CRs - beam halo
  selPhotons.fillBranch("passNPhoGe1SelectionBeamHaloCR", bool(nBHphotons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionBeamHaloCR", bool(nBHphotons == 1));
  selPhotons.fillBranch("passNPhoEq2SelectionBeamHaloCR", bool(nBHphotons == 2));
  //beam halo early CR
  selPhotons.fillBranch("passNPhoGe1SelectionEarlyBeamHaloCR", bool(nBHphotons > 0 && tagged_lead_timesig < earlyTimeCut));
  selPhotons.fillBranch("passNPhoEq1SelectionEarlyBeamHaloCR", bool(nBHphotons == 1 && tagged_lead_timesig < earlyTimeCut));
  selPhotons.fillBranch("passNPhoEq2SelectionEarlyBeamHaloCR", bool(nBHphotons == 2 && tagged_lead_timesig < earlyTimeCut));
  //beam halo late CR
  selPhotons.fillBranch("passNPhoGe1SelectionLateBeamHaloCR", bool(nBHphotons > 0 && tagged_lead_timesig > lateTimeCut));
  selPhotons.fillBranch("passNPhoEq1SelectionLateBeamHaloCR", bool(nBHphotons == 1 && tagged_lead_timesig > lateTimeCut));
  selPhotons.fillBranch("passNPhoEq2SelectionLateBeamHaloCR", bool(nBHphotons == 2 && tagged_lead_timesig > lateTimeCut));
  selPhotons.fillBranch("tagged_lead_timesig", tagged_lead_timesig);
  //not beam halo early CR - FOR VALIDATION ONLY
  selPhotons.fillBranch("passNPhoGe1SelectionNotBHCR", bool(nBHphotons == 0 && nNotBHphotons > 0 && nNonPromptphotons > 0));
  //not beam halo early CR
  selPhotons.fillBranch("passNPhoGe1SelectionEarlyNotBHCR", bool(nBHphotons == 0 && nNotBHphotons > 0 && tagged_lead_timesig < earlyTimeCut && nNonPromptphotons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionEarlyNotBHCR", bool(nBHphotons == 0 && nNotBHphotons == 1 && tagged_lead_timesig < earlyTimeCut && nNonPromptphotons > 0));
  selPhotons.fillBranch("passNPhoEq2SelectionEarlyNotBHCR", bool(nBHphotons == 0 && nNotBHphotons == 2 && tagged_lead_timesig < earlyTimeCut && nNonPromptphotons > 0));

  //delayed SR - not beam halo late SR
  selPhotons.fillBranch("passNPhoGe1SelectionLateNotBHSR", bool(nBHphotons == 0 && nNotBHphotons > 0 && tagged_lead_timesig > lateTimeCut && nNonPromptphotons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionLateNotBHSR", bool(nBHphotons == 0 && nNotBHphotons == 1 && tagged_lead_timesig > lateTimeCut && nNonPromptphotons > 0));
  selPhotons.fillBranch("passNPhoEq2SelectionLateNotBHSR",  bool(nBHphotons == 0 && nNotBHphotons == 2 && tagged_lead_timesig > lateTimeCut && nNonPromptphotons > 0));
  
  ////prompt CRs
  bool ge1notbh0bh = (nNotBHphotons > 0 && nBHphotons == 0 && nNonPromptphotons > 0);
  bool notbhprompt = (tagged_lead_timesig > earlyTimeCut && tagged_lead_timesig < lateTimeCut);
  bool ge1notbh0bhprompt = (ge1notbh0bh && notbhprompt);
  bool eq0nonprompt = (nNonPromptphotons == 0 || ge1notbh0bhprompt); //prompt photon regions are defined by either prompt-prompt OR (at least 1 nonprompt photon but lead notbh photon is prompt (no bh photons) and at least 1 is notbh tagged)
	//debugging logic

  selPhotons.fillBranch("passNPhoGe1SelectionNotBHPrompt", ge1notbh0bhprompt);


  //loose-not-tight iso CR
  selPhotons.fillBranch("passNPhoGe1SelectionPromptLooseNotTightIsoCR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIsophotons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionPromptLooseNotTightIsoCR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIsophotons == 1));
  selPhotons.fillBranch("passNPhoEq2SelectionPromptLooseNotTightIsoCR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIsophotons == 2));

  //validation (split) loose-not-tight iso CR
  selPhotons.fillBranch("passNPhoGe1SelectionPromptLooseNotTightIso1CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons == 0 && nLnTightIso1photons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionPromptLooseNotTightIso1CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons == 0 && nLnTightIso1photons == 1));
  selPhotons.fillBranch("passNPhoEq2SelectionPromptLooseNotTightIso1CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons == 0 && nLnTightIso1photons == 2));
  //validation (split) (loose-not-tight) proxy iso SR (CR)
  selPhotons.fillBranch("passNPhoGe1SelectionPromptLooseNotTightIso2CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionPromptLooseNotTightIso2CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons == 1));
  selPhotons.fillBranch("passNPhoEq2SelectionPromptLooseNotTightIso2CR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 0 && nLnTightIso2photons == 2));


  //prompt SRs
  //tight iso SR
  selPhotons.fillBranch("passNPhoGe1SelectionPromptTightIsoSR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons > 0));
  selPhotons.fillBranch("passNPhoEq1SelectionPromptTightIsoSR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 1));
  selPhotons.fillBranch("passNPhoEq2SelectionPromptTightIsoSR", bool(nBHphotons == 0 && eq0nonprompt && nTightIsophotons == 2));




  //---------------------------------------------------------------------------
  //  saving photon information for rjr
  //---------------------------------------------------------------------------

  float lRJRPhoIndx = -1;
  float slRJRPhoIndx = -1;
  std::vector<float> rjrpho_pt;
  std::vector<float> rjrpho_eta;
  std::vector<float> rjrpho_phi;
  std::vector<float> rjrpho_Mx;
  std::vector<float> rjrpho_My;
  for( int spidx = 0; spidx < nBaseLinePhotons; spidx++ ){

	//if( phoOrderId[spidx] == objID::Base ) continue;
	uInt pIdx = phoOrderIndx[spidx]; 
	uInt exIdx = phoExcIndx[spidx];
	bool isOOT = (*Photon_isOot)[pIdx];
	float lPhoPt = (*Photon_pt)[pIdx];
	float lPhoPhi = (*Photon_phi)[pIdx];
	float lPhoEta = (*Photon_eta)[pIdx];
	float ePhoMx = 0;
    float ePhoMy = 0;
	if( isOOT ){ 
		ePhoMx -= lPhoPt*std::cos(lPhoPhi);
		ePhoMy -= lPhoPt*std::sin(lPhoPhi);
		if( exIdx > -1 ){ 
    		float ePhoPt = (*Photon_pt)[exIdx];
    		float ePhoPhi = (*Photon_phi)[exIdx];
			ePhoMx += ePhoPt*std::cos(ePhoPhi);
			ePhoMy += ePhoPt*std::sin(ePhoPhi);
		}//<<>>if( exIdx > -1 )
	}//<<>>if( isOOT && exIdx > -1 )
	bool isSig = isRJRPho[spidx];
	if( isSig ){ if( lRJRPhoIndx < 0 ) lRJRPhoIndx = spidx; else if( slRJRPhoIndx < 0 ) slRJRPhoIndx = spidx; }
	//std::cout << " -- is sig pho " << isSig << " index " << spidx << std::endl;
	rjrpho_pt.push_back(lPhoPt);
    rjrpho_eta.push_back(lPhoEta);
    rjrpho_phi.push_back(lPhoPhi);
    rjrpho_Mx.push_back(ePhoMx);
    rjrpho_My.push_back(ePhoMy);

  }//<<>>if(  > 0 )
  geVects.set( "selPhoPt", rjrpho_pt ); // selected for RJR 
  geVects.set( "selPhoEta", rjrpho_eta );
  geVects.set( "selPhoPhi", rjrpho_phi );
  geVects.set( "selPhoEMx", rjrpho_Mx );
  geVects.set( "selPhoEMy", rjrpho_My );
  //std::cout << " -- lead : " << lRJRPhoIndx << " sub " << slRJRPhoIndx << " --------------------- " << std::endl;
  geVars.set( "lSigPhoIndx", lRJRPhoIndx );
  geVars.set( "slSigPhoIndx", slRJRPhoIndx );
  //std::cout << " -- is sig pho done --------------------- " << std::endl;

}//<<>>void KUCMSAodSkimmer::processPhoton(){

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setPhotonBranches( TTree* fOutTree ){

  selPhotons.makeBranch( "nPhotons", UINT );
  selPhotons.makeBranch( "nBaseLinePhotons", UINT ); 
  selPhotons.makeBranch( "nblBHPhotons", UINT );
  selPhotons.makeBranch( "nblNotBHPhotons", UINT );
  selPhotons.makeBranch( "nblTightIsoPhotons", UINT );
  selPhotons.makeBranch( "nblLooseNotTightIsoPhotons", UINT );
  selPhotons.makeBranch( "nblNonpromptPhotons", UINT );
  selPhotons.makeBranch( "nblPromptPhotons", UINT );
  selPhotons.makeBranch( "baseLinePhoton_EleVeto", UINT );
  selPhotons.makeBranch( "baseLinePhoton_OOT", VBOOL );
  selPhotons.makeBranch( "baseLinePhoton_SusyId", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Quality", VINT ); 
  selPhotons.makeBranch( "baseLinePhoton_Time", VFLOAT ); 
 
  selPhotons.makeBranch( "baseLinePhoton_RawTime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_WTimeSig", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_LTimeSig", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_STimeSig", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_WTimeSigOld", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_WTime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_WTimeOld", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_LTime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_STime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_WTRes", VFLOAT ); 
  selPhotons.makeBranch( "baseLinePhoton_WTResOld", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_LTRes", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_STRes", VFLOAT );

  selPhotons.makeBranch( "baseLinePhoton_ShasGS", VBOOL );
  selPhotons.makeBranch( "baseLinePhoton_LSCross", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SSCross", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_LWeird", VBOOL );
  selPhotons.makeBranch( "baseLinePhoton_SWeird", VBOOL );

  selPhotons.makeBranch( "baseLinePhoton_Energy", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Eta", VFLOAT ); 
  selPhotons.makeBranch( "baseLinePhoton_Phi", VFLOAT );     
  selPhotons.makeBranch( "baseLinePhoton_Pt", VFLOAT ); 
  selPhotons.makeBranch( "baseLinePhoton_SCx", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SCy", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SCz", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_CorEnergy", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_CorPt", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_MBetaEql", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_MBetaPmt", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_PVtoRH", VFLOAT );

  selPhotons.makeBranch( "baseLinePhoton_HcalTowerSumEtBcConeDR04", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_TrkSumPtHollowConeDR03", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_TrkSumPtHollowConeDR04", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_TrkSumPtSolidConeDR04", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_PixelSeed", VBOOL );
  selPhotons.makeBranch( "baseLinePhoton_EcalRHSumEtConeDR04", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_HadTowOverEM", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_HadOverEM", VFLOAT );

  selPhotons.makeBranch( "baseLinePhoton_HcalTowerSumEtConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_NTrkHollowConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_NTrkSolidConeDR04", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_EcalPFClusterIso", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_HasConversionTracks", VBOOL );   //!
  selPhotons.makeBranch( "baseLinePhoton_HcalPFClusterIso", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_Hoe_PUcorr", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfChargedIso", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfChargedIsoPFPV", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfChargedIsoWorstVtx", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfPhoIso03", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfRelIso03_all_quadratic", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_PfRelIso03_chg_quadratic", VFLOAT );   //!

  selPhotons.makeBranch( "baseLinePhoton_SMaj", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SMin", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SAlp", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_ClstrRn", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_R9", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Nrh", VUINT );
  selPhotons.makeBranch( "baseLinePhoton_CovEtaEta", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_CovEtaPhi", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_CovPhiPhi", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Sieie", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Sieip", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_Sipip", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_EtaWidth", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_PhiWidth", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_S4", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_SigmaIEtaIEta", VFLOAT );   //!
	 
  selPhotons.makeBranch("baseLinePhoton_isoANNScore", VFLOAT);
  selPhotons.makeBranch("baseLinePhoton_nonIsoANNScore", VFLOAT);
  selPhotons.makeBranch("baseLinePhoton_physBkgCNNScore", VFLOAT);
  selPhotons.makeBranch("baseLinePhoton_beamHaloCNNScore", VFLOAT);
  selPhotons.makeBranch("baseLinePhoton_beamHaloID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_notBeamHaloID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_tightIsoID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_looseNotTightIsoID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_looseNotTightIso1ID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_looseNotTightIso2ID",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_beamHaloCR",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_spikeCR",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_GJetsCR",VBOOL);
  selPhotons.makeBranch("baseLinePhoton_DiJetsCR",VBOOL);

  selPhotons.makeBranch("passNPhoGe1SelectionBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionEarlyBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionEarlyBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionEarlyBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionLateBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionLateBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionLateBeamHaloCR", BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionNotBHPrompt", BOOL);
  selPhotons.makeBranch("tagged_lead_timesig", FLOAT);
  selPhotons.makeBranch("passNPhoGe1SelectionNotBHCR", BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionEarlyNotBHCR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionEarlyNotBHCR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionEarlyNotBHCR", BOOL);
  
  selPhotons.makeBranch("passNPhoGe1SelectionLateNotBHSR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionLateNotBHSR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionLateNotBHSR", BOOL);
  ////prompt CRs
  selPhotons.makeBranch("passNPhoGe1SelectionPromptLooseNotTightIsoCR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionPromptLooseNotTightIsoCR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionPromptLooseNotTightIsoCR", BOOL);

  //prompt SRs
  selPhotons.makeBranch("passNPhoGe1SelectionPromptTightIsoSR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionPromptTightIsoSR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionPromptTightIsoSR", BOOL);
  //validation regions (split lnt iso CR into two)
  //validation (split) loose-not-tight iso CR
  selPhotons.makeBranch("passNPhoGe1SelectionPromptLooseNotTightIso1CR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionPromptLooseNotTightIso1CR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionPromptLooseNotTightIso1CR", BOOL);
  //validation (split) (loose-not-tight) proxy iso SR (CR)
  selPhotons.makeBranch("passNPhoGe1SelectionPromptLooseNotTightIso2CR", BOOL);
  selPhotons.makeBranch("passNPhoEq1SelectionPromptLooseNotTightIso2CR", BOOL);
  selPhotons.makeBranch("passNPhoEq2SelectionPromptLooseNotTightIso2CR", BOOL);
  
  selPhotons.makeBranch( "baseLinePhoton_GenPt", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_PhoIsoDr", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_GenIdx", VINT );

  selPhotons.makeBranch( "baseLinePhoton_GenLabTime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_GenLabTimeSig", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_GenTime", VFLOAT );
  selPhotons.makeBranch( "baseLinePhoton_GenTimeSig", VFLOAT );

  selPhotons.makeBranch( "baseLinePhoton_GenSigMomEnergy", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomEta", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomMass", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomPhi", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomPt", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomPx", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomPy", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomPz", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomVx", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomVy", VFLOAT );   //!
  selPhotons.makeBranch( "baseLinePhoton_GenSigMomVz", VFLOAT );   //!

  selPhotons.makeBranch( "endcap_photon_baseline", VBOOL );   //!
  selPhotons.makeBranch( "endcap_photon_energy", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_eta", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_phi", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_pt", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_tspscdr4", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_erhsecdr4", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_htoem", VFLOAT );   //!
  selPhotons.makeBranch( "endcap_photon_isoANNScore", VFLOAT );   //!

  selPhotons.makeBranch( "barrel_photon_baseline", VBOOL );   //!
  selPhotons.makeBranch( "barrel_photon_energy", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_eta", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_phi", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_pt", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_tspscdr4", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_erhsecdr4", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_htoem", VFLOAT );   //!
  selPhotons.makeBranch( "barrel_photon_isoANNScore", VFLOAT );   //!

  selPhotons.makeBranch( "photon_baseline", VBOOL );   //!
  selPhotons.makeBranch( "photon_LeadTime", VFLOAT );
  selPhotons.makeBranch( "photon_LeadTimeSig", VFLOAT );
  selPhotons.makeBranch( "photon_WTime", VFLOAT );
  selPhotons.makeBranch( "photon_WTimeSig", VFLOAT );
  selPhotons.makeBranch( "photon_Eta", VFLOAT );
  selPhotons.makeBranch( "photon_Phi", VFLOAT );
  selPhotons.makeBranch( "photon_Pt", VFLOAT );
  selPhotons.makeBranch( "photon_E", VFLOAT );
  selPhotons.makeBranch( "photon_PixSeed", VBOOL );
  selPhotons.makeBranch( "photon_isoANNScore", VFLOAT );
  selPhotons.makeBranch( "photon_beamHaloCNNScore", VFLOAT );

  // add new photon branches above
  selPhotons.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::TrackMatched(int scidx, double& bestdr){

    //do track matching
    bestdr = 999;
    double bestTrackDr = 999;
    //double maxTrackDr;
    double dr, teta, tphi, de;
    unsigned int detid;

    double dphi = -999;
    //double bestde_dr;
    int trackidx;
    float ec = (*SuperCluster_eta)[scidx];
    float pc = (*SuperCluster_phi)[scidx];

    double track_phi = -999;
    double track_eta = -999;

    //loop through tracks to get best match to this subcluster (tracks are matched to superclusters, if not idx < 0)
    //Track_scIndexs[i][j] is for track i that matched to supercluster j (tracks can match to multiple SCs)
    int nTracks = (*Track_scIndexs).size();
    for(int t = 0; t < nTracks; t++){
            if((*Track_scIndexs)[t][0] < 0) continue; //not matched to any SC
            int nSCs = (*Track_scIndexs)[t].size();
            //cout << "track #" << t << " matched to " << nSCs << " superclusters" << endl;
            for(int sc = 0; sc < nSCs; sc++){
                    int sc_idx = (*Track_scIndexs)[t].at(sc);
                    //get eta, phi of supercluster sc that track is matched to
                    double sc_eta = (*SuperCluster_eta)[sc_idx];
                    double sc_phi = (*SuperCluster_phi)[sc_idx];
                    double sc_phi_02pi = sc_phi;
                    if(sc_phi_02pi < 0) sc_phi_02pi += 2*acos(-1);
                    else if(sc_phi_02pi > 2*acos(-1)) sc_phi_02pi -= 2*acos(-1);
                    else sc_phi_02pi = sc_phi;
                    //cout << "   track #" << t << " matched to sc " << sc << " with eta " << sc_eta << " and phi " << sc_phi << endl;

                    dphi = fabs(pc - sc_phi_02pi);
                    dphi = acos(cos(dphi));
                    dr = sqrt((sc_eta - ec)*(sc_eta - ec) + dphi*dphi);
                    //E = p for photons
                    if(dr < bestTrackDr){
                            bestTrackDr = dr;
                            track_phi = sc_phi;
                            track_eta = sc_eta;
                    }//<<>>if(dr < bestTrackDr)
            }//<<>>for(int sc = 0; sc < nSCs; sc++)
            //cout << "best dr " << bestdr << " " << bestTrackDr << " track eta " << track_eta; 
            //cout << " SC eta " << ec << " track phi " << track_phi << " SC phi " << pc << endl;
            bestdr = bestTrackDr;
    }//<<>>for(int t = 0; t < nTracks; t++)

}//<<>>void KUCMSAodSkimmer::TrackMatched(int scidx, double& bestdr)

//functions needed: dR track matching
int KUCMSAodSkimmer::GetDetBkgCR(int phoidx){

    int scidx = (*Photon_scIndex)[phoidx];
    //tc = seed time
    float tc = (*SuperCluster_SeedTime)[scidx];
    //pc = phi center
    float pc = (*SuperCluster_phi)[scidx];
    //bestdr = dr of best match to track
    double bestdr;
    TrackMatched(scidx, bestdr);
    //tsig = weighted time significance
    float timesignum, timesigdenom;
    float tsig = getTimeSig( scidx, timesignum, timesigdenom );
    //isolation criteria
    float trksum =     (*Photon_trkSumPtSolidConeDR04)[phoidx] < 6.0;
    float ecalrhsum =  (*Photon_ecalRHSumEtConeDR04)[phoidx] < 10.0;
    float htowoverem = (*Photon_hadTowOverEM)[phoidx] < 0.02;
    bool iso = trksum && ecalrhsum && htowoverem;
    bool pcFilter = (pc < 0.1 || (acos(-1) - 0.1 < pc && pc < acos(-1) + 0.1) || 2*acos(-1) - 0.1 < pc );
    bool pcFilter_wide = (pc < 0.3 || (acos(-1) - 0.3 < pc && pc < acos(-1) + 0.3) || 2*acos(-1) - 0.3 < pc );
    bool spikeTime = (tc <= -10 && tc != -999);
    bool bhTime = (-7 < tc && tc <= -2);
    //bool physBkgTime = (-0.5 < tc && tc <= 0.5);
    bool detBkgTimeSig = (tsig < -3);
    bool physBkgTimeSig = (-1 < tsig && tsig < 1);
    bool spikeTrackMatch = (bestdr <= 0.02);
    bool notSpikeTrackVeto = (bestdr > 0.03);
    //spike
    if(spikeTrackMatch && spikeTime && !pcFilter_wide && detBkgTimeSig && iso) return 3;
    //BH
    else if(bhTime && pcFilter && iso && notSpikeTrackVeto && detBkgTimeSig) return 2;
    else return -1;

}//<<>>int KUCMSAodSkimmer::GetDetBkgCR(int phoidx)

//functions needed: dR calculation
bool KUCMSAodSkimmer::GetGJetsCR(int phoidx){

    //event filters applied
    //not in BH sideband (equivalent of BH filter from noise filters)

    float pho_pt = (*Photon_pt)[phoidx];
    float pho_phi = (*Photon_phi)[phoidx];
    float pho_eta = (*Photon_eta)[phoidx];
    //GJets
    //# jets > 1
    int nJets = Jet_energy->size();
    if(nJets < 1) return false;
    //ht > 50
    double ht = 0;
    double j_px = 0;
    double j_py = 0;
    double j_pz = 0;
	float minda = 100;
	int thejet = -1;
    for(int j = 0; j < nJets; j++){
        ht += (*Jet_pt)[j];
        double pt = (*Jet_pt)[j];
        double phi = (*Jet_phi)[j];
        double eta = (*Jet_eta)[j];
        j_px += pt*cos(phi);
        j_py += pt*sin(phi);
        j_pz += pt*sinh(eta);
    	//dphi_jets_photon > pi-0.3
    	double dphi_phojet = pho_phi - phi;
    	dphi_phojet = acos(cos(dphi_phojet));
		float da = std::abs(M_PI - dphi_phojet);
    	if( da < minda ){ minda = da; thejet = j; };
    }//<<>>for(int j = 0; j < nJets; j++)
    double jet_sys_phi = atan2(j_py, j_px);
    if(ht < 50) return false;
    //met < 100
    if(Met_pt > 150) return false; //pass ntuple selection of SVIPM100, originally using MET < 100
    //dphi_jets_photon > pi-0.3
    double dphi_objjet = pho_phi - jet_sys_phi;
    dphi_objjet = acos(cos(dphi_objjet));
    if(dphi_objjet < 4*atan(1) - 0.3) return false;
    float jet_sys_pt = sqrt(j_px*j_px + j_py*j_py);
    float ptasym = std::min(pho_pt, jet_sys_pt) / std::max(pho_pt, jet_sys_pt);
    if(ptasym < 0.6) return false;
	if( gammaJetIndex[0] < 0 ){ gammaJetIndex[0] = phoidx; gammaJetIndex[1] = thejet; }
    return true;

}//<<>>bool KUCMSAodSkimmer::GetGJetsCR(int phoidx)

bool KUCMSAodSkimmer::GetDiJetsCR(int phoidx){

    //event filters applied
    //not in BH sideband (equivalent of BH filter from noise filters)
    //DiJets
    float pho_phi = (*Photon_phi)[phoidx];
    float pho_eta = (*Photon_eta)[phoidx];
    //# jets > 1
    int nJets = Jet_energy->size();
    if(nJets < 2) return false;
    //ht > 50
    double ht = 0;
    double j_px = 0;
    double j_py = 0;
    double j_pz = 0;
    //assuming jets are pt sorted in ntuple - take phi of leading jets
    float j_phi1 = (*Jet_phi)[0];
    float j_phi2 = (*Jet_phi)[1];
    float j_eta1 = (*Jet_eta)[0];
    float j_eta2 = (*Jet_eta)[1];
    for(int j = 0; j < nJets; j++){
        ht += (*Jet_pt)[j];
        double pt = (*Jet_pt)[j];
        double phi = (*Jet_phi)[j];
        double eta = (*Jet_eta)[j];
        j_px += pt*cos(phi);
        j_py += pt*sin(phi);
        j_pz += pt*sinh(eta);
    }//<<>>for(int j = 0; j < nJets; j++)
    if(ht < 50) return false;
    //met < 100
    if(Met_pt > 150) return false; //pass ntuple selection of SVIPM100, originally using MET < 100
    //dphi_jets > pi-0.3
    double dphi_jets = j_phi1 - j_phi2;
    dphi_jets = acos(cos(dphi_jets));
    if(dphi_jets < 4*atan(1) - 0.3) return false;

    double ptasym = (*Jet_pt)[1] / (*Jet_pt)[0];
    if(ptasym < 0.6) return false;
    //min(pho_pt, jet_sys_pt) / max(pho_pt, jet_sys_pt) > 0.6
    //dr_pho_jet > 0.5 for all photons in event
    double pho_jet1 = dR1(pho_eta, pho_phi, j_eta1, j_phi1);
    double pho_jet2 = dR1(pho_eta, pho_phi, j_eta2, j_phi2);
    //double maxdR = 0.5;
    if(pho_jet1 < 0.5 || pho_jet2 < 0.5) return false;
    if( diJetIndex[0] < 0 ){ diJetIndex[0] = 0; diJetIndex[1] = 1; }
    return true;

}//<<>>bool KUCMSAodSkimmer::GetDiJetsCR(int phoidx)

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



