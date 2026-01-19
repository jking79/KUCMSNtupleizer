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
	//cout << "   track #" << t << " matched to supercluster " << sc << " with eta " << sc_eta << " and phi " << sc_phi << endl;
	
	                dphi = fabs(pc - sc_phi_02pi);
	                dphi = acos(cos(dphi));
	
	                dr = sqrt((sc_eta - ec)*(sc_eta - ec) + dphi*dphi);
	
	
	                //E = p for photons
	                if(dr < bestTrackDr){
	                        bestTrackDr = dr;
	                        track_phi = sc_phi;
	                        track_eta = sc_eta;
	                }
	
	        }
	        //cout << "best dr " << bestdr << " " << bestTrackDr << " track eta " << track_eta << " SC eta " << ec << " track phi " << track_phi << " SC phi " << pc << endl;
	        bestdr = bestTrackDr;
	}

}


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
	if(spikeTrackMatch && spikeTime && !pcFilter_wide && detBkgTimeSig && iso)
		return 3;
	//BH
	else if(bhTime && pcFilter && iso && notSpikeTrackVeto && detBkgTimeSig)
		return 2;
	else return -1;
}


//functions needed: dR calculation
bool KUCMSAodSkimmer::GetGJetsCR(int phoidx){
	//event filters applied
	//not in BH sideband (equivalent of BH filter from noise filters)

	float pho_pt = (*Photon_pt)[phoidx];
	float pho_phi = (*Photon_phi)[phoidx];
	//GJets
	//# jets > 1
	int nJets = Jet_energy->size();
	if(nJets < 1) return false;
	//ht > 50
	double ht = 0;
	double j_px = 0;
	double j_py = 0;
	double j_pz = 0;
	for(int j = 0; j < nJets; j++){
		ht += (*Jet_pt)[j];
		double pt = (*Jet_pt)[j];
	       	double phi = (*Jet_phi)[j];
		double eta = (*Jet_eta)[j];
		j_px += pt*cos(phi);
		j_py += pt*sin(phi);
		j_pz += pt*sinh(eta);
	}
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
	return true;
}

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
	}
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
        double maxdR = 0.5;
        if(pho_jet1 < 0.5 || pho_jet2 < 0.5) return false;
	return true;
}


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
  //std::vector<int> phoOrderId;
  std::vector<int> phoExcIndx;
  uInt nSelPhotons = 0; 
  uInt nBasePhos = 0;
  uInt nLoosePhos = 0;
  uInt nTightPhos = 0;
  uInt nEBPhos = 0;
  uInt nPBPhos = 0;
  uInt nLBPhos = 0;
  uInt nELPhos = 0;
  uInt nPLPhos = 0;
  uInt nLLPhos = 0;
  uInt nETPhos = 0;
  uInt nPTPhos = 0;
  uInt nLTPhos = 0;
  uInt nPhotons = Photon_excluded->size();
  //Final State Counters
  int nPhoBH(0), nPhoEarlyBH(0), nPhoLateBH(0), nPhoEarlyPB(0), nPhoPB(0), nPhoLatePB(0);
  int nPhoEndcapNonIso(0);
  //Discriminant Score Cut Values
  float bhCutVal = 0.917252;
  float pbCutVal = 0.81476355;
  float nonIsoCutVal = 0.99142313;
  //Time window cut values
  float earlyTimeCut = -2;
  float lateTimeCut = 2;
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

    auto energy = (*Photon_energy)[it];
    auto pt = (*Photon_pt)[it];
    bool underMinPtEB = pt < 30;
    bool underMinPtEE = pt < 30;
    auto eta = (*Photon_eta)[it];
    auto overMaxEta = std::abs(eta) > 1.479;
    auto phi = (*Photon_phi)[it];

    if( DEBUG ) std::cout << " -- looping photons : getting pho iso " << std::endl;
    auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
    bool passHcalSum = true;
    auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
    bool passTrkSum = tspscdr4/pt < 0.12; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
    auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
    bool passsEcalRhSum = erhsecdr4/pt < 0.2;
    auto htoem = (*Photon_hadTowOverEM)[it];
    bool passHOE = htoem < 0.02;
    bool isoPho = passHOE && passsEcalRhSum && passTrkSum && passHcalSum;
    bool failPhoIsoEB = not isoPho;
	bool failPhoIsoEE = failPhoIsoEB; //not ( tspscdr4 < 7.0 && erhsecdr4 < 12.5 && htoem < 0.05 );

    //  change to skip jets and keep all photons regardless of photon iso with jet
    //bool phoskip = isExcluded || hasPixSeed || overMaxEta || underMinPt || failPhoIso;
	bool phoskip = isExcluded || hasPixSeed;
	//if( not overMaxEta ) phoskip = phoskip || failPhoIsoEB || underMinPtEB;
	//if( overMaxEta ) phoskip = phoskip || failPhoIsoEE || underMinPtEE;
	//removing isolation selection from 'selected' photon requirements to apply an inversion for a sideband
	if( not overMaxEta ) phoskip = phoskip || underMinPtEB;
	if( overMaxEta ) phoskip = phoskip || underMinPtEE;
	//if( overMaxEta ) phoskip = true;

	bool hemEligible1 = pt > 20 && not isExcluded && not hasPixSeed;
    bool hemEligible2 = pt > 30 && not isExcluded && not hasPixSeed; 

	bool isInHemRegion = inHEMRegion( eta, phi );
	hemBits.set( "pho1hvl", isInHemRegion && hemEligible1 );		
    hemBits.set( "pho2hvm", isInHemRegion && hemEligible2 );

	bool isEESig = not isExcluded && not hasPixSeed && overMaxEta;
    if( isEESig ){

		selPhotons.fillBranch( "EESigPho_skipped", phoskip );   //!
    	selPhotons.fillBranch( "EESigPho_energy", (*Photon_energy)[it] );   //!
    	selPhotons.fillBranch( "EESigPho_eta", (*Photon_eta)[it] );   //!
    	selPhotons.fillBranch( "EESigPho_phi", (*Photon_phi)[it] );   //!
    	selPhotons.fillBranch( "EESigPho_tspscdr4", (*Photon_trkSumPtSolidConeDR04)[it] );   //!
        selPhotons.fillBranch( "EESigPho_erhsecdr4", (*Photon_ecalRHSumEtConeDR04)[it] );   //!
        selPhotons.fillBranch( "EESigPho_htoem", (*Photon_hadTowOverEM)[it] );   //!
        selPhotons.fillBranch( "EESigPho_pt", (*Photon_pt)[it] );   //!

    }//if( doGenInfo )

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
    // -------   check isTimeValid on rechit?  -------------------------------------
    for( int sciter = 0; sciter < nrh; sciter++  ){
        uInt pscrhid = rhids[sciter];
        int erhiter = ( rhIDtoIterMap.find(pscrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[pscrhid] : -1;
        if( erhiter != -1 ){
            bool isValid = (*ECALRecHit_isTimeValid)[erhiter];
            float swcrss = (*ECALRecHit_swCross)[erhiter];
            bool isWeird = (*ECALRecHit_isWrd)[erhiter] || (*ECALRecHit_isDiWrd)[erhiter];
            float gainwt = 1;
			float gainwt1 = 1;
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
			if( (*ECALRecHit_hasGS1)[erhiter] ) gainwt1 = 0;
            if( not isValid ){ gainwt = 0; gainwt1 = 0; }
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
                leadTres = ertres; leadSX = swcrss; leadWried = isWeird;
            }//<<>>if( erhe*gainwt > leadE )
            if( erhe > seedE ){
                seedE = erhe; seedTime = ertoftime; seedAres = erampres; seedRHID = pscrhid;
                seedTres = ertres; seedSX = swcrss; seedWried = isWeird; seedGS = hasGainSwitch;
            }//<<>>if( erhe > seedE ) 
			bool hasBadTime = hasGainSwitch || !(*ECALRecHit_isTimeValid)[erhiter];
			//cout << "adding rh with e " << erhe << endl;
			_ca.AddRecHit(erx, ery, erz, erhe, erhct, pscrhid, hasBadTime);
        }//<<>>if( scrhid == rhid )
    }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
    if( sumw == 0 ){ sumw = 1; sumtw = -100; sumtrw = -1000; }
    float phoWTime = sumtw/sumw;
    float phoWRes = sumtrw/sumw;
    float phoWTime1 = sumtw1/sumw1;
    float phoWRes1 = sumtrw1/sumw1;
    //float ltimeres = timeCali->getTimeResoltuion( leadAres, leadRHID, Evt_run, tctag, mctype );
    //if( ltimeres != leadTres ) std::cout << " !!!!!!!   lead res mis : " << ltimeres << " v " << leadTres << std::endl;
    //float stimeres = timeCali->getTimeResoltuion( seedAres, seedRHID, Evt_run, tctag, mctype );
    //if( stimeres != seedTres ) std::cout << " !!!!!!!   seed res mis : " << stimeres << " v " << seedTres << std::endl;
    float leadtimesig = leadTime/std::sqrt(leadTres);
    float seedtimesig = seedTime/std::sqrt(seedTres);
    float wttimesig = phoWTime/std::sqrt(phoWRes);
    float wttimesig1 = phoWTime1/std::sqrt(phoWRes1);

    float timesignum, timesigdenom;
	float phoWTimeSig = getTimeSig( scIndx, timesignum, timesigdenom );

	if( isExcluded ) continue;

    //selPhotons.fillBranch( "allPhoWTime1", phoWTime1 );
    selPhotons.fillBranch( "allPhoWTimeSig", phoWTimeSig );
    selPhotons.fillBranch( "allPhoWTime", phoWTime );
    selPhotons.fillBranch( "allPhoWTimeSig1", wttimesig1 );
    selPhotons.fillBranch( "allPhoEta", eta );
    selPhotons.fillBranch( "allPhoPhi", phi );
    selPhotons.fillBranch( "allPhoPt", pt );
    selPhotons.fillBranch( "allPhoE", energy );
    selPhotons.fillBranch( "allPhoPixSeed", hasPixSeed );

    if( ( geVars("genSigPerfect") != 1 ) && phoskip ){ isSelPho.push_back(false); continue; }		
    if( ( geVars("genSigPerfect") == 1 ) &&  ( not isGenSig ) ){ isSelPho.push_back(false);  continue; }
	isSelPho.push_back(true);

    ///////////  pho selection ////////////////////////////////////////////////////////////////////

// for the tight and loose(!tight) selections the discriminant scores are as follows (with the corresponding ROC curve info):
// - old  tight: iso score > 0.337 (5% noniso contamination with ~97% iso efficiency)
// - old loose: iso score > 0.134 (10% noniso contamination with ~98% iso efficiency)
// - old very loose: iso score > 0.0364 (20% noniso contamination with ~99% iso efficiency)
// - tight: 0.875 (1% bkg eff for 90% sig eff)
// - loose: 0.014 (30% bkg eff for 99.8% sig eff)
//

	float isobkg_score = 0;
	float nonisobkg_score = 0;
	float physbkg_score = 0;
	float bh_score = 0;
	bool gjetscr = false;
	bool dijetscr = false;
	bool spikecr = false;
	bool bhcr = false;
    if(_ca.GetNRecHits() > 2){
    	ClusterObj phoobj;
        _ca.NoClusterRhs(phoobj, true);
     	map<string, double> isomap;
        MakePhotonIsoMap(it, isomap);
        vector<float> photonIDscores;
	if(overMaxEta){ //endcap
		phoobj.CalculateEndcapPhotonIDScores(pt, isomap);
	}
	else{ //barrel
		//phoobj.CalculateBarrelPhotonIDScores(pt, isomap);
		//under development - fill with -999
		photonIDscores.push_back(-999);
		photonIDscores.push_back(-999);
	}
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

	gjetscr = GetDiJetsCR(it);
	dijetscr = GetGJetsCR(it);
	int detbkgcr = GetDetBkgCR(it);
	if(detbkgcr == 2)
		bhcr = true;
	if(detbkgcr == 3)
		spikecr = true;

	//values that provide a 0.1% FPR with a >98% efficiency for the relevant class
	//beam halo > 0.917252
	//phys bkg > 0.81476355 (derived from gogoG ROC)
	if(bh_score > bhCutVal){
		nPhoBH++;
		if(phoWTime < earlyTimeCut)
			nPhoEarlyBH++;
		if(phoWTime > lateTimeCut)
			nPhoLateBH++;
	}
	if(physbkg_score > pbCutVal){
		nPhoPB++;
		if(phoWTime < earlyTimeCut)
			nPhoEarlyPB++;
		if(phoWTime > lateTimeCut)
			nPhoLatePB++;
	}
	//endcap noniso selection
	if(overMaxEta && nonisobkg_score > nonIsoCutVal){
		nPhoEndcapNonIso++;
	}

    	if(DEBUG) cout << "photon " << it << " iso score " << isobkg_score << " noniso score " << nonisobkg_score << endl;
	}//<<>>if(_ca.GetNRecHits() > 2)
    //do photon id score
 
  	_ca.ClearRecHitList();

	bool isTight = isobkg_score > 0.875;
	bool isLoose = isobkg_score > 0.014;
    int phoTSigId = ( phoWTimeSig <= -2.5 ) ? tSigID::Early : ( phoWTimeSig < 2.5 ) ? tSigID::Prompt : tSigID::Late;
    int phoObjId =  ( isTight ) ? objID::Tight : ( isLoose ) ? objID::Loose : objID::Base;

	if( phoObjId == objID::Base ){ 
		nBasePhos++;
        if( phoTSigId == tSigID::Early ) nEBPhos++;
        else if( phoTSigId == tSigID::Prompt ) nPBPhos++;
        else if( phoTSigId == tSigID::Late ) nLBPhos++;
    }//<<>>else if( phoObjId == objID::Base )
	else if( phoObjId == objID::Loose ){
		nLoosePhos++;
    	if( phoTSigId == tSigID::Early ) nELPhos++;
    	else if( phoTSigId == tSigID::Prompt ) nPLPhos++;
    	else if( phoTSigId == tSigID::Late ) nLLPhos++;
	}//<<>>else if( phoObjId == objID::Loose )
	else if( phoObjId == objID::Tight ){ 
		nTightPhos++;
        if( phoTSigId == tSigID::Early ) nETPhos++;
        else if( phoTSigId == tSigID::Prompt ) nPTPhos++;
        else if( phoTSigId == tSigID::Late ) nLTPhos++;
    }//<<>>else if( phoObjId == objID::Tight )

    selPhotons.fillBranch( "selPhoTSigId", phoTSigId );	
    selPhotons.fillBranch( "selPhoObjId", phoObjId );

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

	for( auto as : erhamps ){ hist1d[20]->Fill(as/sumerha); }

    //selPhotons.fillBranch( "selPhoWTime1", phoWTime1 );
    selPhotons.fillBranch( "selPhoWTimeSig", phoWTimeSig );
    selPhotons.fillBranch( "selPhoLTimeSig", leadtimesig );
    selPhotons.fillBranch( "selPhoSTimeSig", seedtimesig );
    selPhotons.fillBranch( "selPhoWTimeSigOld", wttimesig );
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

    selPhotons.fillBranch("selPho_isoANNScore",isobkg_score);
    selPhotons.fillBranch("selPho_nonIsoANNScore",nonisobkg_score);
    selPhotons.fillBranch("selPho_physBkgCNNScore",physbkg_score);
    selPhotons.fillBranch("selPho_beamHaloCNNScore",bh_score);
    selPhotons.fillBranch("selPho_beamHaloCR",bhcr);
    selPhotons.fillBranch("selPho_spikeCR",spikecr);
    selPhotons.fillBranch("selPho_GJetsCR",gjetscr);
    selPhotons.fillBranch("selPho_DiJetsCR",dijetscr);



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
	int phoIsoIndx = -1;
    for( uInt it2 = 0; it2 < nPhotons; it2++ ){
      if( it == it2 ) continue;
      float eta2 = (*Photon_eta)[it2];
      float phi2 = (*Photon_phi)[it2];
      auto dphi12 = dPhi( phi, phi2 );
      auto pho2dr = hypo( eta-eta2, dphi12 ); 
      //std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << pho2dr << std::endl;
      if( pho2dr < phoPhoIsoDr ){ phoPhoIsoDr = pho2dr; phoIsoIndx = it2; }
    }//for( uInt it2 = it+1; it2 < nPhotons; it2++ )
	//bool matchPhoExc = ( phoIsoIndx > -1 ) ? (*Photon_excluded)[it] : false;
	if( phoPhoIsoDr < 0.3 ) phoExcIndx.push_back( phoIsoIndx ); else phoExcIndx.push_back( -1 );

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
    //phoOrderId.push_back(phoObjId);
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
  if( nPhotons == 0 ){ hemBits.set( "pho1hvl", false ); hemBits.set( "pho2hvm", false ); }
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

  bool has1ELPho = ( nELPhos == 1 && nBasePhos < 1 && nTightPhos < 1 ) ? 1 : 0;
  bool has1ETPho = ( nETPhos == 1 && nBasePhos < 1 && nLoosePhos < 1 ) ? 1 : 0;
  bool has1PLPho = ( nPLPhos == 1 && nBasePhos < 1 && nTightPhos < 1 ) ? 1 : 0;
  bool has1PTPho = ( nPTPhos == 1 && nBasePhos < 1 && nLoosePhos < 1 ) ? 1 : 0;
  bool has1LLPho = ( nLLPhos == 1 && nBasePhos < 1 && nTightPhos < 1 ) ? 1 : 0;
  bool has1LTPho = ( nLTPhos == 1 && nBasePhos < 1 && nLoosePhos < 1 ) ? 1 : 0;
  bool has1BasePho = ( nBasePhos > 0 && nLoosePhos < 1 && nTightPhos < 1 ) ? 1 : 0;

  selPhotons.fillBranch( "has1ELPho", has1ELPho );
  selPhotons.fillBranch( "has1ETPho", has1ETPho );
  selPhotons.fillBranch( "has1PLPho", has1PLPho );
  selPhotons.fillBranch( "has1PTPho", has1PTPho );
  selPhotons.fillBranch( "has1LLPho", has1LLPho );
  selPhotons.fillBranch( "has1LTPho", has1LTPho );
  selPhotons.fillBranch( "has1BasePho", has1BasePho );

  selPhotons.fillBranch( "nBasePhotons", nBasePhos );
  selPhotons.fillBranch( "nEBPhotons", nEBPhos );
  selPhotons.fillBranch( "nPBPhotons", nPBPhos );
  selPhotons.fillBranch( "nLBPhotons", nLBPhos );
  selPhotons.fillBranch( "nLoosePhotons", nLoosePhos );
  selPhotons.fillBranch( "nELPhotons", nELPhos );
  selPhotons.fillBranch( "nPLPhotons", nPLPhos );
  selPhotons.fillBranch( "nLLPhotons", nLLPhos );
  selPhotons.fillBranch( "nTightPhotons", nTightPhos );
  selPhotons.fillBranch( "nETPhotons", nETPhos );
  selPhotons.fillBranch( "nPTPhotons", nPTPhos );
  selPhotons.fillBranch( "nLTPhotons", nLTPhos );

  geCnts.set( "nSelPhotons", nSelPhotons );	
  selPhotons.fillBranch( "nSelPhotons",  nSelPhotons );
  selPhotons.fillBranch( "nPhotons", nPhotons );

  bool anyPB = (nPhoPB > 0);
  bool anyBH = (nPhoBH > 0);
  bool anyLatePB = (nPhoLatePB > 0);
  bool anyLateBH = (nPhoLateBH > 0);
  bool anyEarlyBH = (nPhoEarlyBH > 0);
  //do bh CR selection
  selPhotons.fillBranch("passNPhoGe1SelectionBeamHaloCR", bool(nPhoBH>=1 && !anyPB));
  selPhotons.fillBranch("passNPhoGe1SelectionEarlyBeamHaloCR", bool(nPhoEarlyBH>=1 && !anyPB && !anyLateBH));
  selPhotons.fillBranch("passNPhoGe1SelectionLateBeamHaloCR", bool(nPhoLateBH>=1 && !anyPB && !anyEarlyBH));
  selPhotons.fillBranch("passNPhoGe1SelectionEarlyPhysBkgCR", bool(nPhoEarlyPB>=1 && !anyBH && !anyLatePB));
  //do endcap noniso CR selection
  selPhotons.fillBranch("passNPhoGe1SelectionEndcapNonIsoCR", bool(nPhoEndcapNonIso>=1));

  std::vector<float> selpho_pt;
  std::vector<float> selpho_eta;
  std::vector<float> selpho_phi;
  std::vector<float> selpho_Mx;
  std::vector<float> selpho_My;
  for( int spidx = 0; spidx < nSelPhotons; spidx++ ){
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
	selpho_pt.push_back(lPhoPt);
    selpho_eta.push_back(lPhoEta);
    selpho_phi.push_back(lPhoPhi);
    selpho_Mx.push_back(ePhoMx);
    selpho_My.push_back(ePhoMy);
  }//<<>>if( nSelPhotons > 0 )
  geVects.set( "selPhoPt", selpho_pt );
  geVects.set( "selPhoEta", selpho_eta );
  geVects.set( "selPhoPhi", selpho_phi );
  geVects.set( "selPhoEMx", selpho_Mx );
  geVects.set( "selPhoEMy", selpho_My );

}//<<>>void KUCMSAodSkimmer::processPhoton(){

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setPhotonBranches( TTree* fOutTree ){

  selPhotons.makeBranch( "has1ELPho", BOOL );
  selPhotons.makeBranch( "has1ETPho", BOOL );
  selPhotons.makeBranch( "has1PLPho", BOOL );
  selPhotons.makeBranch( "has1PTPho", BOOL );
  selPhotons.makeBranch( "has1LLPho", BOOL );
  selPhotons.makeBranch( "has1LTPho", BOOL );
  selPhotons.makeBranch( "has1BasePho", BOOL );

  selPhotons.makeBranch( "nBasePhotons", UINT );
  selPhotons.makeBranch( "nEBPhotons", UINT );
  selPhotons.makeBranch( "nPBPhotons", UINT );
  selPhotons.makeBranch( "nLBPhotons", UINT );
  selPhotons.makeBranch( "nLoosePhotons", UINT );
  selPhotons.makeBranch( "nELPhotons", UINT );
  selPhotons.makeBranch( "nPLPhotons", UINT );
  selPhotons.makeBranch( "nLLPhotons", UINT );
  selPhotons.makeBranch( "nTightPhotons", UINT );
  selPhotons.makeBranch( "nETPhotons", UINT );
  selPhotons.makeBranch( "nPTPhotons", UINT );
  selPhotons.makeBranch( "nLTPhotons", UINT );

  selPhotons.makeBranch( "nPhotons", UINT );
  selPhotons.makeBranch( "nSelPhotons", UINT ); 
  selPhotons.makeBranch( "selPhoEleVeto", UINT );
  selPhotons.makeBranch( "selPhoOOT", VBOOL );
  selPhotons.makeBranch( "selPhoSusyId", VFLOAT );
  selPhotons.makeBranch( "selPhoQuality", VINT ); 
  selPhotons.makeBranch( "selPhoTime", VFLOAT ); 
 
  selPhotons.makeBranch( "selPhoTSigId", VINT );     
  selPhotons.makeBranch( "selPhoObjId", VINT );

  //selPhotons.makeBranch( "selPhoWTime1", VFLOAT );
  selPhotons.makeBranch( "selPhoWTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoLTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoSTimeSig", VFLOAT );
  selPhotons.makeBranch( "selPhoWTimeSigOld", VFLOAT );
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
	 
  selPhotons.makeBranch("selPho_isoANNScore", VFLOAT);
  selPhotons.makeBranch("selPho_nonIsoANNScore", VFLOAT);
  selPhotons.makeBranch("selPho_physBkgCNNScore", VFLOAT);
  selPhotons.makeBranch("selPho_beamHaloCNNScore", VFLOAT);
  selPhotons.makeBranch("selPho_beamHaloCR",VBOOL);
  selPhotons.makeBranch("selPho_spikeCR",VBOOL);
  selPhotons.makeBranch("selPho_GJetsCR",VBOOL);
  selPhotons.makeBranch("selPho_DiJetsCR",VBOOL);

  selPhotons.makeBranch("passNPhoGe1SelectionBeamHaloCR",BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionEarlyBeamHaloCR",BOOL);
  selPhotons.makeBranch("passNPhoGe1SelectionLateBeamHaloCR",BOOL); 
  selPhotons.makeBranch("passNPhoGe1SelectionEarlyPhysBkgCR",BOOL); 
  selPhotons.makeBranch("passNPhoGe1SelectionEndcapNonIsoCR",BOOL); 

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

  selPhotons.makeBranch( "EESigPho_skipped", VBOOL );   //!
  selPhotons.makeBranch( "EESigPho_energy", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_eta", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_phi", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_pt", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_tspscdr4", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_erhsecdr4", VFLOAT );   //!
  selPhotons.makeBranch( "EESigPho_htoem", VFLOAT );   //!

  selPhotons.makeBranch( "allPhoWTime1", VFLOAT );
  selPhotons.makeBranch( "allPhoWTimeSig1", VFLOAT );
  selPhotons.makeBranch( "allPhoWTime", VFLOAT );
  selPhotons.makeBranch( "allPhoWTimeSig", VFLOAT );
  selPhotons.makeBranch( "allPhoEta", VFLOAT );
  selPhotons.makeBranch( "allPhoPhi", VFLOAT );
  selPhotons.makeBranch( "allPhoPt", VFLOAT );
  selPhotons.makeBranch( "allPhoE", VFLOAT );
  selPhotons.makeBranch( "allPhoPixSeed", VBOOL );

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



