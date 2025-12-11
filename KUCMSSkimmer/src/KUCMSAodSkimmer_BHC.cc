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

//#define true true
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


void KUCMSAodSkimmer::MakePhotonIsoMap(int phoidx, map<string, double>& isomap){
	isomap.clear();
	isomap["hcalTowerSumEtConeDR04"] = (*Photon_hcalTowerSumEtConeDR04)[phoidx];
	isomap["trkSumPtSolidConeDR04"] = (*Photon_trkSumPtSolidConeDR04)[phoidx];
	isomap["trkSumPtHollowConeDR04"] = (*Photon_trkSumPtHollowConeDR04)[phoidx];
	isomap["hadTowOverEM"] = (*Photon_hadTowOverEM)[phoidx];
	isomap["ecalRHSumEtConeDR04"] = (*Photon_ecalRHSumEtConeDR04)[phoidx];
}

void KUCMSAodSkimmer::processBHCPhotons(){

	if( DEBUG ) std::cout << "Finding BC information for photons" << std::endl;
	BHCPhoInfo.clearBranches();

    int selPhoIndex = -1;

	uInt nPhotons = Photon_excluded->size();
	if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons to get BC info" << std::endl;
	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	_ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	_ca.SetTransferFactor(1/30.); //set to 1/min photon pt
	int selPhoBHCIdx = 0;
	vector<ClusterObj> phoobjs;
	float timesig_num, timesig_denom;
	map<int, int> selPhoIdxToIt;
	for( uInt it = 0; it < nPhotons; it++ ){
		//cout << "nphoton #" << it << endl;
		if( not isSelPho[it] ) continue;
		selPhoIndex++;

		std::map< unsigned int, float > rhtresmap;
		std::map< double, unsigned int > energyToId;
        auto scIndx = (*Photon_scIndex)[it];
        auto rhids = (*SuperCluster_rhIds)[scIndx];
		int nSCRecHits = rhids.size();
        auto nRecHits = ECALRecHit_ID->size();
        for( int sciter = 0; sciter < nSCRecHits; sciter++  ){
			auto scrhid = rhids[sciter];
			int erhiter = ( rhIDtoIterMap.find(scrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[scrhid] : -1;
			if( erhiter != -1 ){
                double erhe = (*ECALRecHit_energy)[erhiter];
                bool hasBadTime = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter] || !(*ECALRecHit_isTimeValid)[erhiter];
				bool isEE = fabs((*ECALRecHit_eta)[erhiter]) > 1.479;
				bool badEETime = (*ECALRecHit_timeError)[erhiter] < 1.005;
				if( isEE && badEETime ) hasBadTime = true;
				float erht = erh_corTime[erhiter];
				float rhx = (*ECALRecHit_rhx)[erhiter];
                float rhy = (*ECALRecHit_rhy)[erhiter];
                float rhz = (*ECALRecHit_rhz)[erhiter];
				rhtresmap[scrhid] = erh_timeRes[erhiter];
				energyToId[erhe] = scrhid;
                if( DEBUG ) std::cout << " " << erht << " " << rhx << " " << rhy << " " << rhz << std::endl;
				//skip endcap rechits
				//if(fabs((*ECALRecHit_eta)[erhiter]) > 1.479) continue;
				//cout << "rh e " <<  erhe << " rh t " << erht << " rhid " << scrhid << " eta " << (*ECALRecHit_eta)[erhiter] << " phi " << (*ECALRecHit_phi)[erhiter] << endl;
				//cout << "rh e " <<  erhe << " rh t " << erht << " x " << rhx << " y " << rhy << " z " << rhz << " rhid " << scrhid << endl;
				_ca.AddRecHit(rhx, rhy, rhz, erhe, erht, scrhid, hasBadTime);
				//cout << "processphotons - adding rh to clustering with id " << scrhid << endl;
			}//<<>>if( ecalrhiter != -1 )
		}//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
	
		//if not enough rechits to cluster, skip this object
		if(_ca.GetNRecHits() < 2) continue;

		//
		//  call to BayesianClustering framework to process photon information 
		//  runs BHC (NlnN) algorithm to cluster photon with subcluster regularization
		//  returns lead (highest energy) cluster found

		//put safety in so that if there are < 2 rhs to cluster, skip

		ClusterObj phoobj;
		//cout << "running photon clustering" << endl;
		int clusterret = _ca.RunClustering(phoobj, true);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		//cout << "selPhoIndex " << selPhoIndex << " it " << it << endl;
		selPhoIdxToIt[selPhoIndex] = it;
		phoobj.SetUserIndex(selPhoIndex);
		phoobjs.push_back(phoobj);
	}
	sort(phoobjs.begin(), phoobjs.end(), ptsort);
	BHCPhoInfo.fillBranch("nSelBHCPhos",(int)phoobjs.size());
	for(int p = 0; p < phoobjs.size(); p++){
		ClusterObj phoobj = phoobjs[p];
		//cout << "phoobj #" << p << " usridx " << phoobj.GetUserIndex() << endl; 
		int phoidx = selPhoIdxToIt.at(phoobj.GetUserIndex());
		//cout << "phoidx " << phoidx << endl;
		//center
		float phoeta = phoobj.GetEtaCenter();
		BHCPhoInfo.fillBranch( "selPhoBHC_eta", phoeta);
		float phophi = phoobj.GetPhiCenter();
		BHCPhoInfo.fillBranch( "selPhoBHC_phi", phophi);
		//time at detector face
		float photime = phoobj.GetObjTime_Det();
		BHCPhoInfo.fillBranch( "selPhoBHC_DetTime", photime);
		//time at PV
		float photime_pv = phoobj.GetObjTime_PV();
		BHCPhoInfo.fillBranch( "selPhoBHC_PVTime", photime_pv);
		//cout << "photon - det time " << photime << " pv time " << photime_pv << endl;	
		//covariance
		float etavar = phoobj.GetEtaVar();
		float phivar = phoobj.GetPhiVar();
		float timevar = phoobj.GetTimeVar();
		float etaphicov = phoobj.GetEtaPhiCov();
		double majlen, minlen;
		majlen = -1;
		minlen = -1;
		phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCPhoInfo.fillBranch( "selPhoBHC_selPhoIndex", phoobj.GetUserIndex() );
		BHCPhoInfo.fillBranch( "selPhoBHC_etaVar", (float)etavar);
		BHCPhoInfo.fillBranch( "selPhoBHC_phiVar", (float)phivar);
		BHCPhoInfo.fillBranch( "selPhoBHC_timeVar", (float)timevar);
		BHCPhoInfo.fillBranch( "selPhoBHC_etaPhiCov", (float)etaphicov);
		BHCPhoInfo.fillBranch( "selPhoBHC_majlen", (float)majlen);
		BHCPhoInfo.fillBranch( "selPhoBHC_minlen", (float)minlen);
        	auto scIndx = (*Photon_scIndex)[phoidx];
		map<unsigned int, float> weights;
		//cout << "getting rh weights - bhc" << endl;
		phoobj.GetRecHitWeights(weights); ///these weights are the PU projected out weights if PU cleaning has been applied (otherwise sum to 1)
		//cout << "doing time sig" << endl;
		float timeSignificance = getTimeSig(scIndx, timesig_num, timesig_denom, weights);
		//cout << "did time sig" << endl;
		//cout << "time sig for photon # " << selPhoBHCIdx << ": " << timeSignificance << endl;
		BHCPhoInfo.fillBranch( "selPhoBHC_timeSignficance", timeSignificance);
		BHCPhoInfo.fillBranch( "selPhoBHC_timeSignficanceNum", timesig_num);
		BHCPhoInfo.fillBranch( "selPhoBHC_timeSignficanceDenom", timesig_denom);
		BHCPhoInfo.fillBranch( "selPhoBHC_energy", (float)phoobj.GetEnergy());
		BHCPhoInfo.fillBranch( "selPhoBHC_pt", (float)phoobj.GetPt());

		//calculate PU scores
		phoobj.CalculatePUScores();
		vector<bool> puscores;
		phoobj.GetPUScores(puscores);
		int nsubclusters = phoobj.GetNSubclusters();
		//cout << "nsubclusters " << nsubclusters << endl;
		BHCPhoInfo.fillBranch("selPhoBHC_nSubclusters",nsubclusters);
		for(int k = 0; k < nsubclusters; k++){
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_etaVar",phoobj.GetSubclusterEtaVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_phiVar",phoobj.GetSubclusterPhiVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_timeVar",phoobj.GetSubclusterTimeVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_energy",phoobj.GetSubclusterEnergy(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_selPhotonBHCIndex",selPhoBHCIdx);
			//cout << "subcluster #" << k << " pu score " << puscores[k] << " selPhoBHCIndex " << selPhoBHCIdx << endl;
			if(puscores.size() == 0){
				if(true) cout << "Run CalculatePUScores() BEFORE this loop" << endl;
				continue;
			}
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_puScore",(int)puscores[k]);
		}


		phoobj.CleanOutPU();
		//do PU-cleaned observables
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_selPhoIndex", selPhoIndex );
		//center
		phoeta = phoobj.GetEtaCenter();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_eta", phoeta);
		phophi = phoobj.GetPhiCenter();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_phi", phophi);
		//time at detector face
		photime = phoobj.GetObjTime_Det();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_DetTime", photime);
		//time at PV
		photime_pv = phoobj.GetObjTime_PV();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_PVTime", photime_pv);
		nsubclusters = phoobj.GetNSubclusters();	
		BHCPhoInfo.fillBranch("selPhoBHCPUCleaned_nSubclusters",nsubclusters);
		BHCPhoInfo.fillBranch("selPhoBHCPUCleaned_energy", (float)phoobj.GetEnergy());
		BHCPhoInfo.fillBranch("selPhoBHCPUCleaned_pt", (float)phoobj.GetPt());
	
		//covariance
		etavar = phoobj.GetEtaVar();
		phivar = phoobj.GetPhiVar();
		etaphicov = phoobj.GetEtaPhiCov();
		phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_etaVar", etavar);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_phiVar", phivar);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_etaPhiCov", etaphicov);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_majlen", (float)majlen);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_minlen", (float)minlen);
		//calculate time significance
		//cout << "getting rh weights" << endl;
		phoobj.GetRecHitWeights(weights); ///these weights are the PU projected out weights if PU cleaning has been applied (otherwise sum to 1)
		//cout << "doing cleaned time sig" << endl;
		timeSignificance = getTimeSig(scIndx, timesig_num, timesig_denom, weights);
		//cout << "did cleaned time sig" << endl;
		//cout << "time sig for photon with pu cleaning # " << selPhoBHCIdx << ": " << timeSignificance << endl;
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficance", timeSignificance);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficanceNum", timesig_num);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficanceDenom", timesig_denom);
		vector<vector<float>> detbkgscores;
		phoobj.CalculateDetBkgScores(true);
		//cout << "calculated det bkg scores" << endl;
		phoobj.GetDetBkgScores(detbkgscores);
		//cout << "got detbkg scores" << endl;
		float physbkg_score = detbkgscores[0][0];
		float bh_score = detbkgscores[0][1];
		float spike_score = detbkgscores[0][2];
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_physBkgCNNScore", physbkg_score);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_beamHaloCNNScore", bh_score);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_spikeCNNScore", spike_score);
		
		//need to set isolation information for the original photon
		map<string, double> isomap;
		//cout << "making isomap" << endl;
		MakePhotonIsoMap(phoidx, isomap);
		//cout << "made photon iso map " << endl;
		phoobj.CalculatePhotonIDScores((*Photon_pt)[phoidx], isomap);
		//cout << "calculated photon id scores" << endl;
		vector<float> photonIDscores;
		phoobj.GetPhotonIDScores(photonIDscores);
		float isobkg_score = photonIDscores[0];
		float nonisobkg_score = photonIDscores[1];
		BHCPhoInfo.fillBranch("selPhoBHCPUCleaned_isoBkgANNScore",isobkg_score);
		BHCPhoInfo.fillBranch("selPhoBHCPUCleaned_nonIsoBkgANNScore",nonisobkg_score);


		selPhoBHCIdx++;
		_ca.ClearRecHitList();
        // branches created strating on 2546 in  void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree )
        // on line 2709 and below make a BHCPhoInfo branch : selPhotons.makeBranch( "yourvarname", VFLOAT );		

	}//<<>>for( uInt it = 0; it < nPhotons; it++ )
}//<<>>void KUCMSAodSkimmer::processMLPhotons()

void KUCMSAodSkimmer::processBHCJets(){

    if( DEBUG ) std::cout << "Finding BC information for jets" << std::endl;
    std::cout << "Finding BC information for jets" << std::endl;
	BHCJetInfo.clearBranches();

	int selJetIndex = -1;

	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	_ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	_ca.SetTransferFactor(1/30.); //set to 1/min jet pt
	int nbhcjets = 0;
    std::vector<Jet> jetrhs;
  	uInt nJets = Jet_energy->size();
	float pvtime = 0;
	float norm = 0;
	float pvtime_PUcleaned = 0;
	float norm_PUcleaned = 0;
	float pvtime_PUdetBkgCleaned = 0;
	float norm_PUdetBkgCleaned = 0;
	int selJetBHCIdx = 0;
	std::vector<ClusterObj, Eigen::aligned_allocator<ClusterObj>> jetobjs;
    if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets to get BC info" << std::endl;
  	for( uInt it = 0; it < nJets; it++ ){
		cout << "jet it " << it << endl;
        if( not isSelJet[it] ) continue;
		selJetIndex++;

    	std::map< unsigned int, float > rhtresmap;
    	auto rhids = (*Jet_drRhIds)[it];
        int nJetRecHits = rhids.size();
	double totE = 0;
        for( int jiter = 0; jiter < nJetRecHits; jiter++  ){
            auto jrhid = rhids[jiter];
            int erhiter = ( rhIDtoIterMap.find(jrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[jrhid] : -1;
            if( erhiter != -1 ){
                float erhe = (*ECALRecHit_energy)[erhiter];
                bool hasBadTime = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter] || !(*ECALRecHit_isTimeValid)[erhiter];
                float erht = erh_corTime[erhiter];
                float rhx = (*ECALRecHit_rhx)[erhiter];
                float rhy = (*ECALRecHit_rhy)[erhiter];
                float rhz = (*ECALRecHit_rhz)[erhiter];
                rhtresmap[jrhid] = erh_timeRes[erhiter];
    
            	if( DEBUG ) std::cout << erhe << " " << erht << " " << hasBadTime;
                if( DEBUG ) std::cout << " " << erht << " " << rhx << " " << rhy << " " << rhz << std::endl;
		//skip endcap rechits
		//if(fabs((*ECALRecHit_eta)[erhiter]) > 1.479) continue;
		//cout << "rh e " <<  erhe << " rh t " << erht << " rhid " << jrhid << " eta " << (*ECALRecHit_eta)[erhiter] << endl;
		totE += erhe;
		_ca.AddRecHit(rhx, rhy, rhz, erhe, erht, jrhid, hasBadTime);


            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

		//if not enough rechits to cluster, skip this object
		if(_ca.GetNRecHits() < 2) continue;
		ClusterObj jetobj;
		//cout << "running jet clustering" << endl;
		int clusterret = _ca.RunClustering(jetobj, false);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		//cout << "seljetindx " << selJetIndex << endl;
		jetobj.SetUserIndex(selJetIndex);
		jetobjs.push_back(jetobj);
		cout << "it " << it << " selJet energy " << (*Jet_energy)[it] << " BHC jet energy " << jetobj.GetEnergy() << " tot rhE " << totE << endl;
		_ca.ClearRecHitList();
    }//<<>>for( uInt it = 0; it < nJets; it++ )

	//sort jet objs by energy
	sort(jetobjs.begin(), jetobjs.end(), ptsort);
	for(int j = 0; j < jetobjs.size(); j++){
		//cout << "jetobj #" << j << endl;
		ClusterObj jetobj = jetobjs[j];
		jetobj.CalculatePUScores();
		BHCJetInfo.fillBranch( "selJetBHC_selJetIndex", jetobj.GetUserIndex() );
		nbhcjets++;
		//debug
		//do subcluster observables
		int nk = jetobj.GetNSubclusters();
		BHCJetInfo.fillBranch( "selJetBHC_nSubclusters", nk);
		vector<bool> puscores;
		jetobj.GetPUScores(puscores);
		vector<vector<float>> detbkgscores;
		jetobj.CalculateDetBkgScores(false);
		jetobj.GetDetBkgScores(detbkgscores);
		for(int k = 0; k < nk; k++){
			BHCJetInfo.fillBranch( "selJetBHCSubcl_selJetBHCIndex", j);	
			BHCJetInfo.fillBranch( "selJetBHCSubcl_eta", jetobj.GetSubclusterEtaCenter(k));	
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phi", jetobj.GetSubclusterPhiCenter(k));
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaVar", jetobj.GetSubclusterEtaVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phiVar", jetobj.GetSubclusterPhiVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_timeVar", jetobj.GetSubclusterTimeVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_energy", jetobj.GetSubclusterEnergy(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaPhiCov", jetobj.GetSubclusterEtaPhiCov(k));
			if(puscores.size() == 0){
				if(true) cout << "Run CalculatePUScores() BEFORE this loop" << endl;
				continue;
			}
			//pu score
			BHCJetInfo.fillBranch( "selJetBHCSubcl_puScore", (int)puscores[k]);
		       	//det bkg score
			if(puscores[k] == 0) continue;
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_physBkgCNNScore", detbkgscores[k][0]);	
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_beamHaloCNNScore", detbkgscores[k][1]);	
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_spikeCNNScore", detbkgscores[k][2]);	
    			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_selJetBHCIndex", j);
		}
		//center
		float jeteta = jetobj.GetEtaCenter();
		BHCJetInfo.fillBranch( "selJetBHC_eta", jeteta);
		float jetphi = jetobj.GetPhiCenter();
		BHCJetInfo.fillBranch( "selJetBHC_phi", jetphi);
		//time at detector face
		float jettime = jetobj.GetObjTime_Det();
		BHCJetInfo.fillBranch( "selJetBHC_DetTime", jettime);
		//time at PV
		float jettime_pv = jetobj.GetObjTime_PV();
		BHCJetInfo.fillBranch( "selJetBHC_PVTime", jettime_pv);
		BHCJetInfo.fillBranch("selJetBHC_energy", (float)jetobj.GetEnergy());
		BHCJetInfo.fillBranch( "selJetBHC_pt", (float)jetobj.GetPt());			
		pvtime += jettime_pv*jetobj.GetEnergy();
		norm += jetobj.GetEnergy();

		//covariance
		float etavar = jetobj.GetEtaVar();
		float phivar = jetobj.GetPhiVar();
		float etaphicov = jetobj.GetEtaPhiCov();
		double majlen, minlen;
		majlen = -1;
		minlen = -1;
		jetobj.GetMajMinLengths(majlen, minlen); //space only
		BHCJetInfo.fillBranch( "selJetBHC_etaVar", etavar);
		BHCJetInfo.fillBranch( "selJetBHC_phiVar", phivar);
		BHCJetInfo.fillBranch( "selJetBHC_etaPhiCov", etaphicov);
		BHCJetInfo.fillBranch( "selJetBHC_majlen", (float)majlen);
		BHCJetInfo.fillBranch( "selJetBHC_minlen", (float)minlen);
		//calculate time significance
		float timeSignificance = -999;//waiting for centralized time significance code
		BHCJetInfo.fillBranch( "selJetBHC_timeSignficance", timeSignificance);
		//cout << "cleaning out jet pu " << endl;
		/////////CLEAN OUT PU/////////
		jetobj.CleanOutPU();
		//do PU-cleaned observables
		//center
		jeteta = jetobj.GetEtaCenter();
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_eta", jeteta);
		jetphi = jetobj.GetPhiCenter();
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_phi", jetphi);
		//time at detector face
		jettime = jetobj.GetObjTime_Det();
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_DetTime", jettime);
		//time at PV
		jettime_pv = jetobj.GetObjTime_PV();
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_PVTime", jettime_pv);
		BHCJetInfo.fillBranch("selJetBHCPUCleaned_energy", (float)jetobj.GetEnergy());
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_pt", (float)jetobj.GetPt());			
		
		pvtime_PUcleaned += jettime_pv*jetobj.GetEnergy();
		norm_PUcleaned += jetobj.GetEnergy();
		//covariance
		etavar = jetobj.GetEtaVar();
		phivar = jetobj.GetPhiVar();
		etaphicov = jetobj.GetEtaPhiCov();
		jetobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_etaVar", etavar);
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_phiVar", phivar);
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_etaPhiCov", etaphicov);
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_majlen", (float)majlen);
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_minlen", (float)minlen);
		//calculate time significance
		timeSignificance = -999;//waiting for centralized time significance code
		//vector<double> weights; jetobj.GetRecHitWeights(weights); ///these weights are the PU projected out weights if PU cleaning has been applied via jetobj.CleanOutPU() (otherwise sum to 1)
		BHCJetInfo.fillBranch( "selJetBHCPUCleaned_timeSignficance", timeSignificance);
		//


		/////////CLEAN OUT DET BKG/////////
		/*
		double minscore = 0.9;
		jetobj.CleanOutDetBkg(minscore);
		//do PU-cleaned && det bkg-cleaned observables
		//center
		jeteta = jetobj.GetEtaCenter();
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_eta", jeteta);
		jetphi = jetobj.GetPhiCenter();
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_phi", jetphi);
		//time at detector face
		jettime = jetobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_DetTime", jettime);
		//time at PV
		jettime_pv = jetobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_PVTime", jettime_pv);
		
		pvtime_PUdetBkgCleaned += jettime_pv*jetobj.GetEnergy();
		norm_PUdetBkgCleaned += jetobj.GetEnergy();
	
		//covariance
		etavar = jetobj.GetEtaVar();
		phivar = jetobj.GetPhiVar();
		etaphicov = jetobj.GetEtaPhiCov();
		//jetobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_etaVar", etavar);
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_phiVar", phivar);
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_majlen", (float)majlen);
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_minlen", (float)minlen);
		//calculate time significance
		timeSignificance = jetobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selJetBHCPUDetBkgCleaned_timeSignficance", timeSignificance);
		*/

		


    }//<<>>for( uInt it = 0; it < nPhotons; it++ )
	BHCJetInfo.fillBranch("nSelBHCJets",(int)jetobjs.size());
	pvtime /= norm;
	pvtime_PUcleaned /= norm_PUcleaned;
	pvtime_PUdetBkgCleaned /= norm_PUdetBkgCleaned;

  	BHCJetInfo.fillBranch( "PVtime_BHC", pvtime );
  	BHCJetInfo.fillBranch( "PVtime_BHCPUcleaned", pvtime_PUcleaned );
  	BHCJetInfo.fillBranch( "PVtime_BHCPUdetBkgCleaned", pvtime_PUdetBkgCleaned );

}//<<>>void KUCMSAodSkimmer::processMLJets()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setBCBranches( TTree* fOutTree ){
/*
  for(const string cleanedType : {"BHC", "BHCPUCleaned"}) {
    	BHCPhoInfo.makeBranch("selPho"+cleanedType+"_selPhoIndex", VINT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_eta", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phi", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_DetTime", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_PVTime", VFLOAT);

	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_energy", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_pt", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phiVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_majlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_minlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeSignficance", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeSignficanceNum", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeSignficanceDenom", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_nSubclusters", VINT);
  }
    	BHCPhoInfo.makeBranch("selPhoBHCSubcl_selPhotonBHCIndex", VINT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_physBkgCNNScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_beamHaloCNNScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_spikeCNNScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_isoBkgANNScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_nonIsoBkgANNScore", VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_puScore", VINT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_etaVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_phiVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_timeVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_energy",VFLOAT);
  	BHCPhoInfo.makeBranch("nSelBHCPhos",INT);
		

  BHCPhoInfo.attachBranches( fOutTree );
*/
    BHCJetInfo.makeBranch( "selJetBHC_selJetIndex", VINT);
	BHCJetInfo.makeBranch( "selJetBHC_nSubclusters", VINT);
	BHCJetInfo.makeBranch( "selJetBHCSubcl_time", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_eta", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_phi", VFLOAT);

	BHCJetInfo.makeBranch( "selJetBHCSubcl_etaVar", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_phiVar", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_timeVar", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_energy", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	BHCJetInfo.makeBranch( "selJetBHCSubcl_puScore", VINT);
	//det bkg score
	BHCJetInfo.makeBranch( "selJetBHC_selJetBHCIndex", VINT);	
	BHCJetInfo.makeBranch( "selJetBHCSubcl_selJetBHCIndex", VINT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_physBkgCNNScore", VFLOAT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_beamHaloCNNScore",VFLOAT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_spikeCNNScore", VFLOAT);	
    	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_selJetBHCIndex", VINT);
  for(const string cleanedType : {"BHC", "BHCPUCleaned", "BHCPUDetBkgCleaned"}) {
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_eta", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_phi", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_DetTime", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_PVTime", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_energy", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_pt", VFLOAT);			

	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_etaVar", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_phiVar", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_majlen", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_minlen", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_timeSignficance", VFLOAT);
  }
  BHCJetInfo.makeBranch("nSelBHCJets",INT);
  BHCJetInfo.makeBranch( "PVtime_BHC", FLOAT );
  BHCJetInfo.makeBranch( "PVtime_BHCPUcleaned", FLOAT );
  BHCJetInfo.makeBranch( "PVtime_BHCPUdetBkgCleaned", FLOAT );

  BHCJetInfo.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

