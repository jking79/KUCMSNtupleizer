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
                bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];
				float erht = erh_corTime[erhiter];
				float rhx = (*ECALRecHit_rhx)[erhiter];
                float rhy = (*ECALRecHit_rhy)[erhiter];
                float rhz = (*ECALRecHit_rhz)[erhiter];
				rhtresmap[scrhid] = erh_timeRes[erhiter];
				energyToId[erhe] = scrhid;
                if( DEBUG ) std::cout << " " << erht << " " << rhx << " " << rhy << " " << rhz << std::endl;
				//skip endcap rechits
				if(fabs((*ECALRecHit_eta)[erhiter]) > 1.479) continue;
				//cout << "rh e " <<  erhe << " rh t " << erht << " rhid " << scrhid << " eta " << (*ECALRecHit_eta)[erhiter] << " phi " << (*ECALRecHit_phi)[erhiter] << endl;
				//cout << "rh e " <<  erhe << " rh t " << erht << " x " << rhx << " y " << rhy << " z " << rhz << " rhid " << scrhid << endl;
				_ca.AddRecHit(rhx, rhy, rhz, erhe, erht, scrhid, hasGainSwitch);
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
		int clusterret = _ca.RunClustering(phoobj);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		phoobj.CalculateObjTimes();
		phoobj.CalculateObjTimeSig(rhtresmap);
		phoobj.CalculatePUScores();

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
		
		BHCPhoInfo.fillBranch( "selPhoBHC_selPhoIndex", selPhoIndex );
		BHCPhoInfo.fillBranch( "selPhoBHC_etaVar", (float)etavar);
		BHCPhoInfo.fillBranch( "selPhoBHC_phiVar", (float)phivar);
		BHCPhoInfo.fillBranch( "selPhoBHC_timeVar", (float)timevar);
		BHCPhoInfo.fillBranch( "selPhoBHC_etaPhiCov", (float)etaphicov);
		BHCPhoInfo.fillBranch( "selPhoBHC_majlen", (float)majlen);
		BHCPhoInfo.fillBranch( "selPhoBHC_minlen", (float)minlen);
		//calculate time significance
		//CHECK
		float timeSignificance = phoobj.GetObjTimeSig();
		BHCPhoInfo.fillBranch( "selPhoBHC_timeSignficance", timeSignificance);

		//for studying PU cleaning
		vector<bool> puscores;
		phoobj.GetPUScores(puscores);
		int nsubclusters = phoobj.GetNSubclusters();
		cout << "nsubclusters " << nsubclusters << endl;
		for(int k = 0; k < nsubclusters; k++){
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_etaVar",phoobj.GetSubclusterEtaVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_phiVar",phoobj.GetSubclusterPhiVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_timeVar",phoobj.GetSubclusterTimeVar(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_energy",phoobj.GetSubclusterEnergy(k));
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_selPhotonBHCIndex",selPhoBHCIdx);
			cout << "subcluster #" << k << " pu score " << puscores[k] << " selPhoBHCIndex " << selPhoBHCIdx << " it " << it << endl;
			BHCPhoInfo.fillBranch("selPhoBHCSubcl_puScore",(int)puscores[k]);
		}


		phoobj.CleanOutPU();
		phoobj.CalculateObjTimes();
		phoobj.CalculateObjTimeSig(rhtresmap);
		//do PU-cleaned observables
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
		timeSignificance = phoobj.GetObjTimeSig();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficance", timeSignificance);
		//get detector bkg score for pu-cleaned photon
		phoobj.CalculateDetBkgScores(true);
		vector<vector<float>> detbkgscores;
		phoobj.GetDetBkgScores(detbkgscores);
		cout << "# of detbkg scores " << detbkgscores.size() << endl;
		float physbkg_score = detbkgscores[0][0];
		float bh_score = detbkgscores[0][1];
		float spike_score = detbkgscores[0][2];
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_physBkgMVAScore", physbkg_score);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_beamHaloMVAScore", bh_score);
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_spikeMVAScore", spike_score);
		selPhoBHCIdx++;
		_ca.ClearRecHitList();
        // branches created strating on 2546 in  void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree )
        // on line 2709 and below make a BHCPhoInfo branch : selPhotons.makeBranch( "yourvarname", VFLOAT );		

	}//<<>>for( uInt it = 0; it < nPhotons; it++ )

}//<<>>void KUCMSAodSkimmer::processMLPhotons()

void KUCMSAodSkimmer::processBHCJets(){

    if( DEBUG ) std::cout << "Finding BC information for jets" << std::endl;
	BHCJetInfo.clearBranches();

	int selJetIndex = -1;

	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	_ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	_ca.SetTransferFactor(1/30.); //set to 1/min jet pt

    std::vector<Jet> jetrhs;
  	uInt nJets = Jet_energy->size();
	float pvtime = 0;
	float norm = 0;
	float pvtime_PUcleaned = 0;
	float norm_PUcleaned = 0;
	float pvtime_PUdetBkgCleaned = 0;
	float norm_PUdetBkgCleaned = 0;
	int selJetBHCIdx = 0;
    if( DEBUG ) std::cout << " - Looping over for " << nJets << " photons to get BC info" << std::endl;
  	for( uInt it = 0; it < nJets; it++ ){
        if( not isSelJet[it] ) continue;
		selJetIndex++;

    	std::map< unsigned int, float > rhtresmap;
    	auto rhids = (*Jet_drRhIds)[it];
        int nJetRecHits = rhids.size();
        for( int jiter = 0; jiter < nJetRecHits; jiter++  ){
            auto jrhid = rhids[jiter];
            int erhiter = ( rhIDtoIterMap.find(jrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[jrhid] : -1;
            if( erhiter != -1 ){
                float erhe = (*ECALRecHit_energy)[erhiter];
                bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];
                float erht = erh_corTime[erhiter];
                float rhx = (*ECALRecHit_rhx)[erhiter];
                float rhy = (*ECALRecHit_rhy)[erhiter];
                float rhz = (*ECALRecHit_rhz)[erhiter];
                rhtresmap[jrhid] = erh_timeRes[erhiter];
    
            	if( DEBUG ) std::cout << erhe << " " << erht << " " << hasGainSwitch;
                if( DEBUG ) std::cout << " " << erht << " " << rhx << " " << rhy << " " << rhz << std::endl;
		//skip endcap rechits
		if(fabs((*ECALRecHit_eta)[erhiter]) > 1.479) continue;
		//cout << "rh e " <<  erhe << " rh t " << erht << " rhid " << jrhid << " eta " << (*ECALRecHit_eta)[erhiter] << endl;

		_ca.AddRecHit(rhx, rhy, rhz, erhe, erht, jrhid, hasGainSwitch);


            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

		//if not enough rechits to cluster, skip this object
		if(_ca.GetNRecHits() < 2) continue;

		ClusterObj jetobj;
		int clusterret = _ca.RunClustering(jetobj);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		jetobj.CalculateObjTimes();
		jetobj.CalculatePUScores();
		jetobj.CalculateDetBkgScores(false);
		jetobj.CalculateObjTimeSig(rhtresmap);

		BHCJetInfo.fillBranch( "selJetBHC_selJetIndex", selJetIndex );

		//do subcluster observables
		int nk = jetobj.GetNSubclusters();
		BHCJetInfo.fillBranch( "selJetBHC_nSubclusters", nk);

		vector<bool> puscores;
		jetobj.GetPUScores(puscores);
		vector<vector<float>> detbkgscores;
		jetobj.GetDetBkgScores(detbkgscores);
		for(int k = 0; k < nk; k++){
			BHCJetInfo.fillBranch( "selJetBHCSubcl_selJetBHCIndex", it);	
			BHCJetInfo.fillBranch( "selJetBHCSubcl_eta", jetobj.GetSubclusterEtaCenter(k));	
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phi", jetobj.GetSubclusterPhiCenter(k));
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaVar", jetobj.GetSubclusterEtaVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phiVar", jetobj.GetSubclusterPhiVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_timeVar", jetobj.GetSubclusterTimeVar(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_energy", jetobj.GetSubclusterEnergy(k));			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaPhiCov", jetobj.GetSubclusterEtaPhiCov(k));
			//pu score
			BHCJetInfo.fillBranch( "selJetBHCSubcl_puScore", (int)puscores[k]);
		       	//det bkg score
			if(puscores[k] == 0) continue;
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_jetIndex", it);	
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_physBkgMVAScore", detbkgscores[k][0]);	
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_beamHaloMVAScore", detbkgscores[k][1]);	
			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_spikeMVAScore", detbkgscores[k][2]);	
    			BHCJetInfo.fillBranch( "selJetBHCPUCleanedsubcl_selJetBHCIndex", selJetBHCIdx);
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

		pvtime += jettime_pv*jetobj.GetEnergy();
		norm += jetobj.GetEnergy();

		//covariance
		float etavar = jetobj.GetEtaVar();
		float phivar = jetobj.GetPhiVar();
		float etaphicov = jetobj.GetEtaPhiCov();
		double majlen, minlen;
		majlen = -1;
		minlen = -1;
		//jetobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCJetInfo.fillBranch( "selJetBHC_etaVar", etavar);
		BHCJetInfo.fillBranch( "selJetBHC_phiVar", phivar);
		BHCJetInfo.fillBranch( "selJetBHC_etaPhiCov", etaphicov);
		BHCJetInfo.fillBranch( "selJetBHC_majlen", (float)majlen);
		BHCJetInfo.fillBranch( "selJetBHC_minlen", (float)minlen);
		//calculate time significance
		float timeSignificance = jetobj.GetObjTimeSig();
		BHCJetInfo.fillBranch( "selJetBHC_timeSignficance", timeSignificance);

		//TODO - CHECK PU cleaning and det bkg cleaning
		/*
		/////////CLEAN OUT PU/////////
		jetobj.CleanOutPU();
		//do PU-cleaned observables
		//center
		jeteta = jetobj.GetEtaCenter();
		BHCInfo.fillBranch( "selJetBHCPUCleaned_eta", jeteta);
		jetphi = jetobj.GetPhiCenter();
		BHCInfo.fillBranch( "selJetBHCPUCleaned_phi", jetphi);
		//time at detector face
		jettime = jetobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selJetBHCPUCleaned_DetTime", jettime);
		//time at PV
		jettime_pv = jetobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selJetBHCPUCleaned_PVTime", jettime_pv);
		
		pvtime_PUcleaned += jettime_pv*jetobj.GetEnergy();
		norm_PUcleaned += jetobj.GetEnergy();
	
		//covariance
		etavar = jetobj.GetEtaVar();
		phivar = jetobj.GetPhiVar();
		etaphicov = jetobj.GetEtaPhiCov();
		//jetobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCInfo.fillBranch( "selJetBHCPUCleaned_etaVar", etavar);
		BHCInfo.fillBranch( "selJetBHCPUCleaned_phiVar", phivar);
		BHCInfo.fillBranch( "selJetBHCPUCleaned_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selJetBHCPUCleaned_majlen", (float)majlen);
		BHCInfo.fillBranch( "selJetBHCPUCleaned_minlen", (float)minlen);
		//calculate time significance
		timeSignificance = jetobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selJetBHCPUCleaned_timeSignficance", timeSignificance);
		//
		double minscore = 0.9;
		jetobj.CleanOutDetBkg(minscore);


		/////////CLEAN OUT DET BKG AND PU/////////
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

		selJetBHCIdx++;
		_ca.ClearRecHitList();
		


    }//<<>>for( uInt it = 0; it < nPhotons; it++ )
	pvtime /= norm;
	pvtime_PUcleaned /= norm_PUcleaned;
	pvtime_PUdetBkgCleaned /= norm_PUdetBkgCleaned;

  	BHCJetInfo.fillBranch( "PVtime", pvtime );
  	BHCJetInfo.fillBranch( "PVtime_PUcleaned", pvtime_PUcleaned );
  	BHCJetInfo.fillBranch( "PV_PUdetBkgCleaned", pvtime_PUdetBkgCleaned );

}//<<>>void KUCMSAodSkimmer::processMLJets()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setBCBranches( TTree* fOutTree ){

		/*
	BHCPhoInfo.makeBranch( "selPhoBHC_nSubclusters", VINT);
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_time", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_eta", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_phi", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_etaVar", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_phiVar", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	*/	
    	BHCPhoInfo.makeBranch("selPhoBHC_selPhoIndex", VINT);
  for(const string cleanedType : {"BHC", "BHCPUCleaned"}) {
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_eta", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phi", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_DetTime", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_PVTime", VFLOAT);

	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phiVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_majlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_minlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeSignficance", VFLOAT);
  }
    	BHCPhoInfo.makeBranch("selPhoBHCSubcl_selPhotonBHCIndex", VINT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_physBkgMVAScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_beamHaloMVAScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCPUCleaned_spikeMVAScore", VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_puScore", VINT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_etaVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_phiVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_timeVar",VFLOAT);
	BHCPhoInfo.makeBranch("selPhoBHCSubcl_energy",VFLOAT);
		

  BHCPhoInfo.attachBranches( fOutTree );

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
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_jetIndex", VINT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_physBkgMVAScore", VFLOAT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_beamHaloMVAScore",VFLOAT);	
	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_spikeMVAScore", VFLOAT);	
    	BHCJetInfo.makeBranch( "selJetBHCPUCleanedsubcl_selJetBHCIndex", VINT);
  for(const string cleanedType : {"BHC", "BHCPUCleaned", "BHCPUDetBkgCleaned"}) {
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_eta", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_phi", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_DetTime", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_PVTime", VFLOAT);

	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_etaVar", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_phiVar", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_majlen", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_minlen", VFLOAT);
	BHCJetInfo.makeBranch( "selJet"+cleanedType+"_timeSignficance", VFLOAT);
  }
  BHCJetInfo.makeBranch( "PVtime", FLOAT );
  BHCJetInfo.makeBranch( "PVtime_PUcleaned", FLOAT );
  BHCJetInfo.makeBranch( "PV_PUdetBkgCleaned", FLOAT );

  BHCJetInfo.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

