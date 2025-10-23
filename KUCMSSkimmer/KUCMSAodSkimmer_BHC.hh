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
	ClusterAnalyzer ca;
	//ca.SetVerbosity(0); //can turn on to see low-level warnings
	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	ca.SetTransferFactor(1/30.); //set to 1/min photon pt
	//ca.SetDetIDsEB(_detidmap);
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
				ca.AddRecHit(rhx, rhy, rhz, erhe, erht, scrhid, hasGainSwitch);
			}//<<>>if( ecalrhiter != -1 )
		}//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
	
		//if not enough rechits to cluster, skip this object
		if(ca.GetNRecHits() < 2) continue;

		//
		//  call to BayesianClustering framework to process photon information 
		//  runs BHC (NlnN) algorithm to cluster photon with subcluster regularization
		//  returns lead (highest energy) cluster found

		//put safety in so that if there are < 2 rhs to cluster, skip

		ClusterObj phoobj;
		int clusterret = ca.RunClustering(phoobj);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		phoobj.CalculateObjTimes();
		phoobj.CleanOutPU();
		//phoobj.CalculatePUScores();
		phoobj.CalculateDetBkgScores();
		phoobj.CalculateObjTimeSig(rhtresmap);

		//do NOT do subcluster observables
		/*
		BHCPhoInfo.fillBranch( "selPhoBHC_selPhotonIndex", selPhoIndex );

		//do subcluster observables
		int nk = phoobj.GetNSubclusters();
		BHCPhoInfo.fillBranch( "selPhoBHC_nSubclusters", nk);
		vector<Jet> subcls;
		phoobj.GetSubclusters(subcls);

		vector<bool> puscores;
		phoobj.GetPUScores(puscores);
		vector<pair<int, double>> detbkgscores;
		phoobj.GetDetBkgScores(detbkgscores);
		for(int k = 0; k < nk; k++){
			//cout << "cluster #" << k << endl;
			float subcltime = subcls[k].t();
			float subcleta = subcls[k].eta();
			float subclphi = subcls[k].phi();
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_time", subcltime);			
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_eta", subcleta);			
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_phi", subclphi);
			Matrix subcl_cov = subcls[k].GetCovariance();
			float etavar = subcl_cov.at(0,0);
			float phivar = subcl_cov.at(1,1);
			float etaphicov = subcl_cov.at(0,1);
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_etaVar", etavar);			
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_phiVar", phivar);			
			BHCPhoInfo.fillBranch( "selPhoBHCSubcl_etaPhiCov", etaphicov);
			//pu score
			//CHECK
			//cout << "PU subcl score " << puscores[k] << endl;
			//BHCPhoInfo.fillBranch( "selPhoBHCsubcl_puScore", (float)puscores[k]);
		       	////det bkg score 
			//TODO - for photons, det bkg score is calculated over full PU-cleaned object
			//BHCPhoInfo.fillBranch( "selPhoBHCsubcl_detBkgMaxClassScore", (float)detbkgscores[k].second);
			//BHCPhoInfo.fillBranch( "selPhoBHCsubcl_detBkgMaxClass", detbkgscores[k].first);	
			////CHECK
			BHCPhoInfo.fillBranch( "selPhoBHCsubcl_photonIndex", (int)it);	
		}
		*/
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
		float etaphicov = phoobj.GetEtaPhiCov();
		double majlen, minlen;
		majlen = -1;
		minlen = -1;
		phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCPhoInfo.fillBranch( "selPhoBHC_etaVar", etavar);
		BHCPhoInfo.fillBranch( "selPhoBHC_phiVar", phivar);
		BHCPhoInfo.fillBranch( "selPhoBHC_etaPhiCov", etaphicov);
		BHCPhoInfo.fillBranch( "selPhoBHC_majlen", (float)majlen);
		BHCPhoInfo.fillBranch( "selPhoBHC_minlen", (float)minlen);
		//calculate time significance
		//CHECK
		float timeSignificance = phoobj.GetObjTimeSig();
		BHCPhoInfo.fillBranch( "selPhoBHC_timeSignficance", timeSignificance);


		//TODO - CHECK PU cleaning and det bkg cleaning
		/*
		phoobj.CleanOutPU();
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
		//CHECK
		timeSignificance = phoobj.GetObjTimeSig();
		BHCPhoInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficance", timeSignificance);
		//
		double minscore = 0.9;
		phoobj.CleanOutDetBkg(minscore);
		//do PU-cleaned && det bkg-cleaned observables
		//center
		phoeta = phoobj.GetEtaCenter();
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_eta", phoeta);
		phophi = phoobj.GetPhiCenter();
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_phi", phophi);
		//time at detector face
		photime = phoobj.GetObjTime_Det();
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_DetTime", photime);
		//time at PV
		photime_pv = phoobj.GetObjTime_PV();
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_PVTime", photime_pv);
	
		//covariance
		etavar = phoobj.GetEtaVar();
		phivar = phoobj.GetPhiVar();
		etaphicov = phoobj.GetEtaPhiCov();
		//phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_etaVar", etavar);
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_phiVar", phivar);
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_etaPhiCov", etaphicov);
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_majlen", (float)majlen);
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_minlen", (float)minlen);
		//calculate time significance
		timeSignificance = phoobj.GetObjTimeSig();
		BHCPhoInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_timeSignficance", timeSignificance);
		*/

		ca.ClearRecHitList();
        // branches created strating on 2546 in  void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree )
        // on line 2709 and below make a BHCPhoInfo branch : selPhotons.makeBranch( "yourvarname", VFLOAT );		

	}//<<>>for( uInt it = 0; it < nPhotons; it++ )

}//<<>>void KUCMSAodSkimmer::processMLPhotons()

void KUCMSAodSkimmer::processBHCJets(){

    if( DEBUG ) std::cout << "Finding BC information for jets" << std::endl;
	BHCJetInfo.clearBranches();

	int selJetIndex = -1;

	ClusterAnalyzer ca;
	//ca.SetVerbosity(0); //can turn on to see low-level warnings
	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	ca.SetTransferFactor(1/30.); //set to 1/min jet pt
	//ca.SetDetIDsEB(_detidmap);

    std::vector<Jet> jetrhs;
  	uInt nJets = Jet_energy->size();
	float pvtime = 0;
	float norm = 0;
	float pvtime_PUcleaned = 0;
	float norm_PUcleaned = 0;
	float pvtime_PUdetBkgCleaned = 0;
	float norm_PUdetBkgCleaned = 0;
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

		ca.AddRecHit(rhx, rhy, rhz, erhe, erht, jrhid, hasGainSwitch);


            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

		//if not enough rechits to cluster, skip this object
		if(ca.GetNRecHits() < 2) continue;

		ClusterObj jetobj;
		int clusterret = ca.RunClustering(jetobj);
		if(clusterret == -1) continue; //bad clustering (not enough points, not able to find any clusters, etc)
		jetobj.CalculateObjTimes();
		jetobj.CalculatePUScores();
		jetobj.CalculateDetBkgScores();
		jetobj.CalculateObjTimeSig(rhtresmap);

		BHCJetInfo.fillBranch( "selJetBHC_selJetIndex", selJetIndex );

		//do subcluster observables
		int nk = jetobj.GetNSubclusters();
		BHCJetInfo.fillBranch( "selJetBHC_nSubclusters", nk);
		vector<Jet> subcls;
		jetobj.GetSubclusters(subcls);

		vector<bool> puscores;
		jetobj.GetPUScores(puscores);
		vector<pair<int, double>> detbkgscores;
		jetobj.GetDetBkgScores(detbkgscores);
		for(int k = 0; k < nk; k++){
			float subcltime = subcls[k].t();
			float subcleta = subcls[k].eta();
			float subclphi = subcls[k].phi();
			BHCJetInfo.fillBranch( "selJetBHCSubcl_time", subcltime);			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_eta", subcleta);			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phi", subclphi);
			Matrix subcl_cov = subcls[k].GetCovariance();
			float etavar = subcl_cov.at(0,0);
			float phivar = subcl_cov.at(1,1);
			float etaphicov = subcl_cov.at(0,1);
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaVar", etavar);			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_phiVar", phivar);			
			BHCJetInfo.fillBranch( "selJetBHCSubcl_etaPhiCov", etaphicov);
			//pu score
			BHCJetInfo.fillBranch( "selJetBHCsubcl_puScore", puscores[k]);
		       	//det bkg score
			BHCJetInfo.fillBranch( "selJetBHCsubcl_detBkgMaxClassScore", (float)detbkgscores[k].second);
			BHCJetInfo.fillBranch( "selJetBHCsubcl_detBkgMaxClass", detbkgscores[k].first);	
			BHCJetInfo.fillBranch( "selJetBHCsubcl_jetIndex", it);	
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

		ca.ClearRecHitList();
		


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

    BHCPhoInfo.makeBranch( "selPhoBHC_selPhotonIndex", VINT);
	BHCPhoInfo.makeBranch( "selPhoBHC_nSubclusters", VINT);
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_time", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_eta", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_phi", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_etaVar", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_phiVar", VFLOAT);			
	BHCPhoInfo.makeBranch( "selPhoBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	BHCPhoInfo.makeBranch( "selPhoBHCsubcl_puScore", VFLOAT);
	//det bkg score
	BHCPhoInfo.makeBranch( "selPhoBHCsubcl_detBkgMaxClassScore", VFLOAT);
	BHCPhoInfo.makeBranch( "selPhoBHCsubcl_detBkgMaxClass", VINT);	
	BHCPhoInfo.makeBranch( "selPhoBHCsubcl_photonIndex", VINT);	
  for(const string cleanedType : {"BHC", "BHCPUCleaned", "BHCPUDetBkgCleaned"}) {
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_eta", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phi", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_DetTime", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_PVTime", VFLOAT);

	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_phiVar", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_majlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_minlen", VFLOAT);
	BHCPhoInfo.makeBranch( "selPho"+cleanedType+"_timeSignficance", VFLOAT);
  }

  BHCPhoInfo.attachBranches( fOutTree );

    BHCJetInfo.makeBranch( "selJetBHC_selJetIndex", VINT);
	BHCJetInfo.makeBranch( "selJetBHC_nSubclusters", VINT);
	BHCJetInfo.makeBranch( "selJetBHCSubcl_time", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_eta", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_phi", VFLOAT);
	BHCJetInfo.makeBranch( "selJetBHCSubcl_etaVar", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_phiVar", VFLOAT);			
	BHCJetInfo.makeBranch( "selJetBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	BHCJetInfo.makeBranch( "selJetBHCsubcl_puScore", VFLOAT);
	//det bkg score
	BHCJetInfo.makeBranch( "selJetBHCsubcl_detBkgMaxClassScore", VFLOAT);
	BHCJetInfo.makeBranch( "selJetBHCsubcl_detBkgMaxClass", VINT);	
	BHCJetInfo.makeBranch( "selJetBHCsubcl_jetIndex", VINT);	
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

