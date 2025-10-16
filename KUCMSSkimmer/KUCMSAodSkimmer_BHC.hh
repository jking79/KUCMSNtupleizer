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

    int selPhoIndex = -1;

	uInt nPhotons = Photon_excluded->size();
	if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons to get BC info" << std::endl;
	ClusterAnalyzer ca;
	//ca.SetVerbosity(0); //can turn on to see low-level warnings
	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	ca.SetTransferFactor(1/30.); //set to 1/min photon pt
	ca.SetupDetIDsEB();
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
				//cout << "rh e " <<  erhe << " rh t " << erht << " rhid " << scrhid << " eta " << (*ECALRecHit_eta)[erhiter] << endl;
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
		phoobj.CalculatePUScores();
		phoobj.CalculateDetBkgScores();
		phoobj.CalculateObjTimeSig(rhtresmap);

		BHCInfo.fillBranch( "selPhoBHC_selPhotonIndex", selPhoIndex );

		//do subcluster observables
		int nk = phoobj.GetNSubclusters();
		BHCInfo.fillBranch( "selPhoBHC_nSubclusters", nk);
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
			BHCInfo.fillBranch( "selPhoBHCSubcl_time", subcltime);			
			BHCInfo.fillBranch( "selPhoBHCSubcl_eta", subcleta);			
			BHCInfo.fillBranch( "selPhoBHCSubcl_phi", subclphi);
			Matrix subcl_cov = subcls[k].GetCovariance();
			float etavar = subcl_cov.at(0,0);
			float phivar = subcl_cov.at(1,1);
			float etaphicov = subcl_cov.at(0,1);
			BHCInfo.fillBranch( "selPhoBHCSubcl_etaVar", etavar);			
			BHCInfo.fillBranch( "selPhoBHCSubcl_phiVar", phivar);			
			BHCInfo.fillBranch( "selPhoBHCSubcl_etaPhiCov", etaphicov);
			//pu score
			//CHECK
			//cout << "PU subcl score " << puscores[k] << endl;
			//BHCInfo.fillBranch( "selPhoBHCsubcl_puScore", (float)puscores[k]);
		       	////det bkg score 
			//TODO - for photons, det bkg score is calculated over full PU-cleaned object
			//BHCInfo.fillBranch( "selPhoBHCsubcl_detBkgMaxClassScore", (float)detbkgscores[k].second);
			//BHCInfo.fillBranch( "selPhoBHCsubcl_detBkgMaxClass", detbkgscores[k].first);	
			////CHECK
			BHCInfo.fillBranch( "selPhoBHCsubcl_photonIndex", (int)it);	
		}
		//center
		float phoeta = phoobj.GetEtaCenter();
		BHCInfo.fillBranch( "selPhoBHC_eta", phoeta);
		float phophi = phoobj.GetPhiCenter();
		BHCInfo.fillBranch( "selPhoBHC_phi", phophi);
		//time at detector face
		float photime = phoobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selPhoBHC_DetTime", photime);
		//time at PV
		float photime_pv = phoobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selPhoBHC_PVTime", photime_pv);
	
		//covariance
		float etavar = phoobj.GetEtaVar();
		float phivar = phoobj.GetPhiVar();
		float etaphicov = phoobj.GetEtaPhiCov();
		double majlen, minlen;
		majlen = -1;
		minlen = -1;
		phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCInfo.fillBranch( "selPhoBHC_etaVar", etavar);
		BHCInfo.fillBranch( "selPhoBHC_phiVar", phivar);
		BHCInfo.fillBranch( "selPhoBHC_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selPhoBHC_majlen", (float)majlen);
		BHCInfo.fillBranch( "selPhoBHC_minlen", (float)minlen);
		//calculate time significance
		//CHECK
		float timeSignificance = phoobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selPhoBHC_timeSignficance", timeSignificance);


		//TODO - CHECK PU cleaning and det bkg cleaning
		/*
		phoobj.CleanOutPU();
		//do PU-cleaned observables
		//center
		phoeta = phoobj.GetEtaCenter();
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_eta", phoeta);
		phophi = phoobj.GetPhiCenter();
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_phi", phophi);
		//time at detector face
		photime = phoobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_DetTime", photime);
		//time at PV
		photime_pv = phoobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_PVTime", photime_pv);
	
		//covariance
		etavar = phoobj.GetEtaVar();
		phivar = phoobj.GetPhiVar();
		etaphicov = phoobj.GetEtaPhiCov();
		phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_etaVar", etavar);
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_phiVar", phivar);
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_majlen", (float)majlen);
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_minlen", (float)minlen);
		//calculate time significance
		//CHECK
		timeSignificance = phoobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selPhoBHCPUCleaned_timeSignficance", timeSignificance);
		//
		double minscore = 0.9;
		phoobj.CleanOutDetBkg(minscore);
		//do PU-cleaned && det bkg-cleaned observables
		//center
		phoeta = phoobj.GetEtaCenter();
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_eta", phoeta);
		phophi = phoobj.GetPhiCenter();
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_phi", phophi);
		//time at detector face
		photime = phoobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_DetTime", photime);
		//time at PV
		photime_pv = phoobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_PVTime", photime_pv);
	
		//covariance
		etavar = phoobj.GetEtaVar();
		phivar = phoobj.GetPhiVar();
		etaphicov = phoobj.GetEtaPhiCov();
		//phoobj.GetMajMinLengths(majlen, minlen); //space only
		
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_etaVar", etavar);
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_phiVar", phivar);
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_majlen", (float)majlen);
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_minlen", (float)minlen);
		//calculate time significance
		timeSignificance = phoobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selPhoBHCPUDetBkgCleaned_timeSignficance", timeSignificance);
		*/

		ca.ClearRecHitList();
        // branches created strating on 2546 in  void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree )
        // on line 2709 and below make a BHCInfo branch : selPhotons.makeBranch( "yourvarname", VFLOAT );		

	}//<<>>for( uInt it = 0; it < nPhotons; it++ )

}//<<>>void KUCMSAodSkimmer::processMLPhotons()

void KUCMSAodSkimmer::processBHCJets(){

    if( DEBUG ) std::cout << "Finding BC information for jets" << std::endl;

	int selJetIndex = -1;

	ClusterAnalyzer ca;
	//ca.SetVerbosity(0); //can turn on to see low-level warnings
	//ca.SetDetectorCenter(x, y, z); //set to beamspot center, defaults to (0,0,0)
	ca.SetPV(PV_x, PV_y, PV_z); //set to PV coords (assuming these are global vars)
	ca.SetTransferFactor(1/30.); //set to 1/min jet pt
	ca.SetupDetIDsEB();

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

		BHCInfo.fillBranch( "selJetBHC_selJetIndex", selJetIndex );

		//do subcluster observables
		int nk = jetobj.GetNSubclusters();
		BHCInfo.fillBranch( "selJetBHC_nSubclusters", nk);
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
			BHCInfo.fillBranch( "selJetBHCSubcl_time", subcltime);			
			BHCInfo.fillBranch( "selJetBHCSubcl_eta", subcleta);			
			BHCInfo.fillBranch( "selJetBHCSubcl_phi", subclphi);
			Matrix subcl_cov = subcls[k].GetCovariance();
			float etavar = subcl_cov.at(0,0);
			float phivar = subcl_cov.at(1,1);
			float etaphicov = subcl_cov.at(0,1);
			BHCInfo.fillBranch( "selJetBHCSubcl_etaVar", etavar);			
			BHCInfo.fillBranch( "selJetBHCSubcl_phiVar", phivar);			
			BHCInfo.fillBranch( "selJetBHCSubcl_etaPhiCov", etaphicov);
			//pu score
			BHCInfo.fillBranch( "selJetBHCsubcl_puScore", puscores[k]);
		       	//det bkg score
			BHCInfo.fillBranch( "selJetBHCsubcl_detBkgMaxClassScore", (float)detbkgscores[k].second);
			BHCInfo.fillBranch( "selJetBHCsubcl_detBkgMaxClass", detbkgscores[k].first);	
			BHCInfo.fillBranch( "selJetBHCsubcl_jetIndex", it);	
		}
		//center
		float jeteta = jetobj.GetEtaCenter();
		BHCInfo.fillBranch( "selJetBHC_eta", jeteta);
		float jetphi = jetobj.GetPhiCenter();
		BHCInfo.fillBranch( "selJetBHC_phi", jetphi);
		//time at detector face
		float jettime = jetobj.GetObjTime_Det();
		BHCInfo.fillBranch( "selJetBHC_DetTime", jettime);
		//time at PV
		float jettime_pv = jetobj.GetObjTime_PV();
		BHCInfo.fillBranch( "selJetBHC_PVTime", jettime_pv);

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
		
		BHCInfo.fillBranch( "selJetBHC_etaVar", etavar);
		BHCInfo.fillBranch( "selJetBHC_phiVar", phivar);
		BHCInfo.fillBranch( "selJetBHC_etaPhiCov", etaphicov);
		BHCInfo.fillBranch( "selJetBHC_majlen", (float)majlen);
		BHCInfo.fillBranch( "selJetBHC_minlen", (float)minlen);
		//calculate time significance
		float timeSignificance = jetobj.GetObjTimeSig();
		BHCInfo.fillBranch( "selJetBHC_timeSignficance", timeSignificance);

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

  	BHCInfo.fillBranch( "PVtime", pvtime );
  	BHCInfo.fillBranch( "PVtime_PUcleaned", pvtime_PUcleaned );
  	BHCInfo.fillBranch( "PV_PUdetBkgCleaned", pvtime_PUdetBkgCleaned );

}//<<>>void KUCMSAodSkimmer::processMLJets()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setBCBranches( TTree* fOutTree ){

    BHCInfo.makeBranch( "selPhoBHC_selPhotonIndex", VINT);
	BHCInfo.makeBranch( "selPhoBHC_nSubclusters", VINT);
	BHCInfo.makeBranch( "selPhoBHCSubcl_time", VFLOAT);			
	BHCInfo.makeBranch( "selPhoBHCSubcl_eta", VFLOAT);			
	BHCInfo.makeBranch( "selPhoBHCSubcl_phi", VFLOAT);
	BHCInfo.makeBranch( "selPhoBHCSubcl_etaVar", VFLOAT);			
	BHCInfo.makeBranch( "selPhoBHCSubcl_phiVar", VFLOAT);			
	BHCInfo.makeBranch( "selPhoBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	BHCInfo.makeBranch( "selPhoBHCsubcl_puScore", VFLOAT);
	//det bkg score
	BHCInfo.makeBranch( "selPhoBHCsubcl_detBkgMaxClassScore", VFLOAT);
	BHCInfo.makeBranch( "selPhoBHCsubcl_detBkgMaxClass", VINT);	
	BHCInfo.makeBranch( "selPhoBHCsubcl_photonIndex", VINT);	
  for(const string cleanedType : {"BHC", "BHCPUCleaned", "BHCPUDetBkgCleaned"}) {
	BHCInfo.makeBranch( "selPho"+cleanedType+"_eta", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_phi", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_DetTime", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_PVTime", VFLOAT);

	BHCInfo.makeBranch( "selPho"+cleanedType+"_etaVar", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_phiVar", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_majlen", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_minlen", VFLOAT);
	BHCInfo.makeBranch( "selPho"+cleanedType+"_timeSignficance", VFLOAT);
  }

    BHCInfo.makeBranch( "selJetBHC_selJetIndex", VINT);
	BHCInfo.makeBranch( "selJetBHC_nSubclusters", VINT);
	BHCInfo.makeBranch( "selJetBHCSubcl_time", VFLOAT);			
	BHCInfo.makeBranch( "selJetBHCSubcl_eta", VFLOAT);			
	BHCInfo.makeBranch( "selJetBHCSubcl_phi", VFLOAT);
	BHCInfo.makeBranch( "selJetBHCSubcl_etaVar", VFLOAT);			
	BHCInfo.makeBranch( "selJetBHCSubcl_phiVar", VFLOAT);			
	BHCInfo.makeBranch( "selJetBHCSubcl_etaPhiCov", VFLOAT);
	//pu score
	BHCInfo.makeBranch( "selJetBHCsubcl_puScore", VFLOAT);
	//det bkg score
	BHCInfo.makeBranch( "selJetBHCsubcl_detBkgMaxClassScore", VFLOAT);
	BHCInfo.makeBranch( "selJetBHCsubcl_detBkgMaxClass", VINT);	
	BHCInfo.makeBranch( "selJetBHCsubcl_jetIndex", VINT);	
  for(const string cleanedType : {"BHC", "BHCPUCleaned", "BHCPUDetBkgCleaned"}) {
	BHCInfo.makeBranch( "selJet"+cleanedType+"_eta", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_phi", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_DetTime", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_PVTime", VFLOAT);

	BHCInfo.makeBranch( "selJet"+cleanedType+"_etaVar", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_phiVar", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_etaPhiCov", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_majlen", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_minlen", VFLOAT);
	BHCInfo.makeBranch( "selJet"+cleanedType+"_timeSignficance", VFLOAT);
  }
  BHCInfo.makeBranch( "PVtime", FLOAT );
  BHCInfo.makeBranch( "PVtime_PUcleaned", FLOAT );
  BHCInfo.makeBranch( "PV_PUdetBkgCleaned", FLOAT );

  BHCInfo.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

