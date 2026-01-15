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
		cout << "running jet clustering" << endl;
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

