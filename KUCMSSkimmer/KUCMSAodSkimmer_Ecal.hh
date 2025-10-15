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

void KUCMSAodSkimmer::processRechits(){

	// initilize
    //selECALRecHit.clearBranches(); // <<<<<<<   must do ( but there are no rechit branches right now )
	bool checkRHs = false;
    //bool checkRHs = true;

	// calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

	rhIDtoIterMap.clear();
	erh_corTime.clear();
	erh_timeRes.clear();
    auto nRecHits = ECALRecHit_ID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

		auto rhid = (*ECALRecHit_ID)[it];
		rhIDtoIterMap[rhid] = it;
		auto idinfo = timeCali->getDetIdInfo(rhid);
    	if( DEBUG ) std::cout << " -- TimeCali DetIDInfo : " << idinfo.i2 << " " << idinfo.i1 << " " << idinfo.ecal << std::endl;
		float rht = (*ECALRecHit_time)[it];
		float rhe = (*ECALRecHit_energy)[it];
		float rha = (*ECALRecHit_ampres)[it]; // instead of energy ? ( this is the ADC amplitude in units of the pedistal rms for this rechit )
		float corrht = timeCali->getCorrectedTime( rht, rha, rhid, Evt_run, tctag, mctype );
		float rhtres = timeCali->getTimeResoltuion( rha, rhid, Evt_run, tctag, mctype );
		//if( rhid < 840000000 ){
		//	std::cout << " TC check : inputs : rhid " << rhid << " run " << Evt_run << " tag "  << tctag << " type " << mctype << std::endl;
		//	std::cout << " TC check : time : " << rht << " -> " << corrht << " res : " << rha << " -> " << rhtres << std::endl;
		//}//<<>>if( rhid < 840000000 )
		erh_corTime.push_back( corrht );
        erh_timeRes.push_back( rhtres );
        if( true ){

			hist1d[0]->Fill( rhe, 1 );
            hist1d[1]->Fill( rhe, 1 );
            hist1d[2]->Fill( rhe, 1 );
            hist1d[5]->Fill( rht, 1 );
			//auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] );

		}//<<>>if( (*rhSubdet)[it] == 0 )

		/*  ECAL Rechit information --  not saved currenty
 
		float rhe = (*ECALRecHit_energy)[it];
		float rhx = (*ECALRecHit_rhx)[it];
        float rhy = (*ECALRecHit_rhy)[it];
        float rhz = (*ECALRecHit_rhz)[it];
		float rht = (*ECALRecHit_time)[it];
		float rheta = (*ECALRecHit_eta)[it];
        float rhphi = (*ECALRecHit_phi)[it];
		uInt rhid = (*ECALRecHit_ID)[it];

		*/ // ECAL Rechit information --  not saved currenty

	}//<<>>for( int it = 0; it < nRecHits; it++ )
	if( DEBUG ) std::cout << "Finished w/ rechit loop" << std::endl;
	//  ---  End of RecHit Loop ----------------------------------------------------------------------------------------

    int nMissingRechits = 0;
	int nRechitsPerSC = 0;
	auto nSCs = SuperCluster_nSuperCluster;
	for( int it = 0; it < nSCs; it++ ){

		int nMissingRhInSC = 0;
		bool oot = (*SuperCluster_isOot)[it];
		bool exc = (*SuperCluster_excluded)[it];
		bool orig = (*SuperCluster_original)[it];
		float sce = (*SuperCluster_energyRaw)[it];
		float sceta = (*SuperCluster_eta)[it];
		hist1d[6]->Fill( (*SuperCluster_energyRaw)[it], 1 );
    	//hist1d[3] = new TH1D("sctype","SC !Orig, Orig, OOT, Excl",4,0,4);
    	//hist1d[4] = new TH1D("scorigtype","Orig+OOT, Orig+!OOT, Orig+Exc, !Orig+OOT, !Orig+Exc, !Orig+!OOT",6,0,6);
		if( orig ) hist1d[3]->Fill(1); else hist1d[3]->Fill(0);
        if( oot ) hist1d[3]->Fill(2);
        if( exc ) hist1d[3]->Fill(3);
		if( orig && oot ) hist1d[4]->Fill(0);
        if( orig && not oot ) hist1d[4]->Fill(1);
        if( orig && exc ) hist1d[4]->Fill(2);
        if( not orig && oot ) hist1d[4]->Fill(3);
        if( not orig && not oot ) hist1d[4]->Fill(5);
        if( not orig && exc ) hist1d[4]->Fill(4);
		int xfill = orig ? 1 : 0;
		int yfill = oot ? 1 : exc ? 2 : 0;
		hist2d[1]->Fill(xfill,yfill);

        std::vector<unsigned int> scrhids = (*SuperCluster_rhIds)[it];
		int nSCRhids = scrhids.size();
		if( checkRHs ){
		bool found = false;
		nRechitsPerSC += nSCRhids;
		std::cout << " -- Checking : " << it << " of " << nSCs << " SCs with " << nSCRhids << " rechits. " << std::endl;
		for( int sciter = 0; sciter < nSCRhids; sciter++  ){ 
			auto scrhid = scrhids[sciter];
			for( int iter = 0; iter < nRecHits; iter++ ){
				auto rhid = (*ECALRecHit_ID)[iter];
				//std::cout << " -- Checking : " << rhid << " == " << scrhid << std::endl;
				if( scrhid == rhid ){ found = true; break; } 
			}//<<>>for( int iter = 0; iter < nRecHits; it++ )
			if( not found ){ nMissingRechits++; nMissingRhInSC++; }// std::cout << " ---- !!!!!!  Rechit in SC and not in Rechits !!!!!!! " << std::endl;
		}//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
		}//<<>>if( checkRHs )		
		if( nMissingRhInSC > 0 ){
			std::cout << " -- !! Missing " << nMissingRhInSC << " of " << nSCRhids;
		    std::cout << " rechits in SC " << it << " w/ OOT" << oot << " && ORG " << orig << " rawE " << sce  << std::endl;
		}//<<>>if( nMissingRhInSC > 0 )
		hist1d[8]->Fill(nSCRhids);
        hist1d[9]->Fill(sceta);

	}//<<>>for( int it = 0; it < nSCs; it++ )
	if( checkRHs ){ 
		float rhpsc = ( nSCs > 0 ) ? nRechitsPerSC/nSCs : 0 ;
		std::cout << " Ave# rechits/SC : " << rhpsc << std::endl; 
	}//<<>>if( checkRHs )
	if( nMissingRechits > 0 ) std::cout << " !! ------ !! Missing " << nMissingRechits << " rechits in " << nSCs << " SCs " << std::endl;
	// fill
    auto nETs = Track_pt->size();
    for( int it = 0; it < nETs; it++ ){ hist1d[7]->Fill( (*Track_pt)[it], 1 ); }

}//<<>>void KUCMSAodSkimmer::processRechits()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setEcalBranches( TTree* fOutTree ){

	// nada as of yet

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

