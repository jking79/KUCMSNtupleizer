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


void KUCMSAodSkimmer::processTimeSig(){

	if( DEBUG ) std::cout << "Finding BC information for photons" << std::endl;
	TSPhoInfo.clearBranches();

    int selPhoIndex = -1;

	uInt nPhotons = Photon_excluded->size();
	if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons to get BC info" << std::endl;
	int selPhoBHCIdx = 0;
	for( uInt it = 0; it < nPhotons; it++ ){

		if( not isBaseLinePho[it] ) continue;
		selPhoIndex++;

        int scIndx = (*Photon_scIndex)[it]; // vector<int>
		float timesignum, timesigdenom;
		float wttimesig = getTimeSig( scIndx, timesignum, timesigdenom );

        TSPhoInfo.fillBranch("selPho_selPhoWTimeSig", wttimesig );

	}//<<>>for( uInt it = 0; it < nPhotons; it++ )

}//<<>>void KUCMSAodSkimmer::processMLPhotons()

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setTSBranches( TTree* fOutTree ){

	std::cout << " - Making Branches for TSig." << std::endl;
    TSPhoInfo.makeBranch("selPho_selPhoWTimeSig", VFLOAT);

    TSPhoInfo.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

//--------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------
// Helper functions to calculate time significance
// -----------------------------------------------------------------------------------------------------------

//auto scIndx = (*Photon_scIndex)[it]; // vector<int> 
//auto rhids = (*SuperCluster_rhIds)[scIndx]; // vector<vector<unsigned int> > 
//float KUCMSAodSkimmer::getTimeSig( int scIndex ){
//
//	vector<vector<unsigned int>> rhindexs( (*SuperCluster_rhIds)[scIndex] );
//	return getTimeSig( rhindexs );
//
//}//<<>>void KUCMSAodSkimmer::getTimeSig

//float KUCMSAodSkimmer::getTimeSig( vector<vector<unsigned int>> rhids ){
float KUCMSAodSkimmer::getTimeSig( int scIndex, float& num, float& denom, const map<unsigned int, float>& rhIdToBHCw){

		auto rhids = (*SuperCluster_rhIds)[scIndex];

        double sumtw = 0;
        double sumw = 0;
		double sumrh = 0;
        int nSCRecHits = rhids.size();
        for( int sciter = 0; sciter < nSCRecHits; sciter++  ){
            auto scrhid = rhids[sciter];
            int erhiter = ( rhIDtoIterMap.find(scrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[scrhid] : -1;
            if( erhiter != -1 ){

                double erhct = erh_corTime[erhiter];
                float erx = (*ECALRecHit_rhx)[erhiter];
                float ery = (*ECALRecHit_rhy)[erhiter];
                float erz = (*ECALRecHit_rhz)[erhiter];
                double ertres = erh_timeRes[erhiter];
				float terror = (*ECALRecHit_timeError)[erhiter];
				float eta = (*ECALRecHit_eta)[erhiter];
				bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];

                double cor_cms000 = hypo(erx,ery,erz)/SOL;
                double cor_tofPVtoRH = hypo(erx-PV_x,ery-PV_y,erz-PV_z)/SOL;
                double ertoftime = erhct - cor_cms000 + cor_tofPVtoRH;

				bool isEE = fabs((*ECALRecHit_eta)[erhiter]) > 1.479;
				bool isValid = true;
				if( terror > 900 ) isValid = false;
				if( hasGainSwitch ){
					if( isEE ){ if( terror < 1.0 ) isValid = false; }
					else { if( terror < 0.6 ) isValid = false; } 
				} else {
                	if( isEE ){ if( terror < 1.05 ) isValid = false; }
                    else { if( terror < 0.65 ) isValid = false; }
				}//<<>>if( hasGainSwitch )

                float gainwt = 0;
                if( isValid ) gainwt = 1;

                double invertres = ( ertres > 0 ) ? 1/ertres : -1;
                double erhar = ( invertres > 0 ) ? invertres*gainwt : 0;
				if(rhIdToBHCw.size() > 0){
					double bhcw = 0;
					if(rhIdToBHCw.find(scrhid) != rhIdToBHCw.end()) bhcw = rhIdToBHCw.at(scrhid);
					erhar *= bhcw;
				}//<<>>if(rhIdToBHCw.size() > 0) 
                sumtw += erhar*ertoftime;
				sumrh += gainwt;
                sumw += erhar;

            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

        if( sumw == 0 ){ sumw = 1; sumtw = -70.711; }
        float phoWTime = sumtw/sumw;
        float phoWVar = sumrh/sumw;
        float wttimesig = phoWTime/std::sqrt(phoWVar);
	    num = phoWTime;
	    denom = std::sqrt(phoWVar);
		return wttimesig;

}//<<>>void KUCMSAodSkimmer::getTimeSig

/*     
 *  how to calculate time significance
 *  --------------------------------------------------------
 *
        auto rhids = (*SuperCluster_rhIds)[scIndx]; // vector<vector<unsigned int> > 

        float sumtw = 0;
        float sumw = 0;
        float sumrh = 0;
        int nSCRecHits = rhids.size();
        for( int sciter = 0; sciter < nSCRecHits; sciter++  ){
            auto scrhid = rhids[sciter];
            int erhiter = ( rhIDtoIterMap.find(scrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[scrhid] : -1;
            if( erhiter != -1 ){

                if( DEBUG ) std::cout << " " << erht << " " << rhx << " " << rhy << " " << rhz << std::endl;

                float erhct = erh_corTime[erhiter];
                float erx = (*ECALRecHit_rhx)[erhiter];
                float ery = (*ECALRecHit_rhy)[erhiter];
                float erz = (*ECALRecHit_rhz)[erhiter];
                float ertres = erh_timeRes[erhiter];
                float cor_cms000 = hypo(erx,ery,erz)/SOL;
                float cor_tofPVtoRH = hypo(erx-PV_x,ery-PV_y,erz-PV_z)/SOL;
                float ertoftime = erhct - cor_cms000 + cor_tofPVtoRH;
                bool isValid = (*ECALRecHit_isTimeValid)[erhiter];
                bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];
                float gainwt = 1;
                if( hasGainSwitch ) gainwt = 0;
                if( not isValid ) gainwt = 0;
                float invertres = 1/ertres;
                float erhar = invertres*gainwt;
                sumtw += erhar*ertoftime;
                sumw += erhar;
                sumrh += gainwt;

            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
            
        if( sumw == 0 ){ sumw = 1; sumtw = -70.711; }
        float phoWTime = sumtw/sumw;
        float phoWVar = sumrh/sumw;
        float wttimesig = phoWTime/std::sqrt(phoWVar);
        //float phoWRes = std::sqrt(phoWVar/2);
        //float wttimesig = phoWTime/phoWRes;

*
*   --------------------------------------------------------------
*/

