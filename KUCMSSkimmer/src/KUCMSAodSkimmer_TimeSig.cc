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
                double ertoftime = erhct + cor_cms000 - cor_tofPVtoRH;

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

        if( sumw == 0 ){ sumw = 1; sumtw = -25; sumrh = 0; }
        float phoWTime = sumtw/sumw;
        float phoWVar = 1/sumw;
    	float phoWHVar = sumrh/sumw;
		float sqrtPhoWVar = std::sqrt(phoWVar);
        float wttimesig = phoWTime/sqrtPhoWVar;
	    num = phoWTime;
	    denom = sqrtPhoWVar;
		return wttimesig;

}//<<>>void KUCMSAodSkimmer::getTimeSig

float KUCMSAodSkimmer::getTimeSig( const std::vector<int>& scIndexs, float& num, float& denom, const std::vector<float>& eledelay ){

    double sumtw = 0;
    double sumw = 0;
    double sumrh = 0;
	std::vector<int> usedrechits;
	int nSCs = scIndexs.size();
	if( nSCs == 0 ){ num = -41; denom = -2; return -41; }
	int nRHUsed = 0;
	for( int it = 0; it < nSCs; it++  ){

		int scIndex = scIndexs[it]; 
        auto rhids = (*SuperCluster_rhIds)[scIndex];
        int nSCRecHits = rhids.size();
		float adjust = ( eledelay.size() > 0 ) ? eledelay[it] : 0;
		if( nSCRecHits < 5 ) continue;
        for( int sciter = 0; sciter < nSCRecHits; sciter++  ){
            auto scrhid = rhids[sciter];
            int erhiter = ( rhIDtoIterMap.find(scrhid) != rhIDtoIterMap.end() ) ? rhIDtoIterMap[scrhid] : -1;
            if( erhiter != -1 ){

				// check if rechit is exclusive ?
				if( std::count( usedrechits.begin(), usedrechits.end(), erhiter ) > 0 ) continue;
				usedrechits.push_back( erhiter );

                double erhct = erh_corTime[erhiter];
				if( std::abs(erhct) > 7.5 ) continue; // aviod spike region w/o over weighted delayed tail 
                float erx = (*ECALRecHit_rhx)[erhiter];
                float ery = (*ECALRecHit_rhy)[erhiter];
                float erz = (*ECALRecHit_rhz)[erhiter];
                double ertres = erh_timeRes[erhiter];
                float terror = (*ECALRecHit_timeError)[erhiter];
                float eta = (*ECALRecHit_eta)[erhiter];
                bool hasGainSwitch = (*ECALRecHit_hasGS1)[erhiter] || (*ECALRecHit_hasGS6)[erhiter];

                double cor_cms000 = hypo(erx,ery,erz)/SOL;
                double cor_tofPVtoRH = hypo(erx-PV_x,ery-PV_y,erz-PV_z)/SOL;
                double ertoftime = erhct + cor_cms000 - cor_tofPVtoRH - adjust;

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
				nRHUsed += gainwt;

                double invertres = ( ertres > 0 ) ? 1/ertres : -1;
                double erhar = ( invertres > 0 ) ? invertres*gainwt : 0;
                sumtw += erhar*ertoftime;
                sumrh += gainwt;
                sumw += erhar;

            }//<<>>if( ecalrhiter != -1 )
        }//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )

	}//<<>>for( int scIndex = 0; scIndex < nSCs; scIndex++  )

	if( nRHUsed < 5 ){ num = -43; denom = -1; return -43; }
    if( sumw == 0 ){ num = -42; denom = -1; return -42; }
    float phoWTime = sumtw/sumw;
    float phoWHVar = sumrh/sumw;
    float phoWVar = 1/sumw;
    float sqrtPhoWVar = std::sqrt(phoWVar);
    float wttimesig = phoWTime/sqrtPhoWVar;
    num = phoWTime;
    denom = sqrtPhoWVar;
    return wttimesig;

}//<<>>void KUCMSAodSkimmer::getTimeSig

/*
//should be deltaT = gam - pv
double KUCMSAodSkimmer::CalcGenDeltaT(const Jet& pho){

		//cout << "CalcGenDeltaT - start" << endl;
		//calc times differently for !sig and sig photons
		double dpho = -999;
		int genidx, phoidx, phoid;
		//gen photon coordinates
		double genx, geny, genz, gentheta, geneta, genphi, vx, vy, vz;

		//if no match
		phoidx = pho.GetUserIdx();
		genidx = _base->Photon_genIdx->at(phoidx);
		//cout << "phoidx " << phoidx << " genidx " << genidx << endl;
		if(genidx == -1) phoid = -1;
		else phoid = _base->Gen_susId->at(genidx);
		if(phoid == -1) return dpho;

		//photon production vertex coordinates
		vx = _base->Gen_vx->at(genidx);
		vy = _base->Gen_vy->at(genidx);
		vz = _base->Gen_vz->at(genidx);		
		//momentum components
		double px = _base->Gen_px->at(genidx);
		double py = _base->Gen_py->at(genidx);
		double pz = _base->Gen_pz->at(genidx);		

		//energy
		double e = _base->Gen_energy->at(genidx);

		//find path length from photon production vertex to detector
		double R = 129; //radius of ECAL
		double vt = sqrt(vx*vx + vy*vy); //(l) distance to prod vertex in transverse plane from (0, 0)
		double pt = _base->Gen_pt->at(genidx); //momentum in transverse plane (for betaT)
		double dot = px*vx + py*vy; //pT \dot vT
		double L = sqrt(R*R - vt*vt + (dot/pt)*(dot/pt)) - dot/pt; //flight path from prod vertex to detector
		
		double betaT = pt/e; //beta in transverse plane (would be 1 for photons with z component of momentum)
		double tof = L/(_c*betaT); //time of flight to detector along calculated path length

		double genx_ECAL = vx + (px/e)*_c*tof;
		double geny_ECAL = vy + (py/e)*_c*tof;
		double genz_ECAL = vz + (pz/e)*_c*tof; 

		double ftheta = atan2(sqrt(genx_ECAL*genx_ECAL + geny_ECAL*geny_ECAL),genz_ECAL);
		double feta = -log(tan(ftheta/2.));
		//calculate TOF for that path (transverse plane only)
		//use TOF to propagate x, y, z from gen vertex position (vi) to detector (geni_ECAL)


		double pvx = _base->PV_x;
		double pvy = _base->PV_y;
		double pvz = _base->PV_z;

		double beta;
		if(phoid == 22){
			int momidx = _base->Photon_genSigXMomId->at(phoidx);
			//cout << "signal - phoid " << phoid << " phoidx " << phoidx << " genidx " << genidx << " momidx " << momidx << endl;
			//TODO: remove when ntuples are fixed (10/9/24)
			if(momidx < 0) return dpho;
			//check gen pdgids
			//want production vertex of photon (where LLP -> photon)
			//not production vertex of mother (close to PV)
			vx = _base->Gen_vx->at(genidx);
			vy = _base->Gen_vy->at(genidx);
			vz = _base->Gen_vz->at(genidx);
		
			//production vertex of mother particle (should be close to PV)
			double momvx, momvy, momvz;
			momvx = _base->Gen_vx->at(momidx);
			momvy = _base->Gen_vy->at(momidx);
			momvz = _base->Gen_vz->at(momidx);

			double mompx, mompy, mompz, momE;			

			mompx = _base->Gen_px->at(momidx);			
			mompy = _base->Gen_py->at(momidx);			
			mompz = _base->Gen_pz->at(momidx);			
			momE = _base->Gen_energy->at(momidx);			

			double p = sqrt(mompx*mompx + mompy*mompy + mompz*mompz);
			double m = _base->Gen_mass->at(momidx);
			double gam1 = sqrt(1 + (p/(m))*(p/(m)));
			double beta1 = sqrt(1 - 1/(gam1*gam1));


			//beta/c = p/E = v
			beta = sqrt(mompx*mompx + mompy*mompy + mompz*mompz)/momE;
			//cout << "beta: " << beta << " vel: " << beta*_c << endl;
			//check gen photon energy	
			//cout << "photon energy: " << _base->Photon_energy->at(phoidx) << endl;
			//distance bw photon and production point (where LLP decays to photon)
			dpho = sqrt( (genx_ECAL - vx)*(genx_ECAL - vx) + (geny_ECAL - vy)*(geny_ECAL - vy) + (genz_ECAL - vz)*(genz_ECAL - vz) )/_c;
		
			//distance bw LLP decay point and PV
			//LLP is produced close to PV (should take into account?)
			//LLPdecay - LLPprod?
			dpho += sqrt( (vx - momvx)*(vx - momvx) + (vy - momvy)*(vy - momvy) + (vz - momvz)*(vz - momvz) )/(_c*beta);	

		}
		//assume prompt production
		else{
			//cout << "prompt" << endl;
			//distance bw photon and production point (PV)
			dpho = sqrt( (genx_ECAL - pvx)*(genx_ECAL - pvx) + (geny_ECAL - pvy)*(geny_ECAL - pvy) + (genz_ECAL - pvz)*(genz_ECAL - pvz) )/_c;
		}
			
		//cout << "CalcGenDeltaT - end" << endl;
		return dpho;

}//<<>>double KUCMSAodSkimmer::CalcGenDeltaT(const Jet& pho)
*/

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

