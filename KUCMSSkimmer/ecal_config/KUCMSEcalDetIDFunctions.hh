// -*- C++ -*-
//
//
// Author:  Jack W King III 
// Created:  Wed, 18 Sept 2024 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// system include files
#include <memory>

// basic C++ types
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <tuple>
#include <random>
#include <sys/stat.h>
#include <cmath>

#ifndef KUCMSEcalDetIDFunctions
#define KUCMSEcalDetIDFunctions

enum ECAL {EB, EM, EP, NONE};

struct DetIDStruct {

    DetIDStruct() {}
    DetIDStruct( const int ni1, const int ni2, const int nTT, const int & necal, const float nphi, const float neta ) 
		: i1(ni1), i2(ni2), TT(nTT), ecal(necal), phi(nphi), eta(neta) {}

    int i1; // EB: iphi, EE: ix
    int i2; // EB: ieta, EE: iy
    int TT; // trigger tower
    int ecal; // EB, EM, EP
	float phi; // xtal phi
	float eta; // xtal eta

};//<<>>struct DetIDStruct

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap ){

	const std::string inDir("/uscms/home/jaking/nobackup/el8/ecaltiming/CMSSW_14_0_11/src/GammaResTool/GammaResTool/macros/ecal_config/");
    const std::string detIDConfigEB("fullinfo_v2_detids_EB.txt");
    std::ifstream infile( inDir+detIDConfigEB, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
	float phi, eta;
    std::string pos;

    while( infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos 
				>> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM >> phi >> eta ){ 

		DetIDMap[cmsswId] = {iphi,ieta,TT25,0,phi,eta}; 

	}//<<>>while (infile >>

}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap ){

    const std::string inDir("/uscms/home/jaking/nobackup/el8/ecaltiming/CMSSW_14_0_11/src/GammaResTool/GammaResTool/macros/ecal_config/");
    const std::string detIDConfigEE("fullinfo_v2_detids_EE.txt");
    std::ifstream infile( inDir+detIDConfigEE, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
	float phi, eta;
    std::string EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC 
				>> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant >> phi >> eta){

        int ec = 1;
        if( side > 0 ) ec = 2;
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec,phi,eta};

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )

#endif
//-------------------------------------------------------------------------------------------------------------------

