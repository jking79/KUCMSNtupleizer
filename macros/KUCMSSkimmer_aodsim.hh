//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

//#include "llpgana_hist_base_v2.hh" 
//#include "llpgana_ntuple_base.hh"
//#include "llpgana_hist_base_v10.hh"
//#include "llpgana_hist_rebase_v11.hh"
#include "llpgana_hist_rebase_v18.hh"

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

#define DEBUG false
//#define DEBUG true

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSSkimmer class -----------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class KUCMSSkimmer : root_rebase {

	public:

	void kucmsSkimmer( std::string indir, std::string infilelist, std::string outfilename, int pct );	
	void initHists();
	void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	

    std::map<UInt_t,DetIDStruct> DetIDMap;

	int getRhIdx( uInt rhDetID );
	uInt getLeadRhID( vector<uInt> recHitIds );
	float clstrR9( vector<uInt> recHitIds );
	vector<float> getLeadTofRhTime( vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
	vector<float> getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids );

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

	int nMaps;
	bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( vector<unsigned int> rhcol );

};

//--------------------------------------------------------------------------------------------------------------------------------
