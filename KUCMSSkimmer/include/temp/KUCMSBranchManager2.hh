// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>

//KUCMS includes
#include "KUCMSItemManager.hh"
#include "KUCMSBranch2.hh"

//Root includes
#include "TTree.h"

#ifndef KUCMSBranchManagerHeader_hh
#define KUCMSBranchManagerHeader_hh

////////////////////////////////////////////////////////////////////////////////////////////////

class KUCMSBranchManager {

    public:

    KUCMSBranchManager(){};
    ~KUCMSBranchManager();

    void makeBranch( std::string key, std::string name, BType type, std::string doc = "" );
    void makeBranch( std::string name, BType type, std::string doc = "" );
    void attachBranches( TTree* fOutTree );
    void clearBranches();

    void fillBranch( std::string key, std::vector<float> val );
    void fillBranch( std::string key, std::vector<uInt> val );
    void fillBranch( std::string key, uInt val );
    void fillBranch( std::string key, int val );
    void fillBranch( std::string key, float val );
    void fillBranch( std::string key, std::string val );
    void fillBranch( std::string key, bool val );

    private:

    std::map< std::string, KUCMSBranchBase* > theBranches;
    bool valid( std::string key );

};//<<>>class KUCMSBranchManager 

#endif

