////////////////////////////////////////////////////////////////////// 
// -*- C++ -*- 
// 
// 
// Original Author:  Jack W King III 
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT 
// 
////////////////////////////////////////////////////////////////////// 
 
#ifndef KUCMSTimeCalibrationStructs 
#define KUCMSTimeCalibrationStructs 
 
//--------------------   hh file ------------------------------------------------------------- 
//--------------------------------------------------------------------------------------------- 

//--------------------------------------------------------------------------------------------------------------------- 

struct lumiRunStruct { 
 
    lumiRunStruct(){} 
    lumiRunStruct( int trun, int tfill, float tlumi ) 
        : run(trun), fill(tfill), lumi(tlumi) {} 
 
    int run; 
    int fill; 
    float lumi; 
 
};//<<>>lumiRunStruct 

//--------------------------------------------------------------------------------------------------------------------- 

struct smearTagStruct { 
 
    smearTagStruct(){} 
    smearTagStruct( float tnoise, float tstoch, float tstant ) 
        : noise(tnoise), stoch(tstoch), stant(tstant) {} 
 
    float noise; 
    float stoch; 
    float stant; 
 
};//<<>>smearTagStruct 

//--------------------------------------------------------------------------------------------------------------------- 

struct caliHistStruct { 
     
    caliHistStruct(){} 
    caliHistStruct( TH2F* thist, std::string histname, bool isnew, bool reshist ) 
        : h2f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h1f = NULL; } 
    caliHistStruct( TH1F* thist, std::string histname, bool isnew, bool reshist ) 
        : h1f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h2f = NULL; } 
 
    TH2F* h2f; 
    TH1F* h1f; 
    std::string histName; 
    bool isNew; 
    bool isResHist; 
 
};//<<>>struct caliHistStruct 

//---------------------------------------------------------------------------------------------------------------------  

// fit result struct 
struct TimeFitResult { 

    float mu; 
    float emu; 
    float chi2ndf; 
    float chi2prob; 
    float sigma; 
    float sigmaHigh; 
    float sigmaLow; 
    float esigma; 
    float occ; 
    float rms; 
    float std; 

};//<<>>struct TimeFitResult 

//---------------------------------------------------------------------------------------------------------------------   

struct SigmaFitResult { 

    float noise;//noise 
    float enoise; 
    float stant;//constant 
    float estant; 
    float stoch;//stochastic 
    float estoch; 

};//<<>>struct SigmaFitResult 

//---------------------------------------------------------------------------------------------------------------------   
#endif

