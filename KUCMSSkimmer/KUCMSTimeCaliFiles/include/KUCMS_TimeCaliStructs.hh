////////////////////////////////////////////////////////////////////// 
// -*- C++ -*- 
// 
// 
// Original Author:  Jack W King III 
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT 
// 
////////////////////////////////////////////////////////////////////// 
 
#ifndef KUCMS_TimeCalibrationStructs 
#define KUCMS_TimeCalibrationStructs 
 
//--------------------   hh file ------------------------------------------------------------- 
//--------------------------------------------------------------------------------------------- 

//--------------------------------------------------------------------------------------------------------------------- 

struct kucms_lumiRunStruct { 
 
    kucms_lumiRunStruct(){} 
    kucms_lumiRunStruct( int trun, int tfill, float tlumi ) 
        : run(trun), fill(tfill), lumi(tlumi) {} 
 
    int run; 
    int fill; 
    float lumi; 
 
};//<<>>lumiRunStruct 

//--------------------------------------------------------------------------------------------------------------------- 

struct kucms_resTagStruct { 
 
    kucms_resTagStruct(){} 
    kucms_resTagStruct( 
		float tebnoise, float tebstoch, float tebstant, 
        float tebg1he1res, float tebg1he2res, float tebg1he3res, 
        float tebhg1res, float tebhg2res, float tebhg3res, 
        float teenoise, float teestoch, float teestant, 
        float teeg1he1res, float teeg1he2res, float teeg1he3res,
        float teehg1res, float teehg2res, float teehg3res )
        : ebnoise(tebnoise), ebstoch(tebstoch), ebstant(tebstant), 
		ebg1he1res(tebg1he1res), ebg1he2res(tebg1he2res), ebg1he3res(ebg1he3res),
		ebhg1res(tebhg1res), ebhg2res(tebhg2res), ebhg3res(tebhg3res),
        eenoise(teenoise), eestoch(teestoch), eestant(teestant), 
        eeg1he1res(teeg1he1res), eeg1he2res(teeg1he2res), eeg1he3res(eeg1he3res),
        eehg1res(teehg1res), eehg2res(teehg2res), eehg3res(teehg3res) {}
 
    float ebnoise; 
    float ebstoch; 
    float ebstant; 
	float ebg1he1res;
    float ebg1he2res;
    float ebg1he3res;
    float ebhg1res;
    float ebhg2res;
    float ebhg3res;

    float eenoise;
    float eestoch;
    float eestant;
    float eeg1he1res;
    float eeg1he2res;
    float eeg1he3res;
    float eehg1res;
    float eehg2res;
    float eehg3res;

};//<<>>resTagStruct 

//--------------------------------------------------------------------------------------------------------------------- 

struct kucms_caliHistStruct { 
     
    kucms_caliHistStruct(){} 
    kucms_caliHistStruct( TH2F* thist, std::string histname, bool isnew, bool reshist ) 
        : h2f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h1f = NULL; } 
    kucms_caliHistStruct( TH1F* thist, std::string histname, bool isnew, bool reshist ) 
        : h1f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h2f = NULL; } 
 
    TH2F* h2f; 
    TH1F* h1f; 
    std::string histName; 
    bool isNew; 
    bool isResHist; 
 
};//<<>>struct caliHistStruct 

//---------------------------------------------------------------------------------------------------------------------  

// fit result struct 
struct kucms_TimeFitResult { 

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

struct kucms_SigmaFitResult { 

    float noise;//noise 
    float enoise; 
    float stant;//constant 
    float estant; 
    float stoch;//stochastic 
    float estoch; 

};//<<>>struct SigmaFitResult 

//---------------------------------------------------------------------------------------------------------------------   
#endif

