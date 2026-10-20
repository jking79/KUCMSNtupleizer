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
		ebg1he1res(tebg1he1res), ebg1he2res(tebg1he2res), ebg1he3res(tebg1he3res),
		ebhg1res(tebhg1res), ebhg2res(tebhg2res), ebhg3res(tebhg3res),
        eenoise(teenoise), eestoch(teestoch), eestant(teestant), 
        eeg1he1res(teeg1he1res), eeg1he2res(teeg1he2res), eeg1he3res(teeg1he3res),
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

    float mu = 0;
    float emu = 0;
    float chi2 = 0;
    float chi2ndf = 0;
    float chi2prob = 0;
    float sigma = 0;
    float sigmaHigh = 0;
    float sigmaLow = 0;
    float esigma = 0;
    float modelError = 0;
    float occ = 0;
    float rms = 0;
    float std = 0;
    int ndf = 0;
    int fitStatus = -1;
    int fitStatusHigh = -1;
    int fitStatusLow = -1;
    int covStatus = -1;
    int covStatusHigh = -1;
    int covStatusLow = -1;
    int modelVariations = 0;
    bool sigmaAtLimit = false;
    bool sigmaHighAtLimit = false;
    bool sigmaLowAtLimit = false;
    bool valid = false;
    bool validHigh = false;
    bool validLow = false;
    bool modelValid = false;

};//<<>>struct TimeFitResult 

//---------------------------------------------------------------------------------------------------------------------   

struct kucms_SigmaFitResult { 

    float noise = 0;//noise
    float enoise = 0;
    float stant = 0;//constant
    float estant = 0;
    float stoch = 0;//stochastic
    float estoch = 0;
    float chi2 = 0;
    float chi2ndf = 0;
    float chi2prob = 0;
    int ndf = 0;
    int fitStatus = -1;
    int covStatus = -1;
    int nFitPoints = 0;
    float fractionalScatter = 0;
    float efractionalScatter = 0;
    float objective = 0;
    bool parameterAtLimit = false;
    bool valid = false;

};//<<>>struct SigmaFitResult 

//---------------------------------------------------------------------------------------------------------------------   
#endif

