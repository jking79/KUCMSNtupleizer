//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

// KUCMS
#include "KUCMSItemManager.hh"
//#include "KUCMSRootHelperBaseClass.hh"

#ifndef KUCMSEventSelectionClass
#define KUCMSEventSelectionClass

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

class KUCMSEventSelection {

    public:

    KUCMSEventSelection(){ eventSelectionTag = "none"; };
    ~KUCMSEventSelection(){};

	private:

	std::string eventSelectionTag;
	std::map< std::string, bool > selectionMap;  
	bool metFilter;
	bool phoFilter;

	public:

	void setEventSelectionTag( std::string estag ){ eventSelectionTag = estag; };	
	const bool selectEvent( ItemManager<float>& geVar );
	const bool trackMetFilter(){ return metFilter; };
    const bool trackPhoFilter(){ return phoFilter; };


};//<<>>class KUCMSEventSelection : KUCMSRootHelperBaseClass

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code  -   yes the stuff you put into BlaClass.cpp files for object libaires
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// make constructor to set eventSelectionTag to none

const bool KUCMSEventSelection::selectEvent( ItemManager<float>& geVar ){


    bool hasMinMet150 = geVar("evtMET") > 150;
    bool hasMinMet125 = geVar("evtMET") > 125;
    bool hasMinMet100 = geVar("evtMET") > 100;	
    bool hasMinMet75 = geVar("evtMET") > 75;
	bool hasAL1SelElectron = geVar("nSelEle") > 0;
	bool hasAL1IsoPhoton = geVar("nIsoPhos") > 0;
	bool hasAL1DisSV = geVar("nDisSVs") > 0;
	bool hasAL1NpSC = geVar("nNpSC") > 0;

	bool hasAl1IsoPhoPlusMinMet100 = hasAL1IsoPhoton && hasMinMet100;
    bool hasDSVAl1IsoPhoPlusMinMet100 = ( hasAL1DisSV || hasAL1IsoPhoton ) && hasMinMet100;
    bool hasDSVAl1IsoPhoPlusMinMet125 = ( hasAL1DisSV || hasAL1IsoPhoton ) && hasMinMet125;

	selectionMap["None"] = true;
	selectionMap["MET150"] = hasMinMet150;
    selectionMap["MET100"] = hasMinMet100;
    selectionMap["MET75"] = hasMinMet75;
    selectionMap["AL1IsoPho"] = hasAL1IsoPhoton;
    selectionMap["IsoPhoMet100"] = hasAl1IsoPhoPlusMinMet100;
    selectionMap["AL1SelEle"] = hasAL1SelElectron;
    selectionMap["SVIPMet100"] = hasDSVAl1IsoPhoPlusMinMet100;
    selectionMap["SVIPMet125"] = hasDSVAl1IsoPhoPlusMinMet125;
    selectionMap["AL1DisSV"] = hasAL1DisSV;
    selectionMap["AL1NpSC"] = hasAL1NpSC;

	metFilter = hasMinMet100;
	phoFilter = hasAL1DisSV || hasAL1DisSV;

	//std::cout << "Event select for ntuple with " << eventSelectionTag << " w/ result " << selectionMap[eventSelectionTag] << std::endl;
	return selectionMap[eventSelectionTag];
	
}//<<>>const bool KUCMSEventSelection::selectEvent( ItemManager<float>& geVar )

#endif
//-------------------------------------------------------------------------------------------------------------------
