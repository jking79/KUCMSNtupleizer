//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "skimHistMaker.hh"

//#define DEBUG true
#define DEBUG false
#define doEBEEmaps false

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, int cut, float va, float vb, float vc, float vd ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

    cutselection = cut;
	cutva = va;
    cutvb = vb;
    cutvc = vc;
    cutvd = vd;
    preCutNPhotons = 0;
    preCut30NPhotons = 0;
    preCut100NPhotons = 0;
    postCutNPhotons = 0;
    postCut30NPhotons = 0;
    postCut100NPhotons = 0;

	scale = 137;

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
    std::ifstream infile(indir+infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    auto fConfigTree = new TChain(configtreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << indir+infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		//std::cout << "--  adding from: " << str << std::endl;
        if( str[0] == '#' ) continue;
        if( str == " " ) continue;
		if( str == "" ) continue;
        auto tfilename = eosdir + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
        fConfigTree->Add(tfilename.c_str());
		cnt++;
    }//<<>>while (std::getline(infile,str))

    std::cout << "Setting up For Main Loop." << std::endl;

	Init(fInTree);
	initHists(htitle);

    std::cout << "Filling Config Map." << std::endl;

    UInt_t          nEvents;
    UInt_t          nSelectedEvents;
    string          *sKey;
    Float_t         sCrossSection;
    Float_t         sGMSBGravMass;
    Float_t         sGMSBChi1Mass;
    Float_t         sMCWgt;
    Int_t           sMCType;
    Float_t         sumEvtGenWgt;

    Int_t           cf_m_gt2jets;
    Int_t           cf_met150;
    Int_t           cf_mj_gt1phos;
    Int_t           cf_mjp_leadPhoPt30;
    Int_t           cf_nFltrdEvts;
    Int_t           cf_nTotEvts;	
	Int_t           cf_sel_j;
   	Int_t           cf_sel_m;
   	Int_t           cf_sel_p;
   	Int_t           cf_sel_ppt;
   	Int_t           nFltrdEvts;
   	Int_t           nTotEvts;
   	Float_t         sumEvtWgt;
   	Float_t         sumFltrdEvtWgt;

    TBranch        *b_nEvents;   //!
    TBranch        *b_nSelectedEvents;   //!
    TBranch        *b_sKey;   //!
    TBranch        *b_sCrossSection;   //!
    TBranch        *b_sGMSBGravMass;   //!
    TBranch        *b_sGMSBChi1Mass;   //!
    TBranch        *b_sMCWgt;   //!
    TBranch        *b_sMCType;   //!
    TBranch        *b_sumEvtGenWgt;

   	TBranch        *b_cf_m_gt2jets;   //!
   	TBranch        *b_cf_met150;   //!
   	TBranch        *b_cf_mj_gt1phos;   //!
   	TBranch        *b_cf_mjp_leadPhoPt30;   //!
   	TBranch        *b_cf_nFltrdEvts;   //!
   	TBranch        *b_cf_nTotEvts;   //!
   	TBranch        *b_cf_sel_j;   //!
   	TBranch        *b_cf_sel_m;   //!
   	TBranch        *b_cf_sel_p;   //!
   	TBranch        *b_cf_sel_ppt;   //!
   	TBranch        *b_nFltrdEvts;   //!
   	TBranch        *b_nTotEvts;   //!
   	TBranch        *b_sumEvtWgt;   //!
   	TBranch        *b_sumFltrdEvtWgt;   //!


    sKey = 0;

    fConfigTree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
    fConfigTree->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
    fConfigTree->SetBranchAddress("sKey", &sKey, &b_sKey);
    fConfigTree->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
    fConfigTree->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
    fConfigTree->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
    fConfigTree->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
    fConfigTree->SetBranchAddress("sMCType", &sMCType, &b_sMCType);
    fConfigTree->SetBranchAddress("sumEvtGenWgt", &sumEvtGenWgt, &b_sumEvtGenWgt);

    fConfigTree->SetBranchAddress("cf_m_gt2jets", &cf_m_gt2jets, &b_cf_m_gt2jets);
    fConfigTree->SetBranchAddress("cf_met150", &cf_met150, &b_cf_met150);
    fConfigTree->SetBranchAddress("cf_mj_gt1phos", &cf_mj_gt1phos, &b_cf_mj_gt1phos);
    fConfigTree->SetBranchAddress("cf_mjp_leadPhoPt30", &cf_mjp_leadPhoPt30, &b_cf_mjp_leadPhoPt30);
    fConfigTree->SetBranchAddress("cf_nFltrdEvts", &cf_nFltrdEvts, &b_cf_nFltrdEvts);
    fConfigTree->SetBranchAddress("cf_nTotEvts", &cf_nTotEvts, &b_cf_nTotEvts);
    fConfigTree->SetBranchAddress("cf_sel_j", &cf_sel_j, &b_cf_sel_j);
    fConfigTree->SetBranchAddress("cf_sel_m", &cf_sel_m, &b_cf_sel_m);
    fConfigTree->SetBranchAddress("cf_sel_p", &cf_sel_p, &b_cf_sel_p);
    fConfigTree->SetBranchAddress("cf_sel_ppt", &cf_sel_ppt, &b_cf_sel_ppt);
    fConfigTree->SetBranchAddress("nFltrdEvts", &nFltrdEvts, &b_nFltrdEvts);
    fConfigTree->SetBranchAddress("nTotEvts", &nTotEvts, &b_nTotEvts);
    fConfigTree->SetBranchAddress("sumEvtWgt", &sumEvtWgt, &b_sumEvtWgt);
    fConfigTree->SetBranchAddress("sumFltrdEvtWgt", &sumFltrdEvtWgt, &b_sumFltrdEvtWgt);

    auto nConfigEntries = fConfigTree->GetEntries();
    std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
	configInfo.clear();
    for (Long64_t centry = 0; centry < nConfigEntries; centry++){
		
		auto entry = fConfigTree->LoadTree(centry);

        if(debug) std::cout << " - Getting Branches. " << std::endl;
    	b_nEvents->GetEntry(entry);   //!
    	b_nSelectedEvents->GetEntry(entry);   //!
    	b_sKey->GetEntry(entry);   //!  
    	b_sCrossSection->GetEntry(entry);   //!
    	b_sGMSBGravMass->GetEntry(entry);   //!
    	b_sGMSBChi1Mass->GetEntry(entry);   //!
    	b_sMCWgt->GetEntry(entry);   //!
    	b_sMCType->GetEntry(entry);   //!
        std::string configKey(*sKey);

        b_sumEvtGenWgt->GetEntry(entry);
        b_cf_m_gt2jets->GetEntry(entry);   //!
        b_cf_met150->GetEntry(entry);   //!
        b_cf_mj_gt1phos->GetEntry(entry);   //!
        b_cf_mjp_leadPhoPt30->GetEntry(entry);   //!
        b_cf_nFltrdEvts->GetEntry(entry);   //!
        b_cf_nTotEvts->GetEntry(entry);   //!
        b_cf_sel_j->GetEntry(entry);   //!
        b_cf_sel_m->GetEntry(entry);   //!
        b_cf_sel_p->GetEntry(entry);   //!
        b_cf_sel_ppt->GetEntry(entry);   //!
        b_nFltrdEvts->GetEntry(entry);   //!
        b_nTotEvts->GetEntry(entry);   //!
        b_sumEvtWgt->GetEntry(entry);   //!
        b_sumFltrdEvtWgt->GetEntry(entry);   //!

		if( not configInfo.count(configKey) ){
        	if(debug) std::cout << " - Filling Cutflow/Wieghts. " << std::endl;
        	//fConfigTree->Add(tfilename.c_str());			
			std::map< std::string, float > configValues;
        	configValues["nEvents"] = nEvents;
        	configValues["nSelectedEvents"] = nSelectedEvents;
            configValues["sumEvtGenWgt"] = sumEvtGenWgt;
            configValues["sumEvtWgt"] = sumEvtWgt;
            if(debug) std::cout << " - Filling configValues. " << std::endl;
        	configValues["sCrossSection"] = sCrossSection;
        	configValues["sGMSBGravMass"] = sGMSBGravMass;
        	configValues["sGMSBChi1Mass"] = sGMSBChi1Mass;
        	configValues["sMCWgt"] = sMCWgt;
        	configValues["sMCType"] = sMCType;
            configValues["nTotEvts"] = nTotEvts;
        	if(debug) std::cout << " - Filling configInfo. " << std::endl;
        	configInfo[configKey] = configValues;
		} else {
			auto & configValues = configInfo[configKey];
			configValues["nEvents"] += nEvents;
            configValues["nTotEvts"] += nTotEvts;
			configValues["nSelectedEvents"] += nSelectedEvents;
            configValues["sumEvtGenWgt"] += sumEvtGenWgt;
            configValues["sumEvtWgt"] += sumEvtWgt;
		}//<<>>if( not configInfo.count(configKey) )

		float fillwt = scale * ( sCrossSection * 1000 ) * ( 1 / sumEvtWgt );
        if( not cutflowInfo.count("nTotEvts") ){
            if(debug) std::cout << " - Filling Cutflow/Wieghts. " << std::endl;
            cutflowInfo["cf_m_gt2jets"] = cf_m_gt2jets*fillwt;
            cutflowInfo["cf_met150"] = cf_met150*fillwt;
            cutflowInfo["cf_mj_gt1phos"] = cf_mj_gt1phos*fillwt;
            cutflowInfo["cf_mjp_leadPhoPt30"] = cf_mjp_leadPhoPt30*fillwt;
            cutflowInfo["cf_nFltrdEvts"] = cf_nFltrdEvts*fillwt;
            cutflowInfo["cf_nTotEvts"] = cf_nTotEvts*fillwt;
            cutflowInfo["cf_sel_j"] = cf_sel_j*fillwt;
            cutflowInfo["cf_sel_m"] = cf_sel_m*fillwt;
            cutflowInfo["cf_sel_p"] = cf_sel_p*fillwt;
            cutflowInfo["cf_sel_ppt"] = cf_sel_ppt*fillwt;
            cutflowInfo["nFltrdEvts"] = nFltrdEvts*fillwt;
            cutflowInfo["nTotEvts"] = nTotEvts*fillwt;
            //cutflowWgts["sumEvtWgt"] = sumEvtWgt;
            //cutflowWgts["sumFltrdEvtWgt"] = sumFltrdEvtWgt;
        } else {
            cutflowInfo["cf_m_gt2jets"] += cf_m_gt2jets*fillwt;
            cutflowInfo["cf_met150"] += cf_met150*fillwt;
            cutflowInfo["cf_mj_gt1phos"] += cf_mj_gt1phos*fillwt;
            cutflowInfo["cf_mjp_leadPhoPt30"] += cf_mjp_leadPhoPt30*fillwt;
            cutflowInfo["cf_nFltrdEvts"] += cf_nFltrdEvts*fillwt;
            cutflowInfo["cf_nTotEvts"] += cf_nTotEvts*fillwt;
            cutflowInfo["cf_sel_j"] += cf_sel_j*fillwt;
            cutflowInfo["cf_sel_m"] += cf_sel_m*fillwt;
            cutflowInfo["cf_sel_p"] += cf_sel_p*fillwt;
            cutflowInfo["cf_sel_ppt"] += cf_sel_ppt*fillwt;
            cutflowInfo["nFltrdEvts"] += nFltrdEvts*fillwt;
            cutflowInfo["nTotEvts"] += nTotEvts*fillwt;
            //cutflowWgts["sumEvtWgt"] += sumEvtWgt;
            //cutflowWgts["sumFltrdEvtWgt"] += sumFltrdEvtWgt;
        }//<<>>if( not configInfo.count(configKey) )

	}//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)

    for( auto item : configInfo ){ 
		std::cout << item.first << " ( ";  
		for( auto line : item.second ){ 
			std::cout << line.first <<  " " << line.second << " ";
		}//<<>>for( auto line : item.second )
		std::cout << ")" << std::endl;
	}//<<>>for( auto item : configInfo )


    std::cout << "<<<<<<<< Processing Event Loop <<<<<<<<<<<<<< " << std::endl;
	int loopCounter(100000);
    auto nEntries = fInTree->GetEntries();
    if(debug){ nEntries = 10; loopCounter = 1; }
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
        if(debug) std::cout << "*****************************************************************************" << std::endl;
        auto entry = fInTree->LoadTree(centry);
		if(debug) std::cout << " - getBranches " << std::endl;
		getBranches(entry);
		if(debug) std::cout << " - eventLoop " << std::endl;
		eventLoop(entry);
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
   
    if(debug) std::cout << " - Creating output file " << std::endl;
    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();
	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

	nMaps = 0;
	if( doEBEEmaps ){ for( int it = 0; it < nEBEEMaps; it++ ){ 
		ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
		ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
		ebeeMapR[it]->Write(); delete ebeeMapR[it];
	}}//<<>>for( int it = 0; it < nEBEEMaps; it++ )

    fOutFile->Close();
    std::cout << "histMaker : Thats all Folks!!" << std::endl;
}//<<>>void kucmsSkimmer
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void HistMaker::eventLoop( Long64_t entry ){

	int cs = ( cutselection == 0 ) ? 0 : 1;

    //auto dskey  = *DataSetKey
    //float evtgwt = evtGenWgt;
    float evtgwt = 1;
    float scale = 137;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    //float segwt = (configInfo[configKey])["nEvents"];
    //float segwt = (configInfo[configKey])["nTotEvts"];
    float segwt = (configInfo[configKey])["sumEvtWgt"];
	////float fillwt = scale * ( sCrossSection * 1000 ) * ( 1 / sumEvtWgt );
    float fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
	//if( cutselection == 1 ) fillwt = 1;
    if( DEBUG ) std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
    //if( configKey != "GMSBL"+lambda+"" ) return;

    hist1d[2]->Fill(7,fillwt);

	for( int i = 0; i < 1; i++ ){ // continue loop - can use continue to skip
	if( rjrNRjrPhotons->size() == 0 ) continue;
	//if( selPhoTime->at(0) < 2.0 ) continue;
    hist1d[2]->Fill(8,fillwt);

	//bool genSigTypeCut( genSigType != 11 && genSigType != 12 && genSigType != 22 );

	bool nRjrPhoCut( rjrNRjrPhotons->at(cs) != cutvb );
    if( cutvb == 0 ) nRjrPhoCut = false;

    bool genSigTypeCut( genSigType != cutva );
    if( cutva == 50 ) genSigTypeCut = genSigType == 11 || genSigType == 12 || genSigType == 22;
	else if( cutva == 30 ) genSigTypeCut = genSigType != 11 && genSigType != 12 && genSigType != 22;
	else if( cutva == 0 ) genSigTypeCut = false;

    if( rjrNRjrPhotons->at(cs) == 1 ) hist1d[2]->Fill(9,fillwt);
    if( rjrNRjrPhotons->at(cs) == 2 ) hist1d[2]->Fill(10,fillwt);
    //if( genSigTypeCut ) continue;
    //hist1d[2]->Fill(ecfbin++,fillwt);
    if( nRjrPhoCut ) continue;
    hist1d[2]->Fill(11,fillwt);
	//if( genSigTypeCut || nRjrPhoCut ) continue;
    //hist1d[2]->Fill(ecfbin++,fillwt);

	//if( rjrSMass->at(0) < 2000 ) continue;
	//if( rjrX2NGMean->at(cs) < cutvc ) continue;
    //if( rjrX2GMean->at(cs) < cutvd ) continue;
	if( DEBUG ) std::cout << "RJR cut at " << cutvc << " with value " << rjrX2GMean->at(0) << std::endl;

    auto metCPt = hypo(selCMetPx,selCMetPy);
    auto metPt = hypo(selMetPx,selMetPy);
	if( metCPt < cutvd ) continue;
    hist1d[2]->Fill(12,fillwt);

	//exploritory cut flows
	if( nSelPhotons < 3 ) hist1d[2]->Fill(13,fillwt);
    if( (*rjrNJetsJa)[cs] > 0 && (*rjrNJetsJb)[cs] > 0 ) hist1d[2]->Fill(14,fillwt);
    if( (*rjrNJetsJa)[cs] > 1 && (*rjrNJetsJb)[cs] > 1 ) hist1d[2]->Fill(15,fillwt);
    if( (*rjrNJetsJa)[cs] > 2 && (*rjrNJetsJb)[cs] > 2 ) hist1d[2]->Fill(16,fillwt);

	if( (*rjrNJetsJa)[cs] < 1 || (*rjrNJetsJb)[cs] < 1 ) continue;
    hist1d[2]->Fill(17,fillwt);

	//final cutflow bin
    hist1d[2]->Fill(25,fillwt);//lastbin

	//var hist fill
	
    hist1d[400]->Fill( selCMet, fillwt );
    hist1d[401]->Fill( metCPt, fillwt );
    hist1d[402]->Fill( rjrMET->at(cs), fillwt );
    hist1d[403]->Fill( metPt, fillwt );

    if( genSGMass->size() > 0 ) hist1d[0]->Fill( (*genSGMass)[0] ); else hist1d[0]->Fill( 0 );
    if( genSQMass->size() > 0 ) hist1d[1]->Fill( (*genSQMass)[0] ); else hist1d[1]->Fill( 0 );

    //float rjrVDiff = (*rjrVDiff)[cs];
    //float rjrVSum = (*rjrVSum)[cs];

	float rjrVDiff = (*rjrMVDiff)[cs]; //((*rjrMVa)[cs]-(*rjrMVb)[cs])/((*rjrMVa)[cs]+(*rjrMVb)[cs]); 
	float rjrVSum = (*rjrMVSum)[cs]; //std::sqrt((sq2((*rjrMVa)[cs])+sq2((*rjrMVb)[cs]))/2); 

	//if( rjrVDiff != rjrVDiffC || rjrVSum != rjrVSumC ) std::cout << " V vars miss match !!!!!!" << std::endl;

	float nmass = (*rjrASMass)[cs];
    //float nmass = (*rjrX2QSum)[cs];
	float NormMBetaEql = (*selPhoMBetaEql)[0]/nmass;
    float NormMBetaPmt = (*selPhoMBetaPmt)[0]/nmass;
	float rjrNVSum = 2*rjrVSum/nmass;


	//if( rjrNVSum < 0.2 ) continue;
	//if( std::abs( rjrVDiff ) > 0.6 ) continue;
	//if( nmass < 1250 ) continue;
    //if( (*rjrNJetsJa)[cs] < 3 ) continue;
	//if( (*rjrNJetsJb)[cs] < 4 ) continue;

	//if( (*rjrSMass)[cs] > 1500 ){
	if( DEBUG ) std::cout << " -- Filling Histograms set 1" << std::endl;
    if( true ){

		hist1d[250]->Fill((*rjrSCosA)[cs], fillwt ); //("SCosA", mkht(ht,"SCosA"), 70, -3.5, 3.5);
	    hist1d[251]->Fill((*rjrSMass)[cs], fillwt ); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
	    hist1d[252]->Fill((*rjrX2aMass)[cs], fillwt ); //("X2aMass", mkht(ht,"X2aMass"), 500, 0, 5000);
	    hist1d[253]->Fill((*rjrX2bMass)[cs], fillwt ); //("X2bMass", mkht(ht,"X2bMass"), 500, 0, 5000);
	
	    hist1d[254]->Fill((*rjrDphiMETV)[cs], fillwt ); 
	    hist1d[255]->Fill((*rjrNJetsJa)[cs], fillwt ); 
	    hist1d[256]->Fill((*rjrNJetsJb)[cs], fillwt ); 
	    hist1d[257]->Fill((*rjrSdphiBV)[cs], fillwt ); 
	    hist1d[258]->Fill((*rjrSdphiDA)[cs], fillwt ); 
	    hist1d[259]->Fill((*rjrX2aPtS)[cs], fillwt ); 
	    hist1d[260]->Fill((*rjrX2bPtS)[cs], fillwt ); 

    if( DEBUG ) std::cout << " -- Filling Histograms set 11" << std::endl;

	    hist1d[261]->Fill((*rjrX2Ave)[cs], fillwt ); 
	    hist1d[262]->Fill((*rjrX2Diff)[cs], fillwt ); 
	    hist1d[263]->Fill((*rjrX2GMean)[cs], fillwt ); 
	    hist1d[264]->Fill((*rjrX2NAve)[cs], fillwt ); 
	    hist1d[265]->Fill((*rjrX2NGMean)[cs], fillwt ); 
	    hist1d[266]->Fill((*rjrX2NQSum)[cs], fillwt ); 
	    hist1d[267]->Fill((*rjrX2QSum)[cs], fillwt ); 

    if( DEBUG ) std::cout << " -- Filling Histograms set 12" << std::endl;

        hist1d[268]->Fill((*selPhoMBetaEql)[0], fillwt );
        hist1d[269]->Fill((*selPhoMBetaPmt)[0], fillwt );
        hist1d[270]->Fill((*rjrMVa)[cs], fillwt );
        hist1d[271]->Fill((*rjrMVb)[cs], fillwt );
        hist1d[272]->Fill(rjrVDiff, fillwt );
        hist1d[273]->Fill(rjrNVSum, fillwt );
        hist1d[274]->Fill(NormMBetaEql, fillwt );
        hist1d[275]->Fill(NormMBetaPmt, fillwt );

    if( DEBUG ) std::cout << " -- Filling Histograms set 13" << std::endl;

		//hist2d[0]->Fill( (*rjrX2Diff)[cs], (*rjrX2NQSum)[cs], fillwt );
        //hist2d[2]->Fill( (*rjrSMass)[cs], (*rjrX2NQSum)[cs], fillwt );

	//}//<<>>if( rjrX2NGMean->at(cs) < cutvc )

    //if( (*rjrAX2GMean)[cs] > 300 && (*rjrAX2NGMean)[cs] > 0.1 ){
    //if( (*rjrASMass)[cs] > 1200 ){
    //if( (*rjrASMass)[cs] > 1500 && rjrVSum > 200 ){
    //if( true ){

    if( DEBUG ) std::cout << " -- Filling Histograms set 2" << std::endl;

	    //hist1d[350]->Fill((*rjrSCosA)[cs], fillwt ); //("SCosA", mkht(ht,"SCosA"), 70, -3.5, 3.5);
	    hist1d[351]->Fill((*rjrASMass)[cs], fillwt ); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
	    hist1d[352]->Fill((*rjrAX2aMass)[cs], fillwt ); //("X2aMass", mkht(ht,"X2aMass"), 500, 0, 5000);
	    hist1d[353]->Fill((*rjrAX2bMass)[cs], fillwt ); //("X2bMass", mkht(ht,"X2bMass"), 500, 0, 5000);
	
	    //hist1d[354]->Fill((*rjrDphiMETV)[cs], fillwt );
	    //hist1d[355]->Fill((*rjrNJetsJa)[cs], fillwt );
	    //hist1d[356]->Fill((*rjrNJetsJb)[cs], fillwt );
	    //hist1d[357]->Fill((*rjrSdphiBV)[cs], fillwt );
	    //hist1d[358]->Fill((*rjrSdphiDA)[cs], fillwt );
	    //hist1d[359]->Fill((*rjrX2aPtS)[cs], fillwt );
	    //hist1d[360]->Fill((*rjrX2bPtS)[cs], fillwt );	

	    hist1d[361]->Fill((*rjrAX2Ave)[cs], fillwt );
	    hist1d[362]->Fill((*rjrAX2Diff)[cs], fillwt );
	    hist1d[363]->Fill((*rjrAX2GMean)[cs], fillwt );
	    hist1d[364]->Fill((*rjrAX2NAve)[cs], fillwt );
	    hist1d[365]->Fill((*rjrAX2NGMean)[cs], fillwt );
	    hist1d[366]->Fill((*rjrAX2NQSum)[cs], fillwt );
	    hist1d[367]->Fill((*rjrAX2QSum)[cs], fillwt );
	
        //hist1d[368]->Fill((*selPhoMBetaEql)[0], fillwt );
        //hist1d[369]->Fill((*selPhoMBetaPmt)[0], fillwt );

        //hist1d[370]->Fill((*rjrMVa)[cs], fillwt );
        //hist1d[371]->Fill((*rjrMVb)[cs], fillwt );
        //hist1d[372]->Fill(rjrVDiff, fillwt );
        //hist1d[373]->Fill(rjrVSum, fillwt );
        //hist1d[374]->Fill(NormMBetaEql, fillwt );
        //hist1d[375]->Fill(NormMBetaPmt, fillwt );

        //hist2d[1]->Fill( (*rjrAX2Diff)[cs], (*rjrAX2GMean)[cs], fillwt );

    if( DEBUG ) std::cout << " -- Filling Histograms set 3" << std::endl;

		hist2d[0]->Fill( (*rjrASMass)[cs], (*rjrAX2NQSum)[cs], fillwt );
    	hist2d[1]->Fill( (*rjrASMass)[cs], (*rjrAX2QSum)[cs], fillwt ); // = new TH2D("ASMvX2ANVSum", mkht(ht,"ASMvX2ANVSum;rjrASMass;rjrX2ANVSum").c_str(), 150, 0, 15000, 50, 0, 1.0 );
    	hist2d[2]->Fill( (*rjrASMass)[cs], rjrNVSum, fillwt ); // = new TH2D("ASMvAx2mass", mkht(ht,"ASMvAX2Mass;rjrASMass;rjrAX2Mass").c_str(), 150, 0, 15000, 90, 0, 3600 );
    	hist2d[3]->Fill( (*rjrASMass)[cs], rjrVDiff, fillwt ); // = new TH2D("ASMvVMass", mkht(ht,"ASMvVMass;rjrASMass;rjrVMass").c_str(), 150, 0, 15000, 90, 0, 3600 );

    	hist2d[4]->Fill( (*rjrAX2NQSum)[cs], rjrNVSum, fillwt ); // = new TH2D("ASMvVDiff", mkht(ht,"ASMvVDiff;rjrASMass;rjrVDiff").c_str(), 150, 0, 15000, 100, -1.0, 1.0 );
    	hist2d[5]->Fill( (*rjrAX2NQSum)[cs], rjrVDiff, fillwt ); // = new TH2D("ARatiovAxmkht(ht,"X2ANVSumvAX2Mass;rjrX2ANVSum;rjrAX2Mass").c_str(), 50, 0, 1.0, 90, 0, 3600 );
    	hist2d[6]->Fill( (*rjrAX2NQSum)[cs], (*rjrAX2QSum)[cs], fillwt ); // = new TH2D("ARatiovVMass", mkht(ht,"X2ANVSumvVSVis;rjrX2ANVSum;rjrVSum").c_str(), 50, 0, 1.0, 90, 0, 3600 );

    	hist2d[7]->Fill( (*rjrAX2QSum)[cs], rjrNVSum, fillwt ); // = new TH2D("ARatiovVDiff", mkht(ht,"X2ANVSumvVDiff;rjrX2ANVSum;rjrVDiff").c_str(), 50, 0, 1.0, 100, -1.0, 1.0 );
    	hist2d[8]->Fill( (*rjrAX2QSum)[cs], rjrVDiff, fillwt ); // = new TH2D("Ax2massvVMass", mkht(ht,"AX2MassvVSVis;rjrAX2Mass;rjrVSum").c_str(), 90, 0, 3600, 90, 0, 3600 );

    	hist2d[9]->Fill( rjrNVSum, rjrVDiff, fillwt ); // = new TH2D("Ax2massvVDiff", mkht(ht,"AX2MassvVDiff;rjrAX2Mass;rjrVDiff").c_str(), 90, 0, 3600, 100, -1.0, 1.0 )
		hist2d[10]->Fill( (*rjrNJetsJa)[cs], (*rjrNJetsJb)[cs], fillwt );


	}//<<>>if( rjrAX2NGMean->at(cs) < cutvc && rjrAX2GMean->at(cs) < cutvd )
	}//<<>>for( int i = 0; i < 1; i++ ) -- continue loop

}//<<>>void HistMaker::eventLoop(` Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

    hist1d[2]->SetBinContent(1,cutflowInfo["cf_nTotEvts"]);
    hist1d[2]->SetBinError(1,std::sqrt(cutflowInfo["cf_nTotEvts"]));
    hist1d[2]->SetBinContent(2,cutflowInfo["cf_nFltrdEvts"]);
    hist1d[2]->SetBinError(2,std::sqrt(cutflowInfo["cf_nFltrdEvts"]));
    hist1d[2]->SetBinContent(3,cutflowInfo["cf_met150"]);
    hist1d[2]->SetBinError(3,std::sqrt(cutflowInfo["cf_met150"]));
    hist1d[2]->SetBinContent(4,cutflowInfo["cf_m_gt2jets"]);
    hist1d[2]->SetBinError(4,std::sqrt(cutflowInfo["cf_m_gt2jets"]));
    hist1d[2]->SetBinContent(5,cutflowInfo["cf_mj_gt1phos"]);
    hist1d[2]->SetBinError(5,std::sqrt(cutflowInfo["cf_mj_gt1phos"]));
    hist1d[2]->SetBinContent(6,cutflowInfo["cf_mjp_leadPhoPt30"]);
    hist1d[2]->SetBinError(6,std::sqrt(cutflowInfo["cf_mjp_leadPhoPt30"]));

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------


	//------- jets 0 - 99
    //UInt_t          nJets;
    //UInt_t          nSelJets;
    //std::vector<float>   *selJetEnergy;
    //std::vector<float>   *selJetEta;
    //std::vector<float>   *selJetMass;
    //std::vector<float>   *selJetPhi;
    //std::vector<float>   *selJetPt;
    //std::vector<int>     *selJetQuality;
    //std::vector<float>   *selJetTime;

    std::cout << " title test : " << addstr(ht,"jetPt") << std::endl;

	//----- photons 100 - 249
    //UInt_t          leadSelPho;
    //UInt_t          nPhotons;
    //UInt_t          nSelPhotons;
    //std::vector<float>   *selPhoClstrRn;
    //std::vector<float>   *selPhoEnergy;
    //std::vector<float>   *selPhoEta;
    //std::vector<float>   *selPhoGeoEgnVal;
    //std::vector<float>   *selPhoGeoSMaj;
    //std::vector<float>   *selPhoGeoSMin;
    //std::vector<unsigned int> *selPhoNrh;
    //std::vector<float>   *selPhoPhi;
    //std::vector<float>   *selPhoPt;
    //std::vector<int>     *selPhoQuality;
    //std::vector<float>   *selPhoR9;
    //std::vector<float>   *selPhoSMaj;
    //std::vector<float>   *selPhoSMin;
    //std::vector<float>   *selPhoSieie;
    //std::vector<float>   *selPhoTime;
    //UInt_t          subLeadSelPho;

	//------  genparticles 200 - 249

	//------  rjr 250 - 299
//rjrDphiMETV
//rjrMET
//rjrNJetsJa
//rjrNJetsJb
//rjrNRjrPhotons
//rjrSdphiBV
//rjrSdphiDA

//rjrX2Ave
//rjrX2Diff
//rjrX2GMean
//rjrX2NAve
//rjrX2NGMean
//rjrX2NQSum
//rjrX2QSum
//
//rjrX2aMass
//rjrX2bMass
//rjrX2aPtS
//rjrX2bPtS
//
//genSigType
//

    hist1d[0] = new TH1D("genSigSGluinoMass", addstr(ht,"genSigSGluinoMass").c_str(), 1200, 0, 6000);
    hist1d[1] = new TH1D("genSigSQuarkMass", addstr(ht,"genSigSQuarkMass").c_str(), 1200, 0, 6000);
    hist1d[2] = new TH1D("cutflow", addstr(ht,"cutFlow").c_str(), 25, 0.5, 25.5);
    hist1d[2]->Sumw2();

    hist1d[250] = new TH1D("SCosA", addstr(ht,"SCosA").c_str(), 70, -3.5, 3.5);
    hist1d[251] = new TH1D("SMass", addstr(ht,"SMass").c_str(), 150, 0, 15000);
    hist1d[252] = new TH1D("X2aMass", addstr(ht,"X2aMass").c_str(), 90, 0, 3600);
    hist1d[253] = new TH1D("X2bMass", addstr(ht,"X2bMass").c_str(), 90, 0, 3600);

    hist1d[254] = new TH1D("rjrDphiMETV", addstr(ht,"rjrDphiMETV").c_str(), 70, -3.5, 3.5);
    hist1d[255] = new TH1D("rjrNJetsJa", addstr(ht,"rjrNJetsJa").c_str(), 20, 0, 20);
    hist1d[256] = new TH1D("rjrNJetsJb", addstr(ht,"rjrNJetsJb").c_str(), 20, 0, 20);
    hist1d[257] = new TH1D("rjrSdphiBV", addstr(ht,"rjrSdphiBV").c_str(), 35, 0, 3.5);
    hist1d[258] = new TH1D("rjrSdphiDA", addstr(ht,"rjrSdphiDA").c_str(), 35, 0, 3.5);
    hist1d[259] = new TH1D("rjrX2aPtS", addstr(ht,"rjrX2aPtS").c_str(), 100, 0, 6000);
    hist1d[260] = new TH1D("rjrX2bPtS", addstr(ht,"rjrX2bPtS").c_str(), 100, 0, 6000);

    hist1d[261] = new TH1D("rjrX2Ave", addstr(ht,"rjrX2Ave").c_str(), 90, 0, 3600);
    hist1d[262] = new TH1D("rjrX2Diff", addstr(ht,"rjrX2Diff").c_str(), 100, -1, 1);
    hist1d[263] = new TH1D("rjrX2GMean", addstr(ht,"rjrX2GMean").c_str(), 90, 0, 3600);
    hist1d[264] = new TH1D("rjrX2NAve", addstr(ht,"rjrX2NAve").c_str(), 50, 0, 1.0);
    hist1d[265] = new TH1D("rjrX2NGMean", addstr(ht,"rjrX2NGMean").c_str(), 50, 0, 1.0);
    hist1d[266] = new TH1D("rjrX2NQSum", addstr(ht,"rjrX2NQSum").c_str(), 50, 0, 1.0);
    hist1d[267] = new TH1D("rjrX2QSum", addstr(ht,"rjrX2QSum").c_str(), 90, 0, 3600);

    hist1d[268] = new TH1D("selPhoMBetaEql", addstr(ht,"selPhoBetaMass").c_str(), 60, 0, 1200);//----------- BetaMass
    hist1d[269] = new TH1D("selPhoMBetaPmt", addstr(ht,"selPhoMBetaPmt").c_str(), 60, 0, 1200);

    hist1d[270] = new TH1D("rjrMVa", addstr(ht,"rjrMVa").c_str(), 90, 0, 3600);
    hist1d[271] = new TH1D("rjrMVb", addstr(ht,"rjrMVb").c_str(), 90, 0, 3600);
    hist1d[272] = new TH1D("rjrVDiff", addstr(ht,"NVisMassDiff").c_str(), 25, -1, 1);//---------------VisMassDiff !!!!!!!!!!
    hist1d[273] = new TH1D("rjrVSum", addstr(ht,"NVisMass").c_str(), 25, 0, 1.0);//--------------VisMassSum !!!!!!!!!!

    hist1d[274] = new TH1D("rjrPhoNMBetaEql", addstr(ht,"rjrPhoNBetaEql").c_str(), 100, 0, 2);//----------- BetaMass
    hist1d[275] = new TH1D("rjrPhoNMBetaPmt", addstr(ht,"rjrPhoNMBetaPmt").c_str(), 100, 0, 2);

    // Alternative approach
    //hist1d[350] = new TH1D("ASCosAc", addstr(ht,"ASCosA").c_str(), 70, -3.5, 3.5);
    hist1d[351] = new TH1D("ASMass", addstr(ht,"ASMass").c_str(), 48, 0, 12000); //----------------------ASMass !!!!!!!!!!!!!!!!!!1
    //hist1d[351] = new TH1D("ASMass", addstr(ht,"SMass").c_str(), 150, 0, 15000);
    hist1d[352] = new TH1D("AX2aMass", addstr(ht,"AX2aMass").c_str(), 90, 0, 3600);
    hist1d[353] = new TH1D("AX2bMass", addstr(ht,"AX2bMass").c_str(), 90, 0, 3600);

    //hist1d[354] = new TH1D("rjrDphiMETVcA", addstr(ht,"rjrDphiMETV A").c_str(), 70, -3.5, 3.5);
    //hist1d[355] = new TH1D("rjrNJetsJacA", addstr(ht,"rjrNJetsJa A").c_str(), 20, 0, 20);
    //hist1d[356] = new TH1D("rjrNJetsJbcA", addstr(ht,"rjrNJetsJb A").c_str(), 20, 0, 20);
    //hist1d[357] = new TH1D("rjrSdphiBVcA", addstr(ht,"rjrSdphiBV A").c_str(), 35, 0, 3.5);
    //hist1d[358] = new TH1D("rjrSdphiDAcA", addstr(ht,"rjrSdphiDA A").c_str(), 35, 0, 3.5);
    //hist1d[359] = new TH1D("rjrX2aPtScA", addstr(ht,"rjrX2aPtS A").c_str(), 100, 0, 6000);
    //hist1d[360] = new TH1D("rjrX2bPtScA", addstr(ht,"rjrX2bPtS A").c_str(), 100, 0, 6000);

    hist1d[361] = new TH1D("rjrAX2Ave", addstr(ht,"rjrAX2Ave").c_str(), 90, 0, 3600);
    hist1d[362] = new TH1D("rjrAX2Diff", addstr(ht,"X2MassDiff").c_str(), 100, -1, 1);//----------------X2MassDiff
    //hist1d[362] = new TH1D("rjrAX2Diff", addstr(ht,"X2MassDiff").c_str(), 100, -1, 1);
    hist1d[363] = new TH1D("rjrAX2GMean", addstr(ht,"rjrAX2GMean").c_str(), 90, 0, 3600);
    hist1d[364] = new TH1D("rjrAX2NAve", addstr(ht,"rjrAX2NAve").c_str(), 50, 0, 1.0);
    hist1d[365] = new TH1D("rjrAX2NGMean", addstr(ht,"rjrAX2NGMean").c_str(), 50, 0, 1.0);
    hist1d[366] = new TH1D("rjrAX2NQSum", addstr(ht,"X2MassRatio").c_str(), 25, 0, 1.0);//----------X2MassRatio !!!!!!!!!!!!!!!!!
    //hist1d[366] = new TH1D("rjrAX2NQSum", addstr(ht,"X2MassRatio").c_str(), 50, 0, 1.0)
    //hist1d[367] = new TH1D("rjrAX2QSum", addstr(ht,"X2Mass").c_str(), 90, 0, 3600);//-----------X2Mass  !!!!!!!!!!!!!!!!!!!1
    hist1d[367] = new TH1D("rjrAX2QSum", addstr(ht,"X2Mass").c_str(), 26, 0, 2600);

    //hist1d[368] = new TH1D("selPhoMBetaEqlcA", addstr(ht,"selPhoMBetaEql A").c_str(), 60, 0, 1200);
    //hist1d[369] = new TH1D("selPhoMBetaPmtcA", addstr(ht,"selPhoMBetaPmt A").c_str(), 60, 0, 1200);

    //hist1d[370] = new TH1D("rjrMVacA", addstr(ht,"rjrMVa A").c_str(), 90, 0, 3600);
    //hist1d[371] = new TH1D("rjrMVbcA", addstr(ht,"rjrMVb A").c_str(), 90, 0, 3600);
    //hist1d[372] = new TH1D("rjrVDiffcA", addstr(ht,"rjrVDiff A").c_str(), 100, -1, 1);
    //hist1d[373] = new TH1D("rjrVSumcA", addstr(ht,"rjrVSum Ac").c_str(), 90, 0, 3600);
    //hist1d[374] = new TH1D("selPhoNMBetaEqlcA", addstr(ht,"selPhoNMBetaEql A").c_str(), 100, 0, 2);
    //hist1d[375] = new TH1D("selPhoNMBetaPmtcA", addstr(ht,"selPhoNMBetaPmt A").c_str(), 100, 0, 2);

    //------ ecal rechits 300 - 349

    //------  electrons 350 - 400

	//-------- event vars 400 - 450
    hist1d[400] = new TH1D("evtMetE", addstr(ht,"evtMetE").c_str(), 150, 0, 1500 );
    hist1d[401] = new TH1D("evtMetPt", addstr(ht,"evtMetPt ").c_str(), 150, 0, 1500 );	
    hist1d[402] = new TH1D("evtRjrMetPt", addstr(ht,"evtRjrMetPt ").c_str(), 100, 0, 1500 );
    hist1d[403] = new TH1D("evtRawMetPt", addstr(ht,"evtRawMetPt ").c_str(), 100, 0, 1500 );

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    hist2d[0] = new TH2D("ASMvAX2NQSum", addstr(ht,"ASMvAX2NQSum;ASMass [GeV];X2ANQSum").c_str(), 48, 0, 12000, 25, 0, 1.0 );
    hist2d[1] = new TH2D("ASMvAX2QSum", addstr(ht,"ASMvAX2QSum;ASMass [GeV];AX2QSum [GeV]").c_str(), 48, 0, 12000, 26, 0, 2600 );
    hist2d[2] = new TH2D("ASMvVMass", addstr(ht,"ASMvNVSum;ASMass [GeV];NVSum").c_str(), 48, 0, 12000, 25, 0, 1.0 );
    hist2d[3] = new TH2D("ASMvVDiff", addstr(ht,"ASMvVDiff;ASMass [GeV];VDiff").c_str(), 48, 0, 12000, 25, -1.0, 1.0 );


    hist2d[4] = new TH2D("AX2NQSumvVMass", addstr(ht,"AX2NQSumvNVSum;AX2NQSum;NVSum").c_str(), 25, 0, 1.0, 25, 0, 1.0 );
    hist2d[5] = new TH2D("AX2NQSumvVDiff", addstr(ht,"AX2NQSumvVDiff;AX2NQSum;VDiff").c_str(), 25, 0, 1.0, 25, -1.0, 1.0 );
    hist2d[6] = new TH2D("AX2NQSumvAX2QSum", addstr(ht,"AX2NQSumvAX2QSum;AX2NQSum;X2AQSum [GeV]").c_str(), 25, 0, 1.0, 26, 0, 2600 );

    hist2d[7] = new TH2D("AX2QSumvVMass", addstr(ht,"AX2QSumvNVSum;AX2QSum [GeV];NVSum").c_str(), 26, 0, 2600, 25, 0, 1.0 );
    hist2d[8] = new TH2D("AX2QSumvVDiff", addstr(ht,"AX2QSumvVDiff;AX2QSum [GeV];VDiff").c_str(), 26, 0, 2600, 25, -1.0, 1.0 );

    hist2d[9] = new TH2D("VMassvVDiff", addstr(ht,"NVSumvVDiff;NVSum;VDiff").c_str(), 25, 0, 1.0, 25, -1.0, 1.0 );
    hist2d[10] = new TH2D("NJetsJavNJetsJb", addstr(ht,"NJetsJavNJetsJb;NJetsJa;NJetsJb").c_str(), 20, 0, 20,20, 0, 20 );

	//------- jets ( time ) 0-49 ------------------------------

	//---jet id stuff 50 - 99 ---------------------------------------------------

	//--- Photons 200 - 349 -------------------------------------------
	

	//60 - 63

	//--- rechit collections 350 - 399 -------------------------------------------------

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------
	nMaps = 0;
    std::string label("baseHists");
	if( doEBEEmaps ){ for(int it=0; it<nEBEEMaps; it++){
		fMap[it] = false;
		std::string label(";iEta;iPhi");
        std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
		std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
	}}//<<>>for(int it=0; it<nEBEEMaps; it++)

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "rjr_skim_files/";
                ////auto infilename = "KUCMS_RJR_Skim_List.txt";
                auto infilename0 = "KUCMS_RJR_GMSB100_ootmet_Skim_List.txt";
                auto infilename1 = "KUCMS_RJR_GMSB150_ootmet_Skim_List.txt";
                auto infilename2 = "KUCMS_RJR_GMSB200_ootmet_Skim_List.txt";
                auto infilename3 = "KUCMS_RJR_GMSB250_ootmet_Skim_List.txt";
                auto infilename4 = "KUCMS_RJR_GMSB300_ootmet_Skim_List.txt";

                auto infilename10 = "KUCMS_RJR_GMSB_ct0p1_ootmet_Skim_List.txt";
                auto infilename11 = "KUCMS_RJR_GMSB_ct10_ootmet_Skim_List.txt";
                auto infilename12 = "KUCMS_RJR_GMSB_ct200_ootmet_Skim_List.txt";
                auto infilename13 = "KUCMS_RJR_GMSB_ct1000_ootmet_Skim_List.txt";
                auto infilename14 = "KUCMS_RJR_GMSB_ct10000_ootmet_Skim_List.txt";
                ////auto infilename5 = "KUCMS_RJR_GMSB350_ootmet_Skim_List.txt";
                ////auto infilename6 = "KUCMS_RJR_GMSB400_ootmet_Skim_List.txt";
				auto infilenameJ = "KUCMS_RJR_GIGI_ootmet_Skim_List.txt";

                //auto infilenameG = "KUCMS_RJR_GJETS_rawmet_Skim_List.txt";
                //auto infilenameQ = "KUCMS_RJR_QCD_rawmet_Skim_List.txt";
                //auto infilenameD = "KUCMS_RJR_DEG_rawmet_Skim_List.txt";
                //auto infilename1 = "KUCMS_RJR_GMSB100_ootmet_Skim_List.txt";
                //auto infilename3 = "KUCMS_RJR_GMSB350_ootmet_Skim_List.txt";
                auto infilenameG = "KUCMS_RJR_GJETS_ootmet_Skim_List.txt";
                auto infilenameQ = "KUCMS_RJR_QCD_ootmet_Skim_List.txt";
                //auto infilenameD = "KUCMS_RJR_DEG_ootmet_Skim_List.txt";
                auto infilenameD = "KUCMS_RJR_MET_ootmet_Skim_List.txt";
                //auto infilenameD = "KUCMS_RJR_JetHT_ootmet_Skim_List.txt";
                auto infilenameWJLNu = "KUCMS_RJR_WJLNu_ootmet_Skim_List.txt";
                //std::list<float> cuts = { 10000000 }; 

                int modtype = 0;
			 	//std::string sigtype = "SqkSqk"; modtype = 11; //11
                //std::string sigtype = "GluGlu"; modtype = 22; //22
                //std::string sigtype = "SqkGlu"; modtype = 12; //12
                //std::string sigtype = "XinoXino"; modtype = 33; //33
                //std::string sigtype = "SleptSlept"; modtype = 44; //44
                //std::string sigtype = "BkGrd"; // modtype 0 - no modtype applied
				std::string sigtype = "gogoG";


				std::string version = "_v23_";

                std::string outfilenamed = "KUCMS_MET_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilenamed = "KUCMS_JetHT_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilenamed = "KUCMS_DEG_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilenamegj = "KUCMS_GJets_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilenameqcd = "KUCMS_QCD_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilenamewjlnu = "KUCMS_WJLNu_"+sigtype+version;
                //std::string outfilename100 = "KUCMS_GMSB_L100_T30_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename300 = "KUCMS_GMSB_L350_T30_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename0 = "KUCMS_GMSB_L100_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename1 = "KUCMS_GMSB_L150_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename2 = "KUCMS_GMSB_L200_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename3 = "KUCMS_GMSB_L250_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename4 = "KUCMS_GMSB_L300_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7

                std::string outfilename10 = "KUCMS_GMSB_ct0p1cm_"+sigtype+version;
                std::string outfilename11 = "KUCMS_GMSB_ct10cm_"+sigtype+version;
                std::string outfilename12 = "KUCMS_GMSB_ct200cm_"+sigtype+version;
                std::string outfilename13 = "KUCMS_GMSB_ct1000cm_"+sigtype+version;
                std::string outfilename14 = "KUCMS_GMSB_ct10000cm_"+sigtype+version;
				std::string outfilenameJ = "KUCMS_mGl-2000_mN2-1900_mN1-1_"+sigtype+version;

                //std::string outfilename5 = "KUCMS_GMSB_L350_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename6 = "KUCMS_GMSB_L400_"+sigtype+version; //iso0_Skim_BaseHists.root"; //7

                //std::string ofnending = "wt2_RjrSkim_multiHists_v2.root";
                //std::string ofnending = "NGM0p2_wt2_RjrSkim_multiHists.root";
       
		        //std::string ofnending = "wt2_RjrSkim_v21_phojet_c1b_multiHists.root"; float jrjtype = 1; // 1 = phojet, 0 = phomet
                //std::string ofnending = "wt2_RjrSkim_v21_phomet_c1a_multiHists.root"; float jrjtype = 0; // 1 = phojet, 0 = phomet
                //std::string ofnending = "wt2_RjrSkim_v22_rawmet_phojet_ztc2_multiHists.root"; float jrjtype = 1; // 1 = phojet, 0 = phomet
                //std::string ofnending = "wt2_RjrSkim_v22_rawmet_phomet_ztc2_multiHists.root"; float jrjtype = 0; // 1 = phojet, 0 = phomet

                ////std::string ofnending = "wt2_RjrSkim_v24_ootmet_phomet_multiHists.root"; float jrjtype = 0; // 1 = phojet, 0 = phomet
                std::string ofnending = "wt2_RjrSkim_v24_ootmet_phojet_multiHists.root"; float jrjtype = 1; // 1 = phojet, 0 = phomet

                ////std::string ofnending = "wt2_RjrSkim_v24_ootmet_phojet_multiHists.root"; float jrjtype = 0; // 0 = phojet, no phomet

                std::string htitled = "DEG_"+sigtype+version;
                std::string htitlegj = "GJets_"+sigtype+version;
                std::string htitleqcd = "QCD_"+sigtype+version;
                std::string htitlewjlnu = "WJetsLNu_"+sigtype+version;
                //std::string htitle100 = "GMSB_L100_v19_T30_";
                //std::string htitle300 = "GMSB_L350_v19_T30_";
                std::string htitle0 = "GMSB_L100_"+sigtype+version;
                std::string htitle1 = "GMSB_L150_"+sigtype+version;
                std::string htitle2 = "GMSB_L200_"+sigtype+version;
                std::string htitle3 = "GMSB_L250_"+sigtype+version;
                std::string htitle4 = "GMSB_L300_"+sigtype+version;

				std::string htitle10 = "GMSB_ct0p1cm_"+sigtype+version;
                std::string htitle11 = "GMSB_ct10cm_"+sigtype+version;
                std::string htitle12 = "GMSB_ct200cm_"+sigtype+version;
                std::string htitle13 = "GMSB_ct1000cm_"+sigtype+version;
                std::string htitle14 = "GMSB_ct10000cm_"+sigtype+version;
                //std::string htitle5 = "GMSB_L350_"+sigtype+version;
                //std::string htitle6 = "GMSB_L400_"+sigtype+version;
				std::string htitleJ = "ct10cm_mGl-2000_mN2-1900_mN1-1_"+version;

                HistMaker base;

				float rjrcut = 0; // var cut - unused
                float rjrcut2 = 150;  // cmet > 150
                //float rjrcut2 = 0;  // cmet > 150 BG only

				//float nphos = 2;
                //std::string isoline = "genSigPerfect_nSigPho2_";
                //std::string isolinebg = "nSigPho2_"; 
                float nphos = 1;
                //std::string isoline = "genSigPerfect_nSigPho1_";
                std::string isoline = "nSigPho1_ab1plus_";
                std::string isolinebg = "nSigPho1_"; 
                //float nphos = 0;
                //std::string isoline = "genSigPerfect_nSigPho0_";
                //std::string isolinebg = "nSigPho0_";			
	
                //float nphos = 1;
                //std::string isoline = "gsp_nsp1_nvsvdnjabcut_";
                //std::string isolinebg = "nsp1_nvsvdnjabcut_";

                std::string outfilenamesd = outfilenamed + isolinebg + ofnending;
                std::string outfilenamesgj = outfilenamegj + isolinebg + ofnending;
                std::string outfilenamesqcd = outfilenameqcd + isolinebg + ofnending;
                std::string outfilenameswjlnu = outfilenamewjlnu + isolinebg + ofnending;
                //std::string outfilenames100 = outfilename100 + isoline + ofnending;
                //std::string outfilenames300 = outfilename300 + isoline + ofnending;
				std::string outfilenames0 = outfilename0 + isoline + ofnending;
                std::string outfilenames1 = outfilename1 + isoline + ofnending;
                std::string outfilenames2 = outfilename2 + isoline + ofnending;
                std::string outfilenames3 = outfilename3 + isoline + ofnending;
                std::string outfilenames4 = outfilename4 + isoline + ofnending;

                std::string outfilenames10 = outfilename10 + isoline + ofnending;
                std::string outfilenames11 = outfilename11 + isoline + ofnending;
                std::string outfilenames12 = outfilename12 + isoline + ofnending;
                std::string outfilenames13 = outfilename13 + isoline + ofnending;
                std::string outfilenames14 = outfilename14 + isoline + ofnending;
                //std::string outfilenames5 = outfilename5 + isoline + ofnending;
                //std::string outfilenames6 = outfilename6 + isoline + ofnending;

				std::string outfilenamesJ = outfilenameJ + isoline + ofnending;

                std::string htitlesd =  htitled + isolinebg;
                std::string htitlesgj =  htitlegj + isolinebg;
                std::string htitlesqcd =  htitleqcd + isolinebg;
                std::string htitleswjlnu =  htitlewjlnu + isolinebg;
                //std::string htitles100 =  htitle100 + isoline;
                //std::string htitles300 =  htitle300 + isoline;
                std::string htitle0s =  htitle0 + isoline;
                std::string htitle1s =  htitle1 + isoline;
                std::string htitle2s =  htitle2 + isoline;
                std::string htitle3s =  htitle3 + isoline;
                std::string htitle4s =  htitle4 + isoline;

                std::string htitle10s =  htitle10 + isoline;
                std::string htitle11s =  htitle11 + isoline;
                std::string htitle12s =  htitle12 + isoline;
                std::string htitle13s =  htitle13 + isoline;
                std::string htitle14s =  htitle14 + isoline;
                //std::string htitle5s =  htitle5 + isoline;
                //std::string htitle6s =  htitle6 + isoline;

				std::string htitleJs =  htitleJ + isoline;

				if( sigtype == "gogoG" ){

					base.histMaker( listdir, infilenameJ, outfilenamesJ, htitleJs, jrjtype, modtype, nphos, rjrcut, rjrcut2 );

				} else if( sigtype != "BkGrd" ){
            		base.histMaker( listdir, infilename0, outfilenames0, htitle0s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	base.histMaker( listdir, infilename1, outfilenames1, htitle1s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	base.histMaker( listdir, infilename2, outfilenames2, htitle2s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	base.histMaker( listdir, infilename3, outfilenames3, htitle3s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	base.histMaker( listdir, infilename4, outfilenames4, htitle4s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	////base.histMaker( listdir, infilename5, outfilenames5, htitle5s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	////base.histMaker( listdir, infilename6, outfilenames6, htitle6s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                    //base.histMaker( listdir, infilename10, outfilenames10, htitle10s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                    //base.histMaker( listdir, infilename11, outfilenames11, htitle11s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                    //base.histMaker( listdir, infilename12, outfilenames12, htitle12s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                    //base.histMaker( listdir, infilename13, outfilenames13, htitle13s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                    //base.histMaker( listdir, infilename14, outfilenames14, htitle14s, jrjtype, modtype, nphos, rjrcut, rjrcut2 );

				} else {
                    base.histMaker( listdir, infilenameWJLNu, outfilenameswjlnu, htitleswjlnu, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	////base.histMaker( listdir, infilename1, outfilenames100, htitles100, jrjtype, 30, nphos, rjrcut, rjrcut2 );
                	////base.histMaker( listdir, infilename3, outfilenames300, htitles300, jrjtype, 30, nphos, rjrcut, rjrcut2 );
					//base.histMaker( listdir, infilenameG, outfilenamesgj, htitlesgj, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	//base.histMaker( listdir, infilenameQ, outfilenamesqcd, htitlesqcd, jrjtype, modtype, nphos, rjrcut, rjrcut2 );
                	//base.histMaker( listdir, infilenameD, outfilenamesd, htitlesd, 0, 0, nphos, rjrcut, rjrcut2 );
				}                

				//base.histMaker( listdir, infilename, outfilenames, htitles, fillwt?, gensigtype?, #phos?, rjrcut, n/a );

    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

