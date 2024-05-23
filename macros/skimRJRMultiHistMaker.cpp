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

    TBranch        *b_nEvents;   //!
    TBranch        *b_nSelectedEvents;   //!
    TBranch        *b_sKey;   //!
    TBranch        *b_sCrossSection;   //!
    TBranch        *b_sGMSBGravMass;   //!
    TBranch        *b_sGMSBChi1Mass;   //!
    TBranch        *b_sMCWgt;   //!
    TBranch        *b_sMCType;   //!
    TBranch        *b_sumEvtGenWgt;

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

		if( not configInfo.count(configKey) ){
        	if(debug) std::cout << " - Filling configValues. " << std::endl;
			std::map< std::string, float > configValues;
        	configValues["nEvents"] = nEvents;
        	configValues["nSelectedEvents"] = nSelectedEvents;
        	configValues["sCrossSection"] = sCrossSection;
        	configValues["sGMSBGravMass"] = sGMSBGravMass;
        	configValues["sGMSBChi1Mass"] = sGMSBChi1Mass;
        	configValues["sMCWgt"] = sMCWgt;
        	configValues["sMCType"] = sMCType;
        	if(debug) std::cout << " - Filling configInfo. " << std::endl;
        	configInfo[configKey] = configValues;
		} else {
			auto & configValues = configInfo[configKey];
			configValues["nEvents"] += nEvents;
			configValues["nSelectedEvents"] += nSelectedEvents;
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
    float evtgwt = evtGenWgt;
    //float evtgwt = 1;
    float scale = 1;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nEvents"];
    float fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
	//if( cutselection == 1 ) fillwt = 1;
    //std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
    //if( configKey != "GMSBL"+lambda+"" ) return;

	for( int i = 0; i < 1; i++ ){
	if( rjrNRjrPhotons->size() == 0 ) continue;

	//bool genSigTypeCut( genSigType != 11 && genSigType != 12 && genSigType != 22 );
    bool genSigTypeCut( genSigType != cutva );
	bool nRjrPhoCut( rjrNRjrPhotons->at(cs) != cutvb );
    if( cutva == 50 ) genSigTypeCut = genSigType == 11 || genSigType == 12 || genSigType == 22;
	if( cutva == 30 ) genSigTypeCut = genSigType != 11 && genSigType != 12 && genSigType != 22;
	if( cutva == 0 ) genSigTypeCut = false;
	if( cutvb == 0 ) nRjrPhoCut = false;
	if( genSigTypeCut || nRjrPhoCut ) continue;

	//if( rjrSMass->at(0) < 2000 ) continue;
	//if( rjrX2NGMean->at(cs) < cutvc ) continue;
    //if( rjrX2GMean->at(cs) < cutvd ) continue;
	//std::cout << "RJR cut at " << cutvc << " with value " << rjrX2GMean->at(0) << std::endl;

    auto metCPt = std::sqrt(rad2(selCMetPx,selCMetPy));
    auto metPt = std::sqrt(rad2(selMetPx,selMetPy));

    hist1d[400]->Fill( selCMet, fillwt );
    hist1d[401]->Fill( metCPt, fillwt );
    hist1d[402]->Fill( rjrMET->at(cs), fillwt );
    hist1d[403]->Fill( metPt, fillwt );

	float rjrVDiff = ((*rjrMVa)[cs]-(*rjrMVb)[cs])/((*rjrMVa)[cs]+(*rjrMVb)[cs]); 
	float rjrVSVis = std::sqrt((sq2((*rjrMVa)[cs])+sq2((*rjrMVb)[cs]))/2); 

	//float nmass = (*rjrSMass)[cs];
    float nmass = (*rjrX2VSum)[cs];
	float NormMBetaEql = (*selPhoMBetaEql)[0]/nmass;
    float NormMBetaPmt = (*selPhoMBetaPmt)[0]/nmass;

	//if( (*rjrSMass)[cs] > 1500 ){
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
	
	    hist1d[261]->Fill((*rjrX2Ave)[cs], fillwt ); 
	    hist1d[262]->Fill((*rjrX2Diff)[cs], fillwt ); 
	    hist1d[263]->Fill((*rjrX2GMean)[cs], fillwt ); 
	    hist1d[264]->Fill((*rjrX2NAve)[cs], fillwt ); 
	    hist1d[265]->Fill((*rjrX2NGMean)[cs], fillwt ); 
	    hist1d[266]->Fill((*rjrX2NVSum)[cs], fillwt ); 
	    hist1d[267]->Fill((*rjrX2VSum)[cs], fillwt ); 

        hist1d[268]->Fill((*selPhoMBetaEql)[0], fillwt );
        hist1d[269]->Fill((*selPhoMBetaPmt)[0], fillwt );
        hist1d[270]->Fill((*rjrMVa)[cs], fillwt );
        hist1d[271]->Fill((*rjrMVb)[cs], fillwt );
        hist1d[272]->Fill(rjrVDiff, fillwt );
        hist1d[273]->Fill(rjrVSVis, fillwt );
        hist1d[274]->Fill(NormMBetaEql, fillwt );
        hist1d[275]->Fill(NormMBetaPmt, fillwt );

		hist2d[0]->Fill( (*rjrX2Diff)[cs], (*rjrX2NVSum)[cs], fillwt );
        hist2d[2]->Fill( (*rjrSMass)[cs], (*rjrX2NVSum)[cs], fillwt );

	}//<<>>if( rjrX2NGMean->at(cs) < cutvc )

    //if( (*rjrAX2GMean)[cs] > 300 && (*rjrAX2NGMean)[cs] > 0.1 ){
    //if( (*rjrASMass)[cs] > 1200 ){
    //if( (*rjrASMass)[cs] > 1500 && rjrVSVis > 200 ){
    if( true ){

	    hist1d[350]->Fill((*rjrSCosA)[cs], fillwt ); //("SCosA", mkht(ht,"SCosA"), 70, -3.5, 3.5);
	    hist1d[351]->Fill((*rjrASMass)[cs], fillwt ); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
	    hist1d[352]->Fill((*rjrAX2aMass)[cs], fillwt ); //("X2aMass", mkht(ht,"X2aMass"), 500, 0, 5000);
	    hist1d[353]->Fill((*rjrAX2bMass)[cs], fillwt ); //("X2bMass", mkht(ht,"X2bMass"), 500, 0, 5000);
	
	    hist1d[354]->Fill((*rjrDphiMETV)[cs], fillwt );
	    hist1d[355]->Fill((*rjrNJetsJa)[cs], fillwt );
	    hist1d[356]->Fill((*rjrNJetsJb)[cs], fillwt );
	    hist1d[357]->Fill((*rjrSdphiBV)[cs], fillwt );
	    hist1d[358]->Fill((*rjrSdphiDA)[cs], fillwt );
	    hist1d[359]->Fill((*rjrX2aPtS)[cs], fillwt );
	    hist1d[360]->Fill((*rjrX2bPtS)[cs], fillwt );
	
	    hist1d[361]->Fill((*rjrAX2Ave)[cs], fillwt );
	    hist1d[362]->Fill((*rjrAX2Diff)[cs], fillwt );
	    hist1d[363]->Fill((*rjrAX2GMean)[cs], fillwt );
	    hist1d[364]->Fill((*rjrAX2NAve)[cs], fillwt );
	    hist1d[365]->Fill((*rjrAX2NGMean)[cs], fillwt );
	    hist1d[366]->Fill((*rjrAX2NVSum)[cs], fillwt );
	    hist1d[367]->Fill((*rjrAX2VSum)[cs], fillwt );
	
        hist1d[368]->Fill((*selPhoMBetaEql)[0], fillwt );
        hist1d[369]->Fill((*selPhoMBetaPmt)[0], fillwt );

        hist1d[370]->Fill((*rjrMVa)[cs], fillwt );
        hist1d[371]->Fill((*rjrMVb)[cs], fillwt );
        hist1d[372]->Fill(rjrVDiff, fillwt );
        hist1d[373]->Fill(rjrVSVis, fillwt );
        hist1d[374]->Fill(NormMBetaEql, fillwt );
        hist1d[375]->Fill(NormMBetaPmt, fillwt );

        hist2d[1]->Fill( (*rjrAX2Diff)[cs], (*rjrAX2NVSum)[cs], fillwt );
        hist2d[3]->Fill( (*rjrASMass)[cs], (*rjrAX2NVSum)[cs], fillwt );
		hist2d[4]->Fill( (*rjrAX2VSum)[cs], rjrVSVis, fillwt );

	}//<<>>if( rjrAX2NGMean->at(cs) < cutvc && rjrAX2GMean->at(cs) < cutvd )
	}//<<>>for( int i = 0; i < 1; i++ )

}//<<>>void HistMaker::eventLoop(` Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

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

    std::cout << " title test : " << mkht(ht,"jetPt") << std::endl;

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
//rjrX2NVSum
//rjrX2VSum
//
//rjrX2aMass
//rjrX2bMass
//rjrX2aPtS
//rjrX2bPtS
//
//genSigType
//

    hist1d[250] = new TH1D("SCosA", mkht(ht,"SCosA").c_str(), 70, -3.5, 3.5);
    hist1d[251] = new TH1D("SMass", mkht(ht,"SMass").c_str(), 150, 0, 15000);
    hist1d[252] = new TH1D("X2aMass", mkht(ht,"X2aMass").c_str(), 90, 0, 3600);
    hist1d[253] = new TH1D("X2bMass", mkht(ht,"X2bMass").c_str(), 90, 0, 3600);

    hist1d[254] = new TH1D("rjrDphiMETV", mkht(ht,"rjrDphiMETV").c_str(), 70, -3.5, 3.5);
    hist1d[255] = new TH1D("rjrNJetsJa", mkht(ht,"rjrNJetsJa").c_str(), 20, 0, 20);
    hist1d[256] = new TH1D("rjrNJetsJb", mkht(ht,"rjrNJetsJb").c_str(), 20, 0, 20);
    hist1d[257] = new TH1D("rjrSdphiBV", mkht(ht,"rjrSdphiBV").c_str(), 35, 0, 3.5);
    hist1d[258] = new TH1D("rjrSdphiDA", mkht(ht,"rjrSdphiDA").c_str(), 35, 0, 3.5);
    hist1d[259] = new TH1D("rjrX2aPtS", mkht(ht,"rjrX2aPtS").c_str(), 100, 0, 6000);
    hist1d[260] = new TH1D("rjrX2bPtS", mkht(ht,"rjrX2bPtS").c_str(), 100, 0, 6000);

    hist1d[261] = new TH1D("rjrX2Ave", mkht(ht,"rjrX2Ave").c_str(), 90, 0, 3600);
    hist1d[262] = new TH1D("rjrX2Diff", mkht(ht,"rjrX2Diff").c_str(), 100, -1, 1);
    hist1d[263] = new TH1D("rjrX2GMean", mkht(ht,"rjrX2GMean").c_str(), 90, 0, 3600);
    hist1d[264] = new TH1D("rjrX2NAve", mkht(ht,"rjrX2NAve").c_str(), 50, 0, 0.5);
    hist1d[265] = new TH1D("rjrX2NGMean", mkht(ht,"rjrX2NGMean").c_str(), 50, 0, 0.5);
    hist1d[266] = new TH1D("rjrX2NVSum", mkht(ht,"rjrX2NVSum").c_str(), 50, 0, 0.5);
    hist1d[267] = new TH1D("rjrX2VSum", mkht(ht,"rjrX2VSum").c_str(), 90, 0, 3600);

    hist1d[268] = new TH1D("selPhoMBetaEql", mkht(ht,"selPhoMBetaEql").c_str(), 60, 0, 1200);
    hist1d[269] = new TH1D("selPhoMBetaPmt", mkht(ht,"selPhoMBetaPmt").c_str(), 60, 0, 1200);

    hist1d[270] = new TH1D("rjrMVa", mkht(ht,"rjrMVa").c_str(), 90, 0, 3600);
    hist1d[271] = new TH1D("rjrMVb", mkht(ht,"rjrMVb").c_str(), 90, 0, 3600);
    hist1d[272] = new TH1D("rjrVDiff", mkht(ht,"rjrVDiff").c_str(), 100, -1, 1);
    hist1d[273] = new TH1D("rjrVSVis", mkht(ht,"rjrVSVis").c_str(), 90, 0, 3600);

    hist1d[274] = new TH1D("selPhoNMBetaEql", mkht(ht,"selPhoNMBetaEql").c_str(), 100, 0, 2);
    hist1d[275] = new TH1D("selPhoNMBetaPmt", mkht(ht,"selPhoNMBetaPmt").c_str(), 100, 0, 2);

    hist1d[350] = new TH1D("SCosAcA", mkht(ht,"SCosA Ac").c_str(), 70, -3.5, 3.5);
    hist1d[351] = new TH1D("ASMass", mkht(ht,"ASMass").c_str(), 150, 0, 15000);
    hist1d[352] = new TH1D("AX2aMass", mkht(ht,"AX2aMass").c_str(), 90, 0, 3600);
    hist1d[353] = new TH1D("AX2bMass", mkht(ht,"AX2bMass").c_str(), 90, 0, 3600);

    hist1d[354] = new TH1D("rjrDphiMETVcA", mkht(ht,"rjrDphiMETV Ac").c_str(), 70, -3.5, 3.5);
    hist1d[355] = new TH1D("rjrNJetsJacA", mkht(ht,"rjrNJetsJa Ac").c_str(), 20, 0, 20);
    hist1d[356] = new TH1D("rjrNJetsJbcA", mkht(ht,"rjrNJetsJb Ac").c_str(), 20, 0, 20);
    hist1d[357] = new TH1D("rjrSdphiBVcA", mkht(ht,"rjrSdphiBV Ac").c_str(), 35, 0, 3.5);
    hist1d[358] = new TH1D("rjrSdphiDAcA", mkht(ht,"rjrSdphiDA Ac").c_str(), 35, 0, 3.5);
    hist1d[359] = new TH1D("rjrX2aPtScA", mkht(ht,"rjrX2aPtS Ac").c_str(), 100, 0, 6000);
    hist1d[360] = new TH1D("rjrX2bPtScA", mkht(ht,"rjrX2bPtS Ac").c_str(), 100, 0, 6000);

    hist1d[361] = new TH1D("rjrAX2Ave", mkht(ht,"rjrAX2Ave").c_str(), 90, 0, 3600);
    hist1d[362] = new TH1D("rjrAX2Diff", mkht(ht,"rjrAX2Diff").c_str(), 100, -1, 1);
    hist1d[363] = new TH1D("rjrAX2GMean", mkht(ht,"rjrAX2GMean").c_str(), 90, 0, 3600);
    hist1d[364] = new TH1D("rjrAX2NAve", mkht(ht,"rjrAX2NAve").c_str(), 50, 0, 0.5);
    hist1d[365] = new TH1D("rjrAX2NGMean", mkht(ht,"rjrAX2NGMean").c_str(), 50, 0, 0.5);
    hist1d[366] = new TH1D("rjrAX2NVSum", mkht(ht,"rjrAX2NVSum").c_str(), 50, 0, 0.5);
    hist1d[367] = new TH1D("rjrAX2VSum", mkht(ht,"rjrAX2VSum").c_str(), 90, 0, 3600);

    hist1d[368] = new TH1D("selPhoMBetaEqlcA", mkht(ht,"selPhoMBetaEql Ac").c_str(), 60, 0, 1200);
    hist1d[369] = new TH1D("selPhoMBetaPmtcA", mkht(ht,"selPhoMBetaPmt Ac").c_str(), 60, 0, 1200);

    hist1d[370] = new TH1D("rjrMVacA", mkht(ht,"rjrMVa Ac").c_str(), 90, 0, 3600);
    hist1d[371] = new TH1D("rjrMVbcA", mkht(ht,"rjrMVb Ac").c_str(), 90, 0, 3600);
    hist1d[372] = new TH1D("rjrVDiffcA", mkht(ht,"rjrVDiff Ac").c_str(), 100, -1, 1);
    hist1d[373] = new TH1D("rjrVSViscA", mkht(ht,"rjrVSVis Ac").c_str(), 90, 0, 3600);
    hist1d[374] = new TH1D("selPhoNMBetaEqlcA", mkht(ht,"selPhoNMBetaEql Ac").c_str(), 100, 0, 2);
    hist1d[375] = new TH1D("selPhoNMBetaPmtcA", mkht(ht,"selPhoNMBetaPmt Ac").c_str(), 100, 0, 2);

    //------ ecal rechits 300 - 349

    //------  electrons 350 - 400

	//-------- event vars 400 - 450
    hist1d[400] = new TH1D("evtMetE", mkht(ht,"evtMetE").c_str(), 150, 0, 1500 );
    hist1d[401] = new TH1D("evtMetPt", mkht(ht,"evtMetPt ").c_str(), 150, 0, 1500 );	
    hist1d[402] = new TH1D("evtRjrMetPt", mkht(ht,"evtRjrMetPt ").c_str(), 100, 0, 1500 );
    hist1d[403] = new TH1D("evtRawMetPt", mkht(ht,"evtRawMetPt ").c_str(), 100, 0, 1500 );

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	hist2d[0] = new TH2D("NGMvDiff", mkht(ht,"NGMvDiff;rjrX2Diff;rjrX2NVSum").c_str(), 50, -1, 1, 25, 0, 0.5 );
    hist2d[1] = new TH2D("ANGMvADiff", mkht(ht,"ANGMvADiff;rjrAX2Diff;rjrX2ANVSum").c_str(), 50, -1, 1, 25, 0, 0.5 );

    hist2d[2] = new TH2D("SMvRatio", mkht(ht,"SMvRatio;rjrSMass;rjrX2NVSum").c_str(), 150, 0, 15000, 25, 0, 0.5 );
    hist2d[3] = new TH2D("ASMvARatio", mkht(ht,"ASMvARatio;rjrASMass;rjrX2ANVSum").c_str(), 150, 0, 15000, 25, 0, 0.5 );

    hist2d[4] = new TH2D("rjrAX2VSumVrjrVSVis", mkht(ht,"QMX2 Mass v QMVis Mass;AltX2QMeanMass; VisQMeanMass").c_str(), 90, 0, 3600, 90, 0, 3600 );

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
                //auto infilename = "KUCMS_RJR_Skim_List.txt";
                auto infilename1 = "KUCMS_RJR_GMSB100_Skim_List.txt";
                auto infilename3 = "KUCMS_RJR_GMSB300_Skim_List.txt";
                auto infilenameG = "KUCMS_RJR_GJETS_Skim_List.txt";
                auto infilenameQ = "KUCMS_RJR_QCD_Skim_List.txt";

                //std::list<float> cuts = { 10000000 }; 

			 	//std::string lambda = "100";

                std::string outfilenamegj = "KUCMS_GJets_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilenameqcd = "KUCMS_QCD_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename100 = "KUCMS_GMSB_L100_T30_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename300 = "KUCMS_GMSB_L350_T30_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename0 = "KUCMS_GMSB_L"+lambda+"_T22_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename1 = "KUCMS_GMSB_L"+lambda+"_T11_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename2 = "KUCMS_GMSB_L"+lambda+"_T44_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename3 = "KUCMS_GMSB_L"+lambda+"_T33_v19_"; //iso0_Skim_BaseHists.root"; //7

                //std::string ofnending = "wt2_RjrSkim_multiHists_v2.root";
                //std::string ofnending = "NGM0p2_wt2_RjrSkim_multiHists.root";
       
		        //std::string ofnending = "wt2_RjrSkim_v21_phojet_c1b_multiHists.root"; float jrjtype = 1; // 1 = phojet, 0 = phomet
                //std::string ofnending = "wt2_RjrSkim_v21_phomet_c1a_multiHists.root"; float jrjtype = 0; // 1 = phojet, 0 = phomet
                std::string ofnending = "wt2_RjrSkim_v21a_phojet_multiHists.root"; float jrjtype = 1; // 1 = phojet, 0 = phomet
                //std::string ofnending = "wt2_RjrSkim_v21a_phomet_multiHists.root"; float jrjtype = 0; // 1 = phojet, 0 = phomet

                std::string htitlegj = "GJets_v19_";
                std::string htitleqcd = "QCD_v19_";
                std::string htitle100 = "GMSB_L100_v19_T30_";
                std::string htitle300 = "GMSB_L350_v19_T30_";
                //std::string htitle0 = "GMSB_L"+lambda+"_v19_T22_";
                //std::string htitle1 = "GMSB_L"+lambda+"_v19_T11_";
                //std::string htitle2 = "GMSB_L"+lambda+"_v19_T33_";
                //std::string htitle3 = "GMSB_L"+lambda+"_v19_T44_";

                HistMaker base;

				float rjrcut = 0;
                float rjrcut2 = 0;

				float nphos = 1;
                std::string isoline = "genSigPerfect_nSigPho1_";
                std::string isolinebg = "BG_nSigPho0_";			
	
                std::string outfilenamesgj = outfilenamegj + ofnending;
                std::string outfilenamesqcd = outfilenameqcd + ofnending;
                std::string outfilenames100 = outfilename100 + isoline + ofnending;
                std::string outfilenames300 = outfilename300 + isoline + ofnending;
				//std::string outfilenames0 = outfilename0 + isoline + ofnending;
                //std::string outfilenames1 = outfilename1 + isoline + ofnending;
                //std::string outfilenames2 = outfilename2 + isoline + ofnending;
                //std::string outfilenames3 = outfilename3 + isoline + ofnending;

                std::string htitlesgj =  htitlegj + isolinebg;
                std::string htitlesqcd =  htitleqcd + isolinebg;
                std::string htitles100 =  htitle100 + isoline;
                std::string htitles300 =  htitle300 + isoline;
                //std::string htitle0s =  htitle0 + isoline;
                //std::string htitle1s =  htitle1 + isoline;
                //std::string htitle2s =  htitle2 + isoline;
                //std::string htitle3s =  htitle3 + isoline;

                //base.histMaker( listdir, infilename, outfilenames0, htitle0s, 0, 22, nphos, rjrcut, 0 );
                //base.histMaker( listdir, infilename, outfilenames1, htitle1s, 0, 11, nphos, rjrcut, 0 );
                //base.histMaker( listdir, infilename, outfilenames2, htitle2s, 0, 33, nphos, rjrcut, 0 );
                //base.histMaker( listdir, infilename, outfilenames3, htitle3s, 0, 44, nphos, rjrcut, 0 );

                base.histMaker( listdir, infilename1, outfilenames100, htitles100, jrjtype, 30, nphos, rjrcut, rjrcut2 );
                base.histMaker( listdir, infilename3, outfilenames300, htitles300, jrjtype, 30, nphos, rjrcut, rjrcut2 );
				base.histMaker( listdir, infilenameG, outfilenamesgj, htitlesgj, jrjtype, 0, 0, rjrcut, rjrcut2 );
                base.histMaker( listdir, infilenameQ, outfilenamesqcd, htitlesqcd, jrjtype, 0, 0, rjrcut, rjrcut2 );

				//base.histMaker( listdir, infilename, outfilenames, htitles, fillwt?, gensigtype?, #phos?, rjrcut, n/a );

    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

