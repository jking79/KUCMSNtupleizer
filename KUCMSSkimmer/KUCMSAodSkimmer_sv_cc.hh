//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSVSkimmer.hh"

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
//// KUCMSAodSkimmer class ----------------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------------------------------------

//#define DEBUG true
#define DEBUG false

//#define CLSTRMAPS true
#define CLSTRMAPS false

KUCMSAodSkimmer::KUCMSAodSkimmer(){

	LAB = new LabRecoFrame("LAB","LAB");
	S   = new DecayRecoFrame("S","#tilde{S}");

	X2a = new DecayRecoFrame("X2a","#tilde{#chi}_{2a}");
	X2b = new DecayRecoFrame("X2b","#tilde{#chi}_{2b}");	

	Ja = new DecayRecoFrame("Ja","jets_{a}");
	Jb = new DecayRecoFrame("Jb","jets_{b}");
	J1a = new VisibleRecoFrame("J1a","jets_{1a}");
	J2a = new VisibleRecoFrame("J2a","jets_{2a}");
	J1b = new VisibleRecoFrame("J1b","jets_{1b}");
    J2b = new VisibleRecoFrame("J2b","jets_{2b}");

	X1a = new InvisibleRecoFrame("X1a","#tilde{#chi}_{1a}");
	X1b = new InvisibleRecoFrame("X1b","#tilde{#chi}_{1b}");
	
	LAB->SetChildFrame(*S);
	S->AddChildFrame(*X2a);
	S->AddChildFrame(*X2b);
	X2a->AddChildFrame(*X1a);
	X2b->AddChildFrame(*X1b);
	X2a->AddChildFrame(*Ja);
	X2b->AddChildFrame(*Jb);
	Ja->AddChildFrame(*J1a);
    Ja->AddChildFrame(*J2a);
    Jb->AddChildFrame(*J1b);
    Jb->AddChildFrame(*J2b);

	if(!LAB->InitializeTree()){ std::cout << "Problem initializing tree" << std::endl; }
	
	INV = new InvisibleGroup("INV","Invisible System");
	INV->AddFrame(*X1a);
	INV->AddFrame(*X1b);
	
	InvM = new SetMassInvJigsaw("InvM", "Set inv. system mass");
	INV->AddJigsaw(*InvM);
	
	InvEta = new SetRapidityInvJigsaw("InvEta", "Set inv. system rapidity");
	INV->AddJigsaw(*InvEta);
	InvEta->AddVisibleFrames(S->GetListVisibleFrames());
	
	InvSplit = new MinMassesSqInvJigsaw("InvSplit", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}", 2);
    //InvSplit = new ContraBoostInvJigsaw("InvSplit", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}");
	INV->AddJigsaw(*InvSplit);
	//InvSplit->AddVisibleFrame(Ja->GetVisiblFrames(), 0);
	//InvSplit->AddVisibleFrame(Jb->GetVisiblFrames(), 1);
    InvSplit->AddVisibleFrame(*Ja, 0);
    InvSplit->AddVisibleFrame(*Jb, 1);
	InvSplit->AddInvisibleFrame(*X1a, 0);
	InvSplit->AddInvisibleFrame(*X1b, 1);
	
	COMB_J =  new CombinatoricGroup("COMB_J", "Combinatoric System of Jets");
	CombSplit_J = new MinMassesSqCombJigsaw("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
    CombSplit_Ja = new MinMassesCombJigsaw("CombSplit_Ja", "Minimize M_{J1a} + M_{J2a} ");
    CombSplit_Jb = new MinMassesCombJigsaw("CombSplit_Jb", "Minimize M_{J1b} + M_{J2b} ");

	COMB_J->AddFrame(*J1a);
    COMB_J->AddFrame(*J2a);
	COMB_J->SetNElementsForFrame(*J1a, 0);
    COMB_J->SetNElementsForFrame(*J2a, 1);
	COMB_J->AddFrame(*J1b);
    COMB_J->AddFrame(*J2b);
	COMB_J->SetNElementsForFrame(*J1b, 0);
    COMB_J->SetNElementsForFrame(*J2b, 1);
	
	COMB_J->AddJigsaw(*CombSplit_J);
	CombSplit_J->AddCombFrames(Ja->GetListVisibleFrames(), 0);
	CombSplit_J->AddCombFrames(Jb->GetListVisibleFrames(), 1);
	CombSplit_J->AddObjectFrame(*X2a, 0);
	CombSplit_J->AddObjectFrame(*X2b, 1);// check syntax from example
	 
	//add bloc do b
    COMB_J->AddJigsaw(*CombSplit_Ja);
    CombSplit_Ja->AddCombFrame(*J1a, 0);
    CombSplit_Ja->AddCombFrame(*J2a, 1);
    CombSplit_Ja->AddObjectFrame(*J1a, 0);
    CombSplit_Ja->AddObjectFrame(*J2a, 1);// check syntax from example
    COMB_J->AddJigsaw(*CombSplit_Jb);
    CombSplit_Jb->AddCombFrame(*J1b, 0);
    CombSplit_Jb->AddCombFrame(*J2b, 1);
    CombSplit_Jb->AddObjectFrame(*J2b, 0);
    CombSplit_Jb->AddObjectFrame(*J2b, 1);// check syntax from example

	if(!LAB->InitializeAnalysis()) std::cout << "Problem initializing analysis tree" << std::endl;

	if( false ){
	
		TreePlot tree_plot("TreePlot","TreePlot");

	  	for(int t = 0; t < 2; t++){
	  		tree_plot.SetTree(*LAB);
	  		tree_plot.Draw("ANA_tree", "Reconstruction Tree");

	  		tree_plot.SetTree(*COMB_J);
	  		tree_plot.Draw("ANA_comb", "Combinatoric Jigsaws for jets/sig photons");

	  		//tree_plot.SetTree(*COMB_L);
	  		//tree_plot.Draw("ANA_comb_L", "Combinatoric Jigsaws for leps");

	  		tree_plot.SetTree(*INV);
	  		tree_plot.Draw("ANA_inv", "Invisible Jigsaws");
 
	  	}//<<>>for(int t = 0; t < 2; t++)

		tree_plot.WriteOutput("trees.root");

	}//<<>>treeplot on off

    // Cali Tags : Tags for calibrations to use 
    //std::string r2EOY( "EG_EOY_MINI" ); 
    //std::string r2Fall17AOD( "RunIIFall17DRPremix" ); 
    //std::string r2Fall17MINIAOD( "RunIIFall17MiniAODv2" ); 
    //std::string r2UL( "UL_R2_MINI" ); 

	//timeCali = new KUCMSTimeCalibration();
	//timeCali->setTag(r2UL);

}//<<>>KUCMSAodSkimmer::KUCMSAodSkimmer()

KUCMSAodSkimmer::~KUCMSAodSkimmer(){

    delete LAB;
    delete S;
    delete X2a;
    delete X2b;
    delete Ja;
    delete Jb;
    delete X1a;
    delete X1b;
    delete J1a;
    delete J2a;
    delete J1b;
    delete J2b;      

    delete INV;
    delete InvM;
    delete InvEta;
    delete InvSplit;
      
    delete COMB_J;
    delete CombSplit_J;
    delete CombSplit_Ja;
    delete CombSplit_Jb;

	//delete timeCali;
      
}//<<>>KUCMSAodSkimmer::~KUCMSAodSkimmer()

void KUCMSAodSkimmer::kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect, bool doSVs, int skipCnt, bool useEvtWgts ){


	useEvtWgt = useEvtWgts;
    doGenInfo = hasGenInfo;
    const std::string disphotreename = "tree/llpgtree";
    const std::string configtreename = "tree/configtree";
    //std::string inpath, infiles, key, 
	std::string masterstr; 
    //int mct;
	//float crossSection, gmsblam, gmsbct, mcw;
	std::cout << "Processing Input Lists for : " << infilelist << std::endl;
	std::ifstream masterInfile(listdir+infilelist);
	//while( masterInfile >> inpath >> infiles >> key >> crossSection >> gmsblam >> gmsbct >> mcwgt >> mctype ){

    while( std::getline( masterInfile, masterstr ) ){

		if( DEBUG ) std:: cout << masterstr << std::endl;
        if( masterstr[0] == '#' ) continue;
		if( masterstr == " " ) continue;
		auto instrs = splitString( masterstr, " " );
		if( DEBUG ) std:: cout << instrs.size() << std::endl;
        if( instrs.size() < 8 ) continue;

        auto inpath = instrs[0];
		auto infiles = instrs[1];
		auto key = instrs[2];
		auto crossSection = std::stof( instrs[3] );
		auto gmsbgm = std::stof( instrs[4] );
        auto gmsbxm = std::stof( instrs[5] );
        auto mcw = std::stof( instrs[6] );
        auto mct = std::stoi( instrs[7] );		
		if( DEBUG ) std:: cout << "InPath: " << inpath << std::endl;
        if( DEBUG ) std:: cout << "InFile: " << infiles << std::endl;
        if( DEBUG ) std:: cout << "Key: " << key << std::endl;
        if( DEBUG ) std:: cout << "XSec: " << crossSection << std::endl;
        if( DEBUG ) std:: cout << "GM: " << gmsbgm << std::endl;
        if( DEBUG ) std:: cout << "XM: " << gmsbxm << std::endl;
        if( DEBUG ) std:: cout << "MCw: " << mcw << std::endl;
        if( DEBUG ) std:: cout << "MCt: " << mct << std::endl;

        dataSetKey = key;
        xsctn = crossSection;
        gmass = gmsbgm; // = 0 if not gmsb
        xmass = gmsbxm; // = 0 if not gmsb
        mcwgt = mcw; // default 1
        mctype = mct; // 0 = fullsim

		std::cout << "Processing Events for : " << infiles << std::endl;
	    TChain* fInTree = new TChain(disphotreename.c_str());
        TChain* fInConfigTree = new TChain(configtreename.c_str());
	    std::cout << "Adding files to TChain." << std::endl;
	    std::cout << " - With : " << listdir+infiles << " >> " << fInTree << std::endl;
	    std::string str;
		if( not DEBUG ) std::cout << "--  adding files";

		int nfiles = 0;	
		if( key !=  "single" ){
		std::ifstream infile(listdir+infiles);
	    while( std::getline( infile, str ) ){
			nfiles++;
			if( skipCnt != 0 && ( nfiles%skipCnt != 0 ) ) continue;
	        auto tfilename = eosdir + inpath + str;
			fInTree->Add(tfilename.c_str());
			fInConfigTree->Add(tfilename.c_str());
	        if(DEBUG) std::cout << "--  adding file: " << tfilename << std::endl; else std::cout << ".";
			//if(DEBUG) break;
	    }//<<>>while (std::getline(infile,str))
		} else { // single infile and not a list of infiles
			auto tfilename = eosdir + inpath + infiles;
			fInTree->Add(tfilename.c_str());
			fInConfigTree->Add(tfilename.c_str());
			nfiles++;			
		}//<<>>if( key !=  "test" )
		if( not DEBUG ) std::cout << std::endl;
		if( nfiles == 0 ){ std::cout << " !!!!! no input files !!!!! " << std::endl; return; }
	
		auto fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
	    auto fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");

        Init( fInTree, doGenInfo );
        initHists();
        setOutputBranches(fOutTree);

        SetupDetIDsEB_(DetIDMap);
        SetupDetIDsEE_(DetIDMap);

        startJobs(); // clear && init count varibles

		// ntuple event counts and weights
		// setup config tree inputs

		int nTotEvts;
		int nFltrdEvts;
		float sumEvtWgt;
		float sumFltrdEvtWgt;
		int nMetFltrdEvts;
		int nPhoFltrdEvts;

        TBranch *b_nTotEvts;
        TBranch *b_nFltrdEvts;
        TBranch *b_sumEvtWgt;
        TBranch *b_sumFltrdEvtWgt;
        TBranch *b_nMetFltrdEvts;
        TBranch *b_nPhoFltrdEvts;

    	fInConfigTree->SetBranchAddress("nTotEvts", &nTotEvts, &b_nTotEvts);
        fInConfigTree->SetBranchAddress("nFltrdEvts", &nFltrdEvts, &b_nFltrdEvts);
        fInConfigTree->SetBranchAddress("sumEvtWgt", &sumEvtWgt, &b_sumEvtWgt);
		fInConfigTree->SetBranchAddress("sumFltrdEvtWgt", &sumFltrdEvtWgt, &b_sumFltrdEvtWgt);
        fInConfigTree->SetBranchAddress("nMetFltrdEvts", &nMetFltrdEvts, &b_nMetFltrdEvts);
        fInConfigTree->SetBranchAddress("nPhoFltrdEvts", &nPhoFltrdEvts, &b_nPhoFltrdEvts);

    	auto nConfigEntries = fInConfigTree->GetEntries();
    	std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
    	configCnts.clear();
		configWgts.clear();
    	for (Long64_t centry = 0; centry < nConfigEntries; centry++){

        	auto entry = fInConfigTree->LoadTree(centry);

        	b_nTotEvts->GetEntry(entry);   //!
        	b_nFltrdEvts->GetEntry(entry);   //!
        	b_sumEvtWgt->GetEntry(entry);   //!  
        	b_sumFltrdEvtWgt->GetEntry(entry);   //!
        	b_nMetFltrdEvts->GetEntry(entry);   //!
        	b_nPhoFltrdEvts->GetEntry(entry);   //!

			cutflow["nTotEvts"] += nTotEvts;
			cutflow["nFltrdEvts"] += nFltrdEvts;

			configCnts["nTotEvts"] += nTotEvts;
            configCnts["nFltrdEvts"] += nFltrdEvts;
            configWgts["sumEvtWgt"] += useEvtWgt ? sumEvtWgt : nTotEvts;
            configWgts["sumFltrdEvtWgt"] += useEvtWgt ? sumFltrdEvtWgt : nFltrdEvts;
            configCnts["nMetFltrdEvts"] += nMetFltrdEvts;
            configCnts["nPhoFltrdEvts"] += nPhoFltrdEvts;


    	}//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)

        std::cout << "configCnts ( ";
    	for( auto item : configCnts ){
         	std::cout << item.first <<  " " << item.second << " ";
    	}//<<>>for( auto item : configInfo )
        std::cout << ")" << std::endl;

		std::cout << "configWgts ( ";
        for( auto item : configWgts ){
            std::cout << item.first <<  " " << item.second << " ";
        }//<<>>for( auto item : configInfo )
        std::cout << ")" << std::endl;

	    std::cout << "Setting up For Main Loop." << std::endl;
		int loopCounter(1000);
	    auto nEntries = fInTree->GetEntries();
	    if(DEBUG){ nEntries = 1000; loopCounter = 100; }
	    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
        nEvents = nEntries;
	    for (Long64_t centry = 0; centry < nEntries; centry++){

	        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
	        auto entry = fInTree->LoadTree(centry);
			if(DEBUG) std::cout << " -- Getting Branches " << std::endl;
			getBranches( entry, doGenInfo );
			geCnts.clear();
            geVars.clear();
            //geVects.clear();
			if( genSigPerfect ) geVars.set( "genSigPerfect", 1 ); else geVars.set( "genSigPerfect", 0 );
            if( doSVs ) geVars.set( "doSVs", 1 ); else geVars.set( "doSVs", 0 );
			if(DEBUG) std::cout << " -- Event Loop " << std::endl;
			auto saveToTree = eventLoop(entry);
            if( saveToTree ){ fOutTree->Fill(); }

	    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
		fillConfigTree( fConfigTree ); 

		endJobs();
	
	    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

		auto ext = splitString( infiles, "." );
		std::string extOutFileName( ext[0] + outfilename );
	    TFile* fOutFile = new TFile( extOutFileName.c_str(), "RECREATE" );
	    fOutFile->cd();

		fOutTree->Write();
		fConfigTree->Write();
	
		for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
	    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
	    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }
	
		if( CLSTRMAPS ){
			nMaps = 0;
			for( int it = 0; it < nEBEEMaps; it++ ){ 
				ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
				ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
				ebeeMapR[it]->Write(); delete ebeeMapR[it];
			}//<<>>for( int it = 0; it < nEBEEMaps; it++ )
		}//<<>>f( clsttrMaps )
	
		std::cout << "Finished processing events for : " << infiles << std::endl;
	
	    fOutFile->Close();
	
		delete fInTree;
		delete fOutTree;
		delete fConfigTree;
		delete fOutFile;

	}//<<>>while (std::getline(infile,str))

    std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

}//<<>>void kucmsSkimmer

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// event loop and startup jobs 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


void KUCMSAodSkimmer::startJobs(){

    std::cout << "Runing StartJobs." << std::endl;

	nEvents = 0;
	nSelectedEvents = 0;

    //sumEvtGenWgt = 0.0;

    configCnts.clear();
    configWgts.clear();
    configCnts["nTotEvts"] = 0;
    configCnts["nFltrdEvts"] = 0;
    configWgts["sumEvtWgt"] = 0;
    configWgts["sumFltrdEvtWgt"] = 0;
    configCnts["nMetFltrdEvts"] = 0;
    configCnts["nPhoFltrdEvts"] = 0;

    cutflow.clear();
    cutflow["nTotEvts"] = 0;
    cutflow["nFltrdEvts"] = 0;
    cutflow["met150"] = 0;
    cutflow["m_gt2jets"] = 0;
    cutflow["mj_gt1phos"] = 0;
    cutflow["mjp_leadPhoPt30"] = 0;
    cutflow["sel_m"] = 0;
    cutflow["sel_j"] = 0;
    cutflow["sel_p"] = 0;
    cutflow["sel_ppt"] = 0;

};//<<>>void KUCMSAodSkimmer::startJobs()

bool KUCMSAodSkimmer::eventLoop( Long64_t entry ){

	// counts events and saves event varibles
	// --------------------------------------
	processEvntVars();	
	processRechits();
	processMet();
	processPhotons();
	processElectrons();
	processMuons();
	processJets();
	if( doGenInfo ){ processGenParticles(); }

	// select events to process and store
	//--------------------------------------
	auto saveToTree = eventSelection();	
	if( saveToTree ){ processRJR(0,true); processRJR(1,false); }
	return saveToTree;

}//<<>>void KUCMSAodSkimmer::eventLoop( Long64_t entry )

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

void KUCMSAodSkimmer::processEvntVars(){

    if( DEBUG ) std::cout << "Finding PV Vars" << std::endl;
    selEvtVars.clearBranches(); // <<<<<<<   must do
	bool doEVSVs = true;
	if( not doEVSVs ) std::cout << " ------ !!!!!!!!!! Fill SV Branches turned off !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

	// calc
    selEvtVars.fillBranch( "PVx", PV_x );
    selEvtVars.fillBranch( "PVy", PV_y );
    selEvtVars.fillBranch( "PVz", PV_z ); 
   
	//fill

    if( DEBUG ) std::cout << "Finding Event wts Vars" << std::endl;
    selEvtVars.fillBranch( "dsKey", dataSetKey );
    float evtGenWgt = 1;
	if( useEvtWgt && doGenInfo ) evtGenWgt = Evt_genWgt;
    selEvtVars.fillBranch( "evtGenWgt", evtGenWgt );
    //sumEvtGenWgt += Evt_genWgt;
    selEvtVars.fillBranch( "evtXSection", xsctn );

	float fillWgt = ( ( xsctn * 1000 ) * evtGenWgt  ) / configWgts["sumEvtWgt"];
	selEvtVars.fillBranch( "evtFillWgt", fillWgt );

	// SVs
	if( DEBUG ) std::cout << "Finding SV Vars" << std::endl;

	int nSVs = 0;
	if( doEVSVs )  nSVs = Vertex_mass->size();
	int nLsv = 0;
	int nHsv = 0;

	if( doEVSVs ){
	for( int svit = 0; svit < nSVs; svit++ ){

		float mass = (*Vertex_mass)[svit];
        uInt ntrack = (*Vertex_nTracks)[svit];
        bool peleid = (*Vertex_passLooseElectronID)[svit];
        bool pmuonid = (*Vertex_passLooseMuonID)[svit];

		if( ( ntrack == 2 ) && ( peleid || pmuonid ) ) nLsv++;
		if( ( ntrack >= 5 ) && ( mass/ntrack > 1 ) ) nHsv++;

	}//<<>>for( svit = 0; svit < nSVs; scit++ )
	}//<<>>if( doSVs )

	selEvtVars.fillBranch( "SV_nLeptonic", nLsv );
	selEvtVars.fillBranch( "SV_nHadronic", nHsv );
	geVars.set("nSVLep", nLsv );
	geVars.set("nSVHad", nHsv );

    selEvtVars.fillBranch( "Flag_BadChargedCandidateFilter", Flag_BadChargedCandidateFilter );//not suggested
    selEvtVars.fillBranch( "Flag_BadPFMuonDzFilter", Flag_BadPFMuonDzFilter );//suggested
    selEvtVars.fillBranch( "Flag_BadPFMuonFilter", Flag_BadPFMuonFilter );//suggested
    selEvtVars.fillBranch( "Flag_EcalDeadCellTriggerPrimitiveFilter", Flag_EcalDeadCellTriggerPrimitiveFilter );//suggested
    selEvtVars.fillBranch( "Flag_HBHENoiseFilter", Flag_HBHENoiseFilter );//suggested
    selEvtVars.fillBranch( "Flag_HBHENoiseIsoFilter", Flag_HBHENoiseIsoFilter );//suggested
    selEvtVars.fillBranch( "Flag_ecalBadCalibFilter", Flag_ecalBadCalibFilter );//suggested
    selEvtVars.fillBranch( "Flag_eeBadScFilter", Flag_eeBadScFilter );//suggested
    selEvtVars.fillBranch( "Flag_globalSuperTightHalo2016Filter", Flag_globalSuperTightHalo2016Filter );//suggested
    selEvtVars.fillBranch( "Flag_goodVertices", Flag_goodVertices );//suggested
    selEvtVars.fillBranch( "Flag_hfNoisyHitsFilter", Flag_hfNoisyHitsFilter );//optional
	bool metfilter = Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_HBHENoiseFilter;
	metfilter = metfilter && Flag_HBHENoiseIsoFilter && Flag_ecalBadCalibFilter && Flag_eeBadScFilter && Flag_goodVertices;
	selEvtVars.fillBranch( "Flag_MetFilter",metfilter);

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

	//intilize
	selMet.clearBranches(); // <<<<<<<   must do

	//calc
	//auto met = std::sqrt(sq2(Met_Cpx)+sq2(Met_Cpy));

    geVars.set("metPx", Met_Cpx );
    geVars.set("metPy", Met_Cpy );
	geVars.set("cmet", Met_CPt );

	//fill
	selMet.fillBranch( "cmet", Met_CPt );
    selMet.fillBranch( "cmetPx", Met_Cpx );
    selMet.fillBranch( "cmetPy", Met_Cpy );
    selMet.fillBranch( "met", Met_pt );
    selMet.fillBranch( "metPx", Met_px );
    selMet.fillBranch( "metPy", Met_py );

}//<<>>void KUCMSAodSkimmer::processMet()

void KUCMSAodSkimmer::processRechits(){

	// initilize
    //selECALRecHit.clearBranches(); // <<<<<<<   must do
	bool checkRHs = false;

	// calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    auto nRecHits = ECALRecHit_ID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

		//auto id = (*ECALRecHit_ID)[it];
		//auto idinfo = DetIDMap[id];
		//if( idinfo.ecal == ECAL::EB ){
        if( true ){

			hist1d[0]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[1]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[2]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[5]->Fill( (*ECALRecHit_time)[it], 1 );
			//auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] );

		}//<<>>if( (*rhSubdet)[it] == 0 )

	}//<<>>for( int it = 0; it < nRecHits; it++ )

    int nMissingRechits = 0;
	auto nSCs = SuperCluster_nSuperCluster;
	for( int it = 0; it < nSCs; it++ ){

		int nMissingRhInSC = 0;
		bool oot = (*SuperCluster_isOot)[it];
		bool exc = (*SuperCluster_excluded)[it];
		bool orig = (*SuperCluster_original)[it];
		hist1d[6]->Fill( (*SuperCluster_energyRaw)[it], 1 );
    	//hist1d[3] = new TH1D("sctype","SC !Orig, Orig, OOT, Excl",4,0,4);
    	//hist1d[4] = new TH1D("scorigtype","Orig+OOT, Orig+!OOT, Orig+Exc, !Orig+OOT, !Orig+Exc, !Orig+!OOT",6,0,6);
		if( orig ) hist1d[3]->Fill(1); else hist1d[3]->Fill(0);
        if( oot ) hist1d[3]->Fill(2);
        if( exc ) hist1d[3]->Fill(3);
		if( orig && oot ) hist1d[4]->Fill(0);
        if( orig && not oot ) hist1d[4]->Fill(1);
        if( orig && exc ) hist1d[4]->Fill(2);
        if( not orig && oot ) hist1d[4]->Fill(3);
        if( not orig && not oot ) hist1d[4]->Fill(5);
        if( not orig && exc ) hist1d[4]->Fill(4);
		int xfill = orig ? 1 : 0;
		int yfill = oot ? 1 : exc ? 2 : 0;
		hist2d[1]->Fill(xfill,yfill);

		if( checkRHs ){
		bool found = false;
		std::vector<unsigned int> scrhids = (*SuperCluster_rhIds)[it];
		int nSCRhids = scrhids.size();
		//std::cout << " -- Checking : " << it << " of " << nSCRhids << " SCs " << std::endl;
		for( int sciter = 0; sciter < nSCRhids; sciter++  ){ 
			auto scrhid = scrhids[sciter];
			for( int iter = 0; iter < nRecHits; iter++ ){
				auto rhid = (*ECALRecHit_ID)[iter];
				//std::cout << " -- Checking : " << rhid << " == " << scrhid << std::endl;
				if( scrhid == rhid ){ found = true; break; } 
			}//<<>>for( int iter = 0; iter < nRecHits; it++ )
			if( not found ){ nMissingRechits++; nMissingRhInSC++; }// std::cout << " ---- !!!!!!  Rechit in SC and not in Rechits !!!!!!! " << std::endl;
		}//<<>>for( auto scrhid : (*SuperCluster_rhIds)[it] )
		}//<<>>if( checkRHs )		
		if( nMissingRhInSC > 0 ) std::cout << " !! Missing " << nMissingRhInSC << " rechits in SC " << it << " w/ OOT" << oot << " && ORG " << orig << std::endl;

	}//<<>>for( int it = 0; it < nSCs; it++ )
	if( nMissingRechits > 0 ) std::cout << " !! ------ !! Missing " << nMissingRechits << " rechits in " << nSCs << " SCs "<< std::endl;

	// fill
    auto nETs = Track_pt->size();
    for( int it = 0; it < nETs; it++ ){ hist1d[7]->Fill( (*Track_pt)[it], 1 ); }

}//<<>>void KUCMSAodSkimmer::processRechits()

void KUCMSAodSkimmer::processGenParticles(){

	// initilize
	selGenPart.clearBranches(); // <<<<<<<   must do

	// calc
    if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
    //------------  genparts ------------------------

	int nLSPfXfSg = 0;
	int nQfSg = 0;
    int nQfSqk = 0;
	int nPHOfX = 0;
	int nZfX = 0;
	int nSGlue = 0;
	int nSQuark = 0;
	int nX234 = 0;
	int nLZX = 0;
	int nQfZ = 0;
	int nN0fsqk = 0;
    int nN0fsg = 0;

	//std::cout << "New Event ----------------------------------" << std::endl;
	int nGenParts = Gen_pdgId->size();
    for( int it = 0; it < nGenParts; it++ ){

		float displacment = (*Gen_momDisplacment)[it];
        float energy = (*Gen_energy)[it];
        float eta = (*Gen_eta)[it];
        float phi = (*Gen_phi)[it];
        float pt = (*Gen_pt)[it];
        uInt  pdgId = (*Gen_pdgId)[it];
		int momIndx = (*Gen_motherIdx)[it];
        int   susId = (*Gen_susId)[it];
        float charge = (*Gen_charge)[it];
        float mass = (*Gen_mass)[it];
        int status = (*Gen_status)[it];
        float vx = (*Gen_vx)[it];
        float vy = (*Gen_vy)[it];
        float vz = (*Gen_vz)[it];
        float px = (*Gen_px)[it];
        float py = (*Gen_py)[it];
        float pz = (*Gen_pz)[it];

		float mommass = ( momIndx > -1 ) ? (*Gen_mass)[momIndx] : -1;
		uInt mompdg = ( momIndx > -1 ) ? (*Gen_pdgId)[momIndx] : 0;
		float mompx = ( momIndx > -1 ) ? (*Gen_px)[momIndx] : -1;
        float mompy = ( momIndx > -1 ) ? (*Gen_py)[momIndx] : -1;
        float mompz = ( momIndx > -1 ) ? (*Gen_pz)[momIndx] : -1;
		float genmomp = hypo( mompx, mompy, mompz );
		//float beta = ( mommass > 0 ) ? genmomp/mommass : -1;
		//float gama = ( beta >= 0 ) ? 1/std::sqrt( 1 - beta*beta ) : -1; 
		float gbeta = ( mommass > 0 ) ? genmomp/mommass : -1;
        //float gbeta = ( gama >= 0 && beta >= 0 ) ? gama*beta : -1;
		float ctau = ( ( gbeta >= 0 ) && ( displacment >= 0 ) ) ? displacment/gbeta : -1;
		if( mompdg == 1000023 ){ selGenPart.fillBranch( "genXMomCTau", ctau ); }
		//if( mompdg != 0 ){
        if( false ){
			std::cout << " ctau for : " << pdgId << " mother: " << mompdg << " with mommass " << mommass;
			std::cout << " genp: " << genmomp << " gbeta: " << gbeta << " dis: " << displacment;
        	std::cout << " ctau: " << ctau << std::endl; 
		}//<<>>if( mompdg != 0 )

		if( ( pdgId > 1000000 ) && ( pdgId < 1000007 ) ){ nSQuark++; selGenPart.fillBranch( "genSQMass", mass ); }
        if( pdgId == 1000021 ){ nSGlue++; selGenPart.fillBranch( "genSGMass", mass ); }
        if( pdgId == 1000023 ){ selGenPart.fillBranch( "genLLPMass", mass ); }
        if( pdgId == 1000022 ){ selGenPart.fillBranch( "genLSPMass", mass ); }
        if( pdgId == 1000039 ){ selGenPart.fillBranch( "genGrvtinoMass", mass ); }

		//if( pdgId < 7 ) continue;
		bool hasMom( ( momIndx < nGenParts ) && ( momIndx > -1 ) );
		int gMomIndx = hasMom ? Gen_motherIdx->at(momIndx) : -1;
		bool hasGrandMom( gMomIndx > -1 );
        bool lsp( ( pdgId == 1000039 ) || ( ( pdgId == 1000022 ) && ( status == 1 ) ) );
		bool hasX234( ( pdgId > 1000022 ) && ( pdgId < 1000038 ) );
		bool fromX( hasMom && ( ( Gen_pdgId->at(momIndx) == 1000022 ) || ( Gen_pdgId->at(momIndx) == 1000023 ) ) );
		bool fromSg( hasMom && ( Gen_pdgId->at(momIndx) == 1000021 ) );
        bool fromSqkL( hasMom && ( Gen_pdgId->at(momIndx) > 1000000 ) && ( Gen_pdgId->at(momIndx) < 1000009 ) );
        bool fromSqkR( hasMom && ( Gen_pdgId->at(momIndx) > 2000000 ) && ( Gen_pdgId->at(momIndx) < 2000009 ) );
		bool fromSqk( fromSqkL || fromSqkR );
		bool fromZ( hasMom && ( Gen_pdgId->at(momIndx) == 23 ) );
		bool momFromX( hasGrandMom && ( ( Gen_pdgId->at(gMomIndx) == 1000022 ) || ( Gen_pdgId->at(gMomIndx) == 1000023 ) ) );
		bool momFromSg( hasGrandMom && ( Gen_pdgId->at(gMomIndx) == 1000021 ) );
        bool momFromSqkL( hasGrandMom && ( Gen_pdgId->at(gMomIndx) > 1000000 ) && ( Gen_pdgId->at(gMomIndx) < 1000009 ) );
        bool momFromSqkR( hasGrandMom && ( Gen_pdgId->at(gMomIndx) > 2000000 ) && ( Gen_pdgId->at(gMomIndx) < 2000009 ) );
        bool momFromSqk( momFromSqkL || momFromSqkR );
		bool quark( pdgId < 9 );
		bool photon( pdgId == 22 ); 
        bool zee( pdgId == 23 );
		bool lept( ( pdgId > 10 ) && ( pdgId < 19 ) );
		bool N0( pdgId == 1000022 );

    	//if( lsp && fromX && momFromSg ) nLSPfXfSg++;
        //if( lsp && fromX ) nLSPfXfSg++;
        if( lsp ) nLSPfXfSg++;
    	if( quark && fromSg ) nQfSg++;
        if( quark && fromSqk ) nQfSqk++;
    	if( photon && fromX ) nPHOfX++;
    	if( zee && fromX ) nZfX++;
		if( hasX234 ) nX234++;
		if( lept && fromZ && momFromX ) nLZX++;
		if( quark && fromSqk ) nQfSqk++;
		if( quark && fromZ && momFromX ) nQfZ++;
        if( N0 && fromSqk ) nN0fsqk++;
        if( N0 && fromSg ) nN0fsg++;

		selGenPart.fillBranch( "genPartEnergy", energy );
        selGenPart.fillBranch( "genPartEta", eta );
        selGenPart.fillBranch( "genPartPhi", phi );
        selGenPart.fillBranch( "genPartPt", pt );
        selGenPart.fillBranch( "genPartPdgId", pdgId );
        selGenPart.fillBranch( "genPartSusId", susId );
        selGenPart.fillBranch( "genPartMomIdx", momIndx );
		selGenPart.fillBranch( "genMomCTau", ctau );
        selGenPart.fillBranch( "genCharge", charge );
        selGenPart.fillBranch( "genMass", mass );
        selGenPart.fillBranch( "genStatus", status );
        selGenPart.fillBranch( "genVx", vx );
        selGenPart.fillBranch( "genVy", vy );
        selGenPart.fillBranch( "genVz", vz );
        selGenPart.fillBranch( "genPx", px );
        selGenPart.fillBranch( "genPy", py );
        selGenPart.fillBranch( "genPz", pz );
        //selGenPart.fillBranch( "genSigType", pz );

    }//<<>>for( int it = 0; it < nGenParts; it++ )

	//bool hasLSP( nLSPfXfSg > 0 );
	bool hasLSP( true );
	bool noX234( nX234 == 0 );

    bool has2Nfsqk( nN0fsqk == 2 );
    bool has2Nfsg( nN0fsg == 2 );
	bool has4QfSg( nQfSg > 3 );
    bool has2QfSqk( nQfSqk > 1 );

	bool has2PfN( nPHOfX == 2 );
	bool has2ZfN( nZfX == 2 );

    bool hasLLfz( nLZX == 4 );
    bool hasLfz( nLZX == 2 );
    bool hasNoLfz( nLZX == 0 );

	//if( nZfX > 1 && nLZX > 0 && nX234 == 0 ){
    //if( nX234 == 0 ){
    //if( nQfSg || nQfSqk ){ 
	if( false ){
    	std::cout << "GenEventType::";
		std::cout << " nPHOfX: " << nPHOfX << " nZfX: " << nZfX;
		std::cout << " nX234: " << nX234 << " nN0fsqk: " << nN0fsqk << " nN0fsg: " << nN0fsg;
		std::cout << " nLZX: " << nLZX  << " nQfZ: " << nQfZ; 
		std::cout << " nQfSg: " << nQfSg << " nQfSqk: " << nQfSqk;
		std::cout << std::endl;
	}//<<>>if(

    bool isSTqp( has2Nfsqk && has2QfSqk && has2PfN && noX234 );
	bool isSTqqp( has2Nfsg && has4QfSg && has2PfN && noX234 );
	bool isSTqqzll( has2Nfsg && has4QfSg && has2ZfN && hasLLfz && noX234 ); 
    bool isSTqqzl( has2Nfsg && has4QfSg && has2ZfN && hasLfz && noX234 );
    bool isSTqqz( has2Nfsg && has4QfSg && has2ZfN && hasNoLfz && noX234 );	

    //if( nZfX > 1 && nLZX > 0 && nX234 == 0 ){
    //if( nX234 == 0 ){ 
    //if( isSTqp || isSTqqp || isSTqqzll || isSTqqzl || isSTqqz ){
    //  std::cout << "GenSTFlag:";
	//  std::cout << " isSTqqp: " << isSTqqp << " isSTqqzll: " << isSTqqzll << " isSTqqzl: " << isSTqqzl;
	//  std::cout << " isSTqqz: " << isSTqqz << " isSTqp: " << isSTqp;
    //  std::cout << std::endl;
    //}//<<>>if(

	
    selGenPart.fillBranch( "genSTFlagQQP", isSTqqp );
	selGenPart.fillBranch( "genSTFlagQQZLL", isSTqqzll );
    selGenPart.fillBranch( "genSTFlagQQZL", isSTqqzl );
    selGenPart.fillBranch( "genSTFlagQQZ", isSTqqz );
    selGenPart.fillBranch( "genSTFlagQP", isSTqp );
	selGenPart.fillBranch( "genSigType", Gen_susEvtType->at(0) );

}//<<>>void KUCMSAodSkimmer::processGenParticles()

void KUCMSAodSkimmer::processCalojets(){

	// initilize
    // selCaloJets.clearBranches(); // <<<<<<<   must do

	/* calc
	int nCaloJets = 0;
	for( int it = 0; it < nCaloJets; it++ ){


	}//<<>>for( int it = 0; it < nCaloJets; it++ )
	*/
	//fill


}//<<>>void KUCMSAodSkimmer::processCalojets()

void KUCMSAodSkimmer::processPhotons(){

    //bool verbose = true;
    bool verbose = false;

	// intilize
	selPhotons.clearBranches(); // <<<<<<<   must do

	// calc
    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

	std::vector<int> phoOrderIndx;
    uInt nSelPhotons = 0;
	//uint nEleVeto = 0;
    uInt nPhotons = Photon_excluded->size();	
    if( DEBUG || verbose ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( uInt it = 0; it < nPhotons; it++ ){

		//---------------------------------------------------
    ///////////  pho selection ////////////////////////////////////////////////////////////////////
        if( DEBUG ) std::cout << " -- looping photons : getting pho isEB, has min pt, not excluded, electron veto " << std::endl;
		auto isExcluded = (*Photon_excluded)[it];
		auto scIndx = (*Photon_scIndex)[it];
        auto isEB = (*SuperCluster_seedIsEB)[scIndx];
		auto hasEleVeto = (*Photon_electronVeto)[it];
        bool hasPixSeed = (*Photon_pixelSeed)[it];

		//if( geVars("genSigPerfect") == 1 ) std::cout << " -- pho sel susid " << (*Gen_susId)[(*Photon_genIdx)[it]] << std::endl;
		bool isGenSig = doGenInfo ? ( (*Gen_susId)[(*Photon_genIdx)[it]] == 22 )  : 0;

        auto pt = (*Photon_pt)[it];
		bool underMinPt = pt < 30;
        auto eta = (*Photon_eta)[it];
		auto overMaxEta = std::abs(eta) > 1.479;
        auto phi = (*Photon_phi)[it];

        if( DEBUG ) std::cout << " -- looping photons : getting pho iso " << std::endl;
        auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
        bool passHcalSum = true;
        auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
        bool passTrkSum = tspscdr4 < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
        auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
        bool passsEcalRhSum = erhsecdr4 < 10.0;
        auto htoem = (*Photon_hadTowOverEM)[it];
        bool passHOE = htoem < 0.02;
		bool isoPho = passHOE && passsEcalRhSum && passTrkSum && passHcalSum;
        bool failPhoIso = not isoPho;

/*      // moved to after photon selection
        if( DEBUG ) std::cout << " -- looping photons : getting phojet iso " << std::endl;
        bool isJetPhoton = false;
		phoJetVeto.clear();
		for( int jit = 0; jit < Jet_energy->size(); jit++ ){ phoJetVeto.push_back(false); }
        for( int jit = 0; jit < Jet_energy->size(); jit++ ){

            bool underMinJPt = (*Jet_pt)[jit] < 30.0;
            bool underMinJQual = getJetQuality(it)  < 2;
			auto jeta = (*Jet_eta)[jit];
            auto jphi = (*Jet_phi)[jit];
			auto overMaxJEta = std::abs(jeta) > 2.4;
            if( underMinJPt || underMinJQual || overMaxJEta ) continue;
			
            float dpjeta = jeta - eta;
            float dpjphi = dPhi( jphi, phi );
            float dr = hypo( dpjeta, dpjphi );
			bool minDr = dr < 0.4;
            if( minDr ) isJetPhoton = true; 
			if( minDr && isoPho ) phoJetVeto[jit] = true;

        } // for( int jit = 0; jit < nSelJets; jit++ )
*/

		//  change to skip jets and keep all photons regardless of photon iso with jet
		bool phoskip = isExcluded || hasPixSeed || overMaxEta || underMinPt || failPhoIso;
		//if( geVars("genSigPerfect") == 1 && isGenSig ) std::cout << " -- pho sel: phoskip " << phoskip <<  " isGenSig " << isGenSig << std::endl;  
		//if( geVars("genSigPerfect") == 1 &&  phoskip && isGenSig  ){ 
		//		std::cout << "   -- xsepji: " << isExcluded  << hasPixSeed << overMaxEta << underMinPt << isJetPhoton << failPhoIso << std::endl; }
		if( ( geVars("genSigPerfect") != 1 ) && phoskip ) continue;		
		if( ( geVars("genSigPerfect") == 1 ) &&  ( not isGenSig ) ) continue;

    ///////////  pho selection ////////////////////////////////////////////////////////////////////

        if( DEBUG ) std::cout << " -- pho pull SC info" << std::endl;
        auto scSize = SuperCluster_seedIsEB->size();
        auto rhids = (*SuperCluster_rhIds)[scIndx];
        uInt nrh = rhids.size();

        if( DEBUG ) std::cout << " -- pho pull info" << std::endl;
		auto isOOT = (*Photon_isOot)[it];
		auto time = (*Photon_seedTOFTime)[it];
		auto smaj = (*SuperCluster_smaj)[scIndx];
        auto smin = (*SuperCluster_smin)[scIndx];
        auto r9 = (*Photon_r9)[it];
        auto sieie = 0.f;//(*Photon_sieie)[it];
        auto energy = (*Photon_energy)[it];
        auto cee = (*SuperCluster_covEtaEta)[scIndx];
        auto cep = (*SuperCluster_covEtaPhi)[scIndx];
        auto cpp = (*SuperCluster_covPhiPhi)[scIndx];
        auto s4 = (*Photon_s4)[it];
        auto salp = (*SuperCluster_salp)[scIndx];
        auto etawide = (*SuperCluster_etaWidth)[scIndx];
        auto phiwide = (*SuperCluster_phiWidth)[scIndx];
        auto sigmaieie = (*Photon_SigmaIEtaIEta)[it];   //!
        auto scx = (*SuperCluster_x_calo)[scIndx];
        auto scy = (*SuperCluster_y_calo)[scIndx];
        auto scz = (*SuperCluster_z_calo)[scIndx];

        int genIdx = -1;
        int momIdx = -1;
        float susId = -1;
		float genpt = -10;
        float momEnergy = -10;   //!
        float momEta = -4;   //!
        float momMass = -10;   //!
        float momPhi = -4;   //!
        float momPt = -10;   //!
        float momPx = -10;   //!
        float momPy = -10;   //!
        float momPz = -10;   //!
        float momVx = -1;   //!
        float momVy = -1;   //!
        float momVz = -100;   //!

		float ctau = -10;

		if( doGenInfo ){
			genIdx = (*Photon_genIdx)[it];
			momIdx = (*Gen_motherIdx)[genIdx];
			//momIdx = (*Photon_genSigXMomId)[it];
        	susId = (*Gen_susId)[genIdx];

			if( momIdx > -1.0 ){
                momEnergy = (*Gen_energy)[momIdx];   //!
                momEta = (*Gen_eta)[momIdx];   //!
                momMass = (*Gen_mass)[momIdx];   //!
                momPhi = (*Gen_phi)[momIdx];   //!
                momPt = (*Gen_pt)[momIdx];   //!
                momPx = (*Gen_px)[momIdx];   //!
                momPy = (*Gen_py)[momIdx];   //!
                momPz = (*Gen_pz)[momIdx];   //!
                momVx = (*Gen_vx)[momIdx];   //!
                momVy = (*Gen_vy)[momIdx];   //!
                momVz = (*Gen_vz)[momIdx];   //!
			}//<<>>if( momIdx > -1.0 )

		}//if( doGenInfo )

		float rtime = time;
        float rt2 = sq2(time);
		float re = energy;
        float s1 = 0.99736-0.0035*rtime-0.00319*rt2;
        float s2 = 1.205-0.079*rtime+0.00093*rt2;//<100, <3 // 
        float s3 = std::exp( -1*sq2(rtime-0.2037)/144.942338);
        //float core = ( rtime < -0.549 ) ? 1 : ( rtime < 3.37 ) ? s1 : s2;
        float core = ( rtime < -0.549 ) ? 1 : s1;
        if( ( re < 325 ) && ( rtime > 7.0 ) ) core = s2;
        if( ( re < 200 ) && ( rtime > 6.0 ) ) core = s2;
        if( ( re < 150 ) && ( rtime > 5.0 ) ) core = s2;
        if( ( re < 100 ) && ( rtime > 4.0 ) ) core = s2;
        if( rtime > 16.0 ) core = s3;
        ////float core = ( rtime < 0.204 ) ? 1 : std::exp(-1*sq2(rtime-0.204)/144.942338);
        float ce = energy/core;
        float cpt = pt/core;		

        float cms000 = hypo(scx,scy,scz);
        float calcor = cms000/SOL;
        float tofPVtoRH = hypo(scx-PV_x,scy-PV_y,scz-PV_z);
        float pvtof = tofPVtoRH/SOL;
        float crtime = rtime + calcor - pvtof;
        float tmeasured = rtime*SOL+tofPVtoRH;
        float m1 = std::sqrt(sq2(tofPVtoRH)+8*sq2(tmeasured));
        float m2 = (tofPVtoRH+m1)*(tofPVtoRH/(4*sq2(tmeasured)));
        float m3 = tofPVtoRH/tmeasured;
        float m2_phys = ( m2 > 1 ) ? 1 : m2;
        float m3_phys = ( m3 > 1 ) ? 1 : m3;
        float MBetaEqual = 2*ce*std::sqrt(1-sq2(m2_phys));
        //float MBetaEqual = 2*ce*m2_phys;
        float MBetaPrompt = 2*ce*std::sqrt(1-sq2(m3_phys));
        //float MBetaPrompt = 2*ce*m3_phys;

        auto htsebcdr4 = (*Photon_hcalTowerSumEtBcConeDR04)[it];
        auto tsphcdr3 = (*Photon_trkSumPtHollowConeDR03)[it];
        auto tsphcdr4 = (*Photon_trkSumPtHollowConeDR04)[it];
        //auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
        //auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
        //auto htoem = (*Photon_hadTowOverEM)[it];
        auto hoem = (*Photon_hadOverEM)[it];
        auto sieip = 0.f;//(*Photon_sieip)[it];
        auto sipip = 0.f;//(*Photon_sipip)[it];
        auto phoOOT = (*Photon_isOot)[it];
        //auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
        auto nthcdr4 = (*Photon_nTrkHollowConeDR04)[it];   //!
        auto ntscdr4 = (*Photon_nTrkSolidConeDR04)[it];   //!
        auto epfciso = (*Photon_ecalPFClusterIso)[it];   //!
        auto hct = (*hasConversionTracks)[it];   //!
        auto hpfciso = (*Photon_hcalPFClusterIso)[it];   //!
        auto hoepuc = (*Photon_Hoe_PUcorr)[it];   //!
        auto pfciso = (*Photon_pfChargedIso)[it];   //!
        auto pscisopfpv = (*Photon_pfChargedIsoPFPV)[it];   //!
        auto pfcisowv = (*Photon_pfChargedIsoWorstVtx)[it];   //!
        auto pfphoisso = (*Photon_pfPhoIso03)[it];   //!
        auto pfriso3aq = (*pfRelIso03_all_quadratic)[it];   //!
        auto pfrtso3cq = (*pfRelIso03_chg_quadratic)[it];   //!

        //--------------------------------------------------------------
        if( DEBUG ) std::cout << " -- pho get calclated values" << std::endl;
        int phoQuality = getPhoQuality(it);
        auto phoClstrR9 = clstrR9( rhids );
        //auto phoEigens2D = getRhGrpEigenFromAngles( rhids );
        //auto evaluegeo = phoEigens2D[2];
		//auto geosmaj = phoEigens2D[3];
        //auto geosmin = phoEigens2D[4];

        float phoPhoIsoDr = 10.0;
        for( uInt it2 = 0; it2 < nPhotons; it2++ ){
            if( it == it2 ) continue;
        	float eta2 = (*Photon_eta)[it2];
        	float phi2 = (*Photon_phi)[it2];
			auto dphi12 = dPhi( phi, phi2 );
			auto pho2dr = hypo( eta-eta2, dphi12 ); 
            //std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << pho2dr << std::endl;
			if( pho2dr < phoPhoIsoDr ) phoPhoIsoDr = pho2dr;
		}//for( uInt it2 = it+1; it2 < nPhotons; it2++ )

/*
        // pho object selection ------------------------------------------
        if( DEBUG ) std::cout << " -- pho obj selection" << std::endl;
		auto isMinMedQuality = phoQuality > 1;
		auto underMaxSMaj = smaj == 1.3;
        auto underMaxSMin = smin <= 0.4;
		auto overMinR9 = r9 >= 0.9;
		auto underMaxSieie = sieie <= 0.014;
        auto overMinRhCnt = nrh >= 20;

		//auto phoSelected = isMinMedQuality;
        auto phoSelected = true;
		if( not phoSelected ) continue;
*/
		//float ordpt = pt;
		//while( phoOrderIndx.count(ordpt) > 0 ){ ordpt -= 0.000001; }
		//if( DEBUG ) std::cout << " -- setting pho index : " << it << " for pt : " << ordpt << std::endl;
		phoOrderIndx.push_back(it);
        //if( DEBUG ) std::cout << " -- pho index set : " << phoOrderIndx[ordpt] << std::endl;

		// fill ( vectors )
		if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
        selPhotons.fillBranch( "selPhoOOT", phoOOT );
        selPhotons.fillBranch( "selPhoEleVeto", hasEleVeto );
		selPhotons.fillBranch( "selPhoQuality", phoQuality );
        selPhotons.fillBranch( "selPhoTime", time );
        selPhotons.fillBranch( "selPhoSCx", scx );
        selPhotons.fillBranch( "selPhoSCy", scy );
        selPhotons.fillBranch( "selPhoSCz", scz );

    	selPhotons.fillBranch( "selPhoCorEnergy", ce );
    	selPhotons.fillBranch( "selPhoCorPt", cpt );
    	selPhotons.fillBranch( "selPhoMBetaEql", MBetaEqual );
	    selPhotons.fillBranch( "selPhoMBetaPmt", MBetaPrompt );
        selPhotons.fillBranch( "selPhoPVtoRH", tofPVtoRH );

		//selPhotons.fillBranch( "selPhoGeoEgnVal", evaluegeo );
        selPhotons.fillBranch( "selPhoEta", eta );
        selPhotons.fillBranch( "selPhoPhi", phi );
        selPhotons.fillBranch( "selPhoPt", pt );
        selPhotons.fillBranch( "selPhoGenPt", genpt );
        selPhotons.fillBranch( "selPhoSMaj", smaj );
        selPhotons.fillBranch( "selPhoSMin", smin );
        selPhotons.fillBranch( "selPhoClstrRn", phoClstrR9 );
        selPhotons.fillBranch( "selPhoR9", r9 );
        selPhotons.fillBranch( "selPhoSieie", sieie );
        selPhotons.fillBranch( "selPhoSieip", sieip );
        selPhotons.fillBranch( "selPhoSipip", sipip );
        selPhotons.fillBranch( "selPhoCovEtaEta", cee );
        selPhotons.fillBranch( "selPhoCovEtaPhi", cep );
        selPhotons.fillBranch( "selPhoCovPhiPhi", cpp );
        selPhotons.fillBranch( "selPhoS4", s4 );
        selPhotons.fillBranch( "selPhoSAlp", salp );
        selPhotons.fillBranch( "selPhoEtaWidth", etawide );
        selPhotons.fillBranch( "selPhoPhiWidth", phiwide );
        selPhotons.fillBranch( "selPhoNrh", nrh );
        selPhotons.fillBranch( "selPhoEnergy" , energy );
        selPhotons.fillBranch( "selPhoSusyId" , susId );
        selPhotons.fillBranch( "selPhoSigmaIEtaIEta" , sigmaieie );
        //selPhotons.fillBranch( "selPhoGeoSMaj", geosmaj );
        //selPhotons.fillBranch( "selPhoGeoSMin", geosmin );
        selPhotons.fillBranch( "selPhoHcalTowerSumEtBcConeDR04", htsebcdr4 );
        selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR03", tsphcdr3 );
        selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR04", tsphcdr4 );
        selPhotons.fillBranch( "selPhoTrkSumPtSolidConeDR04", tspscdr4 );
        selPhotons.fillBranch( "selPhoPixelSeed", hasPixSeed );
        selPhotons.fillBranch( "selPhoEcalRHSumEtConeDR04", erhsecdr4 );
        selPhotons.fillBranch( "selPhoHadTowOverEM", htoem );
        selPhotons.fillBranch( "selPhoHadOverEM", hoem );
        selPhotons.fillBranch( "selPhoHcalTowerSumEtConeDR04", htsecdr4  );   //!
        selPhotons.fillBranch( "selPhoNTrkHollowConeDR04", nthcdr4  );   //!
        selPhotons.fillBranch( "selPhoNTrkSolidConeDR04", ntscdr4  );   //!
        selPhotons.fillBranch( "selPhoEcalPFClusterIso", epfciso  );   //!
        selPhotons.fillBranch( "selPhoHasConversionTracks", hct  );   //!
        selPhotons.fillBranch( "selPhoHcalPFClusterIso", hpfciso  );   //!
        selPhotons.fillBranch( "selPhoHoe_PUcorr", hoepuc  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIso", pfciso  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIsoPFPV", pscisopfpv  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIsoWorstVtx", pfcisowv  );   //!
        selPhotons.fillBranch( "selPhoPfPhoIso03", pfphoisso  );   //!
        selPhotons.fillBranch( "selPhoPfRelIso03_all_quadratic", pfriso3aq  );   //!
        selPhotons.fillBranch( "selPhoPfRelIso03_chg_quadratic", pfrtso3cq  );   //!

        selPhotons.fillBranch( "selPhoPhoIsoDr", phoPhoIsoDr );
        selPhotons.fillBranch( "selPhoGenIdx", genIdx );
        selPhotons.fillBranch( "selPhoGenSigMomEnergy", momEnergy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomEta", momEta  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomMass", momMass  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPhi", momPhi  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPt", momPt  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPx", momPx  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPy", momPy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPz", momPz  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVx", momVx  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVy", momVy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVz", momVz  );   //!
		//if( verbose ) std::cout << " -- selPho Pt: " << pt << " phi: " << phi << " geo: " << evaluegeo << " clrn: " << phoClstrR9;
		if( verbose ) std::cout << " nrh: " << nrh << " quality: " << phoQuality << std::endl;

        nSelPhotons++;
    }//<<>>for( int it = 0; it < nPhotons; it++ )
    if( DEBUG ) std::cout << " -- pho loop finished" << std::endl;
	// fill ( other )

/*
	std::vector<unsigned int> phoPtOrder;
	if( phoOrderIndx.size() > 0 ){ 
		for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ ){ 
			phoPtOrder.push_back(phoptit->second); 
			if( DEBUG || verbose ) std::cout << " phoptit pushback : " << phoptit->second << std::endl;
		}//for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ )
	}//if( phoOrderIndx.size() > 0 )
*/

	int loopEnd = ( nSelPhotons > 0 ) ? ( nSelPhotons < 2  ) ? 1 : 2 : 0;
	phoJetVeto.clear();
	for( int jit = 0; jit < Jet_energy->size(); jit++ ){ phoJetVeto.push_back(false); } 
	for( int pit = 0; pit < loopEnd; pit++ ){ 

		float eta = (*Photon_eta)[phoOrderIndx[pit]];
		float phi = (*Photon_phi)[phoOrderIndx[pit]];

        if( DEBUG ) std::cout << " -- looping photons : getting phojet iso " << std::endl;
        //bool isJetPhoton = false;
        for( int jit = 0; jit < Jet_energy->size(); jit++ ){

            bool underMinJPt = (*Jet_pt)[jit] < 30.0;
            bool underMinJQual = getJetQuality(jit)  < 2;
            auto jeta = (*Jet_eta)[jit];
            auto jphi = (*Jet_phi)[jit];
            auto overMaxJEta = std::abs(jeta) > 2.4;
            if( underMinJPt || underMinJQual || overMaxJEta ) continue;

            float dpjeta = jeta - eta;
            float dpjphi = dPhi( jphi, phi );
            float dr = hypo( dpjeta, dpjphi );
            bool minDr = dr < 0.4;
            //if( minDr ) isJetPhoton = true;
            if( minDr ) phoJetVeto[jit] = true;

        } // for( int jit = 0; jit < nSelJets; jit++ )

	}//<<>>for( int pit = 0; pit < loopEnd; pit++ )

    if( DEBUG || verbose ) std::cout << " - Selected " << nSelPhotons << " photons" << std::endl;
    geCnts.set( "nSelPhotons", nSelPhotons );	
	selPhotons.fillBranch( "nSelPhotons",  nSelPhotons );
    selPhotons.fillBranch( "nPhotons", nPhotons );	
    //selPhotons.fillBranch( "nEleVeto", nEleVeto );
    //selPhotons.fillBranch( "selPhoPtOrder", phoPtOrder );

    //if( DEBUG || verbose ) std::cout << " - Setting Leading/Subleading Photon with " << phoPtOrder.size() << std::endl;
	uInt leadPhoIdx = ( nSelPhotons >= 1 ) ? phoOrderIndx[0]  : 9999;
	uInt subLeadPhoIdx = ( nSelPhotons >= 2 ) ? phoOrderIndx[1] : 9999;
    geCnts.set("leadPho",leadPhoIdx);
    geCnts.set("subLeadPho",subLeadPhoIdx);
	//selPhotons.fillBranch( "leadSelPho", leadPho );
	//selPhotons.fillBranch( "subLeadSelPho", subLeadPho );
	//if( verbose ) std::cout << " -- pho lead & sublead selected : " << leadPho << " - " << subLeadPho << std::endl;
	if( DEBUG ) std::cout << " -- pho lead & sublead idx selected : " << leadPhoIdx << " - " << subLeadPhoIdx << std::endl;	
    float lPhoPt = ( nSelPhotons > 0 ) ? (*Photon_pt)[leadPhoIdx] : 0.f;
    geVars.set( "leadPhoPt", lPhoPt );
    float slPhoPt = ( nSelPhotons > 1 ) ? (*Photon_pt)[subLeadPhoIdx] : 0.f;
    geVars.set( "subLeadPhoPt", slPhoPt);
    float lPhoPhi = ( nSelPhotons > 0 ) ? (*Photon_phi)[leadPhoIdx] : 0.f;
    geVars.set( "leadPhoPhi", lPhoPhi );
    float slPhoPhi = ( nSelPhotons > 1 ) ? (*Photon_phi)[subLeadPhoIdx] : 0.f;
    geVars.set( "subLeadPhoPhi", slPhoPhi );
    float lPhoEta = ( nSelPhotons > 0 ) ? (*Photon_eta)[leadPhoIdx] : 0.f;
    geVars.set( "leadPhoEta", lPhoEta );
    float slPhoEta = ( nSelPhotons > 1 ) ? (*Photon_eta)[subLeadPhoIdx] : 0.f;
    geVars.set( "subLeadPhoEta", slPhoEta );

}//<<>>void KUCMSAodSkimmer::processPhoton(){

void KUCMSAodSkimmer::processElectrons(){

    //if( DEBUG ) std::cout << "Finding electrons" << std::endl;
	//-------- electrons --------------------------------------

    // intilize
    selElectrons.clearBranches(); // <<<<<<<   must do

	uInt nElectrons = Electron_energy->size();
	uInt nLooseEle = 0;
    uInt nSelIsoElectrons = 0;
	float nEpDrMatch = 0;
    float nEpSeedIdMatch = 0;
    uInt nEpsidsolo = 0;
    float elePhoIsoMinDr(10.0); 
    //if( DEBUG ) std::cout << " -- Looping electrons: " << std::endl;
    for( uInt itr = 0; itr < nElectrons; itr++ ){

		if( ( doGenInfo && ( (*Electron_genIdx)[itr] > -1 ) ) || not doGenInfo ){
			hist1d[10]->Fill((*Electron_pt)[itr]);
			if( (*Electron_isLVeto)[itr] ) hist1d[11]->Fill((*Electron_pt)[itr]);
        	if( (*Electron_isVeto)[itr] ) hist1d[12]->Fill((*Electron_pt)[itr]);
        	if( (*Electron_isLoose)[itr] ) hist1d[13]->Fill((*Electron_pt)[itr]);
        	if( (*Electron_isMedium)[itr] ) hist1d[14]->Fill((*Electron_pt)[itr]);
		}//<<>>if( (*Electron_hasGenMatch)[itr] )

    	if( DEBUG ) std::cout << " ---- Processing  electron: " << itr << " of: " << nElectrons << std::endl;
		float eta = (*Electron_eta)[itr];
    	float phi = (*Electron_phi)[itr];
		auto scIdx = (*Electron_scIndex)[itr];
		if( DEBUG ) std::cout << " -- doing sc sstuff " << std::endl;
		auto esid = getLeadRhID( (*SuperCluster_rhIds)[scIdx] );
		//std::cout << " -- esid : " << esid << std::endl;
		auto eletime = (*Electron_seedTOFTime)[itr];
		if( std::abs(eletime) > 30 ) continue;

    	uInt nPhotons = Photon_excluded->size();
		float elePhoIsoDr(10.0);
		bool epSeedIdMatch = false;
		bool epsidsolo = true;
		//if( DEBUG ) std::cout << " -- doing ele pho iso sstuff " << std::endl;
    	for( uInt pit2 = 0; pit2 < nPhotons; pit2++ ){
        	float eta2 = (*Photon_eta)[pit2];
        	float phi2 = (*Photon_phi)[pit2];
			auto scIdx2 = (*Photon_scIndex)[pit2];
			auto dphi12 = dPhi( phi, phi2 );
        	float elephodr = hypo( eta-eta2, dphi12 );
        	//std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << elephodr << std::endl;
        	if( elephodr < elePhoIsoDr ) elePhoIsoDr = elephodr;
            //std::cout << " -- scIdx : " << scIdx2 << std::endl; 
			auto rhIds = (*SuperCluster_rhIds)[scIdx2];
			//std::cout << " -- scRh1 : " << rhIds[0] << std::endl;
			auto psid = getLeadRhID(rhIds);
			//std::cout << " -- psid : " << psid << std::endl;
			if( esid == psid ){ if( epSeedIdMatch == true ) epsidsolo = false ; epSeedIdMatch = true; }
    	}//for( uInt pit2 = 0; pit2 < nPhotons; pit2++ )
		//std::cout << " ---- next electrons: " << std::endl;
		elePhoIsoMinDr = elePhoIsoDr;
		bool epDrMatch = ( elePhoIsoDr < 0.2  ) ? true : false;
		bool eleminpt = (*Electron_pt)[itr] > 10;
		if( not epDrMatch & eleminpt ) nSelIsoElectrons++;

		if( epDrMatch ) nEpDrMatch++;
		if( epSeedIdMatch ) nEpSeedIdMatch++;
		if( epsidsolo ) nEpsidsolo++;
		if( DEBUG ) std::cout << " ---- next electrons: " << std::endl;
        //std::cout << " ---- filling electrons: " << std::endl;
		selElectrons.fillBranch( "elePhoMinDr", elePhoIsoMinDr );
		selElectrons.fillBranch( "eleHasSVMatch", (*Electron_hasSVMatch)[itr] );
		selElectrons.fillBranch( "eleIsLoose", (*Electron_isLoose)[itr] );
		if( (*Electron_isLoose)[itr] ) nLooseEle++;

	}//<<>>for( int itr = 0; itr < nElectrons; itr++ )
    if( DEBUG ) std::cout << " -- Finishd looping electrons " << std::endl;

	geCnts.set("nSelEle",nSelIsoElectrons);// passes min pt, isolated from photons, && good time criteria

    selElectrons.fillBranch( "nElectrons", nElectrons );
	if( nElectrons == 0 ) nElectrons = 1;
    //std::cout << " -- Finish nEle looping electrons " << std::endl;
    selElectrons.fillBranch( "epDrMatch", nEpDrMatch/float(nElectrons) );
    selElectrons.fillBranch( "epSeedIdMatch", nEpSeedIdMatch/float(nElectrons) );
    selElectrons.fillBranch( "epSeedIdUnique", nEpsidsolo/float(nElectrons) );
    //std::cout << " -- Finish New looping electrons " << std::endl;
    selElectrons.fillBranch( "nSelElectrons", Electron_nSelElectrons );
    selElectrons.fillBranch( "nEleSVMatched", Electron_nSVMatched );
    selElectrons.fillBranch( "nSelIsoElectrons", nLooseEle );

}//<<>>void KUCMSAodSkimmer::processElectrons

void KUCMSAodSkimmer::processMuons(){

    if( DEBUG ) std::cout << "Finding muons" << std::endl;
    //-------- muons --------------------------------------
    
    selMuons.clearBranches(); // <<<<<<<   must do

    uInt nMuons = Muon_energy->size();
	uInt nLooseMuons = 0;
    if( DEBUG ) std::cout << " -- Looping muons: " << std::endl;
    for( uInt itr = 0; itr < nMuons; itr++ ){
	
		//float energy = (*Muons_energy)[itr];
    	selMuons.fillBranch( "muonIsLoose", (*Muon_isLoose)[itr] );
		selMuons.fillBranch( "muonIsMedium", (*Muon_isMedium)[itr] );

		if( (*Muon_isLoose)[itr] ) nLooseMuons++;

	}//<<>>for( uInt itr = 0; itr < nMuons; itr++ )
    selMuons.fillBranch( "nMuons", nMuons );
    selMuons.fillBranch( "nLooseMuons", nLooseMuons );
	selMuons.fillBranch( "nMuonSVMatched", Muon_nSVMatched );

}//<<>>void KUCMSAodSkimmer::processMuons

void KUCMSAodSkimmer::processJets(){

	// intilize
	selJets.clearBranches(); // <<<<<<<   must do

	// calc
    //if( DEBUG ) std::cout << "Finding jets" << std::endl;
	//--------- jets --------------------------

	uInt nSelJets = 0;
	uInt nJets = Jet_energy->size();
    std::vector<float> seljetpt;
    std::vector<float> seljeteta;
    std::vector<float> seljetphi;
    std::vector<float> seljetmass;
    if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets" << std::endl;
    for( uInt it = 0; it < nJets; it++ ){

		// pull values ---------------------------------------------------
		auto energy = (*Jet_energy)[it];
		auto mass = (*Jet_mass)[it];
        //if( DEBUG ) std::cout << " - Finding Jet Quality" << std::endl;
		auto quality = getJetQuality(it);
		auto pt = (*Jet_pt)[it];
		auto eta = (*Jet_eta)[it];
        auto phi = (*Jet_phi)[it];
		auto rhids = (*Jet_drRhIds)[it];
		if( DEBUG ) std::cout << " -- jettDrId size: " << rhids.size() << std::endl;		

        auto area = (*Jet_area)[it];
        auto chEmEF = (*Jet_chEmEF)[it];
        auto chHEF = (*Jet_chHEF)[it];
        auto chHM = (*Jet_chHM)[it];
        auto muEF = (*Jet_muEF)[it];
        auto neEmEF = (*Jet_neEmEF)[it];
        auto neHEF = (*Jet_neHEF)[it];
        auto neHM = (*Jet_neHM)[it];

        int qrkllpId = 0; //*Jet_genLlpId)[it];
        int gjllpId = 0; //*Jet_genLlpId)[it];
        float jgdpt = 0; //*Jet_genDptMatch)[it];
        float jgdr = 0; //*Jet_genDrMatch)[it];
        float jge = 0; //*Jet_genEnergy)[it];
        float jgeta = 0; //*Jet_genEta)[it];
        float jgimpang = 0; //*Jet_genImpactAngle)[it];
        float llpdp = 0; //*Jet_genLlpDp)[it];
        float llpdr = 0; //*Jet_genLlpDr)[it];
        float jgpt = 0; //*Jet_genPt)[it];
        float jgtof = 0; //*Jet_genTOF)[it];
        float jgt = 0; //(*Jet_genTime)[it];
        float jgllpt = 0; //*Jet_genTimeLLP)[it];
        if( doGenInfo ){
        	//qrkllpId = int((*Jet_genQrkLlpId)[it]);
            gjllpId = int((*Jet_genJetLlpId)[it]);
        	jgdpt = (*Jet_genDptMatch)[it];
        	jgdr = (*Jet_genDrMatch)[it];
        	jge = (*Jet_genEnergy)[it];
        	jgeta = (*Jet_genEta)[it];
        	jgimpang = (*Jet_genImpactAngle)[it];
        	//llpdp = (*Jet_genQrkLlpDp)[it];
        	//llpdr = (*Jet_genQrkLlpDr)[it];
        	jgpt = (*Jet_genPt)[it];
        	jgtof = (*Jet_genTOF)[it];
        	jgt = (*Jet_genTime)[it];
        	jgllpt = (*Jet_genTimeLLP)[it];
		}//if( doGenInfo )

		// get cacluated values-------------------------------------------
		if( DEBUG ) std::cout << " - Finding Jet Time." << std::endl;
		auto rhenergies = getRhGrpEnergies( rhids );
		auto rhtimes = getRhGrpTimes( rhids );
		auto timedist = getDistStats( rhtimes, rhenergies );
		auto time = timedist[6];

		if( DEBUG ) std::cout << " - Jet Obj selection." << std::endl;
		// jet object selection ------------------------------------------
		auto overMinPt = pt > 30; 
		auto underMaxEta = std::abs(eta) <= 2.4;
		auto isMinQuality = quality > 1; // 2 = "tight" 3 = "tighter"

		bool isNotPhoJet = not phoJetVeto[it];

		auto jetSelected = underMaxEta && isMinQuality && overMinPt && isNotPhoJet;
		if( not jetSelected ) continue;
		nSelJets++;

        seljetpt.push_back(pt);
        seljeteta.push_back(eta);
        seljetphi.push_back(phi);
        seljetmass.push_back(mass);

		// fill vectors
		selJets.fillBranch( "selJetGenLlpId", gjllpId );
        selJets.fillBranch( "selJetQrkLlpId", qrkllpId );
		selJets.fillBranch( "selJetQuality", quality );
        selJets.fillBranch( "selJetPt", pt);
        selJets.fillBranch( "selJetMass", mass);
        selJets.fillBranch( "selJetEnergy", energy);
        selJets.fillBranch( "selJetEta", eta);
        selJets.fillBranch( "selJetPhi", phi);
        selJets.fillBranch( "selJetTime", time);

        selJets.fillBranch( "selJetArea", area ); //*Jet_area)[it]; 
        selJets.fillBranch( "selJetChEmEF", chEmEF ); //*Jet_chEmEF)[it]; 
        selJets.fillBranch( "selJetchHEF", chHEF ); //*Jet_chHEF)[it]; 
        selJets.fillBranch( "selJetChHM", chHM ); //*Jet_chHM)[it]; 
        selJets.fillBranch( "selJetMuEF", muEF ); //*Jet_muEF)[it]; 
        selJets.fillBranch( "selJetNeEmEF", neEmEF ); //*Jet_neEmEF)[it]; 
        selJets.fillBranch( "selJetNeHEF", neHEF ); //*Jet_neHEF)[it]; 
        selJets.fillBranch( "selJetNeHM", neHM ); //*Jet_neHM)[it]; 

        selJets.fillBranch( "selGenJetDpt", jgdpt ); //*Jet_genDptMatch)[it]; 
        selJets.fillBranch( "selGenJetdr", jgdr ); //*Jet_genDrMatch)[it]; 
        selJets.fillBranch( "selGenJetEnergy", jge ); //*Jet_genEnergy)[it]; 
        selJets.fillBranch( "selGenJeteta", jgeta ); //*Jet_genEta)[it]; 
        selJets.fillBranch( "selGenJetImpAng", jgimpang ); //*Jet_genImpactAngle)[it]; 
        selJets.fillBranch( "selJetLlpDp", llpdp ); //*Jet_genLlpDp)[it]; 
        selJets.fillBranch( "selJetLlpDr", llpdr ); //*Jet_genLlpDr)[it]; 
        selJets.fillBranch( "selGenJetPt", jgpt ); //*Jet_genPt)[it]; 
        selJets.fillBranch( "selGenJetTof", jgtof ); //*Jet_genTOF)[it]; 
        selJets.fillBranch( "selGenJetTime", jgt ); // (*Jet_genTime)[it]; 
        selJets.fillBranch( "selGenJetLlpTime", jgllpt ); //*Jet_genTimeLLP)[it]; 

	}//<<>>for( int it = 0; it < nJets; it++ )
    if( DEBUG ) std::cout << " - Finished Jet loop." << std::endl;

	// fill other
	geVects.set( "selJetPt", seljetpt );
    geVects.set( "selJetEta", seljeteta );
    geVects.set( "selJetPhi", seljetphi );
    geVects.set( "selJetMass", seljetmass );
	geCnts.set( "nSelJets", nSelJets );
    selJets.fillBranch( "nJets", nJets);
	selJets.fillBranch( "nSelJets", nSelJets );

}//<<>>void KUCMSAodSkimmer::processJets()

//------------------------------------------------------------------------------------------------------------
// process RJR for event
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::processRJR( int type, bool newEvent ){

	//bool verbose = true;
    bool verbose = false;

	if( DEBUG || verbose ) std::cout << "Processing RJR event varibles." << std::endl;
	// ------------------  Loading &  Processing RJR event varibles  ---------------------------------------------------------

	// intilize out branches
	if( newEvent ) selRjrVars.clearBranches(); // <<<<<<<   must do

	// process event
	LAB->ClearEvent();

    auto leadPhoPt = geVars("leadPhoPt"); //selPhotons.getFLBranchValue( "selPhoPt", leadSelPho );
    auto subLeadPhoPt = geVars("subLeadPhoPt"); //selPhotons.getFLBranchValue( "selPhoPt", subLeadSelPho );
    auto leadPhoPhi = geVars("leadPhoPhi"); //selPhotons.getFLBranchValue( "selPhoPhi", leadSelPho );
    auto subLeadPhoPhi = geVars("subLeadPhoPhi"); //selPhotons.getFLBranchValue( "selPhoPhi", subLeadSelPho );
    auto leadPhoEta = geVars("leadPhoEta"); //selPhotons.getFLBranchValue( "selPhoPhi", leadSelPho );
    auto subLeadPhoEta = geVars("subLeadPhoEta"); //selPhotons.getFLBranchValue( "selPhoPhi", subLeadSelPho );

    auto leadSelPho = geCnts("leadPho"); //selPhotons.getUIBranchValue("leadSelPho");
    auto subLeadSelPho = geCnts("subLeadPho"); //selPhotons.getUIBranchValue("subLeadSelPho");
    auto nSelPhotons = geCnts("nSelPhotons");

    auto phoRMetCPx = geVars("metPx"); //selMet.getFLBranchValue("metPx");
    auto phoRMetCPy = geVars("metPy"); //selMet.getFLBranchValue("metPy");
	float unCorMet = hypo(phoRMetCPx,phoRMetCPy);

	int nSelPhos = 0;
	std::vector<RFKey> jetID;
	bool zsig = ( geVars("doSVs") == 1 ) ? true : false;	
	if( not zsig ){
	if( nSelPhotons != 0 ){
		nSelPhos = 1;
		//std::vector<RFKey> phoJetKey;
		if( DEBUG ) std::cout << " - Loading Lead/SubLead Pho" << std::endl;
		if( type == 0 ){
			phoRMetCPx += leadPhoPt*std::cos(leadPhoPhi); 
   			phoRMetCPy += leadPhoPt*std::sin(leadPhoPhi);
		}//<<>>if( type == 0 )
		else if( type == 1 ){
			TLorentzVector phojet;
			phojet.SetPtEtaPhiM( leadPhoPt, leadPhoEta, leadPhoPhi, 0 );
			//phoJetKey.push_back( COMB_J->AddLabFrameFourVector(phojet) );
        	jetID.push_back( COMB_J->AddLabFrameFourVector(phojet) );
    		//jetID.push_back(phoJetKey[0]);
		}//<<>>if( type == 1 )
		else {  std::cout << " !!!!!!!! Valid RJR Photon Processing Option Not Specified !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; }
	}//<<>>if( nSelPhotons > 0 )
	if( nSelPhotons > 1 ){ 
		nSelPhos = 2;
		if( type == 0 ){
        	phoRMetCPx += subLeadPhoPt*std::cos(subLeadPhoPhi);
        	phoRMetCPy += subLeadPhoPt*std::sin(subLeadPhoPhi);
		}//<<>>if( type == 0 )
		else if( type == 1 ){
			TLorentzVector sphojet;
        	sphojet.SetPtEtaPhiM( subLeadPhoPt, subLeadPhoEta, subLeadPhoPhi, 0 );
			//phoJetKey.push_back( COMB_J->AddLabFrameFourVector(sphojet) );
            jetID.push_back( COMB_J->AddLabFrameFourVector(sphojet) );
			//jetID.push_back(phoJetKey[1]);
		}//<<>>if( type == 1 )
	}//<<>>if( nSelPhotons > 1 )
	}//<<>>if( not zsig )

    if( DEBUG ) std::cout << " - Loading MET." << std::endl;
	if( verbose ){
		std::cout << " - Loading MET lPt: " << leadPhoPt << " lPhi: " << leadPhoPhi << std::endl;
    	std::cout << " - Loading MET slPt: " << subLeadPhoPt << " slPhi: " << subLeadPhoPhi << std::endl; 
		std::cout << " - Loading MET x: " << geVars("metPx") << " -> " << phoRMetCPx;
        std::cout << " y: " << geVars("metPy") << " -> " << phoRMetCPy << std::endl;
	}//<<>>if( verbose )
    TVector3 ETMiss(phoRMetCPx,phoRMetCPy,0);
	INV->SetLabFrameThreeVector(ETMiss);

    auto nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
    auto selJetPt = geVects( "selJetPt");
    auto selJetEta = geVects( "selJetEta");
    auto selJetPhi = geVects( "selJetPhi");
    auto selJetMass = geVects( "selJetMass");
    std::vector<RFKey> leadJetKey;
	if( DEBUG ) std::cout << " - Loading Jets." << std::endl;
	//if( type == 0 && nSelJets > 2*nSelPhos ) nSelJets == 2*nSelPhos;
  	for( uInt it = 0; it < nSelJets; it++ ){
		auto sjetPt = selJetPt[it]; //selJets.getFLBranchValue( "selJetPt", it );
		auto sjetEta = selJetEta[it]; //selJets.getFLBranchValue( "selJetEta", it );
		auto sjetPhi = selJetPhi[it]; //selJets.getFLBranchValue( "selJetPhi", it ); 
        auto sjetMass = ( selJetMass[it] > 0 ) ? selJetMass[it] : 0; //selJets.getFLBranchValue( "selJetMass", it );   
		TLorentzVector jet;
		jet.SetPtEtaPhiM( sjetPt, sjetEta, sjetPhi, sjetMass );
		if( verbose ) std::cout << " - Loading Jet Pt: " << sjetPt << " Eta: " << sjetEta;
		if( verbose ) std::cout << " Phi: " << sjetPhi << " M: " << sjetMass << std::endl;
		if( it == 0 ){ leadJetKey.push_back( COMB_J->AddLabFrameFourVector(jet) ); jetID.push_back(leadJetKey[0]); } 
		else jetID.push_back(COMB_J->AddLabFrameFourVector(jet)); 
	}//<<>>for( int i = 0; i < nSelJets; i++ )

	if( !LAB->AnalyzeEvent() ) std::cout << "Something went wrong with tree event analysis" << std::endl;
	
	if( DEBUG ) std::cout << " - Getting RJR varibles." << std::endl;
	// ---------  Finished Processing RJR varibles --------------------------------------------------------------------

	uInt firstObject = 0;
	if( ( type == 1 ) && ( nSelPhos > 0 ) ) firstObject = ( nSelPhos == 2 ) ? 2 : 1 ;
	int nJetsJa = 0;
	int nJetsJb = 0;
	bool isALeadPhoSide = true;	
	int subPhoLocation = 0;
	int nVisObjects = jetID.size();
	//if( type == 1 ){ if( COMB_J->GetFrame(jetID[it]) == *Jb ) isALeadPhoSide = false; } //else 
	//
	//
	//  !!!!!!!!!!!  fix counting for 1a 2a 2b 2a
	//
	//	!!!!!!!!     save 4 vect info for phi/jets in all frames
	//					then get mometum iwth rjr function for 4 vects
	//
	if( COMB_J->GetFrame(jetID[0]) == *J1b || COMB_J->GetFrame(jetID[0]) == *J2b ) isALeadPhoSide = false;
	//  -- redifine A & B side for jets to be : is jet on lead pho side -> A ; is not on lead pho side -> B 
    for( uInt it = firstObject; it < nVisObjects; it++ ){	

        if( COMB_J->GetFrame(jetID[it]) == *J1a || COMB_J->GetFrame(jetID[it]) == *J2a  ){ isALeadPhoSide ? nJetsJa++ : nJetsJb++; } // one for each frame 
		else { isALeadPhoSide ? nJetsJb++ : nJetsJa++; }
        //if( COMB_J->GetFrame(jetID[it]) == *Jb ){ isALeadPhoSide ? nJetsJb++ : nJetsJa++; } // one for each frame 

	}//<<>>for( int i = 0; i < nSelJets; i++ )
	float abDiffSide = isALeadPhoSide ? 1 : -1;
	if( ( type == 1 ) && ( nSelPhos == 2 ) ) subPhoLocation = ( COMB_J->GetFrame(jetID[1]) == *J1a || COMB_J->GetFrame(jetID[1]) == *J2a ) ? 1 : 2;

    selRjrVars.fillBranch( "rjrABSide", isALeadPhoSide );
    selRjrVars.fillBranch( "rjrNJetsJa", nJetsJa );
    selRjrVars.fillBranch( "rjrNJetsJb", nJetsJb );
    selRjrVars.fillBranch( "rjrNJets", int(geCnts("nSelJets")) );
    selRjrVars.fillBranch( "rjrNPhotons", nSelPhos );
    selRjrVars.fillBranch( "rjrNVisObjects", nVisObjects );
    selRjrVars.fillBranch( "rjrSubPhoLocation", subPhoLocation );

	if( type == 1 ){
    	for( uInt it = 0; it < nVisObjects; it++ ){

			bool onAside = true;
			if( COMB_J->GetFrame(jetID[it]) == *J1b || COMB_J->GetFrame(jetID[it]) == *J2b  ) onAside = false;
			if( it < nSelPhotons ) selRjrVars.fillBranch( "rjrVisPhoSide", onAside );
			else selRjrVars.fillBranch( "rjrVisJetSide", onAside );

		}//<<>>for( uInt it = 0; it < nVisObjects; it++ )
	}//<<>>if( type == 1 )

  	float m_MS = S->GetMass();
  	//float m_PS = S->GetMomentum(*CM);
  	float m_cosS  = S->GetCosDecayAngle();
  	float m_dphiS = S->GetDeltaPhiDecayAngle();
  	float m_dphiSI  = S->GetDeltaPhiBoostVisible();
  	float m_PTS = S->GetFourVector().Pt();
  	float m_PZS = S->GetFourVector().Pz();
  
  	float m_MX2a = X2a->GetMass();
  	float m_cosX2a = X2a->GetCosDecayAngle();
  	float m_MX2b = X2b->GetMass();
  	float m_cosX2b = X2b->GetCosDecayAngle();

    selRjrVars.fillBranch( "rjrX2aMass", m_MX2a );
    selRjrVars.fillBranch( "rjrX2aCosA", m_cosX2a );
    selRjrVars.fillBranch( "rjrX2bMass", m_MX2b );
    selRjrVars.fillBranch( "rjrX2bCosA", m_cosX2b );

    selRjrVars.fillBranch( "rjrSMass", m_MS );
    selRjrVars.fillBranch( "rjrSCosA", m_cosS );

    selRjrVars.fillBranch( "rjrSdphiDA", m_dphiS  );
    selRjrVars.fillBranch( "rjrSdphiBV", m_dphiSI );
    selRjrVars.fillBranch( "rjrPTS", m_PTS );
    selRjrVars.fillBranch( "rjrPZS", m_PZS );


	//float X2QSum = std::sqrt((sq2(m_MX2a)+sq2(m_MX2b))/2);
    //float X2NQSum = X2QSum/m_MS;
    //float X2Ave = (m_MX2a+m_MX2b)/2;
	//float X2NAve = X2Ave/m_MS;
    //float X2GMean = std::sqrt(m_MX2a*m_MX2b);
    //float X2NGMean = X2GMean/m_MS;
    //float X2Diff = abDiffSide*(m_MX2a-m_MX2b)/(m_MX2a+m_MX2b);

    //selRjrVars.fillBranch( "rjrX2NQSum", X2NQSum );
    //selRjrVars.fillBranch( "rjrX2NAve", X2NAve );
    //selRjrVars.fillBranch( "rjrX2NGMean", X2NGMean );
    //selRjrVars.fillBranch( "rjrX2QSum", X2QSum );
    //selRjrVars.fillBranch( "rjrX2Ave", X2Ave );
    //selRjrVars.fillBranch( "rjrX2GMean", X2GMean );
    //selRjrVars.fillBranch( "rjrX2Diff", X2Diff );

	//float m_PX1X2a = X1a->GetFourVector(*X2a).P();  // ?
    //float m_PX1X2b = X1b->GetFourVector(*X2b).P();  //?
    //float m_PJX2a = Ja->GetListVisibleFrames().GetFourVector(*X2a).P();
    //float m_PJX2b = Jb->GetListVisibleFrames().GetFourVector(*X2b).P();
	//float m_PX2Sa = X2a->GetFourVector(*S).P();
    //float m_PX2Sb = X2b->GetFourVector(*S).P();
	//float diff_PJX2a = m_PX1X2a - m_PJX2a;
    //float diff_PJX2b = m_PX1X2b - m_PJX2b;

	std::vector< TLorentzVector > p4;
	p4.push_back(Ja->GetFourVector(*S));
    p4.push_back(Jb->GetFourVector(*S));
    p4.push_back(X1a->GetFourVector(*S));
    p4.push_back(X1b->GetFourVector(*S));

	for( int i = 0; i < 4; i++ ){

		float px = p4[i].Px();
		float py = p4[i].Py(); 
		float pz = p4[i].Pz(); 
		float e = sqrt( px*px + py*py + pz*pz );
		p4[i].SetPxPyPzE(px,py,pz,e);	

	}//<<>>for( int i = 0; i < 4; i++ )

	float a_MX2a = ( p4[0] + p4[2] ).M();
	float a_MX2b = ( p4[1] + p4[3] ).M();
	float a_MS = ( p4[0] + p4[2] + p4[1] + p4[3] ).M();// -- !! Mr	

	//float a_MX2a = 2*m_PX1X2a;
    //float a_MX2b = 2*m_PX1X2b;
	//float a_MS = sqrt(sq2(m_PX2Sa)+sq2(a_MX2a)) + sqrt(sq2(m_PX2Sb)+sq2(a_MX2b));

    float AX2QSum = std::sqrt((sq2(a_MX2a)+sq2(a_MX2b))/2);
    float AX2NQSum = 2*AX2QSum/a_MS;//*2 -- !! R
    //float AX2Ave = (a_MX2a+a_MX2b)/2;
    //float AX2NAve = 2*AX2Ave/a_MS;//*2 sch tath 0-1
    //float AX2GMean = std::sqrt(a_MX2a*a_MX2b);
    //float AX2NGMean = 2*AX2GMean/a_MS;
    float AX2Diff = abDiffSide*(a_MX2a-a_MX2b)/(a_MX2a+a_MX2b);
	//float AX2DiffTest = (a_MX2a-a_MX2b)/(a_MX2a+a_MX2b);
	//std::cout << " - AX2Diff : " << abDiffSide << " " << AX2Diff << " " << AX2DiffTest << std::endl;

    float pf_ptX2a = X2a->GetMomentum(*S);
    float pf_ptX2b = X2b->GetMomentum(*S);

    //selRjrVars.fillBranch( "rjrPX1X2a", m_PX1X2a );
    //selRjrVars.fillBranch( "rjrPX1X2b", m_PX1X2b );
    //selRjrVars.fillBranch( "rjrPJX2a", m_PJX2a );
    //selRjrVars.fillBranch( "rjrPJX2b", m_PJX2b );
    //selRjrVars.fillBranch( "rjrPX2Sa", m_PX2Sa );
    //selRjrVars.fillBranch( "rjrPX2Sb", m_PX2Sb );
    //selRjrVars.fillBranch( "rjrDiffPJX2a", diff_PJX2a );
    //selRjrVars.fillBranch( "rjrDiffPJX2b", diff_PJX2b );

    //selRjrVars.fillBranch( "rjrAX2aMass", a_MX2a );
    //selRjrVars.fillBranch( "rjrAX2bMass", a_MX2b );
    //selRjrVars.fillBranch( "rjrASMass", a_MS );
    selRjrVars.fillBranch( "rjr_Mr", a_MS );

    //selRjrVars.fillBranch( "rjrAX2NQSum", AX2NQSum );
    selRjrVars.fillBranch( "rjr_R", AX2NQSum );
    //selRjrVars.fillBranch( "rjrAX2NAve", AX2NAve );
    //selRjrVars.fillBranch( "rjrAX2NGMean", AX2NGMean );
    //selRjrVars.fillBranch( "rjrAX2QSum", AX2QSum );
    //selRjrVars.fillBranch( "rjrAX2Ave", AX2Ave );
    //selRjrVars.fillBranch( "rjrAX2GMean", AX2GMean );
    selRjrVars.fillBranch( "rjrAX2Diff", AX2Diff );

    selRjrVars.fillBranch( "rjrX2aPtS", pf_ptX2a );
    selRjrVars.fillBranch( "rjrX2bPtS", pf_ptX2b );

  	float m_EVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).E();
  	float m_EVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).E();
  	float m_PVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).P();
  	float m_PVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).P();

  	float m_MX1a = X1a->GetMass();
  	float m_cosX1a = X1a->GetCosDecayAngle();
  	float m_MX1b = X1b->GetMass();
  	float m_cosX1b = X1b->GetCosDecayAngle();  
    float x1a_pt = X1a->GetFourVector().Pt();
    float x1a_eta = X1a->GetFourVector().Eta();
    float x1a_phi = X1a->GetFourVector().Phi();
    float x1b_pt = X1b->GetFourVector().Pt();
    float x1b_eta = X1b->GetFourVector().Eta();
    float x1b_phi = X1b->GetFourVector().Phi();

  	//float m_MV = S->GetListVisibleFrames().GetMass();
  	//float m_PV = S->GetListVisibleFrames().GetFourVector(*S).P();
  	//float m_MVa = X2a->GetListVisibleFrames().GetMass();// sub ja ect
    float m_MVa = Ja->GetMass();
  	//float m_MVb = X2b->GetListVisibleFrames().GetMass();
    float m_MVb = Jb->GetMass();

    //float m_MVDiff = abDiffSide*(m_MVa-m_MVb)/(m_MVa+m_MVb);
    float m_MVSum = std::sqrt((sq2(m_MVa)+sq2(m_MVb))/2);
	float m_MVNSum = 2*m_MVSum/a_MS;// -- !! Rv

    //selRjrVars.fillBranch( "rjrMVDiff", m_MVDiff );
    //selRjrVars.fillBranch( "rjrMVSum", m_MVSum );
    //selRjrVars.fillBranch( "rjrMVNSum", m_MVNSum );
    selRjrVars.fillBranch( "rjr_Rv", m_MVNSum );

  	float m_PV_lab    = S->GetListVisibleFrames().GetFourVector().P();
  	float m_dphiMET_V = S->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

    selRjrVars.fillBranch( "rjrX1aMass", m_MX1a );
    selRjrVars.fillBranch( "rjrX1aCosA", m_cosX1a );
    selRjrVars.fillBranch( "rjrX1bMass", m_MX1b );
    selRjrVars.fillBranch( "rjrX1bCosA", m_cosX1b );
    selRjrVars.fillBranch( "rjrX1a_Pt", x1a_pt );
    selRjrVars.fillBranch( "rjrX1a_Eta", x1a_eta );
    selRjrVars.fillBranch( "rjrX1a_Phi", x1a_phi );
    selRjrVars.fillBranch( "rjrX1b_Pt", x1b_pt );
    selRjrVars.fillBranch( "rjrX1b_Eta", x1b_eta );
    selRjrVars.fillBranch( "rjrX1b_Phi", x1b_phi );

    selRjrVars.fillBranch( "rjrEVa", m_EVa );
    selRjrVars.fillBranch( "rjrEVb", m_EVb );
    selRjrVars.fillBranch( "rjrPVa", m_PVa );
    selRjrVars.fillBranch( "rjrPVb", m_PVb );

	//selRjrVars.fillBranch( "rjrMV", m_MV );
    //selRjrVars.fillBranch( "rjrPV", m_PV );
    //selRjrVars.fillBranch( "rjrMVa", m_MVa );
    //selRjrVars.fillBranch( "rjrMVb", m_MVb );

    selRjrVars.fillBranch( "rjrPVlab", m_PV_lab );
    selRjrVars.fillBranch( "rjrDphiMETV", m_dphiMET_V );

	float rjrMet = hypo(phoRMetCPx,phoRMetCPy);

    selRjrVars.fillBranch( "rjrMET", rjrMet ); 

    TLorentzVector x1alab = X1a->GetFourVector(*LAB);
    TLorentzVector x1blab = X1b->GetFourVector(*LAB);
	float n2apx = x1alab.Px();
    float n2apy = x1alab.Py();
    float n2apz = x1alab.Pz();
	float n2bpx = x1blab.Px();
    float n2bpy = x1blab.Py();
    float n2bpz = x1blab.Pz();
	float n2ap = hypo( n2apx, n2apy, n2apz );
    float n2bp = hypo( n2bpx, n2bpy, n2bpz );
	bool aorb = true;
	//bool aorb = n2ap > n2bp;
	//bool aorb = isALeadPhoSide;
	float n2px = ( aorb ) ? n2apx : n2bpx;
    float n2py = ( aorb ) ? n2apy : n2bpy;
    float n2pz = ( aorb ) ? n2apz : n2bpz;
    selRjrVars.fillBranch( "rjrN2Px", n2px );
    selRjrVars.fillBranch( "rjrN2Py", n2py );
    selRjrVars.fillBranch( "rjrN2Pz", n2pz );

}//<<>>void KUCMSAodSkimmer::processRJR( int type, bool newEvent )

//------------------------------------------------------------------------------------------------------------
// decide which events to save to tree
//------------------------------------------------------------------------------------------------------------

bool KUCMSAodSkimmer::eventSelection(){
// select which events to save and fill output branches

	if( DEBUG ) std::cout << "Event selection." << std::endl;
	// determine if we want to save event
	bool zsig = ( geVars("doSVs") == 1 ) ? true : false;

    float evtMet = geVars("cmet");
	int nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
    float nSelPhotons = geCnts("nSelPhotons"); //selPhotons.getUIBranchValue("nSelPhotons");
	float leadSelPho = geCnts("leadPho"); //selPhotons.getUIBranchValue("leadSelPho");
    float subLeadSelPho = geCnts("subLeadPho"); //selPhotons.getUIBranchValue("subLeadSelPho");
	float leadPhoPt = ( nSelPhotons > 0 ) ? geVars("leadPhoPt") : 0;
    float subLeadPhoPt = ( nSelPhotons > 1 ) ? geVars("subLeadPhoPt") : 0;
	if( DEBUG ) std::cout << " - Lead/Sublead Photons: " << leadSelPho << " - " << subLeadSelPho << std::endl;

	bool hasLepSV = geVars("nSVLep") > 0;
    bool hasHadSV = geVars("nSVHad") > 0;
	bool hasSV = hasLepSV || hasHadSV;

	bool met150 = evtMet >= 150;
    bool gt1phos = nSelPhotons >= 1;
    bool gt2jets = nSelJets >= 2;
    bool gt2phos = nSelPhotons >= 2;
	bool leadPhoPt70 = leadPhoPt >= 70;
    bool leadPhoPt30 = leadPhoPt >= 30;
	bool subLeadPhoPt40 = subLeadPhoPt >= 40; 

    bool basesel = met150 && gt2jets;
    //bool svsel = basesel && hasSV;	
    bool phosel = basesel && ( ( gt1phos && leadPhoPt30 ) || hasSV );
	//auto evtSelected = leadPhoPt70 && subLeadPhoPt40 && gt2jets && gt2phos;

	bool evtSelected = zsig ? basesel : phosel;

    if( met150 ) cutflow["met150"]++;
    if( met150 && gt2jets ) cutflow["m_gt2jets"]++;
    if( met150 && gt2jets && gt1phos ){ cutflow["mj_gt1phos"]++; cutflow["sel_ppt"]++; }
    if( met150 && gt2jets && gt1phos && leadPhoPt30 ) cutflow["mjp_leadPhoPt30"]++;

	if( met150 && gt2jets && leadPhoPt30 ) cutflow["sel_p"]++;
    if( met150 && gt1phos && leadPhoPt30 ) cutflow["sel_j"]++;
    if( gt2jets && gt1phos && leadPhoPt30 ) cutflow["sel_m"]++;

	if( DEBUG ){ 
		if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl; 
		else std::cout << " --------------------- Event Failed." << std::endl;}
	//if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl;
	if( evtSelected ){ nSelectedEvents++; return true; } 
	else return false;

}//<<>>void KUCMSAodSkimmer::eventSelection()

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getPhoQuality( int it ){

    // determine pog quality class
    // -----------------------------------------------------
    if( DEBUG ) std::cout << " -- pho id" << std::endl;
    //bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= ( 0.006*(*Photon_pt)[it] + 4.2 );
    bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= 10.0;
    //bool hcalTowIso = (*Photon_hcalTowerSumEtBcConeDR04)[it] <= ( 0.0025*(*Photon_pt)[it] + 2.2 );
    //if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
    //bool hcTrkIsoL = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 3.5 ); //hallow cone track iso
    //bool hcTrkIsoT = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 2 ); //hallow cone track iso
    bool hcTrkIso = (*Photon_trkSumPtSolidConeDR04)[it] <= 6.0;
    //if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
    //bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.05;
    bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.02;
    //auto sieie = (*Photon_sieie)[it];
    //bool sigmaIeieEE = sieie <= 0.03; // tight only EE
    //bool sigmaIeieEB = sieie <= 0.013; // tight only EB

    //if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
    //auto baseCut = rhIso && hcalTowIso && hadOverE;
    //auto looseCut = baseCut && hcTrkIsoL;
    //auto tightCut = baseCut && hcTrkIsoT;
    //auto tightEB = tightCut && sigmaIeieEB;
    //auto tightEE = tightCut && sigmaIeieEE;
	bool isocut = rhIso && hcTrkIso && hadOverE;

	auto phoClass = isocut ? 1 : 0;
    //auto phoClass = tightCut ? 3 : looseCut ? 2 : baseCut ? 1 : 0;

	return phoClass;

}//<<>>int KUCMSAodSkimmer::getPhoQuality( int iter )

int KUCMSAodSkimmer::getJetQuality( int it ){

    const auto eta  = std::abs((*Jet_eta)[it]);	     
    const auto NHF  = (*Jet_neHEF)[it];
    const auto NEMF = (*Jet_neEmEF)[it];
    const auto CHF  = (*Jet_chHEF)[it];
    const auto CEMF = (*Jet_chEmEF)[it];
    const auto NHM  = (*Jet_neHM)[it];
    const auto CHM  = (*Jet_chHM)[it];
    const auto SHM  = NHM + CHM;
    const auto MUF  = (*Jet_muEF)[it];

    int tighter = 3;
    int tightLepVeto = 0;
    int tight = 2;
    int loose = 1;

    bool nhfup  = NHF  <= 0.90; 	// delpho
    bool nhflw  = NHF  >= 0.2;

    bool nemfup1 = NEMF <= 0.90; // delpho
    bool nemfup2 = NEMF <= 0.99;
    bool nemf80 = NEMF <= 0.80;
    bool nemflw = NEMF >= 0.01;
    bool nemf10 = NEMF >= 0.1;
	
    bool shm1  = SHM  >= 1;
    bool muf8  = MUF  <= 0.80;
    bool chf0  = CHF  >= 0;		// delpho
    bool chf10  = CHF  >= 0.10;
    bool chm0  = CHM  >= 0;		// delpho
    bool cemf8 = CEMF >= 0.80;	
    bool nhm2  = NHM  >= 1;
    bool nhm10 = NHM  >= 10;

    bool eta1 = eta <= 2.6;
    bool eta2 = eta <= 2.7;
    bool eta3 = eta <= 3.0;  

  	if (eta1){
    	if      (nhfup && nemfup1 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
    	else if (nhfup && nemf80 && shm1 && chf10 && chm0) return tighter;
    	else if (nhfup && nemfup1 && shm1 && chf0 && chm0) return tight;
    	else    return loose;
    } else if (!eta1 && eta2 ){ //<<>>if (eta1)
    	if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
    	else if (nhfup && nemf80 && chm0) return tighter;
    	else if (nhfup && nemfup2 && chm0) return tight;
    	else    return loose;
    } else if (!eta2 && eta3){ //<<>>if (eta1) : else if
    	if      (nemf10 && nemf80 && nhm2) return tighter;
    	else if (nemflw && nemfup2 && nhm2) return tight;
    	else    return loose;
    } else { //<<>>if (eta1) : else if : else if
    	if      (nhflw && nemfup1 && nhm10) return tight;
    	else    return loose;
    }//<<>>if (eta1) : else if : else if : else

    return -1; // should not happen

}//<<>>int KUCMSAodSkimmer::getJetQuality( int iter )

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree ){

	//fOutTree->Branch( "RunNumber", &RunNumber );
    selEvtVars.makeBranch( "dsKey", "DataSetKey", STR, "Key for source data set of event" );

    selEvtVars.makeBranch( "evtGenWgt", FLOAT );
    selEvtVars.makeBranch( "evtXSection", FLOAT );
    selEvtVars.makeBranch( "evtFillWgt", FLOAT );

    selEvtVars.makeBranch( "PVx", FLOAT );
    selEvtVars.makeBranch( "PVy", FLOAT );
    selEvtVars.makeBranch( "PVz", FLOAT );

    selEvtVars.makeBranch( "SV_nLeptonic", INT );
    selEvtVars.makeBranch( "SV_nHadronic", INT );

    selEvtVars.makeBranch( "Flag_BadChargedCandidateFilter", BOOL );
    selEvtVars.makeBranch( "Flag_BadPFMuonDzFilter", BOOL );
    selEvtVars.makeBranch( "Flag_BadPFMuonFilter", BOOL );
    selEvtVars.makeBranch( "Flag_EcalDeadCellTriggerPrimitiveFilter", BOOL );
    selEvtVars.makeBranch( "Flag_HBHENoiseFilter", BOOL );
    selEvtVars.makeBranch( "Flag_HBHENoiseIsoFilter", BOOL );
    selEvtVars.makeBranch( "Flag_ecalBadCalibFilter", BOOL );
    selEvtVars.makeBranch( "Flag_eeBadScFilter", BOOL );
    selEvtVars.makeBranch( "Flag_globalSuperTightHalo2016Filter", BOOL );
    selEvtVars.makeBranch( "Flag_goodVertices", BOOL );
    selEvtVars.makeBranch( "Flag_hfNoisyHitsFilter", BOOL );
	selEvtVars.makeBranch( "Flag_MetFilter", BOOL );

    selEvtVars.attachBranches( fOutTree );

	//selMet.makeBranch( "Met", FLOAT );
	selMet.makeBranch( "cmet", "selCMet", FLOAT, "Magnitude of event Met corrected for OOT photons" );
    selMet.makeBranch( "cmetPx", "selCMetPx", FLOAT, "Magnitude of event MetPx corrected for OOT photons" );
    selMet.makeBranch( "cmetPy", "selCMetPy", FLOAT, "Magnitude of event MetPy corrected for OOT photons" );
    selMet.makeBranch( "met", "selMet", FLOAT, "Magnitude of event Met" );
    selMet.makeBranch( "metPx", "selMetPx", FLOAT, "Magnitude of event MetPx" );
    selMet.makeBranch( "metPy", "selMetPy", FLOAT, "Magnitude of event MetPy" );

	selMet.attachBranches( fOutTree );

	selGenPart.makeBranch( "genPartEnergy", VFLOAT );
    selGenPart.makeBranch( "genPartEta", VFLOAT );
    selGenPart.makeBranch( "genPartPhi", VFLOAT );
    selGenPart.makeBranch( "genPartPt", VFLOAT );
    selGenPart.makeBranch( "genPartPdgId", VUINT );
    selGenPart.makeBranch( "genPartSusId", VINT );
    selGenPart.makeBranch( "genPartMomIdx", VINT );
    selGenPart.makeBranch( "genXMomCTau", VFLOAT );
    selGenPart.makeBranch( "genMomCTau", VFLOAT );
    selGenPart.makeBranch( "genCharge", VINT );   //!
    selGenPart.makeBranch( "genMass", VFLOAT );   //!   
    selGenPart.makeBranch( "genStatus", VINT );   //!
    selGenPart.makeBranch( "genVx", VFLOAT );   //!
    selGenPart.makeBranch( "genVy", VFLOAT );   //!
    selGenPart.makeBranch( "genVz", VFLOAT );   //!
    selGenPart.makeBranch( "genPx", VFLOAT );   //!
    selGenPart.makeBranch( "genPy", VFLOAT );   //!
    selGenPart.makeBranch( "genPz", VFLOAT );   //!
    selGenPart.makeBranch( "genSigType", INT );   //!
    selGenPart.makeBranch( "genSTFlagQQP", BOOL );   //!
    selGenPart.makeBranch( "genSTFlagQQZLL", BOOL );   //!
    selGenPart.makeBranch( "genSTFlagQQZL", BOOL );   //!
    selGenPart.makeBranch( "genSTFlagQQZ", BOOL );   //!
    selGenPart.makeBranch( "genSTFlagQP", BOOL );   //!
    selGenPart.makeBranch( "genSQMass", VFLOAT );   //! 
    selGenPart.makeBranch( "genSGMass", VFLOAT );   //! 
    selGenPart.makeBranch( "genLLPMass", VFLOAT );   //! 
    selGenPart.makeBranch( "genLSPMass", VFLOAT );   //! 
    selGenPart.makeBranch( "genGrvtinoMass", VFLOAT );   //! 

	selGenPart.attachBranches( fOutTree );

	selElectrons.makeBranch( "nElectrons", UINT );
    selElectrons.makeBranch( "nSelIsoElectrons", UINT );
    selElectrons.makeBranch( "epDrMatch", VFLOAT );
    selElectrons.makeBranch( "epSeedIdMatch", VFLOAT );
    selElectrons.makeBranch( "epSeedIdUnique", VFLOAT );
    selElectrons.makeBranch( "elePhoMinDr", VFLOAT );
    selElectrons.makeBranch( "eleIsLoose", VBOOL );
    selElectrons.makeBranch( "nSelElectrons", UINT );
    selElectrons.makeBranch( "eleHasSVMatch", VBOOL );
    selElectrons.makeBranch( "nEleSVMatched", INT );

    selElectrons.attachBranches( fOutTree );

    selPhotons.makeBranch( "nPhotons", UINT );
    selPhotons.makeBranch( "nSelPhotons", UINT ); 
    selPhotons.makeBranch( "selPhoEleVeto", UINT );
    selPhotons.makeBranch( "selPhoOOT", VBOOL );
    selPhotons.makeBranch( "selPhoSusyId", VFLOAT );
    selPhotons.makeBranch( "selPhoQuality", VINT ); 
    selPhotons.makeBranch( "selPhoTime", VFLOAT ); 

    selPhotons.makeBranch( "selPhoEnergy", VFLOAT );
    selPhotons.makeBranch( "selPhoEta", VFLOAT ); 
    selPhotons.makeBranch( "selPhoPhi", VFLOAT );     
	selPhotons.makeBranch( "selPhoPt", VFLOAT ); 
    selPhotons.makeBranch( "selPhoSCx", VFLOAT );
    selPhotons.makeBranch( "selPhoSCy", VFLOAT );
    selPhotons.makeBranch( "selPhoSCz", VFLOAT );
    selPhotons.makeBranch( "selPhoCorEnergy", VFLOAT );
    selPhotons.makeBranch( "selPhoCorPt", VFLOAT );
    selPhotons.makeBranch( "selPhoMBetaEql", VFLOAT );
    selPhotons.makeBranch( "selPhoMBetaPmt", VFLOAT );
    selPhotons.makeBranch( "selPhoPVtoRH", VFLOAT );

    selPhotons.makeBranch( "selPhoHcalTowerSumEtBcConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR03", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtSolidConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoPixelSeed", VBOOL );
    selPhotons.makeBranch( "selPhoEcalRHSumEtConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoHadTowOverEM", VFLOAT );
    selPhotons.makeBranch( "selPhoHadOverEM", VFLOAT );

    selPhotons.makeBranch( "selPhoHcalTowerSumEtConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoNTrkHollowConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoNTrkSolidConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoEcalPFClusterIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoHasConversionTracks", VBOOL );   //!
    selPhotons.makeBranch( "selPhoHcalPFClusterIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoHoe_PUcorr", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIsoPFPV", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIsoWorstVtx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfPhoIso03", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfRelIso03_all_quadratic", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfRelIso03_chg_quadratic", VFLOAT );   //!

    selPhotons.makeBranch( "selPhoSMaj", VFLOAT );
    selPhotons.makeBranch( "selPhoSMin", VFLOAT );
    selPhotons.makeBranch( "selPhoSAlp", VFLOAT );
    selPhotons.makeBranch( "selPhoClstrRn", VFLOAT );
    selPhotons.makeBranch( "selPhoR9", VFLOAT );
    selPhotons.makeBranch( "selPhoNrh", VUINT );
    selPhotons.makeBranch( "selPhoCovEtaEta", VFLOAT );
    selPhotons.makeBranch( "selPhoCovEtaPhi", VFLOAT );
    selPhotons.makeBranch( "selPhoCovPhiPhi", VFLOAT );
    selPhotons.makeBranch( "selPhoSieie", VFLOAT );
    selPhotons.makeBranch( "selPhoSieip", VFLOAT );
    selPhotons.makeBranch( "selPhoSipip", VFLOAT );
    selPhotons.makeBranch( "selPhoEtaWidth", VFLOAT );
    selPhotons.makeBranch( "selPhoPhiWidth", VFLOAT );
    selPhotons.makeBranch( "selPhoS4", VFLOAT );
    selPhotons.makeBranch( "selPhoSigmaIEtaIEta", VFLOAT );   //!

    selPhotons.makeBranch( "selPhoGenPt", VFLOAT );
    selPhotons.makeBranch( "selPhoPhoIsoDr", VFLOAT );
    selPhotons.makeBranch( "selPhoGenIdx", VINT );

    selPhotons.makeBranch( "selPhoGenSigMomEnergy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomEta", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomMass", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPhi", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPt", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPz", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVz", VFLOAT );   //!

    selPhotons.attachBranches( fOutTree );

    selMuons.makeBranch( "nLooseMuons", UINT );
    selMuons.makeBranch( "nMuons", UINT );
    selMuons.makeBranch( "muonIsLoose", VBOOL );
    selMuons.makeBranch( "muonIsMedium", VBOOL );
    //selMuons.makeBranch( "muonHasSVMatch", VBOOL );
    selMuons.makeBranch( "nMuonSVMatched", INT );

    selMuons.attachBranches( fOutTree );

    //selJets.makeBranch( "JetHt", &JetHt );
    selJets.makeBranch( "nJets", UINT );
    selJets.makeBranch( "nSelJets", UINT );
    selJets.makeBranch( "selJetGenLlpId", VINT ); 
    selJets.makeBranch( "selJetQrkLlpId", VINT );
    selJets.makeBranch( "selJetQuality", VINT ); 
    selJets.makeBranch( "selJetPt", VFLOAT ); 
    selJets.makeBranch( "selJetEnergy", VFLOAT );
    selJets.makeBranch( "selJetEta", VFLOAT ); 
    selJets.makeBranch( "selJetPhi", VFLOAT );
    selJets.makeBranch( "selJetTime", VFLOAT ); 
    selJets.makeBranch( "selJetMass", VFLOAT );

    selJets.makeBranch( "selJetArea", VFLOAT ); //*Jet_area)[it]; 
    selJets.makeBranch( "selJetChEmEF", VFLOAT ); //*Jet_chEmEF)[it]; 
    selJets.makeBranch( "selJetchHEF", VFLOAT ); //*Jet_chHEF)[it]; 
    selJets.makeBranch( "selJetChHM", VFLOAT ); //*Jet_chHM)[it]; 
    selJets.makeBranch( "selJetMuEF", VFLOAT ); //*Jet_muEF)[it]; 
    selJets.makeBranch( "selJetNeEmEF", VFLOAT ); //*Jet_neEmEF)[it]; 
    selJets.makeBranch( "selJetNeHEF", VFLOAT ); //*Jet_neHEF)[it]; 
    selJets.makeBranch( "selJetNeHM", VFLOAT ); //*Jet_neHM)[it]; 

    selJets.makeBranch( "selGenJetDpt", VFLOAT ); //*Jet_genDptMatch)[it]; 
    selJets.makeBranch( "selGenJetdr", VFLOAT ); //*Jet_genDrMatch)[it]; 
    selJets.makeBranch( "selGenJetEnergy", VFLOAT ); //*Jet_genEnergy)[it]; 
    selJets.makeBranch( "selGenJeteta", VFLOAT ); //*Jet_genEta)[it]; 
    selJets.makeBranch( "selGenJetImpAng", VFLOAT ); //*Jet_genImpactAngle)[it]; 
    selJets.makeBranch( "selJetLlpDp", VFLOAT ); //*Jet_genLlpDp)[it]; 
    selJets.makeBranch( "selJetLlpDr", VFLOAT ); //*Jet_genLlpDr)[it]; 
    selJets.makeBranch( "selGenJetPt", VFLOAT ); //*Jet_genPt)[it]; 
    selJets.makeBranch( "selGenJetTof", VFLOAT ); //*Jet_genTOF)[it]; 
    selJets.makeBranch( "selGenJetTime", VFLOAT ); // (*Jet_genTime)[it]; 
    selJets.makeBranch( "selGenJetLlpTime", VFLOAT ); //*Jet_genTimeLLP)[it]; 

    selJets.attachBranches( fOutTree );

    selRjrVars.makeBranch( "rjrX1aMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX1aCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX1bMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX1bCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Pt", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Eta", VFLOAT );
    selRjrVars.makeBranch( "rjrX1a_Phi", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Pt", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Eta", VFLOAT );
    selRjrVars.makeBranch( "rjrX1b_Phi", VFLOAT );

    selRjrVars.makeBranch( "rjrX2aMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX2aCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX2aPtS", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bMass", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrX2bPtS", VFLOAT );

    selRjrVars.makeBranch( "rjrSMass", VFLOAT );
    selRjrVars.makeBranch( "rjrSCosA", VFLOAT );
    selRjrVars.makeBranch( "rjrSdphiDA", VFLOAT );
    selRjrVars.makeBranch( "rjrSdphiBV", VFLOAT );
    selRjrVars.makeBranch( "rjrPTS", VFLOAT );
    selRjrVars.makeBranch( "rjrPZS", VFLOAT );

    selRjrVars.makeBranch( "rjrEVa", VFLOAT );
    selRjrVars.makeBranch( "rjrEVb", VFLOAT );
    selRjrVars.makeBranch( "rjrPVa", VFLOAT );
    selRjrVars.makeBranch( "rjrPVb", VFLOAT );

    //selRjrVars.makeBranch( "rjrMV", VFLOAT );
    //selRjrVars.makeBranch( "rjrPV", VFLOAT );
    //selRjrVars.makeBranch( "rjrMVa", VFLOAT );
    //selRjrVars.makeBranch( "rjrMVb", VFLOAT );

    selRjrVars.makeBranch( "rjrNPhotons", VINT );
    selRjrVars.makeBranch( "rjrMET", VFLOAT );
    selRjrVars.makeBranch( "rjrNJets", VINT );
    selRjrVars.makeBranch( "rjrSubPhoLocation", VINT);
    selRjrVars.makeBranch( "rjrNVisObjects", VINT );
    selRjrVars.makeBranch( "rjrABSide", VBOOL );

    selRjrVars.makeBranch( "rjrVisPhoSide", VBOOL );
    selRjrVars.makeBranch( "rjrVisJetSide", VBOOL );

    selRjrVars.makeBranch( "rjrPVlab", VFLOAT );
    selRjrVars.makeBranch( "rjrDphiMETV", VFLOAT );

    //selRjrVars.makeBranch( "rjrX2NQSum", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2NAve", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2NGMean", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2QSum", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2Ave", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2GMean", VFLOAT );
    //selRjrVars.makeBranch( "rjrX2Diff", VFLOAT );

    selRjrVars.makeBranch( "rjrNJetsJb", VINT );
    selRjrVars.makeBranch( "rjrNJetsJa", VINT );

    //selRjrVars.makeBranch( "rjrPX1X2a", VFLOAT );
    //selRjrVars.makeBranch( "rjrPX1X2b", VFLOAT );
    //selRjrVars.makeBranch( "rjrPJX2a", VFLOAT );
    //selRjrVars.makeBranch( "rjrPJX2b", VFLOAT );
    //selRjrVars.makeBranch( "rjrPX2Sa", VFLOAT );
    //selRjrVars.makeBranch( "rjrPX2Sb", VFLOAT );
    //selRjrVars.makeBranch( "rjrDiffPJX2a", VFLOAT );
    //selRjrVars.makeBranch( "rjrDiffPJX2b", VFLOAT );

    //selRjrVars.makeBranch( "rjrX1a_Sm", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1a_Sp", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1a", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1a", VFLOAT );

    //selRjrVars.makeBranch( "rjrX1X2a", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1X2a", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1b", VFLOAT );
    //selRjrVars.makeBranch( "rjrX1b", VFLOAT );

    //selRjrVars.makeBranch( "rjrAX2aMass", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2bMass", VFLOAT );
    //selRjrVars.makeBranch( "rjrASMass", VFLOAT );

    //selRjrVars.makeBranch( "rjrAX2NQSum", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2NAve", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2NGMean", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2QSum", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2Ave", VFLOAT );
    //selRjrVars.makeBranch( "rjrAX2GMean", VFLOAT );
    selRjrVars.makeBranch( "rjrAX2Diff","rjr_Rdiff", VFLOAT );

    //selRjrVars.makeBranch( "rjrMVDiff", VFLOAT );
    //selRjrVars.makeBranch( "rjrMVSum", VFLOAT );
    //selRjrVars.makeBranch( "rjrMVNSum", VFLOAT );

    selRjrVars.makeBranch( "rjr_Mr", VFLOAT );
    selRjrVars.makeBranch( "rjr_R", VFLOAT );
    selRjrVars.makeBranch( "rjr_Rv", VFLOAT );

    selRjrVars.makeBranch( "rjrN2Px", VFLOAT );
    selRjrVars.makeBranch( "rjrN2Py", VFLOAT );
    selRjrVars.makeBranch( "rjrN2Pz", VFLOAT );

    selRjrVars.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

void KUCMSAodSkimmer::endJobs(){ 

    std::cout << "Running EndJobs." << std::endl;
	hist1d[11]->Divide(hist1d[10]);
    hist1d[12]->Divide(hist1d[10]);
    hist1d[13]->Divide(hist1d[10]);
    hist1d[14]->Divide(hist1d[10]);

}//void KUCMSAodSkimmer::endJobs()

void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree ){

    std::cout << "Filling ConfigTree." << std::endl;
	TBranch *nEventBranch = fConfigTree->Branch( "nEvents", &nEvents );
	nEventBranch->Fill();
	TBranch *nSelectedEventsBranch = fConfigTree->Branch( "nSelectedEvents", &nSelectedEvents );
	nSelectedEventsBranch->Fill();
    TBranch *sKeyBranch = fConfigTree->Branch( "sKey", &dataSetKey ); 
	sKeyBranch->Fill();
    TBranch *sCrossSectionBranch = fConfigTree->Branch( "sCrossSection", &xsctn );
	sCrossSectionBranch->Fill();
    TBranch *sLambdaBranch = fConfigTree->Branch( "sGMSBGravMass", &gmass );
    sLambdaBranch->Fill();
    TBranch *sCTauBranch = fConfigTree->Branch( "sGMSBChi1Mass", &xmass );
    sCTauBranch->Fill();
    TBranch *sMCWgtBranch = fConfigTree->Branch( "sMCWgt", &mcwgt );
    sMCWgtBranch->Fill();
    TBranch *sMCTypeBranch = fConfigTree->Branch( "sMCType", &mctype );
    sMCTypeBranch->Fill();
    //TBranch *sumEvtGenWgtBranch = fConfigTree->Branch( "sumEvtGenWgt", &sumEvtGenWgt );
    //sumEvtGenWgtBranch->Fill();

    for( auto item : cutflow ){
		std::string bname = "cf_" + item.first;
		TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
		cfBranch->Fill();
    }//<<>>for( auto item : configInfo )

    for( auto item : configCnts ){
        std::string bname = item.first;
        TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
        cfBranch->Fill();
    }//<<>>for( auto item : configInfo )

    for( auto item : configWgts ){
        std::string bname = item.first;
        TBranch *cfBranch = fConfigTree->Branch( bname.c_str(), &item.second );
        cfBranch->Fill();
    }//<<>>for( auto item : configInfo )

	fConfigTree->Fill();

}//<<>>void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree, std::string key )

void KUCMSAodSkimmer::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = new TH1D("ecalrhenergy0", "RecHit Energy;rechit E [GeV]",1000,0,1000);
    hist1d[1] = new TH1D("ecalrhenergy1", "RecHit Energy;rechit E [GeV]",200,0,20);
    hist1d[2] = new TH1D("ecalrhenergy2", "RecHit Energy;rechit E [GeV]",200,0,2);
    hist1d[3] = new TH1D("sctype","SC !Orig, Orig, OOT, Excl",4,0,4);
    hist1d[4] = new TH1D("scorigtype","Orig+OOT, Orig+!OOT, Orig+Exc, !Orig+OOT, !Orig+Exc, !Orig+!OOT",6,0,6);
    hist1d[5] = new TH1D("ecalrhtime", "RecHit Time;rechit t [ns]",200,-10,10);
    hist1d[6] = new TH1D("scenergy","SC rawEnergy;energy [GeV]",200,0,1000);
    hist1d[7] = new TH1D("trckpt","Track pt;Track pt [GeV]",500,0,1000);


	hist1d[10] = new TH1D("elept","Electron Pt [GeV];Electron Pt [GeV];a.u.",100,0,1000);
    hist1d[11] = new TH1D("elept_lv","LooseVeto;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[12] = new TH1D("elept_v","Veto;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[13] = new TH1D("elept_l","Loose;Electron Pt [GeV];Eff",100,0,1000);
    hist1d[14] = new TH1D("elept_m","Medium;Electron Pt [GeV];Eff",100,0,1000);

    ////hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    ////hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[1] = new TH2D("scorigcross", "SC Types;Mian-Orig;Std-OOT-Excluded", 2, 0, 2, 3, 0, 3); 

	//------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

    //hist3d[0] = new TH3D("phoNRH_ClR9_phoID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",200,0,200,100,0,1,8,0,4);

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------

    if( CLSTRMAPS ){
		nMaps = 0;
		for(int it=0; it<nEBEEMaps; it++){
			fMap[it] = false;
			std::string label(";iEta;iPhi");
    	    std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
    	    ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
    	    std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
    	    ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
			std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
    	    ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
		}//<<>>for(int it=0; it<nEBEEMaps; it++)
	}//<<>>if( CLSTRMAPS )

}//<<>>void KUCMSAodSkimmer::initHists()

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//	Helper functions that depended on varibles from GetEntry()
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getRhIdx( uInt rhDetID ){

    for( int idx = 0; idx < ECALRecHit_ID->size(); idx++ ){ if( rhDetID == (*ECALRecHit_ID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int KUCMSAodSkimmer::getRhIdx( int rhDetID )

uInt KUCMSAodSkimmer::getLeadRhID( std::vector<uInt> recHitIds ){

    uInt result;
    float enr(0.0);
	if( recHitIds.size() < 1 || ECALRecHit_ID->size() < 1 ){ return 0; }
    for( auto id : recHitIds ){
        auto rhenr = (*ECALRecHit_energy)[getRhIdx(id)];
        if( rhenr > enr ){ enr = rhenr; result = id; }
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit KUCMSAodSkimmer::getLeadRh( rhGroup recHitsi

float KUCMSAodSkimmer::clstrR9( std::vector<uInt> recHitIds ){

	if( recHitIds.size() < 1 || ECALRecHit_ID->size() < 1 ){ return -2.0; }
    auto leadRhID = getLeadRhID( recHitIds );
    auto leadRhEn = (*ECALRecHit_energy)[getRhIdx(leadRhID)];
    float sumRhEn(0);
    if( recHitIds.size() < 1 ){ return -99.0; }
    for ( auto id : recHitIds ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
    return sumRhEn > 0 ? leadRhEn/sumRhEn  : - 1.0;

}//<<>>float KUCMSAodSkimmer::clstrR9( vector<uInt> recHitIds )

std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids ){

	std::vector<float> result;
    if( rechitids.size() < 1 || ECALRecHit_ID->size() < 1 ){ result.push_back(-9.0); return result; }
	for ( auto id : rechitids ){ result.push_back((*ECALRecHit_energy)[getRhIdx(id)]); }
	return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids ){

    std::vector<float> result;
	if( rechitids.size() < 1 || ECALRecHit_ID->size() < 1 ){ result.push_back(-99.0); return result; }
    for ( auto id : rechitids ){
		auto rhtime = (*ECALRecHit_time)[getRhIdx(id)];
		auto rh0tof = (*ECALRecHit_0TOF)[getRhIdx(id)];
        auto rhpvtof = (*ECALRecHit_pvTOF)[getRhIdx(id)];
		result.push_back(rhtime+rh0tof-rhpvtof); 
	}//<<>>for ( auto id : recHitIds )
    return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpEigenFromAngles( std::vector<uInt> rechitids ){

	//bool verbose = true;
    bool verbose = false;

	if( DEBUG || verbose ) std::cout << " Starting getRhGrpEigen_sph  " << std::endl;

    std::vector<float> emptyReturn(9,-9);
    std::vector<float> egwts;
    std::vector<float> rhetas, rhphis;
    std::vector<float> logwtvec, tresvec;
    auto nRecHits = rechitids.size();
    if( nRecHits < 5 || ECALRecHit_ID->size() < 1 ){ if( verbose ) std::cout << " ----  rechit collection has too few rechits" << std::endl; return emptyReturn; }
    float sumRhEn(0);
    for ( auto id : rechitids ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
	if( verbose ) std::cout << " --- EigenAngles sumRhEn : " << sumRhEn << std::endl;
    if( sumRhEn <= 0 ){ if( verbose ) std::cout << " ----  rechit collection has no energy" << std::endl; return emptyReturn; }
    if( DEBUG ) std::cout << "1a, ";
    std::cout << " getRhGrpEigenFromAngles : nRechits: " << nRecHits << " -------------------------- " << std::endl;
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
        auto idinfo = DetIDMap[rechitids[it]];
        auto isEB = idinfo.ecal == ECAL::EB;
        if( rhIDX == -1 ){ if( verbose ) std::cout << " ---- Bad idx !!!!! -- In getRhGrpEigen ---- " << std::endl; return emptyReturn; }
        if( not isEB ){ if( verbose ) std::cout << " ---- rechit group has EE members " << idinfo.ecal << std::endl; return emptyReturn; }

        std::cout << " rhIDX: " << rhIDX << std::endl;
    	const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
        std::cout << " rhEtaPos: " << rhEtaPos << std::endl;
        rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
        const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
        std::cout << " rhPhiPos: " << rhPhiPos << std::endl;
        rhphis.push_back(rhPhiPos-0.5);
        auto rhenergy = (*ECALRecHit_energy)[rhIDX];
        auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));// cut at rh energy < 1.5% of cluster
        logwtvec.push_back(logwt);

    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    std::vector<float> detas, dphis, angles, redlogwtvec;
    auto meta = mean( rhetas, logwtvec );
    auto mphi = meanIPhi( rhphis, logwtvec );
    for( uInt it(0); it < rhetas.size(); it++ ){

        float deta = rhetas[it]-meta;
        float dphi = dIPhi( rhphis[it], mphi );
		if( dphi < 0.5 && deta < 0.5 ) continue;
        detas.push_back(deta);
        dphis.push_back(dphi);
        float angle = getAngle( deta, dphi );
        angles.push_back(angle);
		redlogwtvec.push_back(logwtvec[it]);

    }//<<>>for( uInt it(0); it < etas.size(); it++ )

    auto eigens = getRhGrpEigen( angles, redlogwtvec );
	if( verbose ) std::cout << " --- EigenAngles VAlue : " << eigens[2] << std::endl;

    auto phiCorrFactor = 0.8;
    std::cout << " mean eta: " << meta << " mean phi: " << mphi << std::endl; 
    auto sxx = var( detas, meta, redlogwtvec );
    auto syy = var( dphis, mphi, redlogwtvec, accum(redlogwtvec)/phiCorrFactor );
    auto sxy = cvar( detas, meta, dphis, mphi, redlogwtvec, accum(redlogwtvec)/std::sqrt(phiCorrFactor) );
    std::cout << " sxx: " << sxx << " syy: " << syy << " sxy: " << sxy << std::endl;
    auto smaj = (sxx + syy + std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto smin = (sxx + syy - std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto sang = std::atan((sxx-syy+std::sqrt(sq2(syy-sxx)+4.*sq2(sxy)))/(2.*sxy));

    //eigens[0] //0 geoeigan x vec
    //eigens[1] //1 geoeigan y vec
    //eigens[2] //2 geoeigan mag vec
    eigens.push_back(smaj);//3
    eigens.push_back(smin);//4
    eigens.push_back(sang);//5

    if( DEBUG ) std::cout << " Done" << std::endl;;
    return eigens;

}//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEigen( std::vector<uInt> rechitids )

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////!!!!!!!!!!!!!!!!!11  need rhPosX, rhPosY, & rhPosZ to be saved for following functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

std::vector<float> KUCMSAodSkimmer::getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ ){

    std::vector<float> result;
    if( ( recHitIds.size() < 1 ) || ( ECALRecHit_ID->size() < 1 ) ){ result.push_back(-99); return result; }
    auto lrhid = getLeadRhID(recHitIds);
    auto lrhidx = getRhIdx(lrhid);
    auto X = (*ECALRecHit_rhx)[lrhidx];
    auto Y = (*ECALRecHit_rhy)[lrhidx];
    auto Z = (*ECALRecHit_rhz)[lrhidx];
    const auto d_rh = hypo( X, Y, Z);
    const auto d_pv = hypo( X-vtxX, Y-vtxY, Z-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for( int idx = 0; idx < ECALRecHit_time->size(); idx++ ){result.push_back((*ECALRecHit_time)[idx]-tof);}
    return result;

}//>>>>vector<float> KUCMSAodSkimmer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

