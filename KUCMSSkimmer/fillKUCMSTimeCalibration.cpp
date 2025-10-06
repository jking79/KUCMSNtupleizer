//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//`
//
// Original Author:  Jack W King III
//         Created:  Fri, 20 Dec 2024
//
//////////////////////////////////////////////////////////////////////

#include "ecal_config/KUCMSTimeCalibration.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    ///////////////////////////////////////////////////////////////////
    // this section should be moved to run file --- DONE --------------
    // no need to "recreate" cali maps every time, they will be saved -
    ///////////////////////////////////////////////////////////////////
    // standard usage for applying calibrations :
    // KUCMSTimeCalibration theCali;
    // theCali.setTag( tag ); 
    // 	-- or to use default skip this step 
    // 	-- only needs to be done once unless you want to change the tag to use
    // time = time - theCali.getCalibration( rechitDetID, run );
    //
    // to get smeared time :
    // theCali.setSmearTag( tag ); 
    //  -- or to use default skip this step
    //  -- only needs to be done once unless you want to change the tag to use
    // time = theCali.getSmearedTime( time, amplitude );
    ///////////////////////////////////////////////////////////////////

	// Cali Tags : Tags for calibrations to use
	std::string r2EOY( "EG_EOY_MINI" );
	std::string r2Fall17AOD( "RunIIFall17DRPremix" );
    std::string r2Fall17MINIAOD( "RunIIFall17MiniAODv2" );
    std::string r2UL( "UL_R2_MINI" );

	// IOV tags :  defines IOV periods to make calibration maps
    // for mc use "mc" for TTIov and XIov

    //-----//////////  IOV tags : making tt cali file :
    //std::string r2ulTag( "r2ul" );
    std::string r2ulTag( "r2ultt" );
    std::string mctag( "mc" );
    std::string xiovtag( "prompt" );
    std::string r3TagTT( "r3tt" );
    std::string r3TagX( "r3x" );

    //std::string inputfilename( "kucmsTimeCaliR18TFile.txt");
    std::string inputfilename( "kucmsTimeCaliR24FCCvRtTFile.txt");
    //std::string inputfilename( "kucmsTimeCaliTestFile.txt" ); // MET_AOD_R17_FULL 
    //std::string inputfilename( "kucmsTimeCaliR17File.txt" ); // EG_EOY_MINI PD
    //std::string inputfilename( "kucmsTimeCaliR17Plot2DFullRunFile.txt" );
    //std::string inputfilename( "kucmsTimeCaliRunIIFall17File.txt" ); // RunIIFall17DRPremix

	std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// input parameter!
    //std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");

	//std::string indir("ecalTiming/gammares_llpana_pd/");// input paramter !  
    //std::string indir("/ecalTiming/gammares_llpana_mc/");// input paramter ! 
    //std::string indir("/kuncali/gammares_cali_mc/");
    //std::string indir("/kuncali/gammares_cali/");
    //std::string indir("/ecalTiming/gammares_cali_mc/");
    //std::string indir("ecalTiming/gammares_ECAL_CC_HCAL_DI-v3/");
    //std::string indir("ecalTiming/gammares_r24f_prompt/");
    //std::string indir("ecalTiming/gammares_r25_val/");
    //std::string indir("ecalTiming/gammares_DPG/");
    //std::string indir("ecalTiming/gammares_ul18/EGamma/");
    std::string indir("ecalTiming/gammares_DPG_25/");
    //std::string indir("ecalTiming/gammares_DPG_24/");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//using ( true, true ) will start a new TClifile  - all histograms will be lost
    //////KUCMSTimeCalibration theCali( true, true ); // make a new TClifile and keep TCalifile open
	//KUCMSTimeCalibration theCali( true, true );
	// if Tcalifile kept open be sure to make and save tcalifile ( saving calihist will close tcalifile )

	////KUCMSTimeCalibration theCali; // open tfile in read only then close after setup
	//KUCMSTimeCalibration theCali( true );
	//theCali.SetEosDir(eosdir);
	//theCali.SetInDir(indir);
	//
    //-----//////////  making tt cali  :
    //theCali.SetupIovMaps();
	// for PD R2UL
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
	// for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    // for PD R3
	//theCali.setTTIov( r3TagTT );
	//theCali.setXIov( xiovtag );
    //theCali.makeTTCaliMapEGR( inputfilename, true ); // true == run only subset of events
    //theCali.makeTTCaliMapEGR( inputfilename );
    //theCali.makeCaliHists();
    //theCali.SaveCaliHists();
    //theCali.SaveTTRunFile();

    //-----//////////  making xtal cali :
    //theCali.SetupIovMaps();
    // for PD
    //theCali.setTTIov( r2ulTag );
    //theCali.setXIov( xiovtag );
    // for MC
    //theCali.setTTIov( mctag );
    //theCali.setXIov( mctag );
    //theCali.makeXCaliMapEGR( inputfilename, true ); // true == run only subset of events
    //theCali.makeXCaliMapEGR( inputfilename );
    //theCali.makeCaliHists();
    //theCali.SaveCaliHists();
    //theCali.SaveCaliRunFile();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// set external calibration :
    //KUCMSTimeCalibration theCali;
	//std::string exttag = "r324fpmt";
    //std::string exttag = "r324fcc";
	//std::string extdir = "/uscms/home/jaking/nobackup/el8/ecaltiming/CMSSW_14_0_11/src/GammaResTool/GammaResTool/macros/cali_root_files/";
    //std::string calihistnamecc = "KURes_14011_v12_DEGPD_MiniAOD_ECAL_CC_HCAL_DI-v3_Cali.root";
    //std::string calihistnamert = "KURes_14011_v12_DEGPD_MiniAOD_Run2024F_Cali.root";
    //std::string mapname = "Run3_2024F_Prompt";
    //std::string mapname = "Run3_2024F_CCVal";
	//int startr = 1;
	//int endr = 999999;
    //std::string calihist = extdir+calihistnamert;
	//std::string calihist = extdir+calihistnamecc;
	//theCali.LoadExtCali( calihist, mapname, exttag, startr, endr );
	//theCali.useExtCali( true );
	//theCali.setTag( exttag );	

    // 2d resolution maps  
	//KUCMSTimeCalibration theCali; 
    //theCali.SetEosDir(eosdir);
    //theCali.SetInDir(indir);
    //  - 2d res plot  -- defaults to : false, true, false, "" ( scale, usecali, smear, name ext ) 
    //theCali.plot2dResolutionEGR( inputfilename );

	// 2d resolution maps smeared/uncali
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_nocali" );// input, scale, usecali
    //theCali.setSmearTag("EG300202_DYF17");
    //theCali.plot2dResolutionEGR( inputfilename, true, true, true, "_smeared" );// input, scale, usecali, smear, name ext

	// make resolution paramters for iov periods
    //theCali.doResTimeFits();
	//theCali.SaveCaliRunFile();	

    // make resolution paramters for smearing
    //std::string sourceTag = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist";
	//std::string run = "304476";
	//std::string destTag = "EG_EOY_MINI_" + run + "_X_ZEE_Data_Hist";
	//std::string smearTag = "EG" + run + "_DYF17";
    //theCali.makeSmearTag( sourceTag, destTag, smearTag );
	//theCali.SaveSmearFile();

	// standard resolution 
	//std::string histName = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist_Smeared_v4";
    //std::string histName = "EG_EOY_MINI_304476_X_SRO_Data_Hist_NoCali_nocali_v4";
    //theCali.load2DResHist( histName );
	//theCali.doResTimeFit( histName );

	// extended range 2D
	//KUCMSTimeCalibration theCali;
    //theCali.SetEosDir(eosdir);
    //theCali.SetInDir(indir);
    //theCali.setLowEnergy( true );
    //theCali.setUseEffEnergy( true );
    //theCali.useGainSwitch( true );
	//theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : xr
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
    //theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 40.0 50.0 70.0 90.0" ); // xer
    //theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    //theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : ecal real
 	//theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1200 1600 2200" ); // : xgs ecal
    //////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj : justin
    //////theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
	//--------------------------------------------------------
	// my cali time profile?
    //theCali.SetYBinStr( "CONSTANT 1800 -9 9" ); // jack  ( 10 ps )
    //theCali.SetYBinStr( "CONSTANT 240 -3 3" ); // ecal ( 25 ps ) DPG
    //theCali.SetYBinStr( "CONSTANT 480 -6 6" ); // ecal ( 25 ps )
    ////theCali.SetYBinStr( "CONSTANT 720 -9 9" ); // ( 25 ps )
    //theCali.SetYBinStr( "CONSTANT 360 -9 9" ); // ( 50 ps  )
    //theCali.SetYBinStr( "CONSTANT 240 -12 12" ); ( 100 ps )
    //--------------------------------------------------------
	// justin profile y bins
    //theCali.SetYBinStr( "CONSTANT 240 -12 12" ); // justin ( 100 ps )
    //theCali.SetYBinStr( "CONSTANT 180 -6 6" ); // justin ( 67 ps )
    //theCali.SetYBinStr( "CONSTANT 180 -3 3" ); // justin ( 33 ps )
    //--------------------------------------------------------
    // special 2d input file
    ////std::string inputfilename( "kucmsTimeCaliR17Plot2DFullRunFile.txt" );
	////std::string inputfilename( "kucmsTimeCaliR24FCCvRtTFile.txt");
    //-------------------------------------------------------
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_DegR17Full_xea_pm9b1800_v417" );// scale, cali, smear
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_egr17ef_xea_pm9b180_v327" );// scale, cali, smear
    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_egr17ef_xea_pm3b180_v417" );// scale, cali, smear
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_xa_pm9b720_r325" );// scale, cali, smear
	//theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_deg24fval_xa_pm9b1800_v501" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg24fprmt_xgs_pm3b240_v508" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg24fcc_xer_pm3b240_v508" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg25bc_xa_pm6b480_v616" );
	//theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_deg18c_xa_pm6b480_v521" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg24d_xa_pm6b480_v810" );

    // extended range resfit
    KUCMSTimeCalibration theCali;
    //std::string histName = "ResMap_303832_304616__DRO_Data_Hist_ef_xea_pm9b180_v327";
    //std::string histName = "ResMap_303832_304616__DRO_Data_Hist_NoCali_ef_xea_pm9b180_v327";
    //std::string histName = "ResMap_305044_305081__DRO_Data_Hist_ef_xea_pm9b180_v327";
    //std::string histName = "ResMap_305044_305081__DRO_Data_Hist_NoCali_ef_xea_pm9b180_v327";
	//std::string histName = "ResMap_297114_305081__SRO_Data_Hist_jr17_xea_pm12b240_v327";
    //std::string histName = "ResMap_297114_305081__ZEE_Data_Hist_NoCali_jr17_xea_pm12b240_v327";
	//std::string histName = "ResMap_297114_305081__SRO_Data_Hist_degr17full_xea_pm9b1800_v327";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_deg24fprmt_xa_pm9b1800_v327";
	//std::string histName = "ResMap_297050_306460__SRO_Data_Hist_NoCali_DegR17Full_xea_pm9b1800_v417";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_deg24fval_xa_pm9b1800_v417";
    //std::string histName = "ResMap_297050_306460__SRO_Data_Hist_egr17ef_xea_pm12b240_v417";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_deg24fprmt_xa_pm3b240_v417";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_deg24fprmt_xer_pm3b240_v508";
	//std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_deg25b_xa_pm6b480_v519";
	//std::string histName = "ResMap_0_999999__ZEE_Data_Hist_deg18c_xa_pm6b480_v521";
	//std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_deg25bc_xa_pm6b480_v602";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_deg25bc_xa_pm6b480_v616";
	//std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_eg24f_ee_xa_pm6b480_v808";
    //std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_eg24fcc_ee_xa_pm6b480_v808";
	//std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_eg24fcc_ee_xa_pm6b480_v808";
    //std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_eg25cv2_xa_pm6b480_v807";
    //std::string histName = "ResMap_0_999999__SRO_Data_Hist_NoCali_eg25cv2_xa_pm6b480_v808";
	//std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_eg24fcc_xa_pm6b480_v808";
	std::string histName = "ResMap_0_999999__ZEE_Data_Hist_NoCali_eg24d_xa_pm6b480_v810";

	//
	//theCali.Set2DResTFileName( "egres_DEGPD_AODSIM_R24_Fprompt_v3_14011_v12_resplots.root" );
    //theCali.Set2DResTFileName( "egres_DEGPD_AODSIM_R24_ECAL_v3_14011_v12_resplots.root" );
	//
	//
	theCali.load2DResHist( histName );
	theCali.setLowEnergy( true );
	theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1000 1200 1800" ); // : ecal
    //theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    //theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
    //theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj
    //theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : xr
    //theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 40.0 50.0 70.0 90.0" ); // xer
    //theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1200 1600 2200" ); // : xgs ecal
	theCali.doResTimeFit( histName );

	// plotting of mean time by run with calibraton : filename, start run, end run, usecali
	//theCali.plotMeanRunTimeEGR( inputfilename, 303838, 304796 );
    //theCali.plotMeanRunTimeEGR( inputfilename, 296399, 307554, true );
	//theCali.makeTTDiffMaps();// make trigger tower diffrence maps
	
	//KUCMSTimeCalibration theCali;
	//theCali.makeTTDriftMaps( "r324fcc", 382008, 383814 );	
    //theCali.makeTTDriftMaps( "r324fpmt", 382008, 383814 );
	//theCali.makeTTDriftMaps( "EG_EOY_MINI", 296399, 306460 );
	//theCali.makeTTDriftMaps( "UL_R2_MINI", 316996, 325200 );

    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
