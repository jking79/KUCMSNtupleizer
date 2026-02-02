//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//`
//
// Original Author:  Jack W King III
//         Created:  Fri, 20 Dec 2024
//
//////////////////////////////////////////////////////////////////////

#include "KUCMS_TimeCalibration.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    ///////////////////////////////////////////////////////////////////
    // this section should be moved to run file --- DONE --------------
    // no need to "recreate" cali maps every time, they will be saved -
    ///////////////////////////////////////////////////////////////////
    // standard usage for applying calibrations :
    // KUCMS_TimeCalibration theCali;
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
	std::string r2EOY( "EOY_R2_MINI" );
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

	//std::string inputfilename( "ecal_config/kucmsTimeCaliTestTFile.txt");
    ////std::string inputfilename( "ecal_config/kucmsTimeCaliR17EOYTFile.txt");

    std::string inputfilename( "ecal_config/kucmsTimeCaliR18MCGjetsTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR18ULTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR17ULTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR16ULTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR22PrmptTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR23PrmptTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR24PrmptTFile.txt");

	//std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");// input parameter!
    std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");

    ////std::string indir("/kuncali/gammares_cali_mc/");
    ////std::string indir("/kuncali/gammares_cali/");
    ////std::string indir("kuncali/gammares_mc18/");
    ////std::string indir("kuncali/gammares_cali/DoubleEG/");

    std::string indir("KUCMSNtuple/gammares_ul18/");
    //std::string indir("KUCMSNtuple/gammares_ul17/");
    //std::string indir("KUCMSNtuple/gammares_ul16/DoubleEG/");
    //std::string indir("KUCMSNtuple/gammares_prmt22/EGamma/");
    //std::string indir("KUCMSNtuple/");

    KUCMS_TimeCalibration theCali;
    theCali.SetEosDir(eosdir);
    theCali.SetInDir(indir);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// set resoltuions
	//theCali.SetResParamters( "default", 26.04, 1.19, 0.1078, 38.69, 0.00, 0.0843 );
	//theCali.SetResParamters( "RunIISummer20UL18RECO", 26.04, 1.19, 0.1078, 38.69, 0.00, 0.0843 );
	//theCali.SetResParamters( "r2_ul18", 26.12, 7.60, 0.1006, 42.03, 0.00, 0.3041 );
    //theCali.SetResParamters( "r2_ul18_mc",  26.23, 1.72, 0.1391, 33.36, 1.72, 0.2039 );
    //theCali.SetResParamters( "r2_ul17", 26.55, 1.51, 0.1447, 33.14, 1.59, 0.2075 );
    //theCali.SetResParamters( "r2_ul16", 26.30, 1.87, 0.1606, 34.48, 1.05, 0.1933 );
    //theCali.SetResParamters( "r2_eoy17", 31.3, 2.9, 0.1608, 31.3, 2.9, 0.1608 );
    //theCali.SetResParamters( "r3_p22", 25.50, 3.09, 0.1728, 43.00, 0.00, 0.2495 );
    //theCali.SetResParamters( "r3_p23", 27.81, 7.36, 0.1338, 41.49, 5.52, 0.3146 );
	//theCali.SetResParamters( "r3_p24", 23.41, 1.79, 0.1707, 35.09, 0.00, 0.2533 );
    //theCali.SaveResFile();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// set external calibration :
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
    ////  - 2d res plot  -- defaults to : false, true, false, "" ( scale, usecali, smear, name ext ) 
    ////theCali.plot2dResolutionEGR( inputfilename );

	// 2d resolution maps smeared/uncali
    ////theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_nocali" );// input, scale, usecali
    ////theCali.setSmearTag("EG300202_DYF17");
    ////theCali.plot2dResolutionEGR( inputfilename, true, true, true, "_smeared" );// input, scale, usecali, smear, name ext

	// make resolution paramters for iov periods
    ////theCali.doResTimeFits();
	////theCali.SaveCaliRunFile();	

    // make resolution paramters for smearing
    ////std::string sourceTag = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist";
	////std::string run = "304476";
	////std::string destTag = "EG_EOY_MINI_" + run + "_X_ZEE_Data_Hist";
	////std::string smearTag = "EG" + run + "_DYF17";
    ////theCali.makeSmearTag( sourceTag, destTag, smearTag );
	////theCali.SaveSmearFile();

	// standard resolution 
	////std::string histName = "DYJetsToLLRunIIFall17AOD_0_X_ZEE_Data_Hist_Smeared_v4";
    ////std::string histName = "EG_EOY_MINI_304476_X_SRO_Data_Hist_NoCali_nocali_v4";
    ////theCali.load2DResHist( histName );
	////theCali.doResTimeFit( histName );

	// extended range 2D
    //theCali.setLowEnergy( false ); // true is default
    ////theCali.setUseEffEnergy( true );
    theCali.useGainSwitch( true, true ); // use gainid ? --  do all gians ? : all Gains (HG)
    //theCali.useGainSwitch( true, false ); // use gainid ? --  do all gians ? : only gainID1 (LG)
	////theCali.setDoEE( false ); depreciated - now do EB & EE at the same time
    ////theCali.setDoEE( true ); depreciated - now do EB & EE at the same time
	////theCali.setUseGainID( 1 ); depreciated - now set by useGainSwitch
	////theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : xr
    ////theCali.SetXBinStr( "VARIABLE 5 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1600" ); // : xlg  llpana LG resfit ( old xa )
    theCali.SetXBinStr( "VARIABLE 5 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1600 3200" ); // : xa  llpana HG resfit default
    ////theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 40.0 50.0 70.0 90.0" ); // xer
    ////theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    ////theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : ecal real
 	//-theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1200 1600 2200" ); // : xgs ecal
    //////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj : justin
    //////theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
	//--------------------------------------------------------
	// my cali time profile?
    ////theCali.SetYBinStr( "CONSTANT 1800 -9 9" ); // jack  ( 10 ps )
    //theCali.SetYBinStr( "CONSTANT 240 -3 3" ); // ecal ( 25 ps ) -- ecalDPG
    ////theCali.SetYBinStr( "CONSTANT 480 -6 6" ); // ecal ( 25 ps )
    ////theCali.SetYBinStr( "CONSTANT 720 -9 9" ); // ( 25 ps )
    ////theCali.SetYBinStr( "CONSTANT 360 -9 9" ); // ( 50 ps  )
    ////theCali.SetYBinStr( "CONSTANT 240 -12 12" ); ( 100 ps )
    //--------------------------------------------------------
	// justin profile y bins
	theCali.SetYBinStr( "CONSTANT 1200 -24 24" ); // full range ( 40 ps ) -- llpana default
    ////theCali.SetYBinStr( "CONSTANT 240 -12 12" ); // justin ( 100 ps )
    ////theCali.SetYBinStr( "CONSTANT 360 -12 12" ); // justin ( 67 ps )
    ////theCali.SetYBinStr( "CONSTANT 180 -6 6" ); // justin ( 67 ps )
    ////theCali.SetYBinStr( "CONSTANT 180 -3 3" ); // justin ( 33 ps )
    //--------------------------------------------------------
    // special 2d input file
    ////std::string inputfilename( "kucmsTimeCaliR17Plot2DFullRunFile.txt" );
	////std::string inputfilename( "kucmsTimeCaliR24FCCvRtTFile.txt");
    //std::string inputfilename( "ecal_config/kucmsTimeCaliR18ULTFile.txt");
   	///////std::string inputfilename( "ecal_config/kucmsTimeCaliR18ULMCTFile.txt");

	//-------------------------------------------------------
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_xa_pm9b720_r325" );// scale, cali, smear

	//theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_deg24fval_xa_pm9b1800_v501" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg24fprmt_xgs_pm3b240_v508" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg24fcc_xer_pm3b240_v508" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_deg25bc_xa_pm6b480_v616" );
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg24d_xa_pm6b480_v810" );
	//

	//theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_18UL_xa_pm24b1200_v0121" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_18UL_xa_pm24b1200_v0121" );

    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_17UL_xa_pm24b1200_v0121" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_17UL_xa_pm24b1200_v0121" );

    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_16UL_xa_pm24b1200_v0127" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_16UL_xa_pm24b1200_v0121" );
 
    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_22p_xa_pm24b1200_v0120" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_22p_xa_pm24b1200_v0120" );

    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_23p_xa_pm24b1200_v0115" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_23p_xa_pm24b1200_v0115" );

    //theCali.plot2dResolutionEGR( inputfilename, true, true, false, "_eg_24p_xa_pm24b1200_v0119" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_24p_xa_pm24b1200_v0119" );

	//theCali.setMCResTag("r2_ul16");
	//theCali.plot2dResolutionEGR( inputfilename, true, false, false, "_eg_18MC_xa_pm24b1200_v0116" );
    //theCali.plot2dResolutionEGR( inputfilename, true, false, true, "_eg_18MC_smUL16_xa_pm24b1200_v0201" );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // extended range resfit

	//std::string histName = "ResMap_315257_325172__DRO_Data_Hist_eg_18UL_HG_EB_xa_pm24b1200_v1209";
	//std::string histName = "ResMap_315257_325172__DRO_Data_Hist_eg_18UL_HG_EE_xa_pm24b1200_v1212";
	//std::string histName = "ResMap_296399_306460__DRO_Data_Hist_eg_17UL_HG_EB_xa_pm24b1200_v1212";
    //std::string histName = "ResMap_296399_306460__DRO_Data_Hist_eg_17UL_HG_EE_xa_pm24b1200_v1212";

	//std::string histName = "ResMap_315257_325172__DRO_Data_Hist_eg_18UL_LG_EB_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_315257_325172__DRO_Data_Hist_eg_18UL_LG_EE_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_296399_306460__DRO_Data_Hist_eg_17UL_LG_EB_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_296399_306460__DRO_Data_Hist_eg_17UL_LG_EE_xa_pm24b1200_v0106";

	//std::string histName = "ResMap_275600_283900__DRO_Data_Hist_eg_16UL_HG_EB_xa_pm24b1200_v0105";
    //std::string histName = "ResMap_275600_283900__DRO_Data_Hist_eg_16UL_LG_EB_xa_pm24b1200_v0105";
	//std::string histName = "ResMap_275600_283900__DRO_Data_Hist_eg_16UL_HG_EE_xa_pm24b1200_v0105";
    //std::string histName = "ResMap_275600_283900__DRO_Data_Hist_eg_16UL_LG_EE_xa_pm24b1200_v0105";

    //std::string histName = "ResMap_315257_325172__SRO_Data_Hist_eg_18UL_LG_EB_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_315257_325172__SRO_Data_Hist_eg_18UL_LG_EE_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_296399_306460__SRO_Data_Hist_eg_17UL_LG_EB_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_296399_306460__SRO_Data_Hist_eg_17UL_LG_EE_xa_pm24b1200_v0106";
    //std::string histName = "ResMap_275600_283900__SRO_Data_Hist_eg_16UL_LG_EB_xa_pm24b1200_v0105";
    //std::string histName = "ResMap_275600_283900__SRO_Data_Hist_eg_16UL_LG_EE_xa_pm24b1200_v0105";

	//std::string histName = "ResMap_315257_325172__DRO_Data_Hist_NoCali_eg_18UL_EB_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_315257_325172__SRO_Data_Hist_NoCali_eg_18UL_EB_xa_pm24b1200_v0108";
	//std::string histName = "ResMap_315257_325172__DRO_Data_Hist_NoCali_eg_18UL_EE_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_315257_325172__SRO_Data_Hist_NoCali_eg_18UL_EE_xa_pm24b1200_v0108";

    //std::string histName = "ResMap_296399_306460__DRO_Data_Hist_NoCali_eg_17UL_EB_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_296399_306460__SRO_Data_Hist_NoCali_eg_17UL_EB_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_296399_306460__DRO_Data_Hist_NoCali_eg_17UL_EE_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_296399_306460__SRO_Data_Hist_NoCali_eg_17UL_EE_xa_pm24b1200_v0108";

    //std::string histName = "ResMap_275600_283900__DRO_Data_Hist_NoCali_eg_16UL_EB_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_275600_283900__SRO_Data_Hist_NoCali_eg_16UL_EB_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_275600_283900__DRO_Data_Hist_NoCali_eg_16UL_EE_xa_pm24b1200_v0108";
    //std::string histName = "ResMap_275600_283900__SRO_Data_Hist_NoCali_eg_16UL_EE_xa_pm24b1200_v0108";

	//------------------------------------------------------------------------------------------------

	//std::string histName1 = "ResMap_315257_325172_DRO_Data_Hist_NoCali_eg_18UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_315257_325172_SRO_Data_Hist_NoCali_eg_18UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_315257_325172_SRO_Data_Hist_NoCali_eg_18UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_315257_325172_DRO_Data_Hist_NoCali_eg_18UL_xa_pm24b1200_v0121_EE";

    //std::string histName1 = "ResMap_315257_325172_DRO_Data_Hist_eg_18UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_315257_325172_SRO_Data_Hist_eg_18UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_315257_325172_SRO_Data_Hist_eg_18UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_315257_325172_DRO_Data_Hist_eg_18UL_xa_pm24b1200_v0121_EE";

    //std::string histName1 = "ResMap_296399_306460_DRO_Data_Hist_NoCali_eg_17UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_296399_306460_SRO_Data_Hist_NoCali_eg_17UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_296399_306460_DRO_Data_Hist_NoCali_eg_17UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_296399_306460_SRO_Data_Hist_NoCali_eg_17UL_xa_pm24b1200_v0121_EE";

    //std::string histName1 = "ResMap_296399_306460_DRO_Data_Hist_eg_17UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_296399_306460_SRO_Data_Hist_eg_17UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_296399_306460_DRO_Data_Hist_eg_17UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_296399_306460_SRO_Data_Hist_eg_17UL_xa_pm24b1200_v0121_EE";

    //std::string histName1 = "ResMap_275600_283900_DRO_Data_Hist_NoCali_eg_16UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_275600_283900_SRO_Data_Hist_NoCali_eg_16UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_275600_283900_DRO_Data_Hist_NoCali_eg_16UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_275600_283900_SRO_Data_Hist_NoCali_eg_16UL_xa_pm24b1200_v0121_EE";

    //std::string histName1 = "ResMap_275600_283900_DRO_Data_Hist_eg_16UL_xa_pm24b1200_v0121_EB";
    //std::string histName2 = "ResMap_275600_283900_SRO_Data_Hist_eg_16UL_xa_pm24b1200_v0121_EB";
    //std::string histName3 = "ResMap_275600_283900_DRO_Data_Hist_eg_16UL_xa_pm24b1200_v0121_EE";
    //std::string histName4 = "ResMap_275600_283900_SRO_Data_Hist_eg_16UL_xa_pm24b1200_v0121_EE";

	//std::string histName1 = "ResMap_355794_362760__DRO_Data_Hist_eg_22p_EB_xa_pm24b1200_v0112";
    //std::string histName2 = "ResMap_355794_362760__SRO_Data_Hist_eg_22p_EB_xa_pm24b1200_v0112";
    //std::string histName3 = "ResMap_355794_362760__DRO_Data_Hist_eg_22p_EE_xa_pm24b1200_v0112";
    //std::string histName4 = "ResMap_355794_362760__SRO_Data_Hist_eg_22p_EE_xa_pm24b1200_v0112";

	//std::string histName1 = "ResMap_355794_362760_DRO_Data_Hist_NoCali_eg_22p_xa_pm24b1200_v0120_EB";
    //std::string histName2 = "ResMap_355794_362760_DRO_Data_Hist_NoCali_eg_22p_xa_pm24b1200_v0120_EE";
    //std::string histName3 = "ResMap_355794_362760_SRO_Data_Hist_NoCali_eg_22p_xa_pm24b1200_v0120_EB";
    //std::string histName4 = "ResMap_355794_362760_SRO_Data_Hist_NoCali_eg_22p_xa_pm24b1200_v0120_EE";

    //std::string histName1 = "ResMap_366365_372415_DRO_Data_Hist_eg_23p_xa_pm24b1200_v0115_EB";
    //std::string histName2 = "ResMap_366365_372415_SRO_Data_Hist_eg_23p_xa_pm24b1200_v0115_EB";
    //std::string histName3 = "ResMap_366365_372415_DRO_Data_Hist_eg_23p_xa_pm24b1200_v0115_EE";
    //std::string histName4 = "ResMap_366365_372415_SRO_Data_Hist_eg_23p_xa_pm24b1200_v0115_EE";

    //std::string histName1 = "ResMap_366365_372415_DRO_Data_Hist_NoCali_eg_23p_xa_pm24b1200_v0115_EB";
    //std::string histName2 = "ResMap_366365_372415_SRO_Data_Hist_NoCali_eg_23p_xa_pm24b1200_v0115_EB";
    //std::string histName3 = "ResMap_366365_372415_DRO_Data_Hist_NoCali_eg_23p_xa_pm24b1200_v0115_EE";
    //std::string histName4 = "ResMap_366365_372415_SRO_Data_Hist_NoCali_eg_23p_xa_pm24b1200_v0115_EE";

    //std::string histName1 = "ResMap_0_999999_DRO_Data_Hist_NoCali_eg_18MC_xa_pm24b1200_v0116_EB";
    //std::string histName2 = "ResMap_0_999999_DRO_Data_Hist_NoCali_eg_18MC_xa_pm24b1200_v0116_EE";
    //std::string histName3 = "ResMap_0_999999_SRO_Data_Hist_NoCali_eg_18MC_xa_pm24b1200_v0116_EB";
    //std::string histName4 = "ResMap_0_999999_SRO_Data_Hist_NoCali_eg_18MC_xa_pm24b1200_v0116_EE";

    //std::string histName1 = "ResMap_378971_387721_DRO_Data_Hist_NoCali_eg_24p_xa_pm24b1200_v0119_EB";
    //std::string histName2 = "ResMap_378971_387721_DRO_Data_Hist_NoCali_eg_24p_xa_pm24b1200_v0119_EE";
    //std::string histName3 = "ResMap_378971_387721_SRO_Data_Hist_NoCali_eg_24p_xa_pm24b1200_v0119_EB";
    //std::string histName4 = "ResMap_378971_387721_SRO_Data_Hist_NoCali_eg_24p_xa_pm24b1200_v0119_EE";

    //std::string histName1 = "ResMap_378971_387721_DRO_Data_Hist_eg_24p_xa_pm24b1200_v0119_EB";
    //std::string histName2 = "ResMap_378971_387721_DRO_Data_Hist_eg_24p_xa_pm24b1200_v0119_EE";
    //std::string histName3 = "ResMap_378971_387721_SRO_Data_Hist_eg_24p_xa_pm24b1200_v0119_EB";
    //std::string histName4 = "ResMap_378971_387721_SRO_Data_Hist_eg_24p_xa_pm24b1200_v0119_EE";

	//
	//theCali.Set2DResTFileName( "egres_DEGPD_AODSIM_R24_Fprompt_v3_14011_v12_resplots.root" );
    //theCali.Set2DResTFileName( "egres_DEGPD_AODSIM_R24_ECAL_v3_14011_v12_resplots.root" );
	//

	//std::string histName = "ResMap_275600_283900_DRO_Data_Hist_eg_16UL_xa_pm24b1200_v0127_EB";
    //std::string histName = "ResMap_0_999999_DRO_Data_Hist_NoCali_Smeared_eg_18MC_xa_pm24b1200_v0130_EB";
    //std::string histName = "ResMap_0_999999_DRO_Data_Hist_NoCali_Smeared_eg_18MC_xa_pm24b1200_v0129_EE";
    //std::string histName = "ResMap_0_999999_DRO_Data_Hist_NoCali_Smeared_eg_18MC_smUL17_xa_pm24b1200_v0201_EE";
    std::string histName = "ResMap_0_999999_DRO_Data_Hist_NoCali_Smeared_eg_18MC_smUL16_xa_pm24b1200_v0201_EE";

    //theCali.setLowEnergy( false ); // true is default

	////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800" ); // : xa
    ////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1000 1200 1800" ); // : ecal
    ////theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 30.0 40.0 60.0 120.0" ); // xea
    ////theCali.SetXBinStr( "VARIABLE 5 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200" ); // xb : * w/LE
    ////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 1800" ); // xj
    ////theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 60 80 100 125 150 175 200 225 250 300 400 600" ); // : xr
    ////theCali.SetXBinStr( "VARIABLE 1.0 2.0 3.0 4.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0 25.0 30.0 40.0 50.0 70.0 90.0" ); // xer
    ////theCali.SetXBinStr( "VARIABLE 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1200 1600 2200" ); // : xgs ecal
	////theCali.SetXBinStr( "VARIABLE 5 10 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 1200 1800 2200" );  
	/////////////////////////////// : xa ( xaext+ 2200 )

	////theCali.SetXBinStr( "VARIABLE 5 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1600" ); // : xa 2018UL DRO
    ////theCali.SetXBinStr( "VARIABLE 5 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1600" ); // : xa llpana LG resfit
	theCali.SetXBinStr( "VARIABLE 5 15 20 25 30 40 50 75 100 125 150 175 200 225 250 300 400 600 800 1600 3200" ); // : xa llpana HG resfit

    theCali.load2DResHist( histName );
	theCali.doResTimeFit( histName );

	/*
    theCali.load2DResHist( histName1 );
    theCali.doResTimeFit( histName1 );
    theCali.load2DResHist( histName2 );
    theCali.doResTimeFit( histName2 );
    theCali.load2DResHist( histName3 );
    theCali.doResTimeFit( histName3 );
    theCali.load2DResHist( histName4 );
    theCali.doResTimeFit( histName4 );
	*/

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// plotting of mean time by run with calibraton : filename, start run, end run, usecali ( defaults to usecali true )
	//theCali.plotMeanRunTimeEGR( inputfilename, 315257, 325172 ); // R2 2018 UL 315257_325172 // R2 2017 UL 296399_306460
    //theCali.plotMeanRunTimeEGR( inputfilename, 315257, 325172, false );
    ////theCali.plotMeanRunTimeEGR( inputfilename, 296399, 307554, true );
	////theCali.makeTTDiffMaps();// make trigger tower diffrence maps
	
	/////////////////////////////////////////////////////////////////////

	//KUCMS_TimeCalibration theCali;

	//theCali.makeTTDriftMaps( "r2_ul18", 315257, 325172 );	
    //theCali.makeTTDriftMaps( "r2_ul18", 315000, 318000 );

    /////////////////////////////////////////////////////////////////////
    //-------------------------------------------------------------------
	
	std::cout << " -- Thats All Folks !!!!!!!!! " << std::endl;

    return 1;

}//<<>>int main ( int argc, char *argv[] )
