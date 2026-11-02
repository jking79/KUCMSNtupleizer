//////////////////////////////////////////////////////////////////////
// Batch driver for fitting every resolution TH2F in a res2dPlots file.
//////////////////////////////////////////////////////////////////////

#include "KUCMS_TimeCalibration.hh"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TKey.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

namespace {

struct FitRecord {
  double noise=NAN, enoise=NAN, constant=NAN, econstant=NAN;
  double stochastic=NAN, estochastic=NAN, fractionalScatter=NAN, efractionalScatter=NAN;
  double objective=NAN, chi2=NAN, ndf=NAN, chi2ndf=NAN, probability=NAN;
  double fitStatus=NAN, covarianceStatus=NAN, fitPoints=NAN, parameterAtLimit=NAN, valid=NAN;
};

std::string csvQuote( const std::string& value ){
  std::string escaped = value;
  for( std::string::size_type position=0; (position=escaped.find('"',position)) != std::string::npos; position += 2 )
    escaped.insert(position,"\"");
  return "\"" + escaped + "\"";
}

std::string csvNumber( double value ){
  if( not std::isfinite(value) ) return "";
  std::ostringstream stream; stream << std::setprecision(10) << value; return stream.str();
}

double qualityValue( TH1* quality, const std::string& label ){
  if( not quality ) return NAN;
  for( int bin=1; bin<=quality->GetNbinsX(); bin++ )
    if( label == quality->GetXaxis()->GetBinLabel(bin) ) return quality->GetBinContent(bin);
  return NAN;
}

FitRecord readFitRecord( TFile* file, const std::string& functionName,
                         const std::string& qualityName, bool scatter ){
  FitRecord record;
  if( not file ) return record;
  auto function = dynamic_cast<TF1*>(file->Get(functionName.c_str()));
  auto quality = dynamic_cast<TH1*>(file->Get(qualityName.c_str()));
  if( function ){
    if( function->GetNpar()>0 ){ record.noise=function->GetParameter(0); record.enoise=function->GetParError(0); }
    if( function->GetNpar()>1 ){ record.constant=function->GetParameter(1); record.econstant=function->GetParError(1); }
    if( function->GetNpar()>2 ){ record.stochastic=function->GetParameter(2); record.estochastic=function->GetParError(2); }
    if( scatter && function->GetNpar()>3 ){ record.fractionalScatter=function->GetParameter(3); record.efractionalScatter=function->GetParError(3); }
  }
  record.fitStatus=qualityValue(quality,"fitStatus"); record.covarianceStatus=qualityValue(quality,"covStatus");
  record.chi2=qualityValue(quality,"chi2"); record.ndf=qualityValue(quality,"ndf");
  record.chi2ndf=qualityValue(quality,"chi2ndf"); record.probability=qualityValue(quality,"prob");
  record.parameterAtLimit=qualityValue(quality,"parameterAtLimit"); record.fitPoints=qualityValue(quality,"nFitPoints");
  record.valid=qualityValue(quality,"valid");
  if( scatter ){
    record.fractionalScatter=qualityValue(quality,"fractionalScatter");
    record.efractionalScatter=qualityValue(quality,"fractionalScatterError");
    record.objective=qualityValue(quality,"neg2LogLikelihood");
  }
  return record;
}

void writeRecordHeader( std::ofstream& output, const std::string& prefix ){
  output << prefix << "_noise," << prefix << "_noise_error," << prefix << "_constant," << prefix << "_constant_error,"
    << prefix << "_stochastic," << prefix << "_stochastic_error," << prefix << "_fractional_scatter,"
    << prefix << "_fractional_scatter_error," << prefix << "_objective," << prefix << "_chi2," << prefix << "_ndf,"
    << prefix << "_chi2_ndf," << prefix << "_probability," << prefix << "_fit_status," << prefix << "_covariance_status,"
    << prefix << "_fit_points," << prefix << "_parameter_at_limit," << prefix << "_valid";
}

void writeRecord( std::ofstream& output, const FitRecord& record ){
  output << csvNumber(record.noise) << ',' << csvNumber(record.enoise) << ',' << csvNumber(record.constant) << ','
    << csvNumber(record.econstant) << ',' << csvNumber(record.stochastic) << ',' << csvNumber(record.estochastic) << ','
    << csvNumber(record.fractionalScatter) << ',' << csvNumber(record.efractionalScatter) << ',' << csvNumber(record.objective) << ','
    << csvNumber(record.chi2) << ',' << csvNumber(record.ndf) << ',' << csvNumber(record.chi2ndf) << ',' << csvNumber(record.probability) << ','
    << csvNumber(record.fitStatus) << ',' << csvNumber(record.covarianceStatus) << ',' << csvNumber(record.fitPoints) << ','
    << csvNumber(record.parameterAtLimit) << ',' << csvNumber(record.valid);
}

struct HistogramMetadata { std::string runStart,runEnd,pairType,datasetTag,gainSelection,calibrationState,region; };

HistogramMetadata parseMetadata( const std::string& name ){
  HistogramMetadata metadata;
  static const std::regex pattern("^ResMap_([0-9]+)_([0-9]+)_(DRO|SRO|ZEE)_Data_Hist_(.*)_(EB|EE)$");
  std::smatch match;
  if( std::regex_match(name,match,pattern) ){
    metadata.runStart=match[1]; metadata.runEnd=match[2]; metadata.pairType=match[3]; metadata.datasetTag=match[4]; metadata.region=match[5];
    for( const std::string marker : {"AllGains_","Gain1Only_","NoCali_"} ){
      if( metadata.datasetTag.rfind(marker,0) == 0 ) metadata.datasetTag.erase(0,marker.size());
    }
    metadata.gainSelection=name.find("Gain1Only")!=std::string::npos ? "Gain1Only" : "AllGains";
    metadata.calibrationState=name.find("NoCali")!=std::string::npos ? "NoCali" : "Cali";
  }
  return metadata;
}

bool writeBinCsv( const std::string& csvPath, const std::string& histName, const std::string& rootPath ){
  TFile* file=TFile::Open(rootPath.c_str(),"READ");
  if( not file || file->IsZombie() ){ if(file){file->Close();delete file;} return false; }
  auto source=dynamic_cast<TH2F*>(file->Get(histName.c_str()));
  if( not source ){ file->Close();delete file;return false; }
  auto get=[&](const std::string& suffix){return dynamic_cast<TH1*>(file->Get((histName+suffix).c_str()));};
  TH1* occ=get("_occ"); TH1* mu=get("_mu"); TH1* rms=get("_rms"); TH1* sigma=get("_sigma");
  TH1* sigmaModelErr=get("_sigmaModelErr"); TH1* sigmaModelVariations=get("_sigmaModelVariations");
  TH1* sigmaRange1=get("_sigmaFitRange1"); TH1* sigmaRange4=get("_sigmaFitRange4");
  TH1* fitStatus=get("_gaussFitStatus"); TH1* covStatus=get("_gaussCovStatus"); TH1* fitValid=get("_gaussFitValid");
  TH1* fitStatusRange1=get("_gaussFitStatusRange1"); TH1* covStatusRange1=get("_gaussCovStatusRange1"); TH1* fitValidRange1=get("_gaussFitValidRange1");
  TH1* fitStatusRange4=get("_gaussFitStatusRange4"); TH1* covStatusRange4=get("_gaussCovStatusRange4"); TH1* fitValidRange4=get("_gaussFitValidRange4");
  TH1* sigmaAtLimit=get("_gaussSigmaAtLimit"); TH1* sigmaNscInput=get("_sigmaNSCFitInput");
  TH1* scatterNscInput=get("_sigmaFractionalScatterNSCFitInput"); TH1* scatterNcInput=get("_sigmaFractionalScatterNCFitInput");
  auto metadata=parseMetadata(histName); std::ofstream output(csvPath);
  if( not output.is_open() ){file->Close();delete file;return false;}
  output << "hist_name,run_start,run_end,pair_type,dataset_tag,gain_selection,calibration_state,region,x_bin,x_low,x_high,x_center,input_x_entries,"
    << "occupancy,occupancy_error,mu,mu_error,rms,sigma,sigma_error,sigma_model_error,sigma_model_variations,sigma_fit_range1,sigma_fit_range4,"
    << "gaussian_fit_status,gaussian_covariance_status,gaussian_fit_valid,gaussian_fit_status_low_range,gaussian_covariance_status_low_range,gaussian_fit_valid_low_range,"
    << "gaussian_fit_status_high_range,gaussian_covariance_status_high_range,gaussian_fit_valid_high_range,sigma_at_limit,"
    << "sigma_statistical,sigma_statistical_error,sigma_model_aware_nsc,sigma_model_aware_nsc_error,sigma_fitted_scatter_nsc,sigma_fitted_scatter_nsc_error,"
    << "sigma_fitted_scatter_nc,sigma_fitted_scatter_nc_error\n";
  auto value=[](TH1* histogram,int bin){return histogram?histogram->GetBinContent(bin):NAN;};
  auto error=[](TH1* histogram,int bin){return histogram?histogram->GetBinError(bin):NAN;};
  for( int bin=1;bin<=source->GetNbinsX();bin++ ){
    double inputEntries=0; for( int ybin=1;ybin<=source->GetNbinsY();ybin++ ) inputEntries+=source->GetBinContent(bin,ybin);
    output << csvQuote(histName) << ',' << metadata.runStart << ',' << metadata.runEnd << ',' << metadata.pairType << ',' << csvQuote(metadata.datasetTag) << ','
      << metadata.gainSelection << ',' << metadata.calibrationState << ',' << metadata.region << ',' << bin << ','
      << csvNumber(source->GetXaxis()->GetBinLowEdge(bin)) << ',' << csvNumber(source->GetXaxis()->GetBinUpEdge(bin)) << ',' << csvNumber(source->GetXaxis()->GetBinCenter(bin)) << ','
      << csvNumber(inputEntries) << ',' << csvNumber(value(occ,bin)) << ',' << csvNumber(error(occ,bin)) << ',' << csvNumber(value(mu,bin)) << ',' << csvNumber(error(mu,bin)) << ','
      << csvNumber(value(rms,bin)) << ',' << csvNumber(value(sigma,bin)) << ',' << csvNumber(error(sigma,bin)) << ',' << csvNumber(value(sigmaModelErr,bin)) << ','
      << csvNumber(value(sigmaModelVariations,bin)) << ',' << csvNumber(value(sigmaRange1,bin)) << ',' << csvNumber(value(sigmaRange4,bin)) << ','
      << csvNumber(value(fitStatus,bin)) << ',' << csvNumber(value(covStatus,bin)) << ',' << csvNumber(value(fitValid,bin)) << ','
      << csvNumber(value(fitStatusRange1,bin)) << ',' << csvNumber(value(covStatusRange1,bin)) << ',' << csvNumber(value(fitValidRange1,bin)) << ','
      << csvNumber(value(fitStatusRange4,bin)) << ',' << csvNumber(value(covStatusRange4,bin)) << ',' << csvNumber(value(fitValidRange4,bin)) << ',' << csvNumber(value(sigmaAtLimit,bin)) << ','
      << csvNumber(value(sigma,bin)) << ',' << csvNumber(error(sigma,bin)) << ',' << csvNumber(value(sigmaNscInput,bin)) << ',' << csvNumber(error(sigmaNscInput,bin)) << ','
      << csvNumber(value(scatterNscInput,bin)) << ',' << csvNumber(error(scatterNscInput,bin)) << ',' << csvNumber(value(scatterNcInput,bin)) << ',' << csvNumber(error(scatterNcInput,bin)) << '\n';
  }
  file->Close();delete file;return true;
}

void writeSummaryHeader( std::ofstream& output ){
  output << "hist_name,output_root_file,bin_csv_file,run_start,run_end,pair_type,dataset_tag,gain_selection,calibration_state,region,processing_status,";
  writeRecordHeader(output,"model_aware_nsc");output << ',';writeRecordHeader(output,"model_aware_nc");output << ',';
  writeRecordHeader(output,"stat_only_nsc");output << ',';writeRecordHeader(output,"stat_only_nc");output << ',';
  writeRecordHeader(output,"fitted_scatter_nsc");output << ',';writeRecordHeader(output,"fitted_scatter_nc");output << '\n';
}

void writeSummaryRow( std::ofstream& output,const std::string& histName,const std::string& rootPath,const std::string& csvPath,
                      const FitRecord& modelNsc,const FitRecord& modelNc,const FitRecord& statNsc,const FitRecord& statNc,const FitRecord& scatterNsc,const FitRecord& scatterNc,bool success ){
  auto metadata=parseMetadata(histName);
  output << csvQuote(histName) << ',' << csvQuote(rootPath) << ',' << csvQuote(csvPath) << ',' << metadata.runStart << ',' << metadata.runEnd << ',' << metadata.pairType << ','
    << csvQuote(metadata.datasetTag) << ',' << metadata.gainSelection << ',' << metadata.calibrationState << ',' << metadata.region << ',' << (success?"success":"failed") << ',';
  writeRecord(output,modelNsc);output << ',';writeRecord(output,modelNc);output << ',';writeRecord(output,statNsc);output << ',';writeRecord(output,statNc);output << ',';writeRecord(output,scatterNsc);output << ',';writeRecord(output,scatterNc);output << '\n';
}

} // namespace

int main( int argc,char* argv[] ){
  const std::string inputPath=argc>1?argv[1]:"../../res2dPlotsTFile.root";
  const std::string outputDir=argc>2?argv[2]:"resfit_outputs";
  const std::string nameFilter=argc>3?argv[3]:"";
  if( argc>4 ){std::cerr<<"Usage: "<<argv[0]<<" [input.root] [output_directory] [name_substring_filter]"<<std::endl;return EXIT_FAILURE;}
  std::error_code error;std::filesystem::create_directories(outputDir,error);
  if(error){std::cerr<<"ERROR: cannot create output directory "<<outputDir<<": "<<error.message()<<std::endl;return EXIT_FAILURE;}
  TFile* input=TFile::Open(inputPath.c_str(),"READ");
  if(not input||input->IsZombie()){std::cerr<<"ERROR: cannot open input ROOT file: "<<inputPath<<std::endl;if(input){input->Close();delete input;}return EXIT_FAILURE;}
  KUCMS_TimeCalibration calibration(false,false,true);calibration.SetResFitOutputDir(outputDir);
  std::ofstream summary((std::filesystem::path(outputDir)/"resolution_fit_summary.csv").string());
  if(not summary.is_open()){std::cerr<<"ERROR: cannot create resolution_fit_summary.csv"<<std::endl;input->Close();delete input;return EXIT_FAILURE;}
  writeSummaryHeader(summary);int discovered=0,processed=0,failed=0,skipped=0;TIter next(input->GetListOfKeys());
  while(auto key=dynamic_cast<TKey*>(next())){
    if(std::string(key->GetClassName())!="TH2F")continue;const std::string histName=key->GetName();
    if(histName.find("Amp_v_E")!=std::string::npos){skipped++;continue;}
    if(histName.find("_DRO_Data_Hist_")==std::string::npos&&histName.find("_SRO_Data_Hist_")==std::string::npos&&histName.find("_ZEE_Data_Hist_")==std::string::npos){skipped++;continue;}
    if(not nameFilter.empty()&&histName.find(nameFilter)==std::string::npos)continue;discovered++;
    auto source=dynamic_cast<TH2F*>(key->ReadObj());auto clone=source?dynamic_cast<TH2F*>(source->Clone(histName.c_str())):nullptr;if(source)delete source;
    if(not clone){failed++;continue;}clone->SetDirectory(nullptr);std::cout<<"["<<discovered<<"] Fitting "<<histName<<std::endl;
    const kucms_SigmaFitResult result=calibration.runTimeFitter(clone);delete clone;
    const std::filesystem::path rootPath=std::filesystem::path(outputDir)/(histName+"_NSC_resfit.root");const std::filesystem::path csvPath=std::filesystem::path(outputDir)/(histName+"_NSC_resfit_resolution_bin_results.csv");
    const bool binCsvOk=writeBinCsv(csvPath.string(),histName,rootPath.string());TFile* fitFile=TFile::Open(rootPath.string().c_str(),"READ");const bool fitFileOk=fitFile&&not fitFile->IsZombie();
    FitRecord modelNsc,modelNc,statNsc,statNc,scatterNsc,scatterNc;
    if(fitFileOk){modelNsc=readFitRecord(fitFile,histName+"_sigmaNSCFitInput_fit",histName+"_NSCFitQuality",false);modelNc=readFitRecord(fitFile,histName+"_sigmaNSCFitInput_fit_NC",histName+"_NCFitQuality",false);statNsc=readFitRecord(fitFile,histName+"_sigma_fit_StatOnly_NSC",histName+"_StatOnlyNSCFitQuality",false);statNc=readFitRecord(fitFile,histName+"_sigma_fit_StatOnly_NC",histName+"_StatOnlyNCFitQuality",false);scatterNsc=readFitRecord(fitFile,histName+"_sigma_fit_FractionalScatter_NSC",histName+"_FractionalScatterNSCFitQuality",true);scatterNc=readFitRecord(fitFile,histName+"_sigma_fit_FractionalScatter_NC",histName+"_FractionalScatterNCFitQuality",true);fitFile->Close();delete fitFile;}else if(fitFile){fitFile->Close();delete fitFile;}
    const bool success=result.valid&&binCsvOk&&fitFileOk;writeSummaryRow(summary,histName,rootPath.string(),csvPath.string(),modelNsc,modelNc,statNsc,statNc,scatterNsc,scatterNc,success);if(success)processed++;else failed++;
  }
  input->Close();delete input;summary.close();std::cout<<"Batch resolution fitting complete: discovered="<<discovered<<" successful="<<processed<<" failed="<<failed<<" skipped="<<skipped<<std::endl;
  return failed==0&&discovered>0?EXIT_SUCCESS:EXIT_FAILURE;
}
