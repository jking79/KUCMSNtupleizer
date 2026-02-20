#MINIAOD sample lookup
mc_samples_2018 = [
"/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
"/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM",
]

data_samples_2018 = [
"/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD",
"/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD",
"/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD",
"/EGamma/Run2018D-12Nov2019_UL2018-v2/MINIAOD",
"/JetHT/Run2018A-15Feb2022_UL2018-v2/MINIAOD",
"/JetHT/Run2018B-15Feb2022_UL2018-v1/MINIAOD",
"/JetHT/Run2018C-15Feb2022_UL2018-v1/MINIAOD",
"/JetHT/Run2018D-15Feb2022_UL2018-v1/MINIAOD"
]


mc_samples_2017 = [

]

data_samples_2017 = [
"/DoubleEG/Run2017C-09Aug2019_UL2017-v1/MINIAOD",    

]

MC_PROCS = ["QCD","GJets","TTXJets","WJets","ZJets","DiPhotonJetBox","TTJets","DYJets"]
DATA_NAMES = ["MET","JetMET","EGamma", "DoubleEG","JetHT"]
GLOBALTAGS_DICT = {"2016":"106X_dataRun2_v27", "2017":'106X_dataRun2_v20',"2018":"106X_dataRun2_v36","2022":"140X_dataRun3_v17","2023":"140X_dataRun3_v17","RunIISummer20UL18RECO":"106X_upgrade2018_realistic_v11_L1v1","RunIIAutumn18DRPremix":"94X_mc2017_realistic_v11"}
JSON_DICT = {"2016" : "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt", "2017":'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',"2018":"Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.json","2022":"Cert_Collisions2022_355100_362760_Golden.json","2023":"Cert_Collisions2023_366442_370790_Golden.json"}
MINIAOD_SAMPLES = {"mc_2018" : mc_samples_2018, "data_2018" : data_samples_2018, "mc_2017" : mc_samples_2017,  "data_2017" : data_samples_2017}
