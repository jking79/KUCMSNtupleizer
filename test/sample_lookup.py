#MINIAOD sample lookup

mc_samples_2018 = [
]

data_samples_2018 = [
"/JetHT/Run2018A-15Feb2022_UL2018-v2/MINIAOD",
"/JetHT/Run2018B-15Feb2022_UL2018-v1/MINIAOD"
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
