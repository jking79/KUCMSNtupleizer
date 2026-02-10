
#generateSubmission.py [-h] [--directory DIRECTORY] [--inputSample {GJets,QCD,MET,gogoG}] [--year {2017,2018,2022}] [--HT HT] [--era ERA] [--mGl MGL] [--mN2 MN2]
#  [--mN1 MN1] [--output OUTPUT] [--split SPLIT] [--maxnevts MAXNEVTS](maxevents per output file, runs default over all files) [--verbosity VERBOSITY] [--genSigPerfect] [--noSVorPho]
#--maxnfiles max amount of files to run over
#{DiPJBox,DTBoson,GJets,TTXJets,QCD,WJets,ZJets,gogoG,gogoZ,sqsqG,MET,EGamma}
#--max_mat maxjobs materialized in schedd at one time
#python3 generateSubmission.py --directory testOutput --inputSample GJets --year 2018 --maxnfiles 1 --max_mat 1 --request_memory 4096 --maxnevts 200
#python3 generateSubmission.py --directory testOutput --inputSample gogoG --year 2022 --maxnevts 200 --max_mat 100 --mGl 2000 --mN2 1500 --mN1 1000 --request_memory 4096

#python3 generateSubmission.py --directory testOutput --inputSample QCD --year 2018 --max_mat 100

#python3 generateSubmission.py --directory testOutput --inputSample WJets --year 2018 --max_mat 100 
#python3 generateSubmission.py --directory testOutput --inputSample ZJets --year 2018 --max_mat 100 

#python3 generateSubmission.py --directory testOutput --inputSample TTXJets --year 2018 --max_mat 100
#python3 generateSubmission.py --directory testOutput --inputSample DiPJBox --year 2018 --max_mat 100
#python3 generateSubmission.py --directory testOutput --inputSample GJets --year 2018 --max_mat 100 

#python3 generateSubmission.py --directory testOutput --inputSample gogoG --year 2022 --max_mat 10 
#python3 generateSubmission.py --directory testOutput --inputSample gogoZ --year 2022 --max_mat 15

#python3 generateSubmission.py --directory testOutput --inputSample gogoGZ --year 2022 --max_mat 10 --mGl 2300 --mN2 2200 --mN1 2100
python3 generateSubmission.py --directory testOutput --inputSample gogoGZ --year 2022 --max_mat 10



#commands for ld .so bhc when doin interactive
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/uscms/home/mlazarov/nobackup/CMSSW_13_3_3/src/BayesianClustering/lib



