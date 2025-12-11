
#PROC=DiPJBox
#PROC=GJets
#PROC=QCD
PROC=SMS
#PROC=TTXJets
#PROC=WJets
#PROC=ZJets

#python3 doHadd_Multi.py --dir ./testOutput --proc ${PROC}
#source haddScripts/doHadd_${PROC}_.sh

python3 doXRDCP.py --dir ./testOutput --proc ${PROC} --targetdir skims_v44 #--dryRun
