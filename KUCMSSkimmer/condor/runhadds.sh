
#PROC=DiPJBox
#PROC=GJets
#PROC=QCD
##PROC=SMS
#PROC=TTXJets
PROC=WJets
#PROC=ZJets
TARGETDIR=skims_v46

python3 doHadd_Multi.py --dir ./testOutput --proc ${PROC}
source haddScripts/doHadd_${PROC}_.sh

python3 doXRDCP.py --dir ./testOutput --proc ${PROC} --targetdir ${TARGETDIR} #--dryRun
