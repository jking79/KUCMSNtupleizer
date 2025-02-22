Command line interface for DAS (dasgoclient)
The command line interface for DAS called dasgoclient is available from the DAS web page
￼
 using the "CLI" link in the upper menu. It is located at /cvmfs/cms.cern.ch/common/dasgoclient and its usage is quite trivial. First, one needs to create a proxy for using it by the command voms-proxy-init -voms cms -rfc. The understand the usage of dasgoclient, one can issue the command dasgoclient -help. The output should look something like:
￼
Hide result...
Usage: dasgoclient [options]
  -aggregate
       aggregate results across all data-services
  -daskeys
       Show supported DAS keys
  -dasmaps string
       Specify location of dasmaps
  -examples
       Show examples of supported DAS queries
  -exitCodes
       Show DAS error codes
  -format string
       Compatibility option with python das_client, use json to get das_client behavior
  -funcProfile string
       Specify location of function profile file
  -host string
       Specify hostname to talk to (default "https://cmsweb.cern.ch")
  -idx int
       Compatibility option with python das_client
  -json
       Return results in JSON data-format
  -limit int
       Compatibility option with python das_client
  -profileMode string
       enable profiling mode, one of [cpu, mem, block]
  -query string
       DAS query to run
  -sep string
       Separator to use (default " ")
  -threshold int
       Compatibility option with python das_client, has no effect
  -timeout int
       Timeout for url call
  -token string
       Specify location of token file
  -unique
       Sort results and return unique list
  -urlQueuelimit int
       url queue limit (number of concurrent calls) (default 100)
  -urlRetry int
       urlRetry for url call (default 3)
  -verbose int
       Verbose level, support 0,1,2
  -version
       Show version
Thus, if you want to run your query you'll type (we'll use the same query examples as shown above):
dasgoclient -query="dataset=/Zee*/*/*"
The output would look something like:
￼
Hide result...
/Zee/StoreResults-7TeV-Ele15-314-SUSYPAT-V00-04-12-v2-c44ea04b24888fabe12d6ce1aba555ef/USER
/Zee/StoreResults-7TeV-Photon20-314-SUSYPAT-V00-04-12-v2-c44ea04b24888fabe12d6ce1aba555ef/USER
/Zee/StoreResults-EWKSKIMEMET-946a431cc6a46f4967712f847d809414/USER


Using DBS python client
DAS is a tool which aggregates data from several sources: DBS, Rucio, ReqMgr, SiteDB etc. But not all details of the information stored in those DB's is available, nor is the query as efficient as asking directly one of those services. Therefore DAS should be the first choice when looking for dataset informations, but sophisticated users that find the details or the performance inadequate to their needs can query DBS directly via its python client API. Instructions, examples and guidelines are in this twiki.

Accessing Remote Samples For interactive testing
The ability to access remote files (i.e. located at some Tier2) of various samples is essential to users for interactive testing and debugging. A remote file can be either copied to a local space (e.g. desktop/laptop) or directly opened inside cmsRun, using the Xrootd Service. Please refer to the dedicated chapter in this workbook: Using Xrootd Service (AAA) for Remote Data Access.

Finding existing MC samples for various physics processes
A list of MC samples requested in the latest production campaign can be found at MC co-ordination twiki. The collision data (MINI)AOD from the last year and all (MINI)AODSIM samples from the past couple of production campaigns will always be available at some disk site. If a sample is popular in CRAB, Dynamic Data Management (DDM) team at CMS distributes replicas of such hot samples via an automatic procedure which removes extra-copies of unused datasets.
As an example, to search for samples corresponding to RunII Summer16 MINIAODSIM production campaign one has to do a dump query on DAS like:
dataset=/*/RunIISummer16DR80X*/MINIAODSIM

Availability of Samples
The recent collision data and MC samples are always on the disk at some site. If something is not found on disk or is already archived to tape, one can file a ticket at the JIRA link

Release Validation (CMS.RelVal) samples
As new releases are integrated, and readied for large scale MC production, or data reprocessing, CMS goes through a process referred to as "Release Validation" (CMS.RelVal). As part of that, the Data Operations team makes a variety of samples with that release at small scale. These CMS.RelVal samples are often your best opportunity to develop analysis code for a new release, as they are the first to appear.

* There are CMS.RelVal samples for all major releases.
* These samples have been produced to validate CMSSW pre-releases and releases and the production workflow.
* In general you should run on these with the release with which they were produced (in particular for K_L_M_preX releases)
* One can easily find CMS.RelVal samples using DAS interface. For example, to find ttbar samples your query on DAS will look like:
dataset=/*RelValTTbar*/*/*

dasgoclient -query=
"dataset=“

globalTag=94X_dataRun2_ReReco_EOY17_v1
"/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"


            ['/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'],
            ['/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v2/AODSIM'],
            ['/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            ['/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],
            ['/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'],

            trial          = 'gammares_cali'

            # MC 2017 aod
            config.JobType.pyCfgParams = ['globalTag=94X_mc2017_realistic_v12', 'outputFileName=output.root','doTwoTier=False','doDiag=True']


[jaking@cmslpc236 macros]$ dasgoclient   -daskeys
DAS keys and associated CMS data-service info
---------------------------------------------
 comes from [] services
block comes from [] services, default is [dbs3] system
block,run,lumi comes from [] services
block,run,lumi,events comes from [] services
child comes from [] services
config comes from [] services, default is [reqmgr2] system
dataset comes from [] services, default is [dbs3] system
datatype comes from [] services
era comes from [] services
file comes from [] services, default is [dbs3 rucio] system
file,lumi comes from [] services
file,lumi,events comes from [] services
file,run comes from [] services
file,run,lumi comes from [] services
file,run,lumi,events comes from [] services
group comes from [] services
jobsummary comes from [] services
lumi comes from [] services
lumi,events comes from [] services
mcm comes from [] services
parent comes from [] services
primary_dataset comes from [] services
release comes from [] services
role comes from [] services
rules comes from [] services
run comes from [] services, default is [runregistry dbs3] system
run,lumi comes from [] services
run,lumi,events comes from [] services
site comes from [] services, default is [dbs combined rucio] system
status comes from [] services
summary comes from [] services
tier comes from [] services
user comes from [] services



