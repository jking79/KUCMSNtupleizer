import ROOT
ROOT.EnableImplicitMT()


#define BG and SMS master lists #run this eventcount from the python directory
bkg_master_list = '../ntuple_master_lists/KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt'
SMS_master_list = '../ntuple_master_lists/KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt'

eos_prefix = 'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/'
list_dir = '../ntuple_master_lists/'
output_dir = '../config/'

def WriteEventCountFile( event_counts ):


    with open(output_dir+'EventCount.txt', 'w') as file:
        for key in event_counts:
            event_item = event_counts[key]
            evt_count_str = event_item[0]+' '+str(event_item[1])+' '+str(event_item[2])+'\n'
            file.write(evt_count_str)    
    print(f"Content written to {output_dir}")
    
def SumEventWeights( file_list ):
    
    tree_name = "tree/configtree"
    print("Processing ", len(file_list), "files")
    df = ROOT.RDataFrame(tree_name, file_list)
    ntot_sum_result = df.Sum("nTotEvts")
    evtwt_sum_result = df.Sum("sumEvtWgt")
    ntot_sum = ntot_sum_result.GetValue()
    evtwt_sum = evtwt_sum_result.GetValue()
   #event_count = df.Count()
    #print(f"Total number of events: {event_count.GetValue()}")
    print(f"Total number of events: {ntot_sum}  Sum of evt weights: {evtwt_sum}")
    return ntot_sum, evtwt_sum

def GetNtupleFileList( base_path, ntuple_list ):
    f_ntuple = open(  list_dir+ntuple_list )
    ntuple_file_list = []
    for line in f_ntuple:
        ntuple_file_list.append( (eos_prefix+base_path+line).rstrip() )

    return ntuple_file_list 

def ProcessMasterList( master_list ):
   
    #dictionary keyed by datasetkey, to send to output file generation method
    event_counts = {}

    #read in file list
    f = open( master_list )
    uncommented_lines = []
    for line in f:
        stripped_line = line.strip()
        if not stripped_line or stripped_line.startswith('#'):
                continue
        processed_line = line.partition('#')[0].rstrip()
        if processed_line:
            uncommented_lines.append( processed_line )
    f.close()

    #loop over each element in master list
    for line in uncommented_lines:
        split_line = line.split()
        print("Processing Element:")
        print(split_line)
        base_path = split_line[0]
        ntuple_list_name = split_line[1]
        data_set_key = split_line[2]

        #get the ntuple file list associated with master list element
        ntuple_file_list = GetNtupleFileList( base_path, ntuple_list_name )
        
        #debugging print the file list
        #for ntuple_file in ntuple_file_list:
        #    print(ntuple_file)
        #print(ntuple_file_list)

        #Loop over the ntuple files and get numbers
        ntot_sum, evtwt_sum = SumEventWeights(ntuple_file_list)
        #load dictionary for IO
        event_counts[ data_set_key ] = [ data_set_key, ntot_sum, evtwt_sum ] 
    return event_counts

bkg_event_counts = ProcessMasterList( bkg_master_list )
sms_event_counts = ProcessMasterList( SMS_master_list )
#WriteEventCountFile( sms_event_counts )
combined_event_counts = bkg_event_counts | sms_event_counts
WriteEventCountFile( combined_event_counts )
