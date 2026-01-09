import ROOT
import pickle
import os
import hashlib
from multiprocessing import Pool, cpu_count
from functools import lru_cache
import time

ROOT.EnableImplicitMT()

ROOT.gEnv.SetValue("TFile.AsyncPrefetching", 1)
ROOT.gEnv.SetValue("TFile.MaxCacheSize", 100000000)
ROOT.gEnv.SetValue("TFile.ReadBufferSize", 1048576)


#define BG and SMS master lists #run this eventcount from the python directory
bkg_master_list = '../ntuple_master_lists/KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt'
SMS_master_list = '../ntuple_master_lists/KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt'

eos_prefix = 'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/'
list_dir = '../ntuple_master_lists/'
output_dir = '../config/'
cache_dir = '../cache/'

if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

def WriteEventCountFile( event_counts ):


    with open(output_dir+'EventCount.txt', 'w') as file:
        for key in event_counts:
            event_item = event_counts[key]
            evt_count_str = event_item[0]+' '+str(event_item[1])+' '+str(event_item[2])+'\n'
            file.write(evt_count_str)    
    print(f"Content written to {output_dir}")
    
def GetFileListHash(file_list):
    file_list_str = '|'.join(sorted(file_list))
    return hashlib.md5(file_list_str.encode()).hexdigest()

def ValidateFiles(file_list, max_failures=5):
    valid_files = []
    failures = 0
    
    for file_path in file_list:
        try:
            f = ROOT.TFile.Open(file_path)
            if f and not f.IsZombie() and f.GetNkeys() > 0:
                valid_files.append(file_path)
                if f:
                    f.Close()
            else:
                failures += 1
                if failures >= max_failures:
                    print(f"Warning: {failures} file validation failures, continuing with remaining files")
                    break
        except Exception as e:
            failures += 1
            if failures >= max_failures:
                break
                
    print(f"Validated {len(valid_files)}/{len(file_list)} files ({failures} failures)")
    return valid_files

def SumEventWeights(file_list):
    file_hash = GetFileListHash(file_list)
    cache_file = os.path.join(cache_dir, f"event_weights_{file_hash}.pkl")
    
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                cached_result = pickle.load(f)
            print(f"Using cached result for {len(file_list)} files: ntot={cached_result[0]}, evtwt={cached_result[1]}")
            return cached_result
        except Exception as e:
            print(f"Cache read error: {e}, recomputing...")
    
    print(f"Processing {len(file_list)} files...")
    valid_files = ValidateFiles(file_list)
    
    if not valid_files:
        print("No valid files found!")
        return 0, 0.0
    
    tree_name = "tree/configtree"
    
    try:
        if len(valid_files) == 1:
            df = ROOT.RDataFrame(tree_name, valid_files[0])
        else:
            chain = ROOT.TChain(tree_name)
            for file_path in valid_files:
                chain.Add(file_path)
            df = ROOT.RDataFrame(chain)
        
        ntot_sum_result = df.Sum("nTotEvts")
        evtwt_sum_result = df.Sum("sumEvtWgt")
        ntot_sum = ntot_sum_result.GetValue()
        evtwt_sum = evtwt_sum_result.GetValue()
        
        result = (ntot_sum, evtwt_sum)
        
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(result, f)
        except Exception as e:
            print(f"Cache write warning: {e}")
        
        print(f"Total events: {ntot_sum}, Sum weights: {evtwt_sum}")
        return result
        
    except Exception as e:
        print(f"Error processing files: {e}")
        return 0, 0.0

def GetNtupleFileList( base_path, ntuple_list ):
    f_ntuple = open(  list_dir+ntuple_list )
    ntuple_file_list = []
    for line in f_ntuple:
        ntuple_file_list.append( (eos_prefix+base_path+line).rstrip() )

    return ntuple_file_list 

def ProcessSingleDataset(args):
    base_path, ntuple_list_name, data_set_key = args
    print(f"Processing dataset: {data_set_key}")
    
    try:
        ntuple_file_list = GetNtupleFileList(base_path, ntuple_list_name)
        ntot_sum, evtwt_sum = SumEventWeights(ntuple_file_list)
        return data_set_key, [data_set_key, ntot_sum, evtwt_sum]
    except Exception as e:
        print(f"Error processing dataset {data_set_key}: {e}")
        return data_set_key, [data_set_key, 0, 0.0]

def ProcessMasterList(master_list):
    start_time = time.time()
    
    with open(master_list) as f:
        uncommented_lines = []
        for line in f:
            stripped_line = line.strip()
            if not stripped_line or stripped_line.startswith('#'):
                continue
            processed_line = line.partition('#')[0].rstrip()
            if processed_line:
                uncommented_lines.append(processed_line)
    
    dataset_args = []
    for line in uncommented_lines:
        split_line = line.split()
        if len(split_line) >= 3:
            dataset_args.append((split_line[0], split_line[1], split_line[2]))
    
    print(f"Processing {len(dataset_args)} datasets using {cpu_count()} CPU cores")
    
    num_processes = min(cpu_count(), len(dataset_args))
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(ProcessSingleDataset, dataset_args)
    
    event_counts = dict(results)
    
    elapsed_time = time.time() - start_time
    print(f"Processed {len(event_counts)} datasets in {elapsed_time:.2f} seconds")
    
    return event_counts

if __name__ == "__main__":
    total_start_time = time.time()
    
    print("Starting optimized event counting...")
    print("===========================================")
    
    print("Processing background datasets...")
    bkg_event_counts = ProcessMasterList(bkg_master_list)
    
    print("\nProcessing signal datasets...")
    sms_event_counts = ProcessMasterList(SMS_master_list)
    
    print("\nCombining results...")
    combined_event_counts = bkg_event_counts | sms_event_counts
    
    print(f"\nWriting results for {len(combined_event_counts)} datasets...")
    WriteEventCountFile(combined_event_counts)
    
    total_elapsed = time.time() - total_start_time
    print(f"\nTotal execution time: {total_elapsed:.2f} seconds")
    print("===========================================")
    print("Event counting completed successfully!")
