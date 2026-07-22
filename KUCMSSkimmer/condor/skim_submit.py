#!/usr/bin/env python3
"""
skim_submit.py

Give it an EOS path. It scans EOS directly with xrdfs, finds all
subfolders (HT slices, run eras, signal mass points), and generates
one condor submit file per subfolder with one job per .root file.
Cross sections, keys, and metadata are hard-coded from the master lists
and auto-matched to each subfolder by name.

Example:
  python3 skim_submit.py --eos-path /eos/uscms/store/group/lpcsusylep/jaking/KUCMSNtuple/kucmsntuple_QCD_R18_SVIPM100_v31 --max-mat 50
  python3 skim_submit.py --eos-path /eos/uscms/.../kucmsntuple_GJets_R18_SVIPM100_v31 --dry-run
  python3 skim_submit.py --eos-path /eos/uscms/.../kucmsntuple_MET_R18_SVIPM100_v31
  python3 skim_submit.py --eos-path /eos/uscms/.../kucmsntuple_SMS_Sig_SVIPM100_v31
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from collections import namedtuple

# ---------------------------------------------------------------------------
# Hard-coded cross sections and metadata from master lists
# Keys are substrings that appear in EOS subfolder names.
# mctype: 0 = MC full, 1 = data, 2 = MC fast
# ---------------------------------------------------------------------------

BG_DB = {
    # --- GJets ---
    'GJets_HT-40To100':      {'xsec': 18650,      'key': 'GJets_40t100_UL18'},
    'GJets_HT-100To200':     {'xsec': 8639,        'key': 'GJets_100t200_UL18'},
    'GJets_HT-200To400':     {'xsec': 2173,        'key': 'GJets_200t400_UL18'},
    'GJets_HT-400To600':     {'xsec': 260.7,       'key': 'GJets_400t600_UL18'},
    'GJets_HT-600ToInf':     {'xsec': 86.55,       'key': 'GJets_600up_UL18'},
    # --- QCD ---
    'QCD_HT50to100':         {'xsec': 187700000,   'key': 'QCD_50t100_UL18'},
    'QCD_HT100to200':        {'xsec': 23500000,    'key': 'QCD_100t200_UL18'},
    'QCD_HT200to300':        {'xsec': 1552000,     'key': 'QCD_200t300_UL18'},
    'QCD_HT300to500':        {'xsec': 321100,      'key': 'QCD_300t500_UL18'},
    'QCD_HT500to700':        {'xsec': 30250,       'key': 'QCD_500t700_UL18'},
    'QCD_HT700to1000':       {'xsec': 6398,        'key': 'QCD_700t1000_UL18'},
    'QCD_HT1000to1500':      {'xsec': 1122,        'key': 'QCD_1000t1500_UL18'},
    'QCD_HT1500to2000':      {'xsec': 109.4,       'key': 'QCD_1500t2000_UL18'},
    'QCD_HT2000toInf':       {'xsec': 21.74,       'key': 'QCD_2000up_UL18'},
    # --- WJets ---
    'WJetsToLNu_HT-70To100':    {'xsec': 1283,    'key': 'WJets_70t100_UL18'},
    'WJetsToLNu_HT-100To200':   {'xsec': 1244,    'key': 'WJets_100t200_UL18'},
    'WJetsToLNu_HT-200To400':   {'xsec': 337.8,   'key': 'WJets_200t400_UL18'},
    'WJetsToLNu_HT-400To600':   {'xsec': 45.25,   'key': 'WJets_400t600_UL18'},
    'WJetsToLNu_HT-600To800':   {'xsec': 11.19,   'key': 'WJets_600t800_UL18'},
    'WJetsToLNu_HT-800To1200':  {'xsec': 4.926,   'key': 'WJets_800t1200_UL18'},
    'WJetsToLNu_HT-1200To2500': {'xsec': 1.16,    'key': 'WJets_1200t2500_UL18'},
    'WJetsToLNu_HT-2500ToInf':  {'xsec': 0.026,   'key': 'WJets_2500up_UL18'},
    # --- ZJets ---
    'ZJetsToNuNu_HT-100To200':   {'xsec': 271.3,    'key': 'ZJets_100To200_UL18'},
    'ZJetsToNuNu_HT-200To400':   {'xsec': 72.69,    'key': 'ZJets_200To400_UL18'},
    'ZJetsToNuNu_HT-400To600':   {'xsec': 13.18,    'key': 'ZJets_400To600_UL18'},
    'ZJetsToNuNu_HT-600To800':   {'xsec': 3.257,    'key': 'ZJets_600To800_UL18'},
    'ZJetsToNuNu_HT-800To1200':  {'xsec': 1.076,    'key': 'ZJets_800To1200_UL18'},
    'ZJetsToNuNu_HT-1200To2500': {'xsec': 0.3419,   'key': 'ZJets_1200To2500_UL18'},
    'ZJetsToNuNu_HT-2500ToInf':  {'xsec': 0.005146, 'key': 'ZJets_2500ToInf_UL18'},
    # --- DiPhoton Box ---
    'DiPhotonJetsBox_M40_80':       {'xsec': 306.8,  'key': 'DiPho_M40t80_UL18'},
    'DiPhotonJetsBox_MGG-0to40':    {'xsec': 754.6,  'key': 'DiPho_M0t40_UL18'},
    'DiPhotonJetsBox_MGG-80toInf':  {'xsec': 86.96,  'key': 'DiPho_M80Inf_UL18'},
    # --- TTXJets ---
    'TGJets':   {'xsec': 3,      'key': 'TGJets_UL18'},
    'TTGJets':  {'xsec': 3.77,   'key': 'TTGJets_UL18'},
    'TTJets':   {'xsec': 471.7,  'key': 'TTJets_UL18'},
    'ttWJets':  {'xsec': 0.54,   'key': 'ttWJets_UL18'},
    'ttZJets':  {'xsec': 0.46,   'key': 'ttZJets_UL18'},
    # --- DiBoson/TriBoson ---
    'WGG_5f':   {'xsec': 1.715,   'key': 'WGG_UL18'},
    'WWG_Tune': {'xsec': 0.2158,  'key': 'WWG_UL18'},
    'WWW_4F':   {'xsec': 0.2147,  'key': 'WWW_UL18'},
    'WWZ_4F':   {'xsec': 0.1676,  'key': 'WWZ_UL18'},
    'WZG_Tune': {'xsec': 0.04345, 'key': 'WZG_UL18'},
    'WZ_Tune':  {'xsec': 23.43,   'key': 'WZ_UL18'},
}

# Signal xsec by (process_prefix, gluino_mass_string)
# process_prefix: 'gogo' matches gogoG, gogoZ, gogoGZ; 'sqsq' matches sqsqG
SIGNAL_XSEC = {
    ('gogo', '1500'): 0.02209,
    ('gogo', '2000'): 0.001684,
    ('gogo', '2300'): 0.0004130,
    ('gogo', '2500'): 0.0001687,
    ('sqsq', '1700'): 0.001248,
    ('sqsq', '1850'): 0.0005679,
    ('sqsq', '2000'): 0.0002654,
    ('sqsq', '2150'): 0.0001267,
}

# Timecali tag by year, for data vs MC
TIMECALI = {
    ('R16', 'data'): 'r2_ul16',
    ('R16', 'mc'):   'r2_ul18_mc',
    ('R17', 'data'): 'r2_ul17',
    ('R17', 'mc'):   'r2_ul18_mc',
    ('R18', 'data'): 'r2_ul18',
    ('R18', 'mc'):   'r2_ul18_mc',
    ('R22', 'data'): 'r3_p22',       # TODO: verify r3_p22 exists in caliRunConfig
    ('R22', 'mc'):   'r2_ul18_mc',   # TODO: confirm R3 MC timecali tag
    ('R23', 'data'): 'r3_p23',
    ('R23', 'mc'):   'r2_ul18_mc',   # TODO: confirm R3 MC timecali tag
    ('R24', 'data'): 'r3_p24',
    ('R24', 'mc'):   'r2_ul18_mc',   # TODO: confirm R3 MC timecali tag
    ('R25', 'data'): 'r3_p25unc',
    ('R25', 'mc'):   'r2_ul18_mc',   # TODO: confirm R3 2025 MC timecali tag
}

DATA_KEYWORDS  = {'MET', 'JetMET', 'JetMET0', 'JetMET1', 'JetHT', 'EGamma', 'EGamma0', 'EGamma1', 'EGamma2', 'DisJet'}
SIGNAL_KEYWORDS = {'SMS', 'gogoG', 'gogoZ', 'gogoGZ', 'sqsqG', 'GlGl'}

SANDBOX_DEFAULT  = "/uscms/home/mlazarov/nobackup/sandboxes/sandbox-CMSSW_13_3_3.tar.bz2"
EOS_SKIMS_BASE   = "/store/group/lpcsusylep/malazaro/KUCMSSkims"
EOS_SERVER       = "root://cmseos.fnal.gov"
SKIMMER_DIR      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONDOR_DIR       = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def detect_sample_type(eos_path):
    """Return 'data', 'signal', or 'bg' based on keywords in the EOS path."""
    for kw in SIGNAL_KEYWORDS:
        if kw in eos_path:
            return 'signal'
    for kw in DATA_KEYWORDS:
        # match as word boundary to avoid partial hits
        if re.search(r'(?<![A-Za-z])' + kw + r'(?![A-Za-z])', eos_path):
            return 'data'
    return 'bg'


def detect_year(eos_path):
    """Return 'R16', 'R17', 'R18', 'R22', 'R23', or ''."""
    for yr in ['R16', 'R17', 'R18', 'R22', 'R23', 'R24', 'R25']:
        if yr in eos_path:
            return yr
    return ''


def get_timecali(year, sample_type):
    key = (year, 'data' if sample_type == 'data' else 'mc')
    return TIMECALI.get(key, '')


def default_event_count_path():
    cmssw_base = os.environ.get('CMSSW_BASE', '').strip()
    if cmssw_base:
        return os.path.join(
            cmssw_base,
            'src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer/config/EventCount.txt',
        )
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'config/EventCount.txt')


def read_event_count_keys(path=None):
    if path is None:
        path = default_event_count_path()
    keys = set()
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.split()
                if parts:
                    keys.add(parts[0])
    except OSError:
        print('WARNING: could not read', path,
              '- MC EventCount key checks will be skipped')
    return keys


def normalize_branch_mask_path(mask_path):
    """Resolve a submit-side branch mask path and return the worker-side path.

    Condor jobs run from scratch after unpacking config.tgz, so skimmer
    arguments must refer to files as config/branch_masks/<name>. Locally the
    user may call this script from KUCMSSkimmer/, KUCMSSkimmer/condor/, or
    elsewhere, so accept cwd-relative and skimmer-relative paths.
    """
    if not mask_path:
        return ''

    candidates = []
    if os.path.isabs(mask_path):
        candidates.append(mask_path)
    else:
        candidates.extend([
            os.path.abspath(mask_path),
            os.path.abspath(os.path.join(CONDOR_DIR, mask_path)),
            os.path.abspath(os.path.join(SKIMMER_DIR, mask_path)),
        ])

    resolved = None
    for candidate in candidates:
        if os.path.isfile(candidate):
            resolved = candidate
            break

    if resolved is None:
        raise ValueError('branch mask file not found: ' + mask_path)

    branch_mask_dir = os.path.join(SKIMMER_DIR, 'config', 'branch_masks')
    if os.path.commonpath([branch_mask_dir, resolved]) != branch_mask_dir:
        raise ValueError(
            'branch mask must be under config/branch_masks so it is available '
            'inside config.tgz: ' + resolved
        )

    rel = os.path.relpath(resolved, SKIMMER_DIR)
    return rel.replace(os.sep, '/')


def lookup_bg_meta(subfolder):
    """Find the best matching BG_DB entry for this subfolder name.
    Returns the dict or None."""
    # sort by key length descending so longer (more specific) patterns win
    for pattern in sorted(BG_DB, key=len, reverse=True):
        if pattern in subfolder:
            return BG_DB[pattern]
    return None


def lookup_signal_xsec(subfolder):
    """Extract process prefix and gluino mass from subfolder, look up xsec."""
    proc = 'gogo' if any(p in subfolder for p in ['gogoG', 'gogoZ', 'gogoGZ', 'GlGl']) else 'sqsq'
    m = re.search(r'mGl-(\d+)', subfolder)
    if not m:
        return None, '0', '0', '0'
    mgl = m.group(1)
    xsec = SIGNAL_XSEC.get((proc, mgl), None)
    mn2 = re.search(r'mN2-(\d+)', subfolder)
    mn1 = re.search(r'mN1-(\d+)', subfolder)
    n2 = mn2.group(1) if mn2 else '0'
    n1 = mn1.group(1) if mn1 else '0'
    return xsec, mgl, n2, n1


def signal_process_key(name):
    if 'SqSq' in name or 'sqsq' in name:
        return 'sqsqG'
    if 'GlGl-GZ' in name or 'GlGlGZ' in name or 'gogoGZ' in name:
        return 'gogoGZ'
    if 'GlGl-Z' in name or 'GlGlZ' in name or 'gogoZ' in name:
        return 'gogoZ'
    if 'GlGl-G' in name or 'GlGlG' in name or 'gogoG' in name:
        return 'gogoG'
    return 'gogoGZ'


def signal_event_count_tier(name):
    if re.search(r'(?:FULLMINI|MINIAOD|MINI|MiniAOD|Mini)', name):
        return 'FULLMINI'
    return 'AODSIM'


def normalize_signal_ctau(raw_ctau, tier):
    value = raw_ctau.lstrip('-')
    if tier == 'AODSIM' and value.startswith('0p') and len(value) == 3:
        return value[-1]
    if tier == 'FULLMINI' and value in {'1', '5'}:
        return '0p' + value
    return value


def signal_event_count_key(name):
    mgl = re.search(r'mGl-(\d+)', name)
    mn2 = re.search(r'mN2-(\d+)', name)
    mn1 = re.search(r'mN1-(\d+)', name)
    ctau = re.search(r'(?:(?:^|[_-])ct|\dctau-)(-?[0-9]+p[0-9]+|-?[0-9]+)', name)

    if not (mgl and mn2 and mn1 and ctau):
        match = re.search(
            r'((?:gogoGZ|gogoG|gogoZ|sqsqG)_(?:AODSIM|FULLMINI|MINIAOD|MINI|AOD)_.*)',
            name,
        )
        if match:
            key = re.sub(r'_(MINIAOD|MINI)_', '_FULLMINI_', match.group(1))
            if '_AODSIM_' in key:
                return re.sub(r'_ct0p([15])$', r'_ct\1', key)
            if '_FULLMINI_' in key:
                return re.sub(r'_ct([15])$', r'_ct0p\1', key)
            return key
        raise ValueError('could not derive signal EventCount key from: ' + name)

    process = signal_process_key(name)
    tier = signal_event_count_tier(name)
    ct_key = normalize_signal_ctau(ctau.group(1), tier)
    return (
        f'{process}_{tier}_mGl-{mgl.group(1)}'
        f'_mN2-{mn2.group(1)}_mN1-{mn1.group(1)}_ct{ct_key}'
    )


def signal_metadata_name(subfolder, eos_path):
    """Choose the best available name for signal xsec/EventCount parsing."""
    candidates = [
        subfolder,
        os.path.basename(eos_path.rstrip('/')),
        os.path.basename(os.path.dirname(eos_path.rstrip('/'))),
    ]
    for name in candidates:
        if (re.search(r'mGl-\d+', name)
                and re.search(r'mN2-\d+', name)
                and re.search(r'mN1-\d+', name)
                and re.search(r'(?:ct|\dctau-)(?:-?[0-9]+p[0-9]+|-?[0-9]+)', name)):
            return name
    return subfolder


def make_data_key(subfolder, sample_type_kw):
    """Build a short dataset key for data from the subfolder name."""
    # e.g. 'MET_AOD_Run2018A-15Feb...' -> 'MET18A'
    year_m = re.search(r'(\d{4})', subfolder)
    era_m  = re.search(r'Run\d{4}([A-Z])', subfolder)
    yr  = year_m.group(1)[2:] if year_m else ''
    era = era_m.group(1)       if era_m  else ''
    return sample_type_kw + yr + era


DatasetInfo = namedtuple('DatasetInfo', [
    'collection_dir',
    'collection_name',
    'task_dir',
    'task_name',
    'timestamp_dir',
    'raw_task_suffix',
    'sample_tag',
    'ntuple_tag',
    'canonical_base',
    'root_files',
])


DatasetPaths = namedtuple('DatasetPaths', [
    'sample_tag',
    'ntuple_tag',
    'canonical_base',
    'work_dir',
    'src_dir',
    'log_dir',
    'out_dir',
    'submit_path',
    'ofilename',
])


SubmittedJob = namedtuple('SubmittedJob', [
    'index',
    'input_file',
    'arg_line',
    'submit_path',
])


QueueEntry = namedtuple('QueueEntry', [
    'global_index',
    'arg_line',
])


SubmissionHistory = namedtuple('SubmissionHistory', [
    'submit_files',
    'jobs',
    'eos_out_dir',
    'header_lines',
    'original_count',
    'original_arg_line',
])


def data_keyword(eos_path):
    return next((kw for kw in sorted(DATA_KEYWORDS, key=len, reverse=True)
                 if kw in eos_path), 'Data')


def sample_context(eos_path):
    sample_type = detect_sample_type(eos_path)
    year        = detect_year(eos_path)
    timecali    = get_timecali(year, sample_type)
    if not timecali and sample_type == 'signal':
        timecali = 'r2_ul18_mc'
    return sample_type, year, timecali, data_keyword(eos_path)


def print_sample_context(eos_path, sample_type, year, timecali):
    print('EOS path:    ', eos_path)
    print('Sample type: ', sample_type)
    print('Year:        ', year)
    print('Timecali:    ', timecali)
    print()


def resolve_dataset_metadata(subfolder, eos_path, sample_type, data_kw):
    xsec       = '1'
    key        = subfolder[:subfolder.index('_')] if '_' in subfolder else subfolder
    gluinomass = '0'
    n2mass     = '0'
    n1mass     = '0'
    mctype     = 0
    mc_wt      = '1'

    if sample_type == 'data':
        mctype = 1
        key    = make_data_key(subfolder, data_kw)
    elif sample_type == 'signal':
        metadata_name = signal_metadata_name(subfolder, eos_path)
        sig_xsec, gluinomass, n2mass, n1mass = lookup_signal_xsec(metadata_name)
        if sig_xsec is not None:
            xsec = str(sig_xsec)
        else:
            print('  WARNING: no xsec found for', metadata_name, '- using 0')
            xsec = '0'
        key = signal_event_count_key(metadata_name)
    else:
        meta = lookup_bg_meta(subfolder)
        if meta:
            xsec = str(meta['xsec'])
            key  = meta['key']
        else:
            print('  WARNING: no xsec found for', subfolder, '- using 0')
            xsec = '0'

    return {
        'xsec': xsec,
        'key': key,
        'gluinomass': gluinomass,
        'n2mass': n2mass,
        'n1mass': n1mass,
        'mctype': mctype,
        'mc_wt': mc_wt,
    }


def build_skimmer_flags(metadata, args, timecali):
    flags = (
        ' --xsec '        + metadata['xsec']
        + ' --dataSetKey ' + metadata['key']
        + ' --gluinoMass ' + metadata['gluinomass']
        + ' --N2Mass '     + metadata['n2mass']
        + ' --timeCaliTag ' + timecali
        + ' --MCweight '   + metadata['mc_wt']
        + ' --MCtype '     + str(metadata['mctype'])
    )
    if metadata['mctype'] == 0:
        flags += ' --hasGenInfo'
    if not args.psiche:
        flags += ' --noBHC'
    if args.no_sv:
        flags += ' --noSV'
    if args.hlt_off:
        flags += ' --HLTPathsOff'
    if args.branch_mask:
        flags += ' --branchMask ' + args.branch_mask
    return flags


def dataset_paths(output_dir, dataset, tag):
    odir     = output_dir.rstrip('/') + '/'
    work_dir = odir + dataset.ntuple_tag + '/' + dataset.sample_tag + '/' + tag
    return DatasetPaths(
        sample_tag=dataset.sample_tag,
        ntuple_tag=dataset.ntuple_tag,
        canonical_base=dataset.canonical_base,
        work_dir=work_dir,
        src_dir=work_dir + '/src',
        log_dir=work_dir + '/log',
        out_dir=work_dir + '/out',
        submit_path=work_dir + '/src/submit.sh',
        ofilename='condor_' + dataset.canonical_base + '__' + tag,
    )


def eos_output_dir(args, paths):
    if not args.eos_out:
        return None
    return (os.path.expandvars(args.eos_out).rstrip('/')
            + '/' + paths.ntuple_tag + '/' + paths.sample_tag
            + '/' + args.tag)


def ensure_submit_dirs(paths, eos_out_dir):
    local_dirs = [paths.src_dir, paths.log_dir] if eos_out_dir else [
        paths.src_dir, paths.log_dir, paths.out_dir,
    ]
    for d in local_dirs:
        os.makedirs(d, exist_ok=True)


def warn_missing_event_count_key(metadata, event_count_keys):
    if metadata['mctype'] != 1 and event_count_keys and metadata['key'] not in event_count_keys:
        print('  WARNING: key not found in config/EventCount.txt:', metadata['key'])
        print('           This job will likely get inf evtFillWgt.')


def weight_record(metadata, timecali, root_files, subfolder):
    return {
        'xsec': metadata['xsec'],
        'mc_wt': metadata['mc_wt'],
        'mctype': metadata['mctype'],
        'timecali': timecali,
        'gluinomass': metadata['gluinomass'],
        'n2mass': metadata['n2mass'],
        'n1mass': metadata['n1mass'],
        'n_files': len(root_files),
        'subfolder': subfolder,
    }


def xrdfs_mkdir(path, verbose=False):
    """Create a directory on EOS via xrdfs, including any missing parents."""
    bare = path
    if bare.startswith('/eos/uscms'):
        bare = bare[len('/eos/uscms'):]
    bare = '/' + bare.lstrip('/')
    cmd = ['xrdfs', 'root://cmseos.fnal.gov', 'mkdir', '-p', bare]
    if verbose:
        print(' ', ' '.join(cmd))
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        print('  WARNING: failed to create EOS directory:', bare)
        return False


def xrdfs_ls(path, verbose=False):
    """List entries under path using xrdfs. Returns list of bare paths."""
    # path can be /eos/uscms/... or root://server//...
    if path.startswith('root://'):
        m = re.match(r'(root://[^/]+/)(/.+)', path)
        server = m.group(1).rstrip('/')
        bare   = m.group(2)
    else:
        server = 'root://cmseos.fnal.gov'
        bare   = '/' + path.lstrip('/')
    cmd = ['xrdfs', server, 'ls', bare]
    if verbose:
        print(' ', ' '.join(cmd))
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode()
        return [line.strip() for line in out.splitlines() if line.strip()]
    except subprocess.CalledProcessError:
        return []


def xrdfs_stat(path, verbose=False):
    """Return True if path exists on EOS."""
    bare = path
    if bare.startswith('/eos/uscms'):
        bare = bare[len('/eos/uscms'):]
    bare = '/' + bare.lstrip('/')
    cmd = ['xrdfs', EOS_SERVER, 'stat', bare]
    if verbose:
        print(' ', ' '.join(cmd))
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False


def normalize_eos_path(fpath):
    """Normalize any EOS path form to a bare /store/... path for comparison.
    Handles: root://server//store/..., root://server/store/..., /eos/uscms/store/..., /store/...
    """
    def clean(path):
        return re.sub(r'/+', '/', path.rstrip('/'))

    if fpath.startswith('root://'):
        m = re.match(r'root://[^/]+/(/.+)', fpath)
        if m:
            return clean(m.group(1))
        m = re.match(r'root://[^/]+/(.+)', fpath)
        if m:
            return clean('/' + m.group(1))
    if fpath.startswith('/eos/uscms'):
        return clean(fpath[len('/eos/uscms'):])
    return clean('/' + fpath.lstrip('/'))


def eos_xrootd_url(path, server=EOS_SERVER):
    """Return an EOS path in canonical XRootD URL form."""
    bare = normalize_eos_path(path)
    return server.rstrip('/') + '//' + bare.lstrip('/')


def task_suffix_from_name(collection_name, task_name):
    collection_core = remove_prefix_once(collection_name, 'kucmsntuple_')
    exact_prefix = 'kucmsntuple_' + collection_core + '_'
    if task_name.startswith(exact_prefix):
        return task_name[len(exact_prefix):]
    return task_name


def make_sample_tag(raw_task_suffix):
    tier_tokens = {'MiniAOD', 'MINIAOD', 'AODSIM', 'MINI', 'MIN'}
    tokens = [token for token in raw_task_suffix.split('_') if token not in tier_tokens]
    sample = '_'.join(tokens)
    if sample.startswith('SMS-'):
        sample = sample[len('SMS-'):]
    elif sample.startswith('SMS_'):
        sample = sample[len('SMS_'):]
    return sample


def make_ntuple_tag(collection_name):
    core = remove_prefix_once(collection_name, 'kucmsntuple_')
    tokens = core.split('_')
    start = None
    for i, token in enumerate(tokens):
        if token.startswith('SV'):
            start = i
            break
    if start is None:
        raise ValueError('Collection has no SV* ntuple token: ' + collection_name)
    kept = [token for token in tokens[start:] if token not in {'MiniAOD', 'MINIAOD'}]
    return '_'.join(kept)


def final_merged_filename(sample_tag, ntuple_tag, skim_tag):
    return sample_tag + '__' + ntuple_tag + '__' + skim_tag + '.root'


def find_enclosing_collection(path):
    norm = normalize_eos_path(path)
    parts = [p for p in norm.split('/') if p]
    prefix = ''
    for part in parts:
        prefix += '/' + part
        if part.startswith('kucmsntuple_'):
            return prefix, part
    raise ValueError('EOS path must identify one kucmsntuple_* collection directory: ' + path)


def reject_shard_level_path(path):
    norm = normalize_eos_path(path)
    parent = os.path.basename(os.path.dirname(norm))
    leaf = os.path.basename(norm)
    if _CRAB_TS_RE.match(parent) and re.match(r'^\d+$', leaf):
        raise ValueError(
            'Shard-level paths are not supported; provide the CRAB task or timestamp directory: '
            + path
        )


def _is_dir_path(path):
    return not path.endswith('.root') and '.' not in os.path.basename(path)


def _root_files_for_timestamp(timestamp_dir, max_files=-1, verbose=False):
    files = find_root_files(timestamp_dir, max_files=max_files, verbose=verbose)
    return tuple(normalize_eos_path(f) for f in files)


def _dataset_info(collection_dir, collection_name, task_dir, timestamp_dir,
                  max_files=-1, verbose=False):
    task_dir = normalize_eos_path(task_dir)
    timestamp_dir = normalize_eos_path(timestamp_dir)
    task_name = os.path.basename(task_dir)
    raw_task_suffix = task_suffix_from_name(collection_name, task_name)
    sample_tag = make_sample_tag(raw_task_suffix)
    ntuple_tag = make_ntuple_tag(collection_name)
    canonical_base = sample_tag + '__' + ntuple_tag
    return DatasetInfo(
        collection_dir=collection_dir,
        collection_name=collection_name,
        task_dir=task_dir,
        task_name=task_name,
        timestamp_dir=timestamp_dir,
        raw_task_suffix=raw_task_suffix,
        sample_tag=sample_tag,
        ntuple_tag=ntuple_tag,
        canonical_base=canonical_base,
        root_files=_root_files_for_timestamp(timestamp_dir, max_files=max_files,
                                             verbose=verbose),
    )


def _task_timestamp_dirs(task_dir, verbose=False):
    return sorted(
        normalize_eos_path(e) for e in xrdfs_ls(task_dir, verbose=verbose)
        if _CRAB_TS_RE.match(os.path.basename(e))
    )


def _discover_task_dirs(path, collection_dir, collection_name, max_files=-1,
                        verbose=False):
    entries = [normalize_eos_path(e) for e in xrdfs_ls(path, verbose=verbose)]
    subdirs = [e for e in entries if _is_dir_path(e)]
    timestamps = [e for e in subdirs if _CRAB_TS_RE.match(os.path.basename(e))]
    if timestamps:
        latest = max(timestamps, key=crab_timestamp_key)
        skipped = sorted(os.path.basename(e) for e in timestamps if e != latest)
        if skipped:
            print('  [crab] skipping:', ', '.join(skipped),
                  '-> using:', os.path.basename(latest))
        return [_dataset_info(collection_dir, collection_name, path, latest,
                              max_files=max_files, verbose=verbose)]

    results = []
    for subdir in subdirs:
        results.extend(_discover_task_dirs(subdir, collection_dir, collection_name,
                                           max_files=max_files, verbose=verbose))
    return results


def validate_unique_canonical_bases(datasets):
    seen = {}
    for dataset in datasets:
        previous = seen.get(dataset.canonical_base)
        if previous and previous.task_dir != dataset.task_dir:
            raise ValueError(
                'Duplicate cleaned skim name "{}" for tasks:\n'
                '  {} ({})\n'
                '  {} ({})'.format(
                    dataset.canonical_base,
                    previous.task_dir,
                    previous.timestamp_dir,
                    dataset.task_dir,
                    dataset.timestamp_dir,
                )
            )
        seen[dataset.canonical_base] = dataset


_CRAB_TS_RE = re.compile(r'^(?:\d{6}|\d{8})_\d{6}$')


def remove_prefix_once(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def crab_timestamp_key(path):
    name = os.path.basename(path)
    date, time = name.split('_', 1)
    if len(date) == 6:
        date = '20' + date
    return date + '_' + time


def find_root_files(path, max_files=-1, verbose=False):
    """Recursively find all .root files under path via xrdfs."""
    results = []
    stack = xrdfs_ls(path, verbose=verbose)
    while stack:
        entry = stack.pop()
        if entry.endswith('.root'):
            results.append(entry)
            if max_files > 0 and len(results) >= max_files:
                break
        elif '.' not in os.path.basename(entry):
            stack.extend(xrdfs_ls(entry, verbose=verbose))
    return results


def find_crab_datasets(path, name='', max_files=-1, peek_only=False, verbose=False):
    """Discover CRAB datasets with names derived from collection and task dirs."""
    del name, peek_only
    norm = normalize_eos_path(path)
    reject_shard_level_path(norm)
    collection_dir, collection_name = find_enclosing_collection(norm)

    leaf = os.path.basename(norm)
    parent = normalize_eos_path(os.path.dirname(norm))
    if _CRAB_TS_RE.match(leaf):
        timestamps = _task_timestamp_dirs(parent, verbose=verbose)
        if norm not in timestamps:
            raise ValueError('Explicit timestamp directory does not exist under task: ' + norm)
        datasets = [_dataset_info(collection_dir, collection_name, parent, norm,
                                  max_files=max_files, verbose=verbose)]
    else:
        datasets = _discover_task_dirs(norm, collection_dir, collection_name,
                                       max_files=max_files, verbose=verbose)

    validate_unique_canonical_bases(datasets)
    return datasets


def write_submit(path, log_dir, out_dir, ofilename, root_files, flags, args, eos_out_dir=None, offset=0):
    outname = ofilename + '.$(Process).root'
    if offset > 0 and root_files:
        log_first = offset
        log_last = offset + len(root_files) - 1
        log_stem = 'job.{}_{}'.format(log_first, log_last)
    else:
        log_stem = 'job'

    if eos_out_dir:
        bare_eos = eos_out_dir
        if bare_eos.startswith('/eos/uscms'):
            bare_eos = bare_eos[len('/eos/uscms'):]
        bare_eos = bare_eos.rstrip('/')
        output_lines = [
            'transfer_output_files = ',
            'when_to_transfer_output = NEVER',
        ]
    else:
        abs_cwd = subprocess.check_output(['pwd']).decode().strip()
        remap_target = os.path.join(abs_cwd, out_dir) + '/' + outname
        output_lines = [
            'transfer_output_files = ' + outname,
            'transfer_output_remaps = "' + outname + '=' + remap_target + '"',
            'when_to_transfer_output = ON_EXIT',
        ]

    lines = [
        'universe = vanilla',
        'executable = execute_script.sh',
        'output = ' + os.path.join(log_dir, log_stem + '.$(Process).out'),
        'error = '  + os.path.join(log_dir, log_stem + '.$(Process).err'),
        'log = '    + os.path.join(log_dir, log_stem + '.log'),
        'transfer_input_files = ' + args.sandbox + ', config.tgz,',
        'should_transfer_files = YES',
    ] + output_lines
    if args.max_mat  != -1: lines.append('max_materialize = ' + str(args.max_mat))
    if args.max_idle != -1: lines.append('max_idle = '        + str(args.max_idle))
    if args.memory   != -1: lines.append('request_memory = '  + str(args.memory))

    lines += ['', '', 'queue Arguments from (']
    for i, fpath in enumerate(root_files):
        full = eos_xrootd_url(fpath)
        eos_arg = ' --eos-out ' + bare_eos if eos_out_dir else ''
        job_num = offset + i
        # When offset > 0 (new-inputs submission), use explicit indices so output
        # filenames don't collide with the original submission's $(Process) range.
        proc = str(job_num) if offset else '$(Process)'
        lines.append('###### job' + str(job_num) + ' ######')
        lines.append('-i ' + full + ' ' + flags + ' -o ' + ofilename + '.' + proc + eos_arg)
    lines.append(')')

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# Job checking / resubmission helpers
# ---------------------------------------------------------------------------

def parse_submit_file(submit_path):
    """Parse a condor submit.sh file.
    Returns (header_lines, job_entries, eos_out_dir).
      header_lines  -- lines up to and including 'queue Arguments from ('
      job_entries   -- list of (job_idx, arg_string) for each -i line
      eos_out_dir   -- bare /store/... path from --eos-out, or None
    """
    with open(submit_path) as f:
        lines = f.read().splitlines()

    header_lines = []
    job_entries  = []
    eos_out_dir  = None
    in_queue     = False
    job_idx      = 0

    for line in lines:
        stripped = line.strip()
        if not in_queue:
            header_lines.append(line)
            if stripped == 'queue Arguments from (':
                in_queue = True
            continue
        if stripped == ')':
            break
        if stripped.startswith('#') or not stripped:
            continue
        if stripped.startswith('-i '):
            if eos_out_dir is None:
                m = re.search(r'--eos-out\s+(\S+)', stripped)
                if m:
                    eos_out_dir = m.group(1)
            job_entries.append((job_idx, stripped))
            job_idx += 1

    return header_lines, job_entries, eos_out_dir


def find_submission_history_files(original_submit_path):
    """Return submit.sh, legacy newfiles.sh, then ranged newfiles batches."""
    if not os.path.exists(original_submit_path):
        raise FileNotFoundError(original_submit_path)

    src_dir = os.path.dirname(original_submit_path)
    history = [original_submit_path]

    legacy = os.path.join(src_dir, 'newfiles.sh')
    if os.path.exists(legacy):
        history.append(legacy)

    ranged = []
    pattern = re.compile(r'^newfiles_(\d+)_(\d+)\.sh$')
    for name in os.listdir(src_dir):
        match = pattern.match(name)
        if match:
            ranged.append((int(match.group(1)), int(match.group(2)),
                           os.path.join(src_dir, name)))
    history.extend(path for _, _, path in sorted(ranged))
    return history


def parse_submit_queue_lines(submit_path):
    """Return header lines and all non-comment queue argument lines."""
    header_lines, queue_entries = parse_submit_queue_entries(submit_path)
    return header_lines, [entry.arg_line for entry in queue_entries]


def parse_submit_queue_entries(submit_path):
    """Return header lines and queue entries for supported queue forms."""
    with open(submit_path) as f:
        lines = f.read().splitlines()

    header_lines = []
    queue_entries = []
    in_queue = False
    two_column_queue = False

    for line in lines:
        stripped = line.strip()
        if not in_queue:
            header_lines.append(line)
            if stripped == 'queue Arguments from (':
                in_queue = True
                two_column_queue = False
            elif stripped == 'queue GlobalJobIndex, Arguments from (':
                in_queue = True
                two_column_queue = True
            continue
        if stripped == ')':
            break
        if stripped.startswith('#') or not stripped:
            continue
        if two_column_queue:
            match = re.match(r'^(\d+)\s*,\s*(.+)$', stripped)
            if not match:
                raise ValueError(
                    '{} has malformed GlobalJobIndex queue row: {}'.format(
                        submit_path, stripped)
                )
            queue_entries.append(QueueEntry(int(match.group(1)), match.group(2)))
        else:
            queue_entries.append(QueueEntry(None, stripped))

    return header_lines, queue_entries


def _arg_value(tokens, option, arg_line, submit_path):
    pos = _arg_position(tokens, option, arg_line, submit_path)
    return tokens[pos + 1]


def _arg_position(tokens, option, arg_line, submit_path):
    try:
        pos = tokens.index(option)
    except ValueError:
        raise ValueError('{} entry is missing {}: {}'.format(
            submit_path, option, arg_line))
    if pos + 1 >= len(tokens):
        raise ValueError('{} entry is missing value for {}: {}'.format(
            submit_path, option, arg_line))
    return pos


def submit_arg_tokens(arg_line, submit_path):
    try:
        return shlex.split(arg_line)
    except ValueError as exc:
        raise ValueError('{} has malformed queue arguments: {}'.format(
            submit_path, exc))


def submit_arg_template_signature(arg_line, submit_path):
    tokens = submit_arg_tokens(arg_line, submit_path)
    skip = set()
    for option in ('-i', '-o'):
        pos = _arg_position(tokens, option, arg_line, submit_path)
        skip.add(pos)
        skip.add(pos + 1)
    return tuple(token for idx, token in enumerate(tokens) if idx not in skip)


def build_arg_line_from_template(template_arg_line, submit_path, input_file, output_expr):
    tokens = submit_arg_tokens(template_arg_line, submit_path)
    tokens[_arg_position(tokens, '-i', template_arg_line, submit_path) + 1] = input_file
    tokens[_arg_position(tokens, '-o', template_arg_line, submit_path) + 1] = output_expr
    return ' '.join(shlex.quote(token) for token in tokens)


def replace_output_arg(arg_line, submit_path, output_expr):
    tokens = submit_arg_tokens(arg_line, submit_path)
    tokens[_arg_position(tokens, '-o', arg_line, submit_path) + 1] = output_expr
    return ' '.join(shlex.quote(token) for token in tokens)


def parse_submit_arg_fields(arg_line, submit_path):
    tokens = submit_arg_tokens(arg_line, submit_path)

    input_file = _arg_value(tokens, '-i', arg_line, submit_path)
    output_prefix = _arg_value(tokens, '-o', arg_line, submit_path)
    eos_out_dir = None
    if '--eos-out' in tokens:
        eos_out_dir = _arg_value(tokens, '--eos-out', arg_line, submit_path)
    return input_file, output_prefix, eos_out_dir


def recover_submission_history(original_submit_path):
    """Recover and validate actual global jobs from submit history files."""
    history_files = find_submission_history_files(original_submit_path)
    original_abs = os.path.abspath(original_submit_path)
    path_info = submit_path_parts(original_submit_path)
    expected_prefix = expected_submit_prefix(path_info)

    header_lines = None
    jobs = []
    original_count = 0
    eos_out_dirs = set()
    original_template_signature = None
    original_template_arg_line = None

    for submit_path in history_files:
        file_header, queue_entries = parse_submit_queue_entries(submit_path)
        if os.path.abspath(submit_path) == original_abs:
            header_lines = file_header
            original_count = len(queue_entries)
        is_original = os.path.abspath(submit_path) == original_abs

        for local_idx, queue_entry in enumerate(queue_entries):
            arg_line = queue_entry.arg_line
            input_file, output_expr, eos_out_dir = parse_submit_arg_fields(
                arg_line, submit_path)
            if eos_out_dir:
                eos_out_dirs.add(normalize_eos_path(eos_out_dir))

            template_signature = submit_arg_template_signature(arg_line, submit_path)
            if is_original:
                if original_template_signature is None:
                    original_template_signature = template_signature
                    original_template_arg_line = arg_line
                elif template_signature != original_template_signature:
                    raise ValueError(
                        'original submit.sh jobs do not share one argument template: '
                        + submit_path
                    )
            elif template_signature != original_template_signature:
                raise ValueError(
                    'submission history entry has processing arguments that differ '
                    'from original submit.sh: ' + submit_path
                )

            process_expr = expected_prefix + '.$(Process)'
            explicit_match = re.match(
                r'^' + re.escape(expected_prefix) + r'\.(\d+)$',
                output_expr,
            )
            if output_expr == process_expr:
                if not is_original:
                    raise ValueError(
                        '$(Process) output entry occurs outside original submit.sh: '
                        + submit_path
                )
                job_index = local_idx
            elif explicit_match:
                job_index = int(explicit_match.group(1))
            else:
                raise ValueError(
                    '{} entry uses unexpected output prefix: {} '
                    '(expected {}.$(Process) or {}.<integer>)'.format(
                        submit_path, output_expr, expected_prefix, expected_prefix)
                )
            if queue_entry.global_index is not None:
                if not explicit_match:
                    raise ValueError(
                        '{} GlobalJobIndex row must use an explicit -o index: {}'.format(
                            submit_path, arg_line)
                    )
                if queue_entry.global_index != job_index:
                    raise ValueError(
                        '{} GlobalJobIndex {} does not match -o index {}: {}'.format(
                            submit_path, queue_entry.global_index, job_index, arg_line)
                    )

            jobs.append(SubmittedJob(
                index=job_index,
                input_file=normalize_eos_path(input_file),
                arg_line=arg_line,
                submit_path=submit_path,
            ))

    if original_count == 0:
        raise ValueError('no jobs found in original submit.sh: ' + original_submit_path)

    by_index = {}
    by_input = {}
    for job in jobs:
        if job.index in by_index:
            raise ValueError(
                'duplicate global job index {} in history files:\n  {}\n  {}'.format(
                    job.index, by_index[job.index].submit_path, job.submit_path)
            )
        by_index[job.index] = job

        if job.input_file in by_input:
            raise ValueError(
                'duplicate input file in submission history: {}\n  {}\n  {}'.format(
                    job.input_file, by_input[job.input_file].submit_path,
                    job.submit_path)
            )
        by_input[job.input_file] = job

    expected_indices = set(range(max(by_index) + 1)) if by_index else set()
    missing = sorted(expected_indices - set(by_index))
    if missing:
        raise ValueError(
            'submission history has missing global job indices: '
            + ', '.join(str(i) for i in missing)
        )

    if len(eos_out_dirs) > 1:
        raise ValueError(
            'submission history has conflicting --eos-out directories: '
            + ', '.join(sorted(eos_out_dirs))
        )

    return SubmissionHistory(
        submit_files=history_files,
        jobs=sorted(jobs, key=lambda job: job.index),
        eos_out_dir=next(iter(eos_out_dirs)) if eos_out_dirs else None,
        header_lines=header_lines or [],
        original_count=original_count,
        original_arg_line=original_template_arg_line,
    )


def recover_submission_history_or_exit(submit_path):
    try:
        return recover_submission_history(submit_path)
    except (OSError, ValueError) as exc:
        print('ERROR: invalid submission history for', submit_path)
        print(' ', exc)
        sys.exit(1)


def history_expected_indices(history):
    return {job.index for job in history.jobs}


def new_input_files(current_eos_inputs, submitted_jobs):
    submitted_inputs = {job.input_file for job in submitted_jobs}
    current_inputs = {normalize_eos_path(fpath) for fpath in current_eos_inputs}
    return sorted(current_inputs - submitted_inputs)


def newfiles_submit_path(src_dir, first_index, last_index):
    return os.path.join(
        src_dir,
        'newfiles_{}_{}.sh'.format(first_index, last_index),
    )


def _header_value(line):
    return line.split('=', 1)[1].strip() if '=' in line else ''


def _local_output_dir_from_header(header_lines):
    for line in header_lines:
        m = re.search(r'transfer_output_remaps\s*=\s*"[^=]+=([^"]+)/[^/"]+\.root"', line)
        if m:
            return m.group(1)
    raise ValueError('original submit.sh local output remap could not be parsed')


def new_input_submit_header(history, expected_prefix, first_index, last_index):
    log_stem = 'job.{}_{}'.format(first_index, last_index)
    local_output = history.eos_out_dir is None
    local_output_file = expected_prefix + '.$(GlobalJobIndex).root'
    local_output_dir = (_local_output_dir_from_header(history.header_lines)
                        if local_output else None)
    header = []

    for line in history.header_lines:
        stripped = line.strip()
        if stripped.startswith('output'):
            old_path = _header_value(line)
            log_dir = os.path.dirname(old_path)
            header.append('output = ' + os.path.join(
                log_dir, log_stem + '.$(Process).out'))
        elif stripped.startswith('error'):
            old_path = _header_value(line)
            log_dir = os.path.dirname(old_path)
            header.append('error = ' + os.path.join(
                log_dir, log_stem + '.$(Process).err'))
        elif stripped.startswith('log'):
            old_path = _header_value(line)
            log_dir = os.path.dirname(old_path)
            header.append('log = ' + os.path.join(log_dir, log_stem + '.log'))
        elif local_output and stripped.startswith('transfer_output_files'):
            header.append('transfer_output_files = ' + local_output_file)
        elif local_output and stripped.startswith('transfer_output_remaps'):
            remap_target = os.path.join(local_output_dir, local_output_file)
            header.append(
                'transfer_output_remaps = "'
                + local_output_file + '=' + remap_target + '"'
            )
        elif stripped in {
                'queue Arguments from (',
                'queue GlobalJobIndex, Arguments from (',
        }:
            if local_output:
                header.append('queue GlobalJobIndex, Arguments from (')
            else:
                header.append('queue Arguments from (')
        else:
            header.append(line)
    return header


def write_inherited_new_inputs_submit(path, history, expected_prefix, new_files,
                                      first_index):
    last_index = first_index + len(new_files) - 1
    header = new_input_submit_header(history, expected_prefix, first_index,
                                     last_index)
    local_output = history.eos_out_dir is None

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for line in header:
            handle.write(line + '\n')
        for offset, input_file in enumerate(new_files):
            job_index = first_index + offset
            output_expr = expected_prefix + '.' + str(job_index)
            input_url = eos_xrootd_url(input_file)
            arg_line = build_arg_line_from_template(
                history.original_arg_line,
                history.submit_files[0],
                input_url,
                output_expr,
            )
            handle.write('###### job' + str(job_index) + ' ######\n')
            if local_output:
                handle.write(str(job_index) + ', ' + arg_line + '\n')
            else:
                handle.write(arg_line + '\n')
        handle.write(')\n')


def find_submit_files(output_dir):
    submit_files = []
    for root, dirs, files in os.walk(output_dir):
        if 'submit.sh' in files:
            submit_files.append(os.path.join(root, 'submit.sh'))
    return submit_files


def submit_path_parts(submit_path):
    parts = submit_path.replace('\\', '/').split('/')
    ntuple_tag = parts[-5] if len(parts) >= 5 else ''
    sample_tag = parts[-4] if len(parts) >= 4 else submit_path
    return {
        'ntuple_tag': ntuple_tag,
        'sample_tag': sample_tag,
        'canonical_base': sample_tag + '__' + ntuple_tag if ntuple_tag else sample_tag,
        'tag': parts[-3] if len(parts) >= 3 else 'rjrskim',
    }


def expected_submit_prefix(path_info):
    return 'condor_' + path_info['canonical_base'] + '__' + path_info['tag']


def write_multi_submit_script(path, submit_files, action='Submitting', blank_lines=True):
    with open(path, 'w') as fh:
        fh.write('#!/bin/bash\n')
        for sp in submit_files:
            label = '/'.join(sp.replace('\\', '/').split('/')[-4:-2])
            fh.write('echo "{} {}..."\n'.format(action, label))
            fh.write('condor_submit ' + sp + '\n')
            if blank_lines:
                fh.write('echo ""\n')
    os.chmod(path, 0o755)


def get_job_num_from_filename(fname):
    """Extract job number from condor_xxx.N.root -> N, or None."""
    base = os.path.basename(fname)
    if base.endswith('.root'):
        stem  = base[:-5]
        dot   = stem.rfind('.')
        if dot != -1 and stem[dot+1:].isdigit():
            return int(stem[dot+1:])
    return None


def write_resub_file(submit_path, header_lines, missing_entries):
    """Write a resubmission condor file for missing_entries.
    $(Process) is replaced with the actual job index in each argument line
    so output filenames stay consistent with the original submission.
    Returns the resub file path.
    """
    resub_path = submit_path.replace('submit.sh', 'resub.sh')
    local_output = not any('--eos-out' in submit_arg_tokens(arg_line, submit_path)
                           for _, arg_line in missing_entries)
    expected_prefix = expected_submit_prefix(submit_path_parts(submit_path))
    local_output_file = expected_prefix + '.$(GlobalJobIndex).root'
    local_output_dir = (_local_output_dir_from_header(header_lines)
                        if local_output else None)

    with open(resub_path, 'w') as fh:
        for line in header_lines:
            stripped = line.strip()
            if local_output and stripped.startswith('transfer_output_files'):
                fh.write('transfer_output_files = ' + local_output_file + '\n')
            elif local_output and stripped.startswith('transfer_output_remaps'):
                remap_target = os.path.join(local_output_dir, local_output_file)
                fh.write(
                    'transfer_output_remaps = "'
                    + local_output_file + '=' + remap_target + '"\n'
                )
            elif local_output and stripped in {
                    'queue Arguments from (',
                    'queue GlobalJobIndex, Arguments from (',
            }:
                fh.write('queue GlobalJobIndex, Arguments from (\n')
            else:
                fh.write(line + '\n')
        for job_idx, arg_line in missing_entries:
            output_expr = expected_prefix + '.' + str(job_idx)
            resub_arg = replace_output_arg(arg_line.replace('$(Process)', str(job_idx)),
                                           submit_path, output_expr)
            fh.write('###### job' + str(job_idx) + ' ######\n')
            if local_output:
                fh.write(str(job_idx) + ', ' + resub_arg + '\n')
            else:
                fh.write(resub_arg + '\n')
        fh.write(')\n')
    return resub_path


def print_test_command(submit_path, job_entries, job_idx):
    """Print a local run command for the given job index."""
    matching = [(i, a) for i, a in job_entries if i == job_idx]
    if not matching:
        print('  Job', job_idx, 'not found in', submit_path)
        return
    _, arg_line = matching[0]
    # Replace $(Process) and strip --eos-out so output stays local
    arg_line = arg_line.replace('$(Process)', str(job_idx))
    arg_line = re.sub(r'\s*--eos-out\s+\S+', '', arg_line)
    print()
    print('=' * 60)
    print('Local test command for job', job_idx, '(run with CMSSW set up):')
    print('  ./config/runKUCMSAodSVSkimmer.obj ' + arg_line)
    print('=' * 60)


def terminal_color(text, color_code):
    if sys.stdout.isatty():
        return color_code + text + '\033[0m'
    return text


def red_text(text):
    return terminal_color(text, '\033[31m')


def blue_text(text):
    return terminal_color(text, '\033[34m')


def check_status_text(status):
    if status == 'INCOMPLETE':
        return red_text(status)
    return status


def branch_mask_display(mask_path):
    text = mask_path if mask_path else '[None]'
    return blue_text(text) if mask_path else red_text(text)


def submit_label_from_parts(path_info):
    return path_info['canonical_base'] + ' [' + path_info['tag'] + ']'


def filter_submit_files_by_label(submit_files, name_filter):
    matches = []
    skipped = []
    for submit_path in sorted(submit_files):
        path_info = submit_path_parts(submit_path)
        label = submit_label_from_parts(path_info)
        if not name_filter or name_filter in label:
            matches.append((submit_path, path_info, label))
        else:
            skipped.append((submit_path, label))
    return matches, skipped


def check_jobs(output_dir, verbose=False, test_job=None, name_filter=None):
    """Scan output_dir for submit.sh files, report completeness, and write
    resubmission files for any subfolders with missing output."""

    submit_files = find_submit_files(output_dir)
    if not submit_files:
        print('No submit.sh files found under', output_dir)
        return

    filtered_submit_files, _name_skipped = filter_submit_files_by_label(
        submit_files, name_filter)

    print('Found', len(submit_files), 'submit file(s)')
    if name_filter:
        print('Name filter:', name_filter)
        if not filtered_submit_files:
            print(red_text('No submit files matched filter: ' + name_filter))
            return
        print('Matched', len(filtered_submit_files), 'submit file(s)')
    print()

    resub_files        = []
    first_missing_path = None
    first_missing_idx  = None
    first_missing_jobs = None

    for submit_path, path_info, label in filtered_submit_files:
        expected_prefix = expected_submit_prefix(path_info)

        history = recover_submission_history_or_exit(submit_path)

        header_lines = history.header_lines
        eos_out_dir = history.eos_out_dir
        expected_indices = history_expected_indices(history)
        n_total = len(expected_indices)

        # List output files from EOS or local
        if eos_out_dir:
            entries      = xrdfs_ls(eos_out_dir, verbose=verbose)
            output_files = [os.path.basename(e) for e in entries if e.endswith('.root')]
        else:
            out_dir = None
            for line in header_lines:
                m = re.search(r'transfer_output_remaps\s*=\s*"[^=]+=([^"]+)/[^/"]+\.root"', line)
                if m:
                    out_dir = m.group(1)
                    break
            if out_dir and os.path.isdir(out_dir):
                output_files = [f for f in os.listdir(out_dir) if f.endswith('.root')]
            else:
                output_files = []

        classified = classify_transfer_outputs(output_files, expected_prefix,
                                               expected_indices)
        completed = set(classified.output_for_index)
        n_done = len(completed)
        incomplete_indices = sorted(expected_indices - completed)
        n_missing = len(incomplete_indices)
        status     = 'OK' if n_missing == 0 else 'INCOMPLETE'

        print('{}:'.format(label))
        print('  Original jobs:', history.original_count)
        print('  Added jobs:   ', n_total - history.original_count)
        print('  Expected jobs:', n_total)
        print('  Missing jobs: ', '{}/{} [{}]'.format(
            n_missing, n_total, check_status_text(status)))
        if classified.unexpected_indices:
            print('  unexpected jobs:', classified.unexpected_indices)
        if classified.duplicate_indices:
            print('  duplicate jobs:', classified.duplicate_indices)

        if n_missing > 0:
            missing_indices = incomplete_indices
            missing_set = set(missing_indices)
            missing_entries = [
                (job.index, job.arg_line)
                for job in history.jobs
                if job.index in missing_set
            ]
            resub_path = write_resub_file(submit_path, header_lines, missing_entries)
            resub_files.append(resub_path)
            print('  resub file:  ', resub_path)
            if first_missing_path is None:
                first_missing_path = submit_path
                first_missing_idx  = missing_indices[0]
                first_missing_jobs = [(job.index, job.arg_line) for job in history.jobs]
        print()

    if not resub_files:
        print('All jobs complete.')
        return

    if test_job is not None and first_missing_path is not None:
        target_idx = test_job if test_job >= 0 else first_missing_idx
        print_test_command(first_missing_path, first_missing_jobs, target_idx)

    multi = os.path.join(output_dir, 'multi_resub.sh')
    write_multi_submit_script(multi, resub_files, action='Resubmitting', blank_lines=False)
    print('Multi-resub script:', multi)
    print('Run with:  source', multi)


# ---------------------------------------------------------------------------
# New-inputs: submit only ntuple files that were not in the original submission
# ---------------------------------------------------------------------------

def new_inputs_mode(args):
    """Scan EOS for files not present in existing submit.sh files and generate
    immutable newfiles_<first>_<last>.sh submit files for each dataset that has
    new inputs. Job indices are assigned after the full submission history."""

    eos_path    = args.eos_path.rstrip('/')
    odir        = args.output.rstrip('/') + '/'
    sample_type, year, timecali, data_kw = sample_context(eos_path)
    print_sample_context(eos_path, sample_type, year, timecali)

    print('Scanning EOS tree...')
    datasets = find_crab_datasets(eos_path, max_files=args.max_files,
                                  peek_only=args.dry_run, verbose=args.verbose)
    if not datasets:
        print('ERROR: nothing found under', eos_path)
        sys.exit(1)

    print('Datasets found:', len(datasets))
    new_submit_files = []

    for dataset in datasets:
        paths = dataset_paths(args.output, dataset, args.tag)
        root_files = dataset.root_files

        print('\nDataset:', dataset.sample_tag)

        if not os.path.exists(paths.submit_path):
            print('  No existing submit.sh found — run without --new-inputs first, skipping.')
            continue

        history = recover_submission_history_or_exit(paths.submit_path)

        submitted_bare = {job.input_file for job in history.jobs}
        new_files = new_input_files(root_files, history.jobs)

        cli_eos_out_dir = eos_output_dir(args, paths)
        if history.eos_out_dir:
            if getattr(args, '_eos_out_supplied', False) and cli_eos_out_dir:
                norm_cli_eos = normalize_eos_path(cli_eos_out_dir)
                if norm_cli_eos != history.eos_out_dir:
                    print('ERROR: --eos-out differs from existing production output.')
                    print('  Existing:', history.eos_out_dir)
                    print('  CLI:     ', norm_cli_eos)
                    sys.exit(1)
            eos_out_dir = history.eos_out_dir
        else:
            eos_out_dir = None

        print('  EOS files now    :', len(root_files))
        print('  Already submitted:', len(submitted_bare))
        print('  New files        :', len(new_files))

        if not new_files:
            print('  Nothing new, skipping.')
            continue

        first_index = max(history_expected_indices(history)) + 1
        last_index = first_index + len(new_files) - 1
        newfiles_path = newfiles_submit_path(paths.src_dir, first_index, last_index)
        expected_prefix = expected_submit_prefix(submit_path_parts(paths.submit_path))

        print('  New job indices  :', '{}-{}'.format(first_index, last_index))
        print('  Newfiles submit  :', newfiles_path)
        print('  Configuration    : inherited from', paths.submit_path)
        if eos_out_dir:
            print('  EOS out          :', eos_out_dir)
            print('  Output mode      : direct EOS')
        else:
            print('  Output mode      : local Condor transfer')

        if args.dry_run:
            print('  [DRY RUN]')
            new_submit_files.append(newfiles_path)
            continue

        if eos_out_dir:
            xrdfs_mkdir(eos_out_dir, verbose=args.verbose)

        ensure_submit_dirs(paths, eos_out_dir)
        if os.path.exists(newfiles_path):
            print('ERROR: history file already exists, refusing to overwrite:',
                  newfiles_path)
            sys.exit(1)

        write_inherited_new_inputs_submit(
            newfiles_path,
            history,
            expected_prefix,
            new_files,
            first_index,
        )
        new_submit_files.append(newfiles_path)

    print()
    if not new_submit_files:
        print('No new input files found across all datasets.')
        return

    print('=' * 60)
    if args.dry_run:
        print('[DRY RUN]', len(new_submit_files),
              'newfiles submit file(s) would be generated')
    elif len(new_submit_files) == 1:
        print('To submit:')
        print('  condor_submit', new_submit_files[0])
    else:
        multi = odir + args.tag + '_NewFiles_MultiSub.sh'
        write_multi_submit_script(multi, new_submit_files)
        print('Multi-submit script:', multi)
        print('Run with:  source', multi)
    print('=' * 60)


# ---------------------------------------------------------------------------
# Transfer: hadd per-job EOS output, xrdcp to final skim destination
# ---------------------------------------------------------------------------

def require_transfer_tools():
    for tool in ('hadd', 'xrdcp'):
        try:
            subprocess.check_call(['which', tool],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('ERROR:', tool, 'not found. Run with CMSSW set up.')
            sys.exit(1)


def ensure_transfer_destination(dest_base, verbose=False):
    if not xrdfs_stat(dest_base, verbose=verbose):
        print('Destination folder not found — creating:', dest_base)
        if not xrdfs_mkdir(dest_base, verbose=verbose):
            print('ERROR: could not create destination folder:', dest_base)
            sys.exit(1)
        print('Created:', dest_base)
    else:
        print('Destination folder exists:', dest_base)
    print()


def merged_output_name(path_info):
    return final_merged_filename(path_info['sample_tag'], path_info['ntuple_tag'],
                                 path_info['tag'])


TransferOutputClassification = namedtuple('TransferOutputClassification', [
    'output_for_index',
    'selected_outputs',
    'missing_indices',
    'unexpected_indices',
    'duplicate_indices',
    'unrelated_root_files',
])


def classify_transfer_outputs(entries, expected_prefix, expected_indices):
    pattern = re.compile(r'^' + re.escape(expected_prefix) + r'\.(\d+)\.root$')
    paths_by_index = {}
    unrelated_root_files = []

    for entry in entries:
        base = os.path.basename(entry)
        match = pattern.match(base)
        if match:
            index = int(match.group(1))
            paths_by_index.setdefault(index, []).append(entry)
        elif base.endswith('.root'):
            unrelated_root_files.append(entry)

    duplicate_indices = sorted(index for index, paths in paths_by_index.items()
                               if len(paths) > 1)
    found_indices = set(paths_by_index)
    missing_indices = sorted(expected_indices - found_indices)
    unexpected_indices = sorted(found_indices - expected_indices)
    duplicate_set = set(duplicate_indices)
    output_for_index = {
        index: paths[0]
        for index, paths in paths_by_index.items()
        if index in expected_indices and index not in duplicate_set and len(paths) == 1
    }
    selected_outputs = [
        output_for_index[index]
        for index in sorted(expected_indices)
        if index in output_for_index
    ]

    return TransferOutputClassification(
        output_for_index=output_for_index,
        selected_outputs=selected_outputs,
        missing_indices=missing_indices,
        unexpected_indices=unexpected_indices,
        duplicate_indices=duplicate_indices,
        unrelated_root_files=unrelated_root_files,
    )


def cleanup_existing(paths):
    for path in paths:
        if os.path.exists(path):
            os.remove(path)


def stage_root_files(root_files, scratch):
    local_inputs = []
    for f in root_files:
        src_url = eos_xrootd_url(f)
        local_f = os.path.join(scratch, os.path.basename(f))
        result  = subprocess.run(['xrdcp', '-s', src_url, local_f])
        if result.returncode != 0:
            print('  WARNING: xrdcp staging failed for', os.path.basename(f))
            cleanup_existing(local_inputs)
            return None
        local_inputs.append(local_f)
    return local_inputs


def hadd_root_files(output_path, input_files, scratch, force=False):
    hadd_flags = ['-f'] if force and os.path.exists(output_path) else []
    result = subprocess.run(['hadd', '-k', '-j', '4', '-v', '0', '-d', scratch]
                            + hadd_flags + [output_path] + input_files)
    return result.returncode == 0


def merge_root_files(root_files, hadded, scratch, merged_name, chunk_size, force=False):
    use_chunks = chunk_size > 0 and len(root_files) > chunk_size

    if os.path.exists(hadded) and not force:
        print('  Local hadded file already exists, skipping hadd (use --force to redo)')
        return True

    if use_chunks:
        chunks   = [root_files[i:i+chunk_size] for i in range(0, len(root_files), chunk_size)]
        partials = []
        print('  Chunked merge: {} files in {} chunk(s) of {}'.format(
              len(root_files), len(chunks), chunk_size))
        for ci, chunk in enumerate(chunks):
            partial = os.path.join(scratch, merged_name.replace('.root', '.part{}.root'.format(ci)))
            partials.append(partial)
            print('  Chunk {}/{}: staging {} file(s)...'.format(ci+1, len(chunks), len(chunk)))
            local_inputs = stage_root_files(chunk, scratch)
            if local_inputs is None:
                cleanup_existing(partials)
                return False
            ok = hadd_root_files(partial, local_inputs, scratch)
            cleanup_existing(local_inputs)
            if not ok:
                print('  WARNING: hadd failed on chunk', ci)
                cleanup_existing(partials)
                return False

        print('  Final merge of', len(partials), 'partial(s) ->', hadded)
        ok = hadd_root_files(hadded, partials, scratch, force=force)
        cleanup_existing(partials)
        if not ok:
            print('  WARNING: final hadd failed -- skipping')
        return ok

    # Stage all files at once then hadd locally
    print('  Staging', len(root_files), 'file(s) to scratch...')
    local_inputs = stage_root_files(root_files, scratch)
    if local_inputs is None:
        return False
    print('  hadd:', len(local_inputs), 'file(s) ->', hadded)
    ok = hadd_root_files(hadded, local_inputs, scratch, force=force)
    cleanup_existing(local_inputs)
    if not ok:
        print('  WARNING: hadd failed -- skipping')
    return ok


def xrdcp_to_eos(local_path, dest_url, force=False):
    xrdcp_flags = ['-f'] if force else []
    xrdcp_cmd   = ['xrdcp', '-s'] + xrdcp_flags + [local_path, dest_url]
    print('  xrdcp ->', dest_url)
    result = subprocess.run(xrdcp_cmd)
    return result.returncode == 0


def transfer_jobs(args, name_filter=None):
    """For each submit.sh, hadd per-job output files from EOS and xrdcp to
    /store/group/lpcsusylep/malazaro/KUCMSSkims/{version}/."""

    dest_base = EOS_SKIMS_BASE.rstrip('/') + '/' + args.version.strip('/')
    dest_xrd  = eos_xrootd_url(dest_base)
    submit_files = find_submit_files(args.output)
    if not submit_files:
        print('No submit.sh files found under', args.output)
        return

    filtered_submit_files, _name_skipped = filter_submit_files_by_label(
        submit_files, name_filter)

    print('Found', len(submit_files), 'submit file(s)')
    print('Destination:', dest_xrd)
    print('Tag filter:', args.tag if args.tag else 'all')
    if name_filter:
        print('Name filter:', name_filter)
        if not filtered_submit_files:
            print(red_text('No submit files matched filter: ' + name_filter))
            return
        print('Matched', len(filtered_submit_files), 'submit file(s)')
    print()

    require_transfer_tools()

    if not args.dry_run:
        ensure_transfer_destination(dest_base, verbose=args.verbose)

    n_ok = 0
    n_skip = 0
    n_fail = 0
    n_tag_skip = 0
    tag_skip_examples = []
    for submit_path, path_info, label in filtered_submit_files:
        tag       = path_info['tag']
        if args.tag and tag != args.tag:
            n_tag_skip += 1
            if len(tag_skip_examples) < 5:
                tag_skip_examples.append(submit_path)
            continue
        history = recover_submission_history_or_exit(submit_path)
        eos_out_dir = history.eos_out_dir
        if not eos_out_dir:
            print('SKIP (no --eos-out found):', submit_path)
            n_skip += 1
            continue
        expected_indices = history_expected_indices(history)

        sample_tag  = path_info['sample_tag']
        merged_name = merged_output_name(path_info)
        dest_path   = dest_base.rstrip('/') + '/' + merged_name
        dest_url    = eos_xrootd_url(dest_path)

        print('Sample:', sample_tag)
        print('  Merged :', merged_name)
        print('  Dest   :', dest_url)

        expected_prefix = expected_submit_prefix(path_info)
        entries = xrdfs_ls(eos_out_dir, verbose=args.verbose)
        classified = classify_transfer_outputs(entries, expected_prefix, expected_indices)
        root_files = [
            classified.output_for_index[index]
            for index in sorted(expected_indices)
            if index in classified.output_for_index
        ]

        print('  Original jobs:', history.original_count)
        print('  Added jobs:   ', len(expected_indices) - history.original_count)
        print('  Expected jobs:', len(expected_indices))
        print('  Completed expected jobs:', len(root_files))
        if classified.missing_indices:
            print('  Missing indices:', classified.missing_indices)
        if classified.unexpected_indices:
            print('  Unexpected indices:', classified.unexpected_indices)
        if classified.duplicate_indices:
            print('  Duplicate indices:', classified.duplicate_indices)
        if classified.unrelated_root_files:
            print('  Ignoring unrelated ROOT files:', len(classified.unrelated_root_files))

        if (not history.jobs or classified.missing_indices
                or classified.unexpected_indices or classified.duplicate_indices):
            print('  ERROR: incomplete or inconsistent per-job output set -- not transferring')
            n_fail += 1
            print()
            continue

        if not args.force and xrdfs_stat(dest_path, verbose=args.verbose):
            print('  Already on EOS, skipping (use --force to overwrite)')
            n_skip += 1
            print()
            continue

        print('  Files  :', len(root_files))

        scratch = os.path.expandvars(args.scratch).rstrip('/')
        hadded  = os.path.join(scratch, merged_name)

        if args.dry_run:
            print('  [dry-run] hadd', len(root_files), 'file(s) ->', hadded)
            print('  [dry-run] xrdcp ->', dest_url)
            print()
            continue

        if not os.path.isdir(scratch):
            os.makedirs(scratch, exist_ok=True)

        if not merge_root_files(root_files, hadded, scratch, merged_name,
                                args.chunk_size, force=args.force):
            n_fail += 1
            print()
            continue

        if not xrdcp_to_eos(hadded, dest_url, force=args.force):
            print('  WARNING: xrdcp failed')
            n_fail += 1
            print()
            continue

        os.remove(hadded)
        print('  [OK] Transfer complete.')
        print('  Removed scratch files.')

        if args.clean:
            print('  Removing', len(root_files), 'per-job EOS files...')
            for f in root_files:
                subprocess.run(['xrdfs', EOS_SERVER, 'rm', f],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print('  Done.')

        n_ok += 1
        print()

    print('=' * 60)
    print('Transfer complete: {} merged, {} skipped, {} failed'.format(n_ok, n_skip, n_fail))
    if n_tag_skip:
        print('Ignored {} submit file(s) that did not match job tag {}'.format(
              n_tag_skip, args.tag))
        for path in tag_skip_examples:
            print('  ignored:', path)
    print('=' * 60)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate HTCondor skimmer submissions from an EOS path.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--new-inputs', dest='new_inputs', action='store_true',
        help='Generate newfiles_<first>_<last>.sh for each dataset containing '
             'only EOS files not present in the submission history '
             '(requires --eos-path)')
    parser.add_argument('--check', nargs='?', const='', default=None, metavar='FILTER',
        help='Check job completeness and generate resub files for missing output; '
             'optionally filter submit files by job-name substring')
    parser.add_argument('--test-job', dest='test_job', type=int, nargs='?', const=-1, default=None,
        metavar='N', help='With --check: print local run command for job N (omit N for first missing)')
    parser.add_argument('--transfer', nargs='?', const='', default=None, metavar='FILTER',
        help='hadd per-job EOS output per sample and xrdcp to final skim destination; '
             'optionally filter submit files by job-name substring')
    parser.add_argument('--version', default=None,
        metavar='N', help='Skim version number; appends _v{N} to the job tag during submission '
                          'and sets the skims_v{N} destination during --transfer')
    parser.add_argument('--scratch', default='/uscmst1b_scratch/lpc1/3DayLifetime/$USER',
        metavar='DIR', help='With --transfer: local scratch directory for temporary hadded files')
    parser.add_argument('--force', action='store_true',
        help='With --transfer: redo hadd/xrdcp even if output already exists')
    parser.add_argument('--clean', action='store_true',
        help='With --transfer: remove per-job EOS files after successful xrdcp')
    parser.add_argument('--chunk-size', dest='chunk_size', type=int, default=0,
        metavar='N', help='With --transfer: merge in chunks of N files (default: 0 = disabled)')
    parser.add_argument('--eos-path', '-e', default=None,
        help='EOS directory to scan (bare /eos/... or root://... URL); required unless --check or --transfer')
    parser.add_argument('--output', '-o', default='Output',
        help='Output directory (default: Output)')
    parser.add_argument('--tag', '-t', default=None,
        help='Job name tag (default: rjrskim for submission; all tags for --transfer)')
    parser.add_argument('--sandbox', default=SANDBOX_DEFAULT,
        help='CMSSW sandbox tarball path on worker node')
    parser.add_argument('--max-mat',  dest='max_mat',  type=int, default=-1,
        help='max_materialize condor option')
    parser.add_argument('--max-idle', dest='max_idle', type=int, default=-1,
        help='max_idle condor option')
    parser.add_argument('--memory',   type=int, default=-1,
        help='request_memory in MB')
    parser.add_argument('--max-files', dest='max_files', type=int, default=-1,
        help='Max .root files per subfolder (useful for testing)')
    parser.add_argument('--no-sv',   action='store_true', help='Pass --noSV to skimmer')
    parser.add_argument('--hlt-off', action='store_true', help='Pass --HLTPathsOff to skimmer')
    parser.add_argument('--branch-mask', dest='branch_mask', default='',
        help='Branch mask file passed to skimmer, e.g. config/branch_masks/sv_analysis_core.txt')
    parser.add_argument('--psiche',  action='store_true', help='Enable PSICHE jets (default: off)')
    parser.add_argument('--eos-out', dest='eos_out', default='/eos/uscms/store/user/$USER/LLPSkims',
        help='Write output directly to this EOS path (default: /eos/uscms/store/user/$USER/LLPSkims)')
    parser.add_argument('--dry-run', action='store_true', help='Print plan, write nothing')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print xrdfs commands')
    args = parser.parse_args()
    args._eos_out_supplied = any(
        item == '--eos-out' or item.startswith('--eos-out=')
        for item in sys.argv[1:]
    )

    if args.branch_mask:
        try:
            original_branch_mask = args.branch_mask
            args.branch_mask = normalize_branch_mask_path(args.branch_mask)
            if args.branch_mask != original_branch_mask:
                print('Branch mask:', original_branch_mask, '->', args.branch_mask)
        except ValueError as exc:
            parser.error(str(exc))

    if args.new_inputs:
        if not args.eos_path:
            parser.error('--new-inputs requires --eos-path')
        if args.tag is None:
            args.tag = 'rjrskim'
        if args.version:
            args.tag = args.tag + '_v' + args.version.lstrip('v')
        new_inputs_mode(args)
        return

    if args.check is not None:
        check_jobs(args.output, verbose=args.verbose, test_job=args.test_job,
                   name_filter=args.check)
        return

    if args.transfer is not None:
        if args.version is None:
            detected = None
            for submit_path in find_submit_files(args.output):
                tag_part = submit_path_parts(submit_path)['tag']
                m = re.search(r'_v(\d+)$', tag_part)
                if m:
                    detected = m.group(1)
                    break
            if detected:
                version_str = 'skims_v' + detected
                ans = input('Auto-detected version {}. Proceed? [Y/n] '.format(version_str)).strip().lower()
            else:
                version_str = 'skims_v49'
                ans = input('No --version specified. Proceed with {} ? [Y/n] '.format(version_str)).strip().lower()
            if ans not in ('', 'y', 'yes'):
                print('Aborted.')
                return
            args.version = version_str
        else:
            args.version = 'skims_v{}'.format(args.version.lstrip('v'))
        transfer_jobs(args, name_filter=args.transfer)
        return

    if not args.eos_path:
        parser.error('--eos-path is required unless --check or --transfer is specified')

    if args.tag is None:
        args.tag = 'rjrskim'

    if args.version:
        args.tag = args.tag + '_v' + args.version.lstrip('v')

    eos_path    = args.eos_path.rstrip('/')
    odir        = args.output.rstrip('/') + '/'
    sample_type, year, timecali, data_kw = sample_context(eos_path)
    print_sample_context(eos_path, sample_type, year, timecali)

    # -----------------------------------------------------------------------
    # Discover datasets
    # -----------------------------------------------------------------------
    print('Scanning EOS tree...')
    datasets = find_crab_datasets(eos_path, max_files=args.max_files,
                                  peek_only=args.dry_run, verbose=args.verbose)
    if not datasets:
        print('ERROR: nothing found under', eos_path)
        print('Check the path and that xrdfs is available on this system.')
        sys.exit(1)

    print('Datasets found:', len(datasets))
    event_count_keys = read_event_count_keys() if sample_type != 'data' else set()

    # -----------------------------------------------------------------------
    # Process each dataset
    # -----------------------------------------------------------------------
    weights      = {}
    submit_files = []

    for dataset in datasets:
        root_files = dataset.root_files
        paths = dataset_paths(args.output, dataset, args.tag)
        eos_out_dir = eos_output_dir(args, paths)
        final_skim = final_merged_filename(dataset.sample_tag, dataset.ntuple_tag,
                                           args.tag)

        print('\nDataset:', dataset.sample_tag)

        if not root_files and not args.dry_run:
            print('  No .root files found, skipping.')
            continue

        dataset_data_kw = data_keyword(dataset.raw_task_suffix) if sample_type == 'data' else data_kw
        metadata = resolve_dataset_metadata(dataset.raw_task_suffix, dataset.task_name,
                                            sample_type, dataset_data_kw)
        flags    = build_skimmer_flags(metadata, args, timecali)

        print('  Ntuple tag:     ', dataset.ntuple_tag)
        print('  Skim tag:       ', args.tag)
        print('  Branch mask:    ', branch_mask_display(args.branch_mask))
        print('  Input files:    ', len(root_files))
        print('  Jobs:           ', len(root_files))
        print()
        print('  Submit file:    ', paths.submit_path)
        if eos_out_dir:
            print('  Job output dir: ', eos_out_dir)
        else:
            print('  Job output dir: ', paths.out_dir)
        print('  Job output name:', paths.ofilename + '.<job>.root')
        print('  Final skim:     ', final_skim)
        print()
        print('  Metadata:       ', 'xsec=' + metadata['xsec'] + ', key=' + metadata['key'])
        warn_missing_event_count_key(metadata, event_count_keys)

        if args.verbose:
            print()
            print('  Collection:     ', dataset.collection_dir)
            print('  Task:           ', dataset.task_dir)
            print('  Timestamp:      ', dataset.timestamp_dir)
            print('  Raw task suffix:', dataset.raw_task_suffix)
            print('  Canonical base: ', dataset.canonical_base)
            print('  Work directory: ', paths.work_dir)

        # --- record weights ---
        weights[metadata['key']] = weight_record(metadata, timecali, root_files,
                                                 dataset.raw_task_suffix)

        submit_files.append(paths.submit_path)

        if args.dry_run:
            print('  [DRY RUN]')
            continue

        if eos_out_dir:
            print('  EOS out:', eos_out_dir)
            xrdfs_mkdir(eos_out_dir, verbose=args.verbose)

        ensure_submit_dirs(paths, eos_out_dir)

        write_submit(paths.submit_path, paths.log_dir, paths.out_dir, paths.ofilename,
                     root_files, flags, args, eos_out_dir=eos_out_dir)

    # -----------------------------------------------------------------------
    # Save weights and print instructions
    # -----------------------------------------------------------------------
    if not args.dry_run and submit_files:
        os.makedirs(odir, exist_ok=True)
        with open(odir + 'weights.json', 'w') as fh:
            json.dump(weights, fh, indent=2)
        with open(odir + 'weights.py', 'w') as fh:
            fh.write('# Auto-generated by skim_submit.py\n')
            fh.write('# xsec in pb. mctype: 0=MC, 1=data, 2=MC fast.\n\n')
            fh.write('WEIGHTS = ' + json.dumps(weights, indent=2) + '\n')
        print('\nWeights saved to:', odir + 'weights.json')

    if not submit_files:
        print('\nNothing generated.')
        return

    print()
    print('=' * 60)
    if len(submit_files) == 1:
        print('To submit:')
        print('  condor_submit', submit_files[0])
    else:
        multi = odir + args.tag + '_MultiSub.sh'
        if not args.dry_run:
            write_multi_submit_script(multi, submit_files)
            print('Multi-submit script:', multi)
            print('Run with:  source', multi)
        else:
            print('[DRY RUN]', len(submit_files), 'submit files would be generated')
    print('=' * 60)


if __name__ == '__main__':
    main()
