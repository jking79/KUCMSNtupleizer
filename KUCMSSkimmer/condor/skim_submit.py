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
import subprocess
import sys

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
}

DATA_KEYWORDS  = {'MET', 'JetMET', 'JetMET0', 'JetMET1', 'JetHT', 'EGamma', 'EGamma0', 'EGamma1', 'EGamma2', 'DisJet'}
SIGNAL_KEYWORDS = {'SMS', 'gogoG', 'gogoZ', 'gogoGZ', 'sqsqG', 'GlGl'}

SANDBOX_DEFAULT  = "/uscms/home/mlazarov/nobackup/sandboxes/sandbox-CMSSW_13_3_3.tar.bz2"
#EOS_SKIMS_BASE   = "/store/group/lpcsusylep/malazaro/KUCMSSkims"
EOS_SKIMS_BASE   = "/store/user/malazaro/LLPSkims"
EOS_SERVER       = "root://cmseos.fnal.gov"

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
    for yr in ['R16', 'R17', 'R18', 'R22', 'R23', 'R24']:
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


def normalize_signal_ctau(raw_ctau):
    if raw_ctau.startswith('-0p') and len(raw_ctau) == 4:
        return raw_ctau[3:]
    value = raw_ctau.lstrip('-')
    if value.startswith('0p') and len(value) == 3:
        return value
    return value


def signal_event_count_key(name):
    mgl = re.search(r'mGl-(\d+)', name)
    mn2 = re.search(r'mN2-(\d+)', name)
    mn1 = re.search(r'mN1-(\d+)', name)
    ctau = re.search(r'(?:^|[_-])ct(-?[0-9]+p[0-9]+|-?[0-9]+)', name)

    if not (mgl and mn2 and mn1 and ctau):
        match = re.search(
            r'((?:gogoGZ|gogoG|gogoZ|sqsqG)_(?:AODSIM|FULLMINI|MINIAOD|MINI|AOD)_.*)',
            name,
        )
        if match:
            key = re.sub(r'_(FULLMINI|MINIAOD|MINI|AOD)_',
                         '_AODSIM_', match.group(1))
            return re.sub(r'_ct0p([15])$', r'_ct\1', key)
        raise ValueError('could not derive signal EventCount key from: ' + name)

    process = signal_process_key(name)
    ct_key = normalize_signal_ctau(ctau.group(1))
    return (
        f'{process}_FULLMINI_mGl-{mgl.group(1)}'
        f'_mN2-{mn2.group(1)}_mN1-{mn1.group(1)}_ct{ct_key}'
    )


def make_data_key(subfolder, sample_type_kw):
    """Build a short dataset key for data from the subfolder name."""
    # e.g. 'MET_AOD_Run2018A-15Feb...' -> 'MET18A'
    year_m = re.search(r'(\d{4})', subfolder)
    era_m  = re.search(r'Run\d{4}([A-Z])', subfolder)
    yr  = year_m.group(1)[2:] if year_m else ''
    era = era_m.group(1)       if era_m  else ''
    return sample_type_kw + yr + era


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
    if fpath.startswith('root://'):
        m = re.match(r'root://[^/]+/(/.+)', fpath)
        if m:
            return m.group(1)
        m = re.match(r'root://[^/]+/(.+)', fpath)
        if m:
            return '/' + m.group(1)
    if fpath.startswith('/eos/uscms'):
        return fpath[len('/eos/uscms'):]
    return '/' + fpath.lstrip('/')


def collection_tag_from_url(url):
    """Extract collection tag from an XRootD input URL.
    e.g. root://...//store/.../kucmsntuple_QCD_R18_SVIPM100_MiniAOD_v34/... -> QCD_R18_SVIPM100_MiniAOD_v34
    Returns empty string if not found."""
    m = re.search(r'kucmsntuple_([^/]+)', url)
    return m.group(1) if m else ''


_CRAB_TS_RE = re.compile(r'^\d{6}_\d{6}$')


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
    """Walk the EOS tree and return [(group_name, [root_files])] per dataset.

    group_name accumulates the basenames of directories traversed below the
    entry point, so it reflects scope regardless of where --eos-path points:
      year level   -> 'JetMET0_kucmsntuple_..._Run2024C-PromptReco-v1'
      stream level -> 'kucmsntuple_..._Run2024C-PromptReco-v1'
      era level    -> 'kucmsntuple_..._Run2024C-PromptReco-v1'

    Dirs whose children are all CRAB timestamps (YYMMDD_HHMMSS) are treated as
    dataset dirs; only the latest timestamp is used.  peek_only skips deep
    recursion for dry-run estimates.
    """
    entries    = xrdfs_ls(path, verbose=verbose)
    root_files = [e for e in entries if e.endswith('.root')]
    subdirs    = [e for e in entries if not e.endswith('.root')
                  and '.' not in os.path.basename(e)]

    group = name or os.path.basename(path)

    # Case 1: .root files live directly here
    if root_files:
        files = root_files[:max_files] if max_files > 0 else root_files
        return [(group, files)]

    if not subdirs:
        return []

    # Case 2: all subdirs are CRAB timestamps -> dataset dir, use only latest
    if all(_CRAB_TS_RE.match(os.path.basename(e)) for e in subdirs):
        latest  = max(subdirs, key=os.path.basename)
        skipped = sorted(os.path.basename(e) for e in subdirs if e != latest)
        if skipped:
            print('  [crab] skipping:', ', '.join(skipped),
                  '-> using:', os.path.basename(latest))
        if peek_only:
            files = [e for e in xrdfs_ls(latest, verbose=verbose) if e.endswith('.root')]
        else:
            files = find_root_files(latest, max_files=max_files, verbose=verbose)
        return [(group, files)]

    # Case 3: intermediate dirs -> recurse, accumulating name
    results = []
    for sd in subdirs:
        child_name = (name + '_' + os.path.basename(sd)) if name else os.path.basename(sd)
        results.extend(find_crab_datasets(sd, name=child_name,
                                          max_files=max_files, peek_only=peek_only,
                                          verbose=verbose))
    return results


def write_submit(path, log_dir, out_dir, ofilename, root_files, flags, args, eos_out_dir=None, offset=0):
    outname = ofilename + '.$(Process).root'

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
        'output = ' + os.path.join(log_dir, 'job.$(Process).out'),
        'error = '  + os.path.join(log_dir, 'job.$(Process).err'),
        'log = '    + os.path.join(log_dir, 'job.log'),
        'transfer_input_files = ' + args.sandbox + ', config.tgz,',
        'should_transfer_files = YES',
    ] + output_lines
    if args.max_mat  != -1: lines.append('max_materialize = ' + str(args.max_mat))
    if args.max_idle != -1: lines.append('max_idle = '        + str(args.max_idle))
    if args.memory   != -1: lines.append('request_memory = '  + str(args.memory))

    lines += ['', '', 'queue Arguments from (']
    for i, fpath in enumerate(root_files):
        # fpath is a bare EOS path; build the xrootd URL without the /eos/uscms
        # prefix so the skimmer's eosdir guard (//store/...) matches correctly
        if fpath.startswith('root://'):
            full = fpath
        else:
            bare = fpath
            if bare.startswith('/eos/uscms'):
                bare = bare[len('/eos/uscms'):]
            full = 'root://cmseos.fnal.gov/' + bare
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
    with open(resub_path, 'w') as fh:
        for line in header_lines:
            fh.write(line + '\n')
        for job_idx, arg_line in missing_entries:
            fh.write('###### job' + str(job_idx) + ' ######\n')
            fh.write(arg_line.replace('$(Process)', str(job_idx)) + '\n')
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


def check_jobs(output_dir, verbose=False, test_job=None):
    """Scan output_dir for submit.sh files, report completeness, and write
    resubmission files for any subfolders with missing output."""

    submit_files = []
    for root, dirs, files in os.walk(output_dir):
        if 'submit.sh' in files:
            submit_files.append(os.path.join(root, 'submit.sh'))

    if not submit_files:
        print('No submit.sh files found under', output_dir)
        return

    print('Found', len(submit_files), 'submit file(s)')
    print()

    resub_files        = []
    first_missing_path = None
    first_missing_idx  = None
    first_missing_jobs = None

    for submit_path in sorted(submit_files):
        header_lines, job_entries, eos_out_dir = parse_submit_file(submit_path)
        n_total = len(job_entries)

        # Sample name is 4 levels up from submit.sh
        # e.g. Output/QCD/QCD_HT100to200_.../rjrskim/src/submit.sh
        parts    = submit_path.replace('\\', '/').split('/')
        subfolder = parts[-4] if len(parts) >= 4 else submit_path

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

        completed  = {get_job_num_from_filename(f) for f in output_files} - {None}
        n_done     = len(completed)
        n_missing  = n_total - n_done
        status     = 'OK' if n_missing == 0 else 'INCOMPLETE'

        print('{}: {}/{} missing  [{}]'.format(subfolder, n_missing, n_total, status))

        if n_missing > 0:
            missing_indices = sorted(set(range(n_total)) - completed)
            missing_entries = [(i, a) for i, a in job_entries if i in set(missing_indices)]
            resub_path = write_resub_file(submit_path, header_lines, missing_entries)
            resub_files.append(resub_path)
            print('  missing jobs:', missing_indices)
            print('  resub file:  ', resub_path)
            if first_missing_path is None:
                first_missing_path = submit_path
                first_missing_idx  = missing_indices[0]
                first_missing_jobs = job_entries

    print()
    if not resub_files:
        print('All jobs complete.')
        return

    if test_job is not None and first_missing_path is not None:
        target_idx = test_job if test_job >= 0 else first_missing_idx
        print_test_command(first_missing_path, first_missing_jobs, target_idx)

    multi = os.path.join(output_dir, 'multi_resub.sh')
    with open(multi, 'w') as fh:
        fh.write('#!/bin/bash\n')
        for rp in resub_files:
            label = '/'.join(rp.replace('\\', '/').split('/')[-4:-2])
            fh.write('echo "Resubmitting {}..."\n'.format(label))
            fh.write('condor_submit {}\n'.format(rp))
    os.chmod(multi, 0o755)
    print('Multi-resub script:', multi)
    print('Run with:  source', multi)


# ---------------------------------------------------------------------------
# New-inputs: submit only ntuple files that were not in the original submission
# ---------------------------------------------------------------------------

def new_inputs_mode(args):
    """Scan EOS for files not present in existing submit.sh files and generate
    newfiles.sh submit files for each dataset that has new inputs.
    Job indices are offset past the original submission so output filenames
    (condor_SAMPLE_tag.N.root) never collide."""

    eos_path    = args.eos_path.rstrip('/')
    odir        = args.output.rstrip('/') + '/'
    sample_type = detect_sample_type(eos_path)
    year        = detect_year(eos_path)
    timecali    = get_timecali(year, sample_type)
    if not timecali and sample_type == 'signal':
        timecali = 'r2_ul18_mc'

    print('EOS path:    ', eos_path)
    print('Sample type: ', sample_type)
    print('Year:        ', year)
    print('Timecali:    ', timecali)
    print()

    data_kw = next((kw for kw in sorted(DATA_KEYWORDS, key=len, reverse=True)
                    if kw in eos_path), 'Data')

    print('Scanning EOS tree...')
    datasets = find_crab_datasets(eos_path, max_files=args.max_files,
                                  peek_only=args.dry_run, verbose=args.verbose)
    if not datasets:
        print('ERROR: nothing found under', eos_path)
        sys.exit(1)

    print('Datasets found:', len(datasets))
    event_count_keys = read_event_count_keys() if sample_type != 'data' else set()

    new_submit_files = []

    for subfolder, root_files in datasets:
        prefix      = subfolder.split('_')[0]
        work_dir    = odir + prefix + '/' + subfolder + '/' + args.tag
        src_dir     = work_dir + '/src'
        log_dir     = work_dir + '/log'
        out_dir     = work_dir + '/out'
        submit_path = src_dir  + '/submit.sh'

        print('\nDataset:', subfolder)

        if not os.path.exists(submit_path):
            print('  No existing submit.sh found — run without --new-inputs first, skipping.')
            continue

        header_lines, job_entries, eos_out_dir_parsed = parse_submit_file(submit_path)

        # Build set of already-submitted files as normalized bare /store/... paths
        submitted_bare = set()
        for _, arg_line in job_entries:
            m = re.search(r'-i\s+(\S+)', arg_line)
            if m:
                submitted_bare.add(normalize_eos_path(m.group(1)))

        new_files = [f for f in root_files
                     if normalize_eos_path(f) not in submitted_bare]

        print('  EOS files now    :', len(root_files))
        print('  Already submitted:', len(submitted_bare))
        print('  New files        :', len(new_files))

        if not new_files:
            print('  Nothing new, skipping.')
            continue

        # New job indices start immediately after the original submission
        offset = len(job_entries)

        # --- metadata lookup (mirrors main loop) ---
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
            sig_xsec, gluinomass, n2mass, n1mass = lookup_signal_xsec(subfolder)
            if sig_xsec is not None:
                xsec = str(sig_xsec)
            else:
                print('  WARNING: no xsec found for', subfolder, '- using 0')
                xsec = '0'
            key = signal_event_count_key(subfolder)
        else:
            meta = lookup_bg_meta(subfolder)
            if meta:
                xsec = str(meta['xsec'])
                key  = meta['key']
            else:
                print('  WARNING: no xsec found for', subfolder, '- using 0')
                xsec = '0'

        flags = (
            ' --xsec '        + xsec
            + ' --dataSetKey ' + key
            + ' --gluinoMass ' + gluinomass
            + ' --N2Mass '     + n2mass
            + ' --timeCaliTag ' + timecali
            + ' --MCweight '   + mc_wt
            + ' --MCtype '     + str(mctype)
        )
        if mctype == 0:
            flags += ' --hasGenInfo'
        if not args.psiche:
            flags += ' --noBHC'
        if args.no_sv:
            flags += ' --noSV'
        if args.hlt_off:
            flags += ' --HLTPathsOff'

        print('  xsec:', xsec, '| key:', key)
        if mctype != 1 and event_count_keys and key not in event_count_keys:
            print('  WARNING: key not found in config/EventCount.txt:', key)
            print('           This job will likely get inf evtFillWgt.')

        # EOS output dir: prefer --eos-out from CLI, fall back to parsed original
        eos_out_dir = None
        if args.eos_out:
            eos_out_dir = (os.path.expandvars(args.eos_out).rstrip('/')
                           + '/' + prefix + '/' + subfolder + '/' + args.tag)
        elif eos_out_dir_parsed:
            eos_out_dir = eos_out_dir_parsed

        newfiles_path = src_dir + '/newfiles.sh'
        ofilename     = 'condor_' + subfolder + '_' + args.tag

        print('  Job offset       :', offset, '(indices start at', offset, ')')
        print('  Newfiles submit  :', newfiles_path)
        if eos_out_dir:
            print('  EOS out          :', eos_out_dir)

        if args.dry_run:
            print('  [DRY RUN]')
            new_submit_files.append(newfiles_path)
            continue

        if eos_out_dir:
            xrdfs_mkdir(eos_out_dir, verbose=args.verbose)

        local_dirs = [src_dir, log_dir] if eos_out_dir else [src_dir, log_dir, out_dir]
        for d in local_dirs:
            os.makedirs(d, exist_ok=True)

        write_submit(newfiles_path, log_dir, out_dir, ofilename,
                     new_files, flags, args, eos_out_dir=eos_out_dir, offset=offset)
        new_submit_files.append(newfiles_path)

    print()
    if not new_submit_files:
        print('No new input files found across all datasets.')
        return

    print('=' * 60)
    if len(new_submit_files) == 1:
        print('To submit:')
        print('  condor_submit', new_submit_files[0])
    else:
        multi = odir + args.tag + '_NewFiles_MultiSub.sh'
        if not args.dry_run:
            with open(multi, 'w') as fh:
                fh.write('#!/bin/bash\n')
                for sp in new_submit_files:
                    label = '/'.join(sp.split('/')[-4:-2])
                    fh.write('echo "Submitting {}..."\n'.format(label))
                    fh.write('condor_submit ' + sp + '\n')
                    fh.write('echo ""\n')
            os.chmod(multi, 0o755)
            print('Multi-submit script:', multi)
            print('Run with:  source', multi)
        else:
            print('[DRY RUN]', len(new_submit_files),
                  'newfiles submit file(s) would be generated')
    print('=' * 60)


# ---------------------------------------------------------------------------
# Transfer: hadd per-job EOS output, xrdcp to final skim destination
# ---------------------------------------------------------------------------

def transfer_jobs(args):
    """For each submit.sh, hadd per-job output files from EOS and xrdcp to
    /store/group/lpcsusylep/malazaro/KUCMSSkims/{version}/."""

    dest_base = EOS_SKIMS_BASE.rstrip('/') + '/' + args.version.strip('/')
    dest_xrd  = EOS_SERVER + '//' + dest_base.lstrip('/')
    version = args.version.split("_")[1]
    submit_files = []
    for root, dirs, files in os.walk(args.output):
        if 'submit.sh' in files:
            submit_files.append(os.path.join(root, 'submit.sh'))
    if not submit_files:
        print('No submit.sh files found under', args.output)
        return

    print('Found', len(submit_files), 'submit file(s)')
    print('Destination:', dest_xrd)
    print()

    for tool in ('hadd', 'xrdcp'):
        try:
            subprocess.check_call(['which', tool],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('ERROR:', tool, 'not found. Run with CMSSW set up.')
            sys.exit(1)

    if not args.dry_run:
        if not xrdfs_stat(dest_base, verbose=args.verbose):
            print('Destination folder not found — creating:', dest_base)
            if not xrdfs_mkdir(dest_base, verbose=args.verbose):
                print('ERROR: could not create destination folder:', dest_base)
                sys.exit(1)
            print('Created:', dest_base)
        else:
            print('Destination folder exists:', dest_base)
        print()

    n_ok = 0
    n_skip = 0
    n_fail = 0
    for submit_path in sorted(submit_files):
        if f"rjrskim_{version}" not in submit_path:
            continue
        header_lines, job_entries, eos_out_dir = parse_submit_file(submit_path)
        if not eos_out_dir:
            print('SKIP (no --eos-out found):', submit_path)
            n_skip += 1
            continue

        parts     = submit_path.replace('\\', '/').split('/')
        subfolder = parts[-4] if len(parts) >= 4 else 'unknown'
        tag       = parts[-3] if len(parts) >= 3 else 'rjrskim'

        coll_tag = ''
        for _, arg_line in job_entries[:1]:
            m = re.search(r'-i\s+(\S+)', arg_line)
            if m:
                coll_tag = collection_tag_from_url(m.group(1))

        # Strip 'kucmsntuple_{coll_tag}_' prefix from subfolder if present
        # (older R22 ntuple dirs were named kucmsntuple_{tag}_{dataset}_...,
        # while R23/R24 dirs are just {dataset}_...; normalize to the clean form)
        clean_subfolder = subfolder
        kucms_prefix = 'kucmsntuple_' + coll_tag + '_'
        if coll_tag and clean_subfolder.startswith(kucms_prefix):
            clean_subfolder = clean_subfolder[len(kucms_prefix):]
        elif clean_subfolder.startswith('kucmsntuple_'):
            clean_subfolder = clean_subfolder[len('kucmsntuple_'):]

        merged_name = (coll_tag + '_' if coll_tag else '') + clean_subfolder + '_' + tag + '.root'
        dest_path   = dest_base.rstrip('/') + '/' + merged_name
        dest_url    = EOS_SERVER + '//' + dest_path.lstrip('/')

        print('Sample:', subfolder)
        print('  Merged :', merged_name)
        print('  Dest   :', dest_url)

        if not args.force and xrdfs_stat(dest_path, verbose=args.verbose):
            print('  Already on EOS, skipping (use --force to overwrite)')
            n_skip += 1
            print()
            continue

        entries    = xrdfs_ls(eos_out_dir, verbose=args.verbose)
        root_files = [e for e in entries if e.endswith('.root')]

        if not root_files:
            print('  WARNING: no output .root files found -- skipping')
            n_skip += 1
            print()
            continue

        print('  Files  :', len(root_files))
        xrd_inputs = [EOS_SERVER + '//' + f.lstrip('/') for f in root_files]

        scratch = os.path.expandvars(args.scratch).rstrip('/')
        hadded  = os.path.join(scratch, merged_name)

        if args.dry_run:
            print('  [dry-run] hadd', len(xrd_inputs), 'file(s) ->', hadded)
            print('  [dry-run] xrdcp ->', dest_url)
            print()
            continue

        if not os.path.isdir(scratch):
            os.makedirs(scratch, exist_ok=True)

        use_chunks = args.chunk_size > 0 and len(root_files) > args.chunk_size

        if os.path.exists(hadded) and not args.force:
            print('  Local hadded file already exists, skipping hadd (use --force to redo)')
        elif use_chunks:
            chunks   = [root_files[i:i+args.chunk_size] for i in range(0, len(root_files), args.chunk_size)]
            partials = []
            chunk_ok = True
            print('  Chunked merge: {} files in {} chunk(s) of {}'.format(
                  len(root_files), len(chunks), args.chunk_size))
            for ci, chunk in enumerate(chunks):
                partial = os.path.join(scratch, merged_name.replace('.root', '.part{}.root'.format(ci)))
                partials.append(partial)
                print('  Chunk {}/{}: staging {} file(s)...'.format(ci+1, len(chunks), len(chunk)))
                local_inputs = []
                stage_ok = True
                for f in chunk:
                    src_url = EOS_SERVER + '//' + f.lstrip('/')
                    local_f = os.path.join(scratch, os.path.basename(f))
                    result  = subprocess.run(['xrdcp', '-s', src_url, local_f])
                    if result.returncode != 0:
                        print('  WARNING: xrdcp staging failed for', os.path.basename(f))
                        stage_ok = False
                        break
                    local_inputs.append(local_f)
                if not stage_ok:
                    for f in local_inputs:
                        if os.path.exists(f): os.remove(f)
                    chunk_ok = False
                    break
                hadd_cmd = ['hadd', '-k', '-j', '4', '-v', '0', '-d', scratch, partial] + local_inputs
                result   = subprocess.run(hadd_cmd)
                for f in local_inputs:
                    if os.path.exists(f): os.remove(f)
                if result.returncode != 0:
                    print('  WARNING: hadd failed on chunk', ci)
                    chunk_ok = False
                    break
            if not chunk_ok:
                for p in partials:
                    if os.path.exists(p): os.remove(p)
                n_fail += 1
                print()
                continue
            print('  Final merge of', len(partials), 'partial(s) ->', hadded)
            hadd_flags = ['-f'] if args.force and os.path.exists(hadded) else []
            result = subprocess.run(['hadd', '-k', '-j', '4', '-v', '0', '-d', scratch]
                                    + hadd_flags + [hadded] + partials)
            for p in partials:
                if os.path.exists(p): os.remove(p)
            if result.returncode != 0:
                print('  WARNING: final hadd failed -- skipping')
                n_fail += 1
                print()
                continue
        else:
            # Stage all files at once then hadd locally
            local_inputs = []
            print('  Staging', len(root_files), 'file(s) to scratch...')
            stage_ok = True
            for f in root_files:
                src_url = EOS_SERVER + '//' + f.lstrip('/')
                local_f = os.path.join(scratch, os.path.basename(f))
                result  = subprocess.run(['xrdcp', '-s', src_url, local_f])
                if result.returncode != 0:
                    print('  WARNING: xrdcp staging failed for', os.path.basename(f))
                    stage_ok = False
                    break
                local_inputs.append(local_f)
            if not stage_ok:
                for f in local_inputs:
                    if os.path.exists(f): os.remove(f)
                n_fail += 1
                print()
                continue
            hadd_flags = ['-f'] if args.force and os.path.exists(hadded) else []
            hadd_cmd   = ['hadd', '-k', '-j', '4', '-v', '0', '-d', scratch] + hadd_flags + [hadded] + local_inputs
            print('  hadd:', len(local_inputs), 'file(s) ->', hadded)
            result = subprocess.run(hadd_cmd)
            for f in local_inputs:
                if os.path.exists(f): os.remove(f)
            if result.returncode != 0:
                print('  WARNING: hadd failed -- skipping')
                n_fail += 1
                print()
                continue

        xrdcp_flags = ['-f'] if args.force else []
        xrdcp_cmd   = ['xrdcp', '-s'] + xrdcp_flags + [hadded, dest_url]
        print('  xrdcp ->', dest_url)
        result = subprocess.run(xrdcp_cmd)
        if result.returncode != 0:
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
        help='Generate newfiles.sh for each dataset containing only EOS files '
             'not present in the existing submit.sh (requires --eos-path)')
    parser.add_argument('--check', action='store_true',
        help='Check job completeness and generate resub files for missing output')
    parser.add_argument('--test-job', dest='test_job', type=int, nargs='?', const=-1, default=None,
        metavar='N', help='With --check: print local run command for job N (omit N for first missing)')
    parser.add_argument('--transfer', action='store_true',
        help='hadd per-job EOS output per sample and xrdcp to final skim destination')
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
    parser.add_argument('--tag', '-t', default='rjrskim',
        help='Job name tag (default: rjrskim)')
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
    parser.add_argument('--psiche',  action='store_true', help='Enable PSICHE jets (default: off)')
    parser.add_argument('--eos-out', dest='eos_out', default='/eos/uscms/store/user/$USER/LLPSkims',
        help='Write output directly to this EOS path (default: /eos/uscms/store/user/$USER/LLPSkims)')
    parser.add_argument('--dry-run', action='store_true', help='Print plan, write nothing')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print xrdfs commands')
    args = parser.parse_args()

    if args.new_inputs:
        if not args.eos_path:
            parser.error('--new-inputs requires --eos-path')
        if args.version:
            args.tag = args.tag + '_v' + args.version.lstrip('v')
        new_inputs_mode(args)
        return

    if args.check:
        check_jobs(args.output, verbose=args.verbose, test_job=args.test_job)
        return

    if args.transfer:
        if args.version is None:
            detected = None
            for root, dirs, files in os.walk(args.output):
                if 'submit.sh' in files:
                    parts = root.replace('\\', '/').split('/')
                    tag_part = parts[-2] if len(parts) >= 2 else ''
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
        transfer_jobs(args)
        return

    if not args.eos_path:
        parser.error('--eos-path is required unless --check or --transfer is specified')

    if args.version:
        args.tag = args.tag + '_v' + args.version.lstrip('v')

    eos_path    = args.eos_path.rstrip('/')
    odir        = args.output.rstrip('/') + '/'
    sample_type = detect_sample_type(eos_path)
    year        = detect_year(eos_path)
    timecali    = get_timecali(year, sample_type)
    if not timecali and sample_type == 'signal':
        timecali = 'r2_ul18_mc'

    print('EOS path:    ', eos_path)
    print('Sample type: ', sample_type)
    print('Year:        ', year)
    print('Timecali:    ', timecali)
    print()

    # detect the dominant keyword for data key generation
    data_kw = next((kw for kw in sorted(DATA_KEYWORDS, key=len, reverse=True) if kw in eos_path), 'Data')

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

    for subfolder, root_files in datasets:
        print('\nDataset:', subfolder)

        if args.dry_run:
            print('  .root files (peek):', len(root_files),
                  '(dry-run; real count may be higher)')
        else:
            print('  Files:', len(root_files))

        if not root_files and not args.dry_run:
            print('  No .root files found, skipping.')
            continue

        # --- metadata lookup ---
        xsec       = '1'
        key        = subfolder[:subfolder.index('_')] if '_' in subfolder else subfolder
        gluinomass = '0'
        n2mass     = '0'
        n1mass     = '0'
        mctype     = 0
        mc_wt      = '1'

        if sample_type == 'data':
            mctype   = 1
            key      = make_data_key(subfolder, data_kw)

        elif sample_type == 'signal':
            sig_xsec, gluinomass, n2mass, n1mass = lookup_signal_xsec(subfolder)
            if sig_xsec is not None:
                xsec = str(sig_xsec)
            else:
                print('  WARNING: no xsec found for', subfolder, '- using 0')
                xsec = '0'
            key = signal_event_count_key(subfolder)

        else:  # bg
            meta = lookup_bg_meta(subfolder)
            if meta:
                xsec = str(meta['xsec'])
                key  = meta['key']
            else:
                print('  WARNING: no xsec found for', subfolder, '- using 0')
                xsec = '0'

        # --- build skimmer flags ---
        flags = (
            ' --xsec '        + xsec
            + ' --dataSetKey ' + key
            + ' --gluinoMass ' + gluinomass
            + ' --N2Mass '     + n2mass
            + ' --timeCaliTag ' + timecali
            + ' --MCweight '   + mc_wt
            + ' --MCtype '     + str(mctype)
        )
        if mctype == 0:
            flags += ' --hasGenInfo'
        if not args.psiche:
            flags += ' --noBHC'
        if args.no_sv:
            flags += ' --noSV'
        if args.hlt_off:
            flags += ' --HLTPathsOff'

        print('  xsec:', xsec, '| key:', key)
        if mctype != 1 and event_count_keys and key not in event_count_keys:
            print('  WARNING: key not found in config/EventCount.txt:', key)
            print('           This job will likely get inf evtFillWgt.')

        # --- output paths ---
        prefix      = subfolder.split('_')[0]
        work_dir    = odir + prefix + '/' + subfolder + '/' + args.tag
        src_dir     = work_dir + '/src'
        log_dir     = work_dir + '/log'
        out_dir     = work_dir + '/out'
        submit_path = src_dir  + '/submit.sh'
        ofilename   = 'condor_' + subfolder + '_' + args.tag

        print('  Submit:', submit_path)

        # --- record weights ---
        weights[key] = {
            'xsec': xsec, 'mc_wt': mc_wt, 'mctype': mctype,
            'timecali': timecali, 'gluinomass': gluinomass,
            'n2mass': n2mass, 'n1mass': n1mass,
            'n_files': len(root_files), 'subfolder': subfolder,
        }

        submit_files.append(submit_path)

        if args.dry_run:
            print('  [DRY RUN]')
            continue

        eos_out_dir = None
        if args.eos_out:
            eos_out_dir = os.path.expandvars(args.eos_out).rstrip('/') + '/' + prefix + '/' + subfolder + '/' + args.tag
            print('  EOS out:', eos_out_dir)
            xrdfs_mkdir(eos_out_dir, verbose=args.verbose)

        local_dirs = [src_dir, log_dir] if eos_out_dir else [src_dir, log_dir, out_dir]
        for d in local_dirs:
            os.makedirs(d, exist_ok=True)

        write_submit(submit_path, log_dir, out_dir, ofilename,
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
            with open(multi, 'w') as fh:
                fh.write('#!/bin/bash\n')
                for sp in submit_files:
                    label = '/'.join(sp.split('/')[-4:-2])
                    fh.write('echo "Submitting ' + label + '..."\n')
                    fh.write('condor_submit ' + sp + '\n')
                    fh.write('echo ""\n')
            os.chmod(multi, 0o755)
            print('Multi-submit script:', multi)
            print('Run with:  source', multi)
        else:
            print('[DRY RUN]', len(submit_files), 'submit files would be generated')
    print('=' * 60)


if __name__ == '__main__':
    main()
