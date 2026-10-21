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
import fnmatch
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile
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

# Run 3 (13.6 TeV) background xsec [pb], keyed by the primary dataset name up
# to its '_Tune...' suffix. Values are XSDB records (LO, 13.6 TeV), from the
# XSDB dump jking79/KUCMSNtupleizer macros/info_XSDB_2025-04-04_22-16.json and,
# for QCD HT-40to70 / HT-70to100 (absent from the dump), the XSDB values quoted
# in cms-flaf/FLAF config/crossSections13p6TeV.yaml (commit d01db63). GJets are
# the exception (see below). 'key' is the
# EventCount.txt key stem; the era ('_R23', ...) is appended at lookup.
BG_DB_RUN3 = {
    # --- QCD ---
    'QCD-4Jets_HT-40to70':      {'xsec': 311400000, 'key': 'QCD_40t70'},
    'QCD-4Jets_HT-70to100':     {'xsec': 58500000,  'key': 'QCD_70t100'},
    'QCD-4Jets_HT-100to200':    {'xsec': 25400000,  'key': 'QCD_100t200'},
    'QCD-4Jets_HT-200to400':    {'xsec': 1961000,   'key': 'QCD_200t400'},
    'QCD-4Jets_HT-400to600':    {'xsec': 95620,     'key': 'QCD_400t600'},
    'QCD-4Jets_HT-600to800':    {'xsec': 13540,     'key': 'QCD_600t800'},
    'QCD-4Jets_HT-800to1000':   {'xsec': 3033,      'key': 'QCD_800t1000'},
    'QCD-4Jets_HT-1000to1200':  {'xsec': 883.7,     'key': 'QCD_1000t1200'},
    'QCD-4Jets_HT-1200to1500':  {'xsec': 383.5,     'key': 'QCD_1200t1500'},
    'QCD-4Jets_HT-1500to2000':  {'xsec': 125.2,     'key': 'QCD_1500t2000'},
    'QCD-4Jets_HT-2000':        {'xsec': 26.49,     'key': 'QCD_2000up'},
    # --- GJets (Sherpa PTG-120 sample: no value yet) ---
    # From toicca/dijet_rdf src/jec4prompt/skim.py (commit 170dad9), listed there
    # for Summer24 with no XSDB reference; not cross-checked against XSDB.
    'GJ-4Jets_HT-40to70':    {'xsec': 15240, 'key': 'GJets_40t70'},
    'GJ-4Jets_HT-70to100':   {'xsec': 8111,  'key': 'GJets_70t100'},
    'GJ-4Jets_HT-100to200':  {'xsec': 7327,  'key': 'GJets_100t200'},
    'GJ-4Jets_HT-200to400':  {'xsec': 1541,  'key': 'GJets_200t400'},
    'GJ-4Jets_HT-400to600':  {'xsec': 167.6, 'key': 'GJets_400t600'},
    'GJ-4Jets_HT-600':       {'xsec': 54.39, 'key': 'GJets_600up'},
    # --- WJets ---
    'WtoLNu-4Jets_MLNu-0to120_HT-40to100':    {'xsec': 4254,   'key': 'WJets_MLNu0t120_40t100'},
    'WtoLNu-4Jets_MLNu-0to120_HT-100to400':   {'xsec': 1626,   'key': 'WJets_MLNu0t120_100t400'},
    'WtoLNu-4Jets_MLNu-0to120_HT-400to800':   {'xsec': 59.99,  'key': 'WJets_MLNu0t120_400t800'},
    'WtoLNu-4Jets_MLNu-0to120_HT-800to1500':  {'xsec': 6.23,   'key': 'WJets_MLNu0t120_800t1500'},
    'WtoLNu-4Jets_MLNu-0to120_HT-1500to2500': {'xsec': 0.4477, 'key': 'WJets_MLNu0t120_1500t2500'},
    'WtoLNu-4Jets_MLNu-0to120_HT-2500':       {'xsec': 0.03075, 'key': 'WJets_MLNu0t120_2500up'},
    'WtoLNu-4Jets_MLNu-120_HT-40to100':     {'xsec': 20.56,     'key': 'WJets_MLNu120up_40t100'},
    'WtoLNu-4Jets_MLNu-120_HT-100to400':    {'xsec': 10.19,     'key': 'WJets_MLNu120up_100t400'},
    'WtoLNu-4Jets_MLNu-120_HT-400to800':    {'xsec': 0.5239,    'key': 'WJets_MLNu120up_400t800'},
    'WtoLNu-4Jets_MLNu-120_HT-800to1500':   {'xsec': 0.06255,   'key': 'WJets_MLNu120up_800t1500'},
    'WtoLNu-4Jets_MLNu-120_HT-1500to2500':  {'xsec': 0.005066,  'key': 'WJets_MLNu120up_1500t2500'},
    'WtoLNu-4Jets_MLNu-120_HT-2500':        {'xsec': 0.0003788, 'key': 'WJets_MLNu120up_2500up'},
    # --- ZJets ---
    'Zto2Nu-4Jets_HT-100to200':   {'xsec': 273.7,    'key': 'ZJets_100t200'},
    'Zto2Nu-4Jets_HT-200to400':   {'xsec': 75.96,    'key': 'ZJets_200t400'},
    'Zto2Nu-4Jets_HT-400to800':   {'xsec': 13.19,    'key': 'ZJets_400t800'},
    'Zto2Nu-4Jets_HT-800to1500':  {'xsec': 1.364,    'key': 'ZJets_800t1500'},
    'Zto2Nu-4Jets_HT-1500to2500': {'xsec': 0.09865,  'key': 'ZJets_1500t2500'},
    'Zto2Nu-4Jets_HT-2500':       {'xsec': 0.006699, 'key': 'ZJets_2500up'},
}
RUN3_YEARS = {'R22', 'R23', 'R24', 'R25'}

# Signal xsec [pb] by (process_prefix, gluino_mass_string)
# process_prefix: 'gogo' matches gogoG, gogoZ, gogoGZ (and T5Zg); 'sqsq' matches sqsqG
# Run 3 (13.6 TeV). Signal samples without an explicit year in their name are Run 3.
SIGNAL_XSEC = {
    ('gogo', '1500'): 0.02209,
    ('gogo', '1800'): 4.524e-03,
    ('gogo', '2000'): 0.001684,
    ('gogo', '2300'): 0.0004130,
    ('gogo', '2500'): 0.0001687,
    ('sqsq', '1700'): 0.001248,
    ('sqsq', '1850'): 0.0005679,
    ('sqsq', '2000'): 0.0002654,
    ('sqsq', '2150'): 0.0001267,
}

# Run 2 (13 TeV) signal xsec [pb], used only for samples explicitly tagged
# with a Run 2 year (e.g. '_2018').
SIGNAL_XSEC_RUN2 = {
    ('gogo', '1800'): 3.349e-03,
    ('gogo', '2300'): 2.871e-04,
}
RUN2_YEARS = {'R16', 'R17', 'R18'}

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
SIGNAL_KEYWORDS = {'SMS', 'gogoG', 'gogoZ', 'gogoGZ', 'sqsqG', 'GlGl', 'T5Zg'}
BRANCH_MASKS = {
    'llpCombineSkim': 'config/branch_masks/llpcombine_analysis_exact.txt',
    'svSkim': 'config/branch_masks/sv_analysis_core.txt',
    'svFullFastValidation': 'config/branch_masks/sv_fullfast_validation.txt',
    'llpStandardPlots': 'config/branch_masks/llpstandardplots_min.txt',
}

SANDBOX_DEFAULT  = "/uscms/home/mlazarov/nobackup/sandboxes/sandbox-CMSSW_13_3_3.tar.bz2"
EOS_SKIMS_BASE   = "/store/group/lpcsusylep/malazaro/KUCMSSkims"
EOS_SERVER       = "root://cmseos.fnal.gov"
EOS_OUT_DEFAULT  = "/eos/uscms/store/user/$USER/LLPSkims"
SKIMMER_DIR      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONDOR_DIR       = os.path.dirname(os.path.abspath(__file__))
CONFIG_TAR       = os.path.join(CONDOR_DIR, 'config.tgz')  # `make configtar` output
EVENT_COUNT_REL  = 'config/EventCount.txt'
DEFAULT_MC_TIMECALI   = 'r2_ul18_mc'
DEFAULT_SKIM_VERSION  = 'skims_v49'  # --transfer fallback when no version is given/detected
BG_DB_YEAR       = 'R18'  # every BG_DB xsec/key is UL18

QUEUE_ONE_COLUMN = 'queue Arguments from ('
QUEUE_TWO_COLUMN = 'queue GlobalJobIndex, Arguments from ('

_CRAB_TS_RE = re.compile(r'^(?:\d{6}|\d{8})_\d{6}$')

# Signal mass-point tokens. Parse names only after canonical_signal_name()
# so 'mGo-' spellings are already 'mGl-'.
_MGL_RE  = re.compile(r'mGl-(\d+)')
_MN2_RE  = re.compile(r'mN2-(\d+)')
_MN1_RE  = re.compile(r'mN1-(\d+)')
_CTAU_RE = re.compile(r'(?:(?:^|[_-])ct|\dctau-)(-?[0-9]+p[0-9]+|-?[0-9]+)')


class EosError(RuntimeError):
    """An xrdfs command failed for a reason other than the path not existing."""

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
    # Always this checkout's copy: it is the one `make configtar` (run in
    # SKIMMER_DIR) tars for the workers and the one the staleness check compares.
    return os.path.join(SKIMMER_DIR, EVENT_COUNT_REL)


def read_event_count_keys(path=None):
    if path is None:
        path = default_event_count_path()
    if not os.path.isfile(path):
        print('WARNING: could not read', path,
              '- MC EventCount key checks will be skipped')
        return set()
    return set(read_event_count_file(path))


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


def branch_mask_mapping_text():
    if not BRANCH_MASKS:
        return 'Available branch mask aliases: [none configured]'
    lines = ['Available branch mask aliases:']
    for alias in sorted(BRANCH_MASKS):
        lines.append('  {} -> {}'.format(alias, BRANCH_MASKS[alias]))
    return '\n'.join(lines)


def resolve_branch_mask(mask):
    """Resolve a branch mask alias or path.

    Returns (worker_side_path, alias). alias is None when the user supplied an
    explicit path instead of a known mask name.
    """
    alias_path = BRANCH_MASKS.get(mask)
    if alias_path:
        return normalize_branch_mask_path(alias_path), mask
    try:
        return normalize_branch_mask_path(mask), None
    except ValueError as exc:
        raise ValueError(str(exc) + '\n' + branch_mask_mapping_text())


def lookup_bg_meta(subfolder, year=''):
    """Find the best matching BG_DB entry for this subfolder name.
    Returns the dict or None. BG_DB only holds UL18 xsecs/keys, so a match on
    a sample tagged with any other year is refused rather than silently
    normalized against the UL18 EventCount entry."""
    # sort by key length descending so longer (more specific) patterns win
    for pattern in sorted(BG_DB, key=len, reverse=True):
        if pattern in subfolder:
            if year and year != BG_DB_YEAR:
                raise ValueError('BG_DB only has {} entries but this {} sample matched '
                                 '{!r}; add {} xsecs/keys before submitting: {}'.format(
                                     BG_DB_YEAR, year, pattern, year, subfolder))
            return BG_DB[pattern]
    return None


def lookup_bg_meta_run3(names, year=''):
    """Run 3 counterpart of lookup_bg_meta. `names` are candidate sample names
    (primary dataset directory first: CRAB task names truncate the HT bin, e.g.
    'HT-1500to25'). Returns {'xsec', 'key'} with the era appended to the key
    (event counts differ per era even though the 13.6 TeV xsec does not), or
    None when no name is a BG_DB_RUN3 sample."""
    for name in names:
        for pattern, meta in BG_DB_RUN3.items():
            if name == pattern or name.startswith(pattern + '_'):
                if year not in RUN3_YEARS:
                    raise ValueError('{!r} is a Run 3 sample but its collection year is '
                                     '{!r}; expected one of {}'.format(
                                         name, year, sorted(RUN3_YEARS)))
                return {'xsec': meta['xsec'], 'key': meta['key'] + '_' + year}
    return None


def signal_mass_point(name):
    """(mGl, mN2, mN1, ctau) strings parsed from a canonical signal name; any
    token that is absent comes back as None."""
    matches = (regex.search(name) for regex in (_MGL_RE, _MN2_RE, _MN1_RE, _CTAU_RE))
    return tuple(m.group(1) if m else None for m in matches)


def lookup_signal_xsec(subfolder, year=''):
    """Extract process prefix and gluino mass from subfolder, look up xsec.
    Samples explicitly tagged with a Run 2 year use the 13 TeV table."""
    mgl, mn2, mn1, _ = signal_mass_point(subfolder)
    if mgl is None:
        return None, '0', '0', '0'
    proc = signal_process_key(subfolder)[:4]  # 'gogo' or 'sqsq'
    table = SIGNAL_XSEC_RUN2 if year in RUN2_YEARS else SIGNAL_XSEC
    xsec = table.get((proc, mgl), None)
    return xsec, mgl, mn2 or '0', mn1 or '0'


def signal_process_key(name):
    if 'SqSq' in name or 'sqsq' in name:
        return 'sqsqG'
    if 'T5Zg' in name:
        return 'gogoGZ'
    if 'GlGl-GZ' in name or 'GlGlGZ' in name or 'gogoGZ' in name:
        return 'gogoGZ'
    if 'GlGl-Z' in name or 'GlGlZ' in name or 'gogoZ' in name:
        return 'gogoZ'
    if 'GlGl-G' in name or 'GlGlG' in name or 'gogoG' in name:
        return 'gogoG'
    # Some naming schemes carry the decay mode as its own '_'-delimited token
    # instead of hyphenated onto the base (e.g. "SMS-GlGl_mGl-.._mN1-.._GZ_ct..").
    # Check GZ before Z/G so a "_GZ_" token isn't mistaken for a bare Z or G one.
    if re.search(r'(?:^|_)GZ(?:_|$)', name):
        return 'gogoGZ'
    if re.search(r'(?:^|_)Z(?:_|$)', name):
        return 'gogoZ'
    if re.search(r'(?:^|_)G(?:_|$)', name):
        return 'gogoG'
    # A gluino sample with no decay-mode token is GZ (same convention as
    # normalize_signal_process_tag); anything else is unrecognized.
    if 'GlGl' in name or 'gogo' in name:
        return 'gogoGZ'
    raise ValueError('could not classify signal process (GlGl/SqSq/gogo*/sqsq*/T5Zg) '
                     'from: ' + name)


def signal_event_count_tier(name):
    # Order matters: check Fast* patterns before the Mini* check, since
    # some FastSim tags (e.g. 'FASTMINI') contain 'MINI' as a substring
    # and would otherwise be misclassified as FullSim.
    if re.search(r'fast', name, re.IGNORECASE):
        return 'FASTSIM'
    if re.search(r'(?:FULLMINI|MINIAOD|MINI|MiniAOD|Mini)', name):
        return 'FULLMINI'
    if 'AODSIM' in name:
        return 'AODSIM'
    # No silent default: FullSim and FastSim samples are separate productions
    # with separate EventCount.txt entries (see EventCount.txt tier tags
    # FULLMINI/FASTSIM/AODSIM). Guessing here previously caused FastSim
    # samples to be silently normalized against an unrelated AODSIM campaign.
    raise ValueError(
        "could not classify signal EventCount tier (matched none of "
        "Fast*/Mini*/AODSIM) from: " + name
    )


def normalize_signal_ctau(raw_ctau, tier):
    value = raw_ctau.lstrip('-')
    if tier in ('AODSIM', 'FASTSIM') and value.startswith('0p') and len(value) == 3:
        return value[-1]
    if tier == 'FULLMINI' and value in {'1', '5'}:
        return '0p' + value
    return value


def signal_tier_from_collection(collection_name):
    """EventCount tier from the kucmsntuple_* collection name (e.g.
    'kucmsntuple_SMS_GZ_SVHPM100_Fast1_v37'). This is the single source of
    truth: the collection is one production, whereas task suffixes carry CRAB
    noise (a stray 'AODSIM' token) that must never decide FastSim vs FullSim.
    Returns None if the collection carries no tier tag."""
    if not collection_name:
        return None
    try:
        return signal_event_count_tier(collection_name)
    except ValueError:
        return None


_SIGNAL_YEAR_RE = re.compile(r'(?:^|_)20(16|17|18|22|23|24|25)(?=_|$)')


def signal_year(*names):
    """Era code ('R22', 'R18', ...) from an explicit '_20YY' token in the signal
    dataset names, or '' if none carries one. Only an explicit year counts:
    pre-existing signal samples are Run 3 with no year in their names and must
    keep their unsuffixed EventCount keys. The token must be '_'-delimited so
    masses like 'mGo-2022' are never mistaken for a year."""
    years = {'R' + yy for name in names if name for yy in _SIGNAL_YEAR_RE.findall(name)}
    if len(years) > 1:
        raise ValueError('conflicting signal years %s in %r' % (sorted(years), names))
    return years.pop() if years else ''


def canonical_signal_name(name):
    """Normalize gluino-mass spellings to the 'mGl-' form every parser expects
    (e.g. T5Zg samples use 'mGo-1800')."""
    return re.sub(r'mGo-(?=\d)', 'mGl-', name)


def signal_event_count_key(name, tier=None, year=''):
    """EventCount key for a signal sample. A non-empty year (e.g. 'R22') is
    appended as '_R22'; samples without an explicit year keep the unsuffixed key."""
    suffix = '_' + year if year else ''
    mgl, mn2, mn1, ctau = signal_mass_point(name)

    if not (mgl and mn2 and mn1 and ctau):
        match = re.search(
            r'((?:gogoGZ|gogoG|gogoZ|sqsqG)_'
            r'(?:AODSIM|FASTSIM|FASTSIMAOD|FASTMINI|FASTAOD|FULLMINI|MINIAOD|MINI|AOD)_.*)',
            name,
        )
        if match:
            key = re.sub(r'_(MINIAOD|MINI)_', '_FULLMINI_', match.group(1))
            key = re.sub(r'_(FASTSIMAOD|FASTMINI|FASTAOD)_', '_FASTSIM_', key)
            if tier is not None and '_' + tier + '_' not in key:
                raise ValueError('EventCount key %s contradicts tier %s' % (key, tier))
            # An already-formatted key may carry its year suffix; strip it so the
            # end-anchored ctau normalization below still applies.
            key_year = re.search(r'_(R\d\d)$', key)
            if key_year:
                if year and key_year.group(1) != year:
                    raise ValueError('EventCount key %s contradicts year %s' % (key, year))
                suffix = '_' + key_year.group(1)
                key = key[:key_year.start()]
            if '_AODSIM_' in key or '_FASTSIM_' in key:
                key = re.sub(r'_ct0p([15])$', r'_ct\1', key)
            elif '_FULLMINI_' in key:
                key = re.sub(r'_ct([15])$', r'_ct0p\1', key)
            return key + suffix
        raise ValueError('could not derive signal EventCount key from: ' + name)

    process = signal_process_key(name)
    if tier is None:
        tier = signal_event_count_tier(name)
    ct_key = normalize_signal_ctau(ctau, tier)
    return f'{process}_{tier}_mGl-{mgl}_mN2-{mn2}_mN1-{mn1}_ct{ct_key}{suffix}'


_SIGNAL_TIER_TAG_RE = re.compile(
    r'(?:FASTSIM|FASTSIMAOD|FASTMINI|FASTAOD|Fast1|FastSim'
    r'|FULLMINI|MINIAOD|MINI|MiniAOD|Mini|AODSIM)'
)


def signal_metadata_name(subfolder, eos_path):
    """Choose the best available name for signal xsec/EventCount parsing.

    Some naming schemes (e.g. ".../kucmsntuple_SMS_GZ_SVHPM100_Fast1_v37/...")
    only encode the production tier (Fast1/AODSIM/MINI/...) in the collection
    name, not in the per-mass-point task suffix -- even though the task suffix
    alone already has mGl/mN2/mN1/ctau and would otherwise be picked first.
    Prefer a candidate that carries both the mass point *and* a tier tag;
    fall back to the first mass-point-complete candidate if none do (tier
    classification will then fail loudly in signal_event_count_tier()
    instead of silently guessing).
    """
    candidates = [
        subfolder,
        os.path.basename(eos_path.rstrip('/')),
        os.path.basename(os.path.dirname(eos_path.rstrip('/'))),
    ]

    def has_mass_point(name):
        return all(signal_mass_point(canonical_signal_name(name)))

    complete = [name for name in candidates if has_mass_point(name)]
    if not complete:
        return subfolder
    for name in complete:
        if _SIGNAL_TIER_TAG_RE.search(name):
            return name
    return complete[0]


def make_data_key(subfolder, sample_type_kw):
    """Build a short dataset key for data from the subfolder name."""
    # e.g. 'MET_AOD_Run2018A-15Feb...' -> 'MET18A'
    year_m = re.search(r'Run(\d{4})', subfolder) or re.search(r'(\d{4})', subfolder)
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
        timecali = DEFAULT_MC_TIMECALI
    return sample_type, year, timecali, data_keyword(eos_path)


def dataset_timecali(metadata, sample_type, timecali):
    """Timecali for one dataset: a year found in the dataset's own name
    (signal only, see signal_year) overrides the collection-level one."""
    if metadata.get('year'):
        return get_timecali(metadata['year'], sample_type) or timecali
    return timecali


def print_sample_context(eos_path, sample_type, year, timecali):
    print('EOS path:    ', eos_path)
    print('Sample type: ', sample_type)
    print('Year:        ', year)
    print('Timecali:    ', timecali)
    print()


def resolve_dataset_metadata(subfolder, eos_path, sample_type, data_kw, collection_name=None,
                             primary_dataset=None):
    xsec       = '1'
    key        = subfolder[:subfolder.index('_')] if '_' in subfolder else subfolder
    gluinomass = '0'
    n2mass     = '0'
    n1mass     = '0'
    mctype     = 0
    mc_wt      = '1'
    year       = ''

    if sample_type == 'data':
        mctype = 1
        key    = make_data_key(subfolder, data_kw)
    elif sample_type == 'signal':
        metadata_name = canonical_signal_name(signal_metadata_name(subfolder, eos_path))
        year = signal_year(subfolder, eos_path, collection_name)
        sig_xsec, gluinomass, n2mass, n1mass = lookup_signal_xsec(metadata_name, year)
        if sig_xsec is None:
            # A zero xsec silently zeroes evtFillWgt for the whole skim.
            table = 'SIGNAL_XSEC_RUN2' if year in RUN2_YEARS else 'SIGNAL_XSEC'
            raise ValueError('no signal xsec for %r (year %r): add its gluino/squark '
                             'mass to %s' % (metadata_name, year or 'Run 3', table))
        xsec = str(sig_xsec)
        tier = signal_tier_from_collection(collection_name)
        if tier is None:
            # No tier in the collection name: fall back to the task names, but
            # refuse to guess if they disagree.
            tiers = set()
            for candidate in (subfolder, eos_path):
                try:
                    tiers.add(signal_event_count_tier(candidate))
                except ValueError:
                    pass
            if len(tiers) > 1:
                raise ValueError('conflicting signal tiers %s in %r / %r and no tier '
                                 'in collection %r' % (sorted(tiers), subfolder,
                                                       eos_path, collection_name))
        key = signal_event_count_key(metadata_name, tier=tier, year=year)
        # dataSetKey's tier is the source of truth for FastSim vs FullSim: it is
        # what EventCount.txt was actually keyed and counted under.
        if '_FASTSIM_' in key:
            mctype = 2
    else:
        bg_year = detect_year(collection_name or eos_path)
        meta = (lookup_bg_meta_run3([n for n in (primary_dataset, subfolder) if n], bg_year)
                or lookup_bg_meta(subfolder, bg_year))
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
        'year': year,
    }


def dataset_metadata(dataset, sample_type, data_kw):
    """The one place datasets are resolved to metadata, so every caller gets
    the collection name (the authoritative source of the signal tier)."""
    dataset_data_kw = data_keyword(dataset.raw_task_suffix) if sample_type == 'data' else data_kw
    return resolve_dataset_metadata(dataset.raw_task_suffix, dataset.task_name,
                                    sample_type, dataset_data_kw,
                                    collection_name=dataset.collection_name,
                                    primary_dataset=os.path.basename(
                                        os.path.dirname(dataset.task_dir)))


def resolve_dxy_scale(metadata, dxy_scale):
    """The --dxy-scale value actually applied to this dataset, mirroring the
    gating in build_skimmer_flags exactly (single source of truth so the
    skimmer command line and the on-disk ntuple_tag naming never disagree).
    Returns (applied, warning): applied is 'off' unless this dataset is
    FastSim with usable mN2/mN1; warning is a printable message or None."""
    if dxy_scale == 'off':
        return 'off', None
    if metadata['mctype'] != 2:
        return 'off', ('  WARNING: --dxy-scale requested but ' + metadata['key']
                        + ' is not a FastSim sample - skimmer dxySig scaling not applied')
    if metadata['n2mass'] == '0' and metadata['n1mass'] == '0':
        return 'off', ('  WARNING: --dxy-scale requested but mN2/mN1 unavailable for '
                        + metadata['key'] + ' - skimmer dxySig scaling not applied')
    return dxy_scale, None


def insert_sxy_scale_tag(ntuple_tag, dxy_scale):
    """Insert the SxyScale{Nominal,Up,Down} variant marker into ntuple_tag
    immediately after its FastSim reco-tier token (Fast1/FASTSIM/FastSim/...),
    so nominal/up/down dxySig-scale submissions of the same FastSim ntuple
    never collide on disk and the variant is visible right next to the tier
    it modifies (Sxy = the dxy significance this scales). No-op for 'off'.
    Only ever called with a non-'off' dxy_scale for FastSim datasets (see
    resolve_dxy_scale), so a Fast* token should always be present; falls back
    to appending at the end if one somehow isn't."""
    if dxy_scale == 'off':
        return ntuple_tag
    suffix = 'SxyScale' + dxy_scale.capitalize()
    tokens = ntuple_tag.split('_')
    for i, token in enumerate(tokens):
        if re.search(r'fast', token, re.IGNORECASE):
            return '_'.join(tokens[:i + 1] + [suffix] + tokens[i + 1:])
    return ntuple_tag + '_' + suffix


def remove_sxy_scale_tag(ntuple_tag):
    """Inverse of insert_sxy_scale_tag. Job output names (-o) are built from
    dataset.canonical_base, which never carries the SxyScale marker, so
    anything that re-derives the job prefix from an on-disk ntuple_tag must
    strip it first. The variant stays in the directory and merged-file names."""
    return '_'.join(token for token in ntuple_tag.split('_')
                    if not re.match(r'^SxyScale(Nominal|Up|Down)$', token))


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
    if metadata['mctype'] != 1:
        flags += ' --hasGenInfo'
    if not args.psiche:
        flags += ' --noBHC'
    if args.no_sv:
        flags += ' --noSV'
    if args.hlt_off:
        flags += ' --HLTPathsOff'
    if args.branch_mask:
        flags += ' --branchMask ' + args.branch_mask
    applied_dxy_scale, warning = resolve_dxy_scale(metadata, args.dxy_scale)
    if warning:
        print(warning)
    elif applied_dxy_scale != 'off':
        delta_m = float(metadata['n2mass']) - float(metadata['n1mass'])
        flags += ' --dxySigScale ' + applied_dxy_scale + ' --svDxyDeltaM ' + str(delta_m)
    return flags, applied_dxy_scale


def dataset_paths(output_dir, dataset, tag, dxy_scale='off'):
    odir       = output_dir.rstrip('/') + '/'
    ntuple_tag = insert_sxy_scale_tag(dataset.ntuple_tag, dxy_scale)
    work_dir   = odir + ntuple_tag + '/' + dataset.sample_tag + '/' + tag
    return DatasetPaths(
        sample_tag=dataset.sample_tag,
        ntuple_tag=ntuple_tag,
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


# ---------------------------------------------------------------------------
# EventCount.txt: reading/writing the full table, and counting straight from
# the EOS ntuples themselves (no ntuple_master_lists dependency). ROOT is only
# imported lazily, on first use, so plain job-generation/--check/--transfer
# runs never need it.
# ---------------------------------------------------------------------------

_ROOT = None
_ROOT_CONFIGURED = False


def _setup_root():
    global _ROOT, _ROOT_CONFIGURED
    if _ROOT is None:
        import ROOT as ROOT_MODULE
        _ROOT = ROOT_MODULE
    if _ROOT_CONFIGURED:
        return
    _ROOT.gEnv.SetValue('TFile.AsyncPrefetching', 1)
    _ROOT.gEnv.SetValue('TFile.MaxCacheSize', 100000000)
    _ROOT.gEnv.SetValue('TFile.ReadBufferSize', 1048576)
    _ROOT_CONFIGURED = True


def sum_event_weights(root_files):
    """Sum nTotEvts/sumEvtWgt out of the configtree of root_files (bare EOS
    paths, as stored on DatasetInfo.root_files) directly from EOS."""
    _setup_root()
    urls = [eos_xrootd_url(path) for path in root_files]
    tree_name = 'tree/configtree'
    if len(urls) == 1:
        df = _ROOT.RDataFrame(tree_name, urls[0])
    else:
        chain = _ROOT.TChain(tree_name)
        for url in urls:
            chain.Add(url)
        df = _ROOT.RDataFrame(chain)
    ntot = df.Sum('nTotEvts').GetValue()
    sumw = df.Sum('sumEvtWgt').GetValue()
    return ntot, sumw


def read_event_count_file(path=None):
    """Full key -> (nTotEvts_str, sumEvtWgt_str) table from EventCount.txt."""
    if path is None:
        path = default_event_count_path()
    data = {}
    try:
        with open(path) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 3:
                    data[parts[0]] = (parts[1], parts[2])
    except OSError:
        pass
    return data


def write_event_count_file(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for key in sorted(data):
            ntot, sumw = data[key]
            handle.write('{} {} {}\n'.format(key, ntot, sumw))


def missing_event_count_entries(datasets, sample_type, data_kw, event_count_keys):
    """(dataset, key) pairs for datasets whose EventCount.txt key doesn't exist
    yet. One entry per distinct missing key (HT slices etc. can share a key)."""
    missing = []
    seen_keys = set()
    for dataset in datasets:
        if not dataset.root_files:
            continue
        metadata = dataset_metadata(dataset, sample_type, data_kw)
        key = metadata['key']
        if metadata['mctype'] == 1 or key in event_count_keys or key in seen_keys:
            continue
        seen_keys.add(key)
        missing.append((dataset, key))
    return missing


def fill_missing_event_counts(datasets, sample_type, data_kw, event_count_keys, args):
    """Find EventCount.txt keys this submission needs but doesn't have, and
    optionally count them directly from the EOS ntuples' configtree instead of
    requiring a hand-maintained master list. Returns the updated key set."""
    missing = missing_event_count_entries(datasets, sample_type, data_kw, event_count_keys)
    if not missing:
        return event_count_keys

    print()
    print('EventCount.txt is missing', len(missing), 'key(s) needed by this submission:')
    for dataset, key in missing:
        print('  -', key, '(' + str(len(dataset.root_files)) + ' files, ' +
              dataset.sample_tag + ')')

    if args.auto_count_events:
        do_count = True
    else:
        ans = input('\nCount these now from the EOS ntuples and add them to '
                     'config/EventCount.txt? [y/N] ').strip().lower()
        do_count = ans in ('y', 'yes')

    if not do_count:
        print('Skipping. These jobs will FAIL at runtime (missing dataSetKey) until '
              'config/EventCount.txt is updated.')
        return event_count_keys

    event_count_path = default_event_count_path()
    counted = read_event_count_file(event_count_path)
    for dataset, key in missing:
        print('Counting', key, '(' + str(len(dataset.root_files)), 'files)...')
        ntot, sumw = sum_event_weights(dataset.root_files)
        print('  nTotEvts =', ntot, ' sumEvtWgt =', sumw)
        counted[key] = (str(ntot), str(sumw))
    write_event_count_file(event_count_path, counted)
    print('Wrote', len(missing), 'new key(s) to', event_count_path)
    # config.tgz is now stale; ensure_configtar_current() handles the rebuild.

    return event_count_keys | {key for _, key in missing}


# ---------------------------------------------------------------------------
# config.tgz: the worker only sees what was tarred by `make configtar`
# ---------------------------------------------------------------------------

def offer_configtar(args):
    """Run `make configtar` if --auto-configtar or the user agrees. Returns
    True if the tarball was rebuilt."""
    if args.auto_configtar:
        do_tar = True
    else:
        ans = input('Run `make configtar` now from ' + SKIMMER_DIR + '? [y/N] ').strip().lower()
        do_tar = ans in ('y', 'yes')
    if not do_tar:
        return False
    print('Running `make configtar` in', SKIMMER_DIR, '...')
    subprocess.run(['make', 'configtar'], cwd=SKIMMER_DIR, check=True)
    print('config.tgz rebuilt.')
    return True


def configtar_member_problem(rel_path, tar_path=None):
    """Why the worker would not see SKIMMER_DIR/rel_path as it is locally, or
    None if config.tgz contains an identical copy. rel_path is the worker-side
    path, e.g. config/branch_masks/sv_analysis_core.txt."""
    if tar_path is None:
        tar_path = CONFIG_TAR
    if not os.path.isfile(tar_path):
        return tar_path + ' does not exist'
    with open(os.path.join(SKIMMER_DIR, rel_path), 'rb') as handle:
        local_bytes = handle.read()
    try:
        with tarfile.open(tar_path, 'r:gz') as tar:
            member = None
            for candidate in tar.getmembers():
                if os.path.normpath(candidate.name) == os.path.normpath(rel_path):
                    member = candidate
            if member is None or not member.isfile():
                return rel_path + ' is not in ' + tar_path
            tar_bytes = tar.extractfile(member).read()
    except (tarfile.TarError, OSError) as exc:
        return 'could not read ' + tar_path + ': ' + str(exc)
    if tar_bytes != local_bytes:
        return rel_path + ' in ' + tar_path + ' differs from the local file'
    return None


def configtar_required_files(args, sample_type):
    """(worker-side path, what goes wrong if config.tgz has a stale copy)."""
    required = []
    if args.branch_mask:
        required.append((args.branch_mask,
                         'jobs would run WITHOUT this mask and write full (unmasked) skims'))
    if sample_type != 'data' and os.path.isfile(os.path.join(SKIMMER_DIR, EVENT_COUNT_REL)):
        required.append((EVENT_COUNT_REL,
                         'jobs would fail or be mis-normalized (missing/old dataSetKey entries)'))
    return required


def configtar_problems(required):
    problems = []
    for rel_path, consequence in required:
        problem = configtar_member_problem(rel_path)
        if problem is not None:
            problems.append((problem, consequence))
    return problems


def ensure_configtar_current(args, sample_type):
    """Stop before writing submit files if the worker's config.tgz is missing,
    or has an old copy of, a file this submission depends on (branch mask,
    EventCount.txt). Offers a single `make configtar` covering all of them."""
    required = configtar_required_files(args, sample_type)
    problems = configtar_problems(required)
    if not problems:
        return

    print()
    print('WARNING: config.tgz is stale for this submission:')
    for problem, consequence in problems:
        print('  ' + problem)
        print('    -> ' + consequence)
    if args.dry_run:
        print('  [DRY RUN] Run `make configtar` from', SKIMMER_DIR, 'before submitting.')
        return

    if offer_configtar(args):
        problems = configtar_problems(required)
        if not problems:
            return
        print('ERROR: config.tgz still does not match after rebuilding:')
        for problem, _ in problems:
            print('  ' + problem)
    sys.exit('Aborting: run `make configtar` from ' + SKIMMER_DIR + ' and resubmit.')


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


def _xrdfs_target(path):
    """(server, bare /store/... path) for any supported EOS path form."""
    m = re.match(r'(root://[^/]+)/', path)
    return (m.group(1) if m else EOS_SERVER), normalize_eos_path(path)


def _xrdfs(subcommand, path, extra=(), verbose=False):
    server, bare = _xrdfs_target(path)
    cmd = ['xrdfs', server, subcommand] + list(extra) + [bare]
    if verbose:
        print(' ', ' '.join(cmd))
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)


def xrdfs_mkdir(path, verbose=False):
    """Create a directory on EOS via xrdfs, including any missing parents."""
    if _xrdfs('mkdir', path, extra=['-p'], verbose=verbose).returncode == 0:
        return True
    print('  WARNING: failed to create EOS directory:', normalize_eos_path(path))
    return False


def xrdfs_ls(path, verbose=False):
    """List entries under path using xrdfs. Returns list of bare paths.
    A path that does not exist lists as empty; any other failure (EOS/xrootd
    trouble) raises EosError instead of being mistaken for an empty directory."""
    result = _xrdfs('ls', path, verbose=verbose)
    if result.returncode == 0:
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not xrdfs_stat(path, verbose=verbose):
        return []
    raise EosError('xrdfs ls failed for {}: {}'.format(
        normalize_eos_path(path), result.stderr.strip()))


def xrdfs_stat(path, verbose=False):
    """Return True if path exists on EOS."""
    return _xrdfs('stat', path, verbose=verbose).returncode == 0


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


def normalize_signal_process_tag(sample):
    """Collapse a 'GlGl'/'SqSq' base plus its decay-mode token (bare 'G',
    'Z', or 'GZ' -- whether hyphenated onto the base or left as its own
    '_'-delimited token further down the name, e.g. from mass-point tokens
    sitting in between) into the single canonical process tag ('gogoG',
    'gogoZ', 'gogoGZ', 'sqsqG') that downstream skim readers key off of.
    Leaves non-signal sample names untouched.
    """
    hyphen_combos = (
        ('GlGl-GZ', 'gogoGZ'), ('GlGl-Z', 'gogoZ'), ('GlGl-G', 'gogoG'),
        ('SqSq-G', 'sqsqG'),
    )
    for combo, canon in hyphen_combos:
        if sample == combo or sample.startswith(combo + '_'):
            return canon + sample[len(combo):]

    tokens = sample.split('_')
    base_map = {'GlGl': 'gogo', 'SqSq': 'sqsq'}
    proc = base_map.get(tokens[0])
    if proc is None:
        return sample

    valid_decays = {'GZ', 'Z', 'G'} if proc == 'gogo' else {'G'}
    for i, tok in enumerate(tokens[1:], start=1):
        if tok in valid_decays:
            remaining = tokens[1:i] + tokens[i + 1:]
            return '_'.join([proc + tok] + remaining)

    default_decay = 'GZ' if proc == 'gogo' else 'G'
    return '_'.join([proc + default_decay] + tokens[1:])


_N2CTAU_TOKEN_RE = re.compile(r'(?<![^_])N2ctau-(\d+p\d+)(?![^_])')


def normalize_signal_ctau_tag(sample):
    """Rewrite the raw CRAB-task lifetime token 'N2ctau-0p5' (meters) to the
    'ct0p5' spelling that LLPCombine's BFTool::NormalizeCtauToken parses
    (ct0p1 -> process key _10, ct0p5 -> _50).  Other names are untouched."""
    return _N2CTAU_TOKEN_RE.sub(r'ct\1', sample)


def make_sample_tag(raw_task_suffix):
    tier_tokens = {'MiniAOD', 'MINIAOD', 'AODSIM', 'MINI', 'MIN'}
    tokens = [token for token in raw_task_suffix.split('_') if token not in tier_tokens]
    sample = '_'.join(tokens)
    if sample.startswith('SMS-'):
        sample = sample[len('SMS-'):]
    elif sample.startswith('SMS_'):
        sample = sample[len('SMS_'):]
    return normalize_signal_ctau_tag(normalize_signal_process_tag(sample))


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


def skim_version_number(version):
    """'50', 'v50', 'skims_v50' (optionally with a trailing '/') -> '50'."""
    version = version.strip().rstrip('/')
    return remove_prefix_once(version, 'skims_v').lstrip('v')


def normalize_skim_version(version):
    if not version:
        return None
    return 'skims_v' + skim_version_number(version)


def tag_version(tag):
    """Version number from a job tag's _vN suffix, or None."""
    match = re.search(r'_v(\d+)$', tag)
    return match.group(1) if match else None


def resolve_job_tag(args):
    """Submission job tag: --tag (default: branch mask alias or rjrskim),
    plus _v<N> when --version is given."""
    tag = args.tag if args.tag is not None else (args.branch_mask_alias or 'rjrskim')
    if args.version:
        tag += '_v' + skim_version_number(args.version)
    return tag


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
        elif _is_dir_path(entry):
            stack.extend(xrdfs_ls(entry, verbose=verbose))
    return results


def find_crab_datasets(path, max_files=-1, verbose=False):
    """Discover CRAB datasets with names derived from collection and task dirs."""
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


def write_submit(path, log_dir, out_dir, ofilename, root_files, flags, args, eos_out_dir=None):
    outname = ofilename + '.$(Process).root'

    if eos_out_dir:
        eos_arg = ' --eos-out ' + normalize_eos_path(eos_out_dir)
        output_lines = [
            'transfer_output_files = ',
            'when_to_transfer_output = NEVER',
        ]
    else:
        eos_arg = ''
        remap_target = os.path.join(os.getcwd(), out_dir) + '/' + outname
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

    lines += ['', '', QUEUE_ONE_COLUMN]
    for i, fpath in enumerate(root_files):
        lines.append('###### job' + str(i) + ' ######')
        lines.append('-i ' + eos_xrootd_url(fpath) + ' ' + flags
                     + ' -o ' + ofilename + '.$(Process)' + eos_arg)
    lines.append(')')

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# Job checking / resubmission helpers
# ---------------------------------------------------------------------------

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
            if stripped == QUEUE_ONE_COLUMN:
                in_queue = True
                two_column_queue = False
            elif stripped == QUEUE_TWO_COLUMN:
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


def rewrite_submit_header(header_lines, job_log_stem, cluster_log, two_column,
                          local_output_file=None):
    """Adapt an original submit.sh header for a follow-up submit file (a
    newfiles batch or a resub). Per-job logs become <job_log_stem>.out/.err and
    the cluster log <cluster_log>, in the original log directory, so follow-ups
    never overwrite the original jobs' logs. local_output_file, when given,
    redirects local Condor output transfer to that name in the original remap
    directory (it must use $(GlobalJobIndex), so two_column must be True)."""
    local_output_dir = (_local_output_dir_from_header(header_lines)
                        if local_output_file else None)
    header = []
    for line in header_lines:
        key = line.split('=', 1)[0].strip() if '=' in line else None
        stripped = line.strip()
        if key in ('output', 'error'):
            log_dir = os.path.dirname(_header_value(line))
            ext = '.out' if key == 'output' else '.err'
            header.append(key + ' = ' + os.path.join(log_dir, job_log_stem + ext))
        elif key == 'log':
            log_dir = os.path.dirname(_header_value(line))
            header.append('log = ' + os.path.join(log_dir, cluster_log))
        elif local_output_file and key == 'transfer_output_files':
            header.append('transfer_output_files = ' + local_output_file)
        elif local_output_file and key == 'transfer_output_remaps':
            remap_target = os.path.join(local_output_dir, local_output_file)
            header.append('transfer_output_remaps = "'
                          + local_output_file + '=' + remap_target + '"')
        elif stripped in (QUEUE_ONE_COLUMN, QUEUE_TWO_COLUMN):
            header.append(QUEUE_TWO_COLUMN if two_column else QUEUE_ONE_COLUMN)
        else:
            header.append(line)
    return header


def write_indexed_submit(path, header, indexed_args, two_column):
    """Write header plus one queue row per (global_index, arg_line)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for line in header:
            handle.write(line + '\n')
        for job_index, arg_line in indexed_args:
            handle.write('###### job' + str(job_index) + ' ######\n')
            if two_column:
                handle.write(str(job_index) + ', ' + arg_line + '\n')
            else:
                handle.write(arg_line + '\n')
        handle.write(')\n')


def local_output_file_for(expected_prefix):
    return expected_prefix + '.$(GlobalJobIndex).root'


def write_inherited_new_inputs_submit(path, history, expected_prefix, new_files,
                                      first_index):
    last_index = first_index + len(new_files) - 1
    local_output = history.eos_out_dir is None
    log_stem = 'job.{}_{}'.format(first_index, last_index)
    header = rewrite_submit_header(
        history.header_lines, log_stem + '.$(Process)', log_stem + '.log',
        two_column=local_output,
        local_output_file=local_output_file_for(expected_prefix) if local_output else None,
    )
    indexed_args = []
    for offset, input_file in enumerate(new_files):
        job_index = first_index + offset
        indexed_args.append((job_index, build_arg_line_from_template(
            history.original_arg_line,
            history.submit_files[0],
            eos_xrootd_url(input_file),
            expected_prefix + '.' + str(job_index),
        )))
    write_indexed_submit(path, header, indexed_args, two_column=local_output)


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
    job_ntuple_tag = remove_sxy_scale_tag(ntuple_tag)
    return {
        'ntuple_tag': ntuple_tag,
        'sample_tag': sample_tag,
        'canonical_base': sample_tag + '__' + ntuple_tag if ntuple_tag else sample_tag,
        'job_base': sample_tag + '__' + job_ntuple_tag if job_ntuple_tag else sample_tag,
        'tag': parts[-3] if len(parts) >= 3 else 'rjrskim',
    }


def expected_submit_prefix(path_info):
    return 'condor_' + path_info['job_base'] + '__' + path_info['tag']


def write_multi_submit_script(path, submit_files, action='Submitting', blank_lines=True):
    # Submit files use paths relative to the directory they were generated from
    # (executable, config.tgz, log dirs), so condor_submit must run there no
    # matter where this script is sourced. Each line is a one-line subshell so
    # the caller's cwd is untouched and the script parses when sourced from
    # tcsh (LPC login shell) as well as bash: tcsh rejects multi-line ( ... ).
    cd_prefix = '( cd ' + shlex.quote(os.getcwd()) + ' && condor_submit '
    with open(path, 'w') as fh:
        fh.write('#!/bin/bash\n')
        for sp in submit_files:
            label = '/'.join(sp.replace('\\', '/').split('/')[-4:-2])
            fh.write('echo "{} {}..."\n'.format(action, label))
            fh.write(cd_prefix + sp + ' )\n')
            if blank_lines:
                fh.write('echo ""\n')
    os.chmod(path, 0o755)


def write_resub_file(submit_path, header_lines, missing_entries):
    """Write a resubmission condor file for missing_entries.
    $(Process) is replaced with the actual job index in each argument line
    so output filenames stay consistent with the original submission. Rows
    carry their GlobalJobIndex so each resubmitted job logs to
    resub.<global index>.out/.err instead of clobbering job.<N>.* logs.
    Returns the resub file path.
    """
    resub_path = os.path.join(os.path.dirname(submit_path), 'resub.sh')
    local_output = not any('--eos-out' in submit_arg_tokens(arg_line, submit_path)
                           for _, arg_line in missing_entries)
    expected_prefix = expected_submit_prefix(submit_path_parts(submit_path))
    header = rewrite_submit_header(
        header_lines, 'resub.$(GlobalJobIndex)', 'resub.log', two_column=True,
        local_output_file=local_output_file_for(expected_prefix) if local_output else None,
    )
    indexed_args = [
        (job_idx, replace_output_arg(arg_line.replace('$(Process)', str(job_idx)),
                                     submit_path, expected_prefix + '.' + str(job_idx)))
        for job_idx, arg_line in missing_entries
    ]
    write_indexed_submit(resub_path, header, indexed_args, two_column=True)
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
    if status == 'OK':
        return blue_text(status)
    return status


def print_job_completion_summary(history, classified, expected_indices):
    n_total = len(expected_indices)
    n_missing = len(classified.missing_indices)
    status = 'OK' if n_missing == 0 else 'INCOMPLETE'

    print('  Original jobs:', history.original_count)
    print('  Added jobs:   ', n_total - history.original_count)
    print('  Expected jobs:', n_total)
    print('  Missing jobs: ', '{}/{} [{}]'.format(
        n_missing, n_total, check_status_text(status)))
    if classified.unexpected_indices:
        print('  unexpected jobs:', classified.unexpected_indices)
    if classified.duplicate_indices:
        print('  duplicate jobs:', classified.duplicate_indices)


def transfer_version_for_check(path_info, version=None):
    return normalize_skim_version(version or tag_version(path_info['tag']))


def print_transfer_check_summary(path_info, version=None, verbose=False):
    skim_version = transfer_version_for_check(path_info, version=version)
    if not skim_version:
        print('  Transfer:      UNKNOWN (no --version and tag has no _vN suffix)')
        return

    merged_name = merged_output_name(path_info)
    dest_path = EOS_SKIMS_BASE.rstrip('/') + '/' + skim_version + '/' + merged_name
    if xrdfs_stat(dest_path, verbose=verbose):
        status = blue_text('TRANSFERRED')
    else:
        status = red_text('NOT TRANSFERRED')
    print('  Transfer:      {} ({})'.format(status, eos_xrootd_url(dest_path)))


def branch_mask_display(mask_path):
    text = mask_path if mask_path else '[None]'
    return blue_text(text) if mask_path else red_text(text)


def submit_label_from_parts(path_info):
    return path_info['canonical_base'] + ' [' + path_info['tag'] + ']'


def label_matches_filter(label, name_filter):
    if not name_filter:
        return True
    if any(char in name_filter for char in '*?[]'):
        return fnmatch.fnmatchcase(label, name_filter)
    return name_filter in label


def filter_submit_files(submit_files, name_filter=None, tag_filter=None):
    matches = []
    name_skipped = []
    tag_skipped = []
    for submit_path in sorted(submit_files):
        path_info = submit_path_parts(submit_path)
        label = submit_label_from_parts(path_info)
        if not label_matches_filter(label, name_filter):
            name_skipped.append((submit_path, label))
            continue
        if tag_filter and path_info['tag'] != tag_filter:
            tag_skipped.append((submit_path, path_info['tag']))
            continue
        matches.append((submit_path, path_info, label))
    return matches, name_skipped, tag_skipped


def print_filter_result(name_filter, tag_filter, matched):
    """Print the name-filter line and match count shared by --check and
    --transfer. Returns False (after saying so) when a filter matched nothing."""
    if name_filter:
        print('Name filter:', name_filter)
    if (name_filter or tag_filter) and not matched:
        if name_filter and not tag_filter:
            filter_text = name_filter
        elif tag_filter and not name_filter:
            filter_text = 'tag ' + tag_filter
        else:
            filter_text = 'name ' + name_filter + ', tag ' + tag_filter
        print(red_text('No submit files matched filter: ' + filter_text))
        return False
    if name_filter or tag_filter:
        print('Matched', len(matched), 'submit file(s)')
    print()
    return True


def submission_dirs(output_dir):
    """Directories condor_submit may have run from. In --eos-out mode a job
    whose worker-side xrdcp failed leaves its output in the scratch dir, and
    the empty transfer_output_files then brings it back to the submit
    directory instead of EOS. Candidates: the cwd, this script's directory,
    and every `cd <dir> && condor_submit` target in output_dir's *.sh scripts."""
    dirs = [os.getcwd(), os.path.dirname(os.path.abspath(__file__))]
    cd_pattern = re.compile(r'\(\s*cd\s+(.+?)\s+&&\s+condor_submit\b')
    for root, _subdirs, files in os.walk(output_dir):
        for name in files:
            if not name.endswith('.sh') or name == 'submit.sh':
                continue
            try:
                with open(os.path.join(root, name)) as fh:
                    for line in fh:
                        m = cd_pattern.search(line)
                        if m:
                            dirs.append(shlex.split(m.group(1))[0])
            except (OSError, ValueError):
                continue
    seen = set()
    unique = []
    for d in dirs:
        d = os.path.realpath(d)
        if d not in seen and os.path.isdir(d):
            seen.add(d)
            unique.append(d)
    return unique


def find_stranded_outputs(search_dirs, expected_prefix, missing_indices):
    """Return {index: local path} for missing jobs whose non-empty output
    landed in one of search_dirs. An index found in more than one directory
    is skipped as ambiguous."""
    pattern = re.compile(r'^' + re.escape(expected_prefix) + r'\.(\d+)\.root$')
    missing = set(missing_indices)
    found = {}
    ambiguous = set()
    for d in search_dirs:
        for name in os.listdir(d):
            m = pattern.match(name)
            if not m or int(m.group(1)) not in missing:
                continue
            path = os.path.join(d, name)
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                continue
            idx = int(m.group(1))
            if idx in found:
                ambiguous.add(idx)
            found[idx] = path
    for idx in sorted(ambiguous):
        print('  WARNING: job {} output found in several submit dirs; '
              'not moving it'.format(idx))
        del found[idx]
    return found


def recover_stranded_outputs(stranded, eos_out_dir, verbose=False, dry_run=False):
    """xrdcp each stranded local output to eos_out_dir and delete the local
    copy once the EOS file is confirmed. Returns the indices now on EOS."""
    recovered = []
    for idx in sorted(stranded):
        local = stranded[idx]
        dest = normalize_eos_path(eos_out_dir) + '/' + os.path.basename(local)
        if dry_run:
            print('  [dry-run] would move', local, '->', eos_xrootd_url(dest))
            continue
        cmd = ['xrdcp', '-f', local, eos_xrootd_url(dest)]
        if verbose:
            print(' ', ' '.join(cmd))
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        if result.returncode != 0 or not xrdfs_stat(dest, verbose=verbose):
            print('  WARNING: failed to move {} to EOS: {}'.format(
                local, result.stderr.strip()))
            continue
        os.remove(local)
        print('  Moved stranded output:', os.path.basename(local), '->', eos_xrootd_url(dest))
        recovered.append(idx)
    return recovered


def check_jobs(output_dir, verbose=False, test_job=None, name_filter=None,
               tag_filter=None, version=None, dry_run=False):
    """Scan output_dir for submit.sh files, report completeness, and write
    resubmission files for any subfolders with missing output."""

    submit_files = find_submit_files(output_dir)
    if not submit_files:
        print('No submit.sh files found under', output_dir)
        return

    filtered_submit_files, _name_skipped, _tag_skipped = filter_submit_files(
        submit_files, name_filter=name_filter, tag_filter=tag_filter)

    print('Found', len(submit_files), 'submit file(s)')
    if tag_filter:
        print('Tag filter:', tag_filter)
    if not print_filter_result(name_filter, tag_filter, filtered_submit_files):
        return

    # SxyScale nominal/up/down submissions share one job output prefix (see
    # remove_sxy_scale_tag), so a stranded file with such a prefix cannot be
    # attributed to a variant. Count over all submit files, not just the
    # filtered ones, so a filtered --check still sees the sibling variants.
    prefix_counts = {}
    for sp in submit_files:
        prefix = expected_submit_prefix(submit_path_parts(sp))
        prefix_counts[prefix] = prefix_counts.get(prefix, 0) + 1

    resub_files        = []
    first_missing_path = None
    first_missing_idx  = None
    first_missing_jobs = None
    search_dirs        = None

    for submit_path, path_info, label in filtered_submit_files:
        history = recover_submission_history_or_exit(submit_path)
        expected_indices = history_expected_indices(history)

        # List output files from EOS or local
        if history.eos_out_dir:
            output_files = [os.path.basename(e)
                            for e in xrdfs_ls(history.eos_out_dir, verbose=verbose)]
        else:
            try:
                out_dir = _local_output_dir_from_header(history.header_lines)
            except ValueError:
                out_dir = None
            output_files = os.listdir(out_dir) if out_dir and os.path.isdir(out_dir) else []

        expected_prefix = expected_submit_prefix(path_info)
        classified = classify_transfer_outputs(output_files, expected_prefix,
                                               expected_indices)

        print('{}:'.format(label))
        # Outputs that never reached EOS (worker xrdcp failed) come back to
        # the submit dir; move them over and count them as done. Duplicated
        # indices are left alone: they need a human look, not another copy.
        if history.eos_out_dir and classified.missing_indices:
            recoverable = sorted(set(classified.missing_indices)
                                 - set(classified.duplicate_indices))
            if search_dirs is None:
                search_dirs = submission_dirs(output_dir)
            stranded = find_stranded_outputs(search_dirs, expected_prefix, recoverable)
            if stranded and prefix_counts.get(expected_prefix, 0) > 1:
                print('  WARNING: {} output(s) found in the submit dir, but {} submissions '
                      '(SxyScale variants) share the prefix {} -- not moving; move them '
                      'by hand to the right variant'.format(
                          len(stranded), prefix_counts[expected_prefix], expected_prefix))
                stranded = {}
            if stranded:
                print('  Found {} output(s) in the submit dir instead of EOS'.format(
                    len(stranded)))
                recovered = recover_stranded_outputs(
                    stranded, history.eos_out_dir, verbose=verbose, dry_run=dry_run)
                output_files += [expected_prefix + '.' + str(idx) + '.root'
                                 for idx in recovered]
                classified = classify_transfer_outputs(output_files, expected_prefix,
                                                       expected_indices)
        missing_indices = classified.missing_indices
        print_job_completion_summary(history, classified, expected_indices)
        if not missing_indices and not classified.unexpected_indices:
            print_transfer_check_summary(path_info, version=version, verbose=verbose)

        if missing_indices:
            missing_set = set(missing_indices)
            missing_entries = [
                (job.index, job.arg_line)
                for job in history.jobs
                if job.index in missing_set
            ]
            resub_path = write_resub_file(submit_path, history.header_lines, missing_entries)
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

def discover_datasets(args):
    """Classify --eos-path and discover its CRAB datasets (shared by normal
    submission and --new-inputs). Exits if nothing is found.
    Returns (sample_type, timecali, data_kw, datasets)."""
    eos_path = args.eos_path.rstrip('/')
    sample_type, year, timecali, data_kw = sample_context(eos_path)
    print_sample_context(eos_path, sample_type, year, timecali)

    print('Scanning EOS tree...')
    datasets = find_crab_datasets(eos_path, max_files=args.max_files, verbose=args.verbose)
    if not datasets:
        print('ERROR: nothing found under', eos_path)
        print('Check the path and that xrdfs is available on this system.')
        sys.exit(1)
    print('Datasets found:', len(datasets))
    return sample_type, timecali, data_kw, datasets


def new_inputs_mode(args):
    """Scan EOS for files not present in existing submit.sh files and generate
    immutable newfiles_<first>_<last>.sh submit files for each dataset that has
    new inputs. Job indices are assigned after the full submission history."""

    odir = args.output.rstrip('/') + '/'
    sample_type, _timecali, data_kw, datasets = discover_datasets(args)
    new_submit_files = []

    for dataset in datasets:
        root_files = dataset.root_files

        print('\nDataset:', dataset.sample_tag)

        metadata = dataset_metadata(dataset, sample_type, data_kw)
        applied_dxy_scale, warning = resolve_dxy_scale(metadata, args.dxy_scale)
        if warning:
            print(warning)
        paths = dataset_paths(args.output, dataset, args.tag, applied_dxy_scale)

        if not os.path.exists(paths.submit_path):
            print('  No existing submit.sh found — run without --new-inputs first, skipping.')
            continue

        history = recover_submission_history_or_exit(paths.submit_path)

        new_files = new_input_files(root_files, history.jobs)

        cli_eos_out_dir = eos_output_dir(args, paths)
        if history.eos_out_dir:
            if args.eos_out_supplied and cli_eos_out_dir:
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
        print('  Already submitted:', len(history.jobs))
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
    unexpected_indices = sorted(found_indices - expected_indices)
    output_for_index = {
        index: paths[0]
        for index, paths in paths_by_index.items()
        if index in expected_indices and len(paths) == 1
    }
    # A duplicated index has no usable output either, so it counts as missing.
    missing_indices = sorted(expected_indices - set(output_for_index))
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


def hadd_root_files(output_path, input_files, scratch):
    # No -k: a corrupt or unreadable per-job file must fail the merge rather
    # than be silently dropped (and then deleted by --clean).
    result = subprocess.run(['hadd', '-f', '-j', '4', '-v', '0', '-d', scratch,
                             output_path] + input_files)
    return result.returncode == 0


def merge_root_files(root_files, hadded, scratch, merged_name, chunk_size):
    """Always merges from scratch: a leftover hadded/partial file from an
    interrupted run may be truncated, so it is never reused."""
    use_chunks = chunk_size > 0 and len(root_files) > chunk_size

    if os.path.exists(hadded):
        print('  Removing leftover local hadded file from a previous run')
        os.remove(hadded)

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
        ok = hadd_root_files(hadded, partials, scratch)
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
    ok = hadd_root_files(hadded, local_inputs, scratch)
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

    filtered_submit_files, _name_skipped, _tag_skipped = filter_submit_files(
        submit_files, name_filter=name_filter, tag_filter=args.tag)

    print('Found', len(submit_files), 'submit file(s)')
    print('Destination:', dest_xrd)
    print('Tag filter:', args.tag if args.tag else 'all')
    if not print_filter_result(name_filter, args.tag, filtered_submit_files):
        return

    require_transfer_tools()

    if not args.dry_run:
        ensure_transfer_destination(dest_base, verbose=args.verbose)

    n_ok = 0
    n_skip = 0
    n_fail = 0
    for submit_path, path_info, label in filtered_submit_files:
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
        root_files = classified.selected_outputs

        print_job_completion_summary(history, classified, expected_indices)
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
                                args.chunk_size):
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
        help='With --transfer: overwrite the merged skim on EOS if it already exists')
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
        help='Branch mask alias or file passed to skimmer, e.g. llpCombineSkim '
             'or config/branch_masks/llpcombine_analysis_exact.txt')
    parser.add_argument('--psiche',  action='store_true', help='Enable PSICHE jets (default: off)')
    parser.add_argument('--dxy-scale', dest='dxy_scale', default='off',
        choices=['off', 'nominal', 'up', 'down'],
        help='FastSim SV dxySig ("Sxy") scale correction (hadronic tailfrac calibration, '
             'propagated in deltaM=mN2-mN1) passed to the skimmer as --dxySigScale/'
             '--svDxyDeltaM. Only applied to signal samples with a FastSim dataSetKey; '
             'ignored (with a warning) for FullSim/data samples. When applied, the '
             'variant (SxyScaleNominal/Up/Down) is inserted into the ntuple_tag right '
             'after the Fast* token, so nominal/up/down submissions never collide on '
             'disk; the skim tag (--branch-mask/--tag) is left untouched. Default: off.')
    parser.add_argument('--eos-out', dest='eos_out', default=None,
        help='Write output directly to this EOS path (default: ' + EOS_OUT_DEFAULT + ')')
    parser.add_argument('--dry-run', action='store_true', help='Print plan, write nothing')
    parser.add_argument('--verbose', '-v', action='store_true', help='Print xrdfs commands')
    parser.add_argument('--auto-count-events', dest='auto_count_events', action='store_true',
        help='If EventCount.txt is missing a needed key, count it from the EOS ntuples '
             'and add it without prompting')
    parser.add_argument('--auto-configtar', dest='auto_configtar', action='store_true',
        help='When config.tgz is stale (EventCount.txt or the --branch-mask file '
             'missing/outdated in it), run `make configtar` without prompting')
    args = parser.parse_args()
    args.eos_out_supplied = args.eos_out is not None
    if not args.eos_out_supplied:
        args.eos_out = EOS_OUT_DEFAULT

    if args.branch_mask:
        try:
            original_branch_mask = args.branch_mask
            args.branch_mask, args.branch_mask_alias = resolve_branch_mask(args.branch_mask)
            if args.branch_mask != original_branch_mask:
                print('Branch mask:', original_branch_mask, '->', args.branch_mask)
        except ValueError as exc:
            parser.error(str(exc))
    else:
        args.branch_mask_alias = None

    if args.new_inputs:
        if not args.eos_path:
            parser.error('--new-inputs requires --eos-path')
        args.tag = resolve_job_tag(args)
        new_inputs_mode(args)
        return

    if args.check is not None:
        check_jobs(args.output, verbose=args.verbose, test_job=args.test_job,
                   name_filter=args.check, tag_filter=args.tag,
                   version=args.version, dry_run=args.dry_run)
        return

    if args.transfer is not None:
        if args.version is None:
            # Only look at the submit files this transfer will actually touch.
            matches, _, _ = filter_submit_files(find_submit_files(args.output),
                                                name_filter=args.transfer,
                                                tag_filter=args.tag)
            detected = sorted({tag_version(path_info['tag']) for _, path_info, _ in matches}
                              - {None}, key=int)
            if len(detected) > 1:
                parser.error('matched submit files span versions {}; pass --version '
                             'or narrow with --tag/filter'.format(', '.join(detected)))
            if detected:
                version_str = 'skims_v' + detected[0]
                ans = input('Auto-detected version {}. Proceed? [Y/n] '.format(version_str)).strip().lower()
            else:
                version_str = DEFAULT_SKIM_VERSION
                ans = input('No --version specified. Proceed with {} ? [Y/n] '.format(version_str)).strip().lower()
            if ans not in ('', 'y', 'yes'):
                print('Aborted.')
                return
            args.version = version_str
        else:
            args.version = normalize_skim_version(args.version)
        transfer_jobs(args, name_filter=args.transfer)
        return

    if not args.eos_path:
        parser.error('--eos-path is required unless --check or --transfer is specified')

    args.tag = resolve_job_tag(args)
    odir = args.output.rstrip('/') + '/'
    sample_type, timecali, data_kw, datasets = discover_datasets(args)
    event_count_keys = read_event_count_keys() if sample_type != 'data' else set()

    if sample_type != 'data' and not args.dry_run:
        event_count_keys = fill_missing_event_counts(datasets, sample_type, data_kw,
                                                      event_count_keys, args)

    # After the EventCount step so a single `make configtar` covers both it
    # and the branch mask.
    ensure_configtar_current(args, sample_type)

    # -----------------------------------------------------------------------
    # Process each dataset
    # -----------------------------------------------------------------------
    weights      = {}
    submit_files = []

    for dataset in datasets:
        root_files = dataset.root_files

        print('\nDataset:', dataset.sample_tag)

        if not root_files and not args.dry_run:
            print('  No .root files found, skipping.')
            continue

        metadata = dataset_metadata(dataset, sample_type, data_kw)
        ds_timecali = dataset_timecali(metadata, sample_type, timecali)
        flags, applied_dxy_scale = build_skimmer_flags(metadata, args, ds_timecali)

        paths = dataset_paths(args.output, dataset, args.tag, applied_dxy_scale)
        eos_out_dir = eos_output_dir(args, paths)
        final_skim = final_merged_filename(dataset.sample_tag, paths.ntuple_tag,
                                           args.tag)

        print('  Ntuple tag:     ', paths.ntuple_tag)
        print('  Skim tag:       ', args.tag)
        print('  Branch mask:    ', branch_mask_display(args.branch_mask))
        if applied_dxy_scale != 'off':
            print('  dxySig scale:   ', applied_dxy_scale)
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
        if metadata.get('year'):
            print('  Year:           ', metadata['year'] + ' (timecali ' + ds_timecali + ')')
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
        weights[metadata['key']] = weight_record(metadata, ds_timecali, root_files,
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
        # Merge into the existing file: one --output dir usually collects
        # several submissions (QCD, GJets, ...), and each run only knows its own.
        try:
            with open(odir + 'weights.json') as fh:
                weights = dict(json.load(fh), **weights)
        except (OSError, ValueError):
            pass
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
    try:
        main()
    except EosError as exc:
        sys.exit('ERROR: ' + str(exc))
