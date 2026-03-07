import os, sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# ---------------------------------------------------------------------------
# Inherit the full base configuration from llpgana_miniaod.py.
#
# We exec() it with a stripped sys.argv so its internal parseArguments()
# uses all defaults and doesn't conflict with our own argument parsing.
# Any future edits to llpgana_miniaod.py propagate here automatically.
# ---------------------------------------------------------------------------
_base_cfg = os.path.join(
    os.environ.get('CMSSW_BASE', ''),
    'src/KUCMSNtupleizer/KUCMSNtupleizer/test/llpgana_miniaod.py'
)
if not os.path.isfile(_base_cfg):
    raise RuntimeError(
        'Base config not found: {}\n'
        'Make sure CMSSW_BASE is set and cmsenv has been run.'.format(_base_cfg)
    )

_saved_argv = sys.argv[:]
sys.argv    = [sys.argv[0]]   # suppress base config's VarParsing.parseArguments()
exec(open(_base_cfg).read())  # populates `process` in this scope
sys.argv    = _saved_argv

# ---------------------------------------------------------------------------
# Parse our own arguments.
# VarParsing('python') pre-registers: inputFiles (list), outputFile (str),
# secondaryInputFiles (list), maxEvents (int, default -1), section, tag.
# ---------------------------------------------------------------------------
_opts = VarParsing('python')
_opts.register('eventFilter', '',
               VarParsing.multiplicity.singleton,
               VarParsing.varType.string,
               'event filter / skim selection (overrides config default when set)')
_opts.register('outputFileName', 'output.root',
               VarParsing.multiplicity.singleton,
               VarParsing.varType.string,
               'output file name (matches llpgana_miniaod.py convention)')
_opts.parseArguments()

# Override input source
if _opts.inputFiles:
    process.source.fileNames = cms.untracked.vstring(_opts.inputFiles)

# Override output file
process.TFileService.fileName = cms.string(_opts.outputFileName)

# Override maxEvents (-1 = all events, the right default for full jobs)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(_opts.maxEvents)
)

# Override event filter if explicitly supplied
if _opts.eventFilter:
    process.tree.fltrSelection = cms.string(_opts.eventFilter)
