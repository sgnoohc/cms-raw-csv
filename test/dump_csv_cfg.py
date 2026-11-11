# Unpack CMS RAW data (RawToDigi) and write events as CSV files, one directory per event with one file per subsystem.
#   cmsRun dump_csv_cfg.py inputFiles=file:RAW.root                            # the first event in the file
#   cmsRun dump_csv_cfg.py inputFiles=file:RAW.root maxEvents=5                # the first 5 events
#   cmsRun dump_csv_cfg.py inputFiles=file:RAW.root events=383631:994851190    # chosen events (run:event, comma-separated)
# More options: outputDir (default event_csv), globalTag (auto:run3_data_prompt), era (Run3_2024), threads (4).
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('events', '', VarParsing.multiplicity.list, VarParsing.varType.string,
                 'events to dump as run:event, comma-separated (default: the first maxEvents events)')
options.register('outputDir', 'event_csv', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'directory that gets one run<run>_evt<event>/ directory per event')
options.register('globalTag', 'auto:run3_data_prompt', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'conditions: cabling maps and geometry')
options.register('era', 'Run3_2024', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 'data-taking era, as in Configuration/Eras/python/Era_<era>_cff.py')
options.register('threads', 4, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'number of threads')
options.setDefault('inputFiles', 'file:RAW.root')
options.setDefault('maxEvents', 1)
options.parseArguments()

era = getattr(importlib.import_module(f'Configuration.Eras.Era_{options.era}_cff'), options.era)
process = cms.Process('CSV', era)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(
    [f if ':' in f else 'file:' + f for f in options.inputFiles]))  # a bare path is a local file
if options.events:
    process.source.eventsToProcess = cms.untracked.VEventRange(options.events)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1 if options.events else options.maxEvents))
process.options.numberOfThreads = options.threads

process.csv = cms.EDAnalyzer('EventCsvDumper', dir=cms.string(options.outputDir))
process.p = cms.Path(process.RawToDigi * process.csv)
