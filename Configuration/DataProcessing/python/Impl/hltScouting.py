#!/usr/bin/env python3
"""
_hltScouting_

Scenario supporting proton collisions with input HLT scouting data

"""

from __future__ import print_function

import os
import sys

from Configuration.DataProcessing.Scenario import *
from Configuration.DataProcessing.Utils import stepALCAPRODUCER,stepSKIMPRODUCER,addMonitoring,dictIO,dqmIOSource,harvestingMode,dqmSeq,nanoFlavours,gtNameAndConnect
import FWCore.ParameterSet.Config as cms
from Configuration.DataProcessing.RecoTLR import customisePrompt,customiseExpress

class hltScouting(Scenario):
    def __init__(self):
        Scenario.__init__(self)
        self.recoSeq=''
        self.cbSc='pp'
        self.isRepacked=False
        self.promptCustoms= [ 'Configuration/DataProcessing/RecoTLR.customisePrompt' ]
        #self.expressCustoms=[ ]
        #self.alcaHarvCustoms=[]
        #self.visCustoms=[ ]
        self.promptModifiers = cms.ModifierChain()
        #self.expressModifiers = modifyExpress
        #self.visModifiers = modifyExpress
    """
    _hltScouting_

    Implement configuration building for data processing for proton
    collision data taking with input HLT scouting data for Era_Run3_2024
    """

    def _checkRepackedFlag(self, options, **args):
        if 'repacked' in args:
            if args['repacked'] == True:
                options.isRepacked = True
            else:
                options.isRepacked = False
        else:
            # not sure what we should do in this case?
            options.isRepacked = False

    def promptReco(self, globalTag, **args):
        """
        _promptReco_

        Proton collision data taking prompt reco with input HLT scouting data

        """
        if not 'skims' in args:
            args['skims']=['@allForPrompt']

        if not 'customs' in args:
            args['customs']= [ ]

        for c in self.promptCustoms:
            args['customs'].append(c)

        #step = stepALCAPRODUCER(args['skims'])
        PhysicsSkimStep = ''
        if ("PhysicsSkims" in args) :
            PhysicsSkimStep = stepSKIMPRODUCER(args['PhysicsSkims'])
        #dqmStep = dqmSeq(args,'')
        options = Options()
        options.__dict__.update(defaultOptions.__dict__)
        options.scenario = self.cbSc
        if ('nThreads' in args) :
            options.nThreads=args['nThreads']

        miniAODStep = ''
        nanoAODStep = ''
        
        if 'outputs' in args:
            print(args['outputs'])
            for a in args['outputs']:
                if a['dataTier'] == 'MINIAOD':
                    miniAODStep = ',PAT'
                if a['dataTier'] in ['NANOAOD', 'NANOEDMAOD']:
                    if "nanoFlavours" in args:
                        nanoAODStep = ',NANO'+nanoFlavours(args['nanoFlavours'])
                    else:
                        nanoAODStep = ',NANO:@PHYS+@L1'
                        
        self._checkRepackedFlag(options, **args)

        options.customisation_file=args['customs']

        options.step = self.recoSeq + PhysicsSkimStep
        options.step += miniAODStep + nanoAODStep
        #options.step += ',DQM' + dqmStep
        options.step += ',ENDJOB'

        dictIO(options,args)
        options.conditions = gtNameAndConnect(globalTag, args)
        
        process = cms.Process('RECO', cms.ModifierChain(self.eras, self.promptModifiers) )
        cb = ConfigBuilder(options, process = process, with_output = True)

        # Input source
        process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring()
        )
        cb.prepare()

        addMonitoring(process)
        
        return process
