#!/usr/bin/env python3
"""
_hltScouting_Run3_2024_

Scenario supporting proton collisions with input HLT scouting data for 2024

"""

import os
import sys

from Configuration.DataProcessing.Reco import Reco
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

from Configuration.DataProcessing.Impl.hltScouting import hltScouting

class hltScouting_Run3_2024(hltScouting):
    def __init__(self):
        hltScouting.__init__(self)
        self.recoSeq=''
        self.cbSc='pp'
        self.eras=Run3_2024
        self.promptCustoms += [ 'Configuration/DataProcessing/RecoTLR.customisePostEra_Run3_2024' ]
        #self.expressCustoms += [ 'Configuration/DataProcessing/RecoTLR.customisePostEra_Run3_2024' ]
        #self.visCustoms += [ 'Configuration/DataProcessing/RecoTLR.customisePostEra_Run3_2024' ]
    """
    _hltScouting_Run3_2024_
    Implement configuration building for data processing for proton
    collision data taking with input HLT scouting data for Era_Run3_2024
    """
