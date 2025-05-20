import FWCore.ParameterSet.Config as cms

def Cut(cut, src=None, valtype=None):
    if   valtype == float: valtype = "float"
    elif valtype == int:   valtype = "int"
    elif valtype == bool:  valtype = "bool"
    if (src): # value map
        return cms.PSet(
                    cut = cms.string(cut),
                    src = src if isinstance(src, cms.InputTag)
                              else cms.InputTag(src),
                    type = cms.string(valtype)
               )
    else:
        return cms.PSet(cut = cms.string(cut))

def Cuts(cuts, mode="", expr=""):
    if (mode and expr):
        raise ValueError("Must specify mode or expr, but not both!")
    if (not mode and not expr): # both are empty so default to mode = "ALL"
        mode = "ALL"
    cut_list = list()
    if (expr): # substitute name of a cut with index
        for cut_idx, (cut_name, cut) in enumerate(cuts.items()):
            expr = expr.replace(cut_name, f"c[{cut_idx}]")
            cut_list.append(cut)
        return cms.PSet(
            expr = cms.string(expr),
            cuts = cms.VPSet(*cut_list)
        )
    else:
        for cut_name, cut in cuts.items():
            cut_list.append(cut)
        return cms.PSet(
            mode = cms.string(mode),
            cuts = cms.VPSet(*cut_list)
        )
