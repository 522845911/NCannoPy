import sys
import os

NOAHCARE=os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), '..','..','..'))
config={
    "trans" : NOAHCARE+"/config/nc_v12/anno/tr_nc.lst",
    "ExAC_POP"           : "EAS",
    "FLANK_LEN"           : 2,
    "intron_edge_remain" : 10,
    "autoInterp_freq_threshold" : {
    "TGP_AF"    : 0.05,
    "ESP6500AF" : 0.05,
    "PanelAF"   : 0.05,
    "ExAC_AF"   : 0.05,
    },
    "kickoutExcel_function" : {
           "."        : 1,
        "promoter" : 1,
        "intron"   : 1,
        "utr-3"    : 1,
        "utr-5"    : 1,
        "unknown" : 1,
    "no-change" : 1,
    "annotation-fail" : 1,
    "unknown-no-call" : 1,
    },

    "possible_deleterious" : 
    "coding-synon missense stop-retained splice splice-3 splice-5 span abnormal-intron".split()
    ,

    "Panel_Control" : NOAHCARE+"/annodb/panelDB/Hcan2_V1.0.HIGHQ.bed.gz",

#   If you want to specify your own source tag
#   you can uncomment GAP_RULE_SOURCE, and delete BGIGaP_rule
#   GAP_RULE_SOURCE : "IAD;BGI;ClinVar;LOVD:EMV;MSV3d",

    "GeneTestCode"  : NOAHCARE+"/config/nc_v12/anno/ht_part.tcd.gz",
    "CoreFuncRegion" : NOAHCARE+"/config/MET/met_exon_skip.bed.gz",
}