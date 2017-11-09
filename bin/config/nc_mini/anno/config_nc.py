import sys
import os
NOAHCARE=os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), '..','..','..'))
config={
    "trans" : NOAHCARE+"/config/nc_mini/anno/tr_nc.lst",
    "ExAC_POP"           : "EAS",
    "FLANK_LEN"	       : 2,
    "intron_edge_remain" : 10,
    "autoInterp_freq_threshold" : {
	"TGP_AF"    : 0.05,
	"ESP6500AF" : 0.05,
	"PVFD_AF"   : 0.05,
	"PanelAF"   : 0.05,
	"ExAC_AF"   : 0.05,
    },
    "kickoutExcel_function" : {
       	"."        : 1,
	"no-change" : 1,
	"annotation-fail" : 1,
	"unknown-no-call" : 1,
    },

    "possible_deleterious" : 
	"coding-synon missense stop-retained splice splice-3 splice-5 span abnormal-intron".split()
    ,

    "Panel_Control" : NOAHCARE+"/db/anno/panelDB/NC_V1.SNV.bed.gz",

}