#!/usr/bin/env python
# coding: utf-8
#$ -S /usr/bin/python
#$ -cwd
# Copyright (c) GenePlus 2015
# Program       : NCanno.py
# Author        : Liu Tao
# Program Date  : 2015-07-08
# Modifier      : 
# Last Modified : 
# Description   : annotate vcf by BedAnno, and database API, give a well formatted serious mutants, 
#                 grouped by transcript id, together with the protein id
# Dependency    : This script depend on the bgzip,tabix and Faidx module
#                 to be available, and the whole genome fasta to be read from 
#                 db/anno/aln_db/hg19/ named "hg19_chM.fa", which change 
#                 the original chrM of hg19 to chrM_NC_012920.1.
#                 Also the BedAnno annotation databases are needed to annotate.

from __future__ import division
from string import maketrans
import multiprocessing
from multiprocessing import Process,Queue,Pool,Manager
from pprint import pprint
import path
import os
import sys
import argparse
import copy
import re
import gzip
import pprint
import time


sys.path.append(os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'BedAnno')))
sys.path.append(os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'NClib')))

import BedAnno
import CheckInOut
import pysam
import argparse
import time
import vcf
import bz2
import GetOmim
import GetHGMD
import GetClinVar
import GetCGDinfo
import GetVcfAF
import logging
import functools
import BedAnnoVar
import GetCytoBand
import GetPfam
import GetPrediction
import CondelPred
import GetCytoBand
import GetRepeatTag
import GetGWAS
import GetPfam
import CondelPred
import GetPhyloP46wayScore
import GetDBSNP
import GetCOSMIC
import GetVcfAF
import GetExAC
import GetGAD
import GetCGpub
def LogExc(f):
    @functools.wraps(f)
    def wrapper(*largs, **kwargs):
        try:
            res=f(*largs, **kwargs)
        except Exception, e:
            multiprocessing.get_logger().exception(e)
            raise e
        return res
    return wrapper
@LogExc
def normalise_null(result):
	def map_func(i):
		if i is None:
			i='.'
		i=re.sub(r'^\\N$|^N\/A$|^\s*$|^null$', '.', str(i), flags=re.IGNORECASE)
		return i
	return map(map_func,result)

# check the predict summary.
@LogExc
def judge_pred(mPred):
	mPred = re.sub(r'\.','', mPred) # remove no prediction record
	total = len(mPred)
	dnum = mPred.count("D")
	if dnum >= total / 2.0:
		return "D"
	else:
		return "P"

# Reverse and complement
@LogExc
def revcom(Seq):
	Seq=Seq[::-1]
	intab =  "ATCGatcgRYKMSWBDHVrykmswbdhv"
	outtab = "TAGCtagcYRMKSWVHDByrmkswvhdb"
	trantab = maketrans(intab, outtab)
	Seq=Seq.translate(trantab)
	return Seq


@LogExc
def read_HomoRef_var(file):
	logger=multiprocessing.get_logger()
	logger.info("113 file=")
	logger.info(file)
	HA_h=None
	if re.search(r'\.gz$',file):
		try:
			HA_h = gzip.open(file, 'rb')
		except Exception, e:
			raise Error("Error: [" + file + "] GunzipError+\n"+e)
	else:
		HA_h = open(file, 'rb')
	homoRef_var = dict()
	for line in HA_h:
		line=re.sub(r'\s+$','',line)
		if re.search(r'^\s*#|^\s*$',line):
			continue
		chr, start, end, ref=re.split(r'\s+',line,4)[:4]
		if not ref:
			raise Exception("Error format ["+file+"]")
		cur_var ={
			"chr":chr,
			"begin":start,
			"end":end,
			"referenceSequence":ref,
			"variantSequence":ref
		}
		if chr not in homoRef_var:
			homoRef_var[chr]=list()
		homoRef_var[chr].append(cur_var)
	HA_h.close()
	def sort_vars_f(a, b):
		return cmp(a['begin'],b['begin']) or cmp(a['end'],b['end'])
	for chrom in homoRef_var.keys():
		sorted_vars = sorted(homoRef_var[chrom],cmp=sort_vars_f)
		homoRef_var[chrom]=sorted_vars
	return homoRef_var

# read sample infos.
@LogExc
def read_sample_list(f):
	sample_info = dict()
	F=None
	try:
		F=open(f,"r")
	except Exception, e:
		raise Error("Error: ["+f+"] \n"+e)
	for line in F:
		line=re.sub(r'\s+$','',line)
		if re.search(r'^\s*#|^\s*$',line):
			continue
		all_ids = re.split(r'\s+',line)
		if all_ids[0] not in sample_info:
			sample_info[all_ids[0]]=dict()
		sample_info[all_ids[0]]['sex']=all_ids[1]
		sample_info[all_ids[0]]['fam']=all_ids[2]
	F.close()
	return sample_info
@LogExc
def read_inmsqc(msqc_result):
	msqc_input = dict()
	MSQCRST=None
	try:
		MSQCRST=open(msqc_result,"r")
	except Exception, e:
		raise Error("Error: ["+msqc_result+"] \n"+e)
	for qcrst in MSQCRST:
		if re.search(r'^#',qcrst):
			continue
		qcrst = qcrst.strip('\r\n')
		itms = re.split(r'\s+',qcrst)
		sampid=itms.pop(0)
		msqc_input[sampid]=[i.upper() for i in itms]
	MSQCRST.close()
	return msqc_input
@LogExc
def read_MSQC_info(msqc_info_file):
	msqc_sites = dict()
	MSQCINFO=None
	try:
		MSQCINFO=open(msqc_info_file,'r')
	except Exception, e:
		raise Error("Error: ["+msqc_info_file+"] \n"+e)
	for qcl in MSQCINFO:
		if re.search(r'^#',qcl):
			continue
		qcl=qcl.strip("\r\n")
		itms = re.split(r'\t',qcl,4)
		sid = itms[0]
		sid =re.sub(r'^.+\.','',sid)
		if re.match(r'^\d+$',sid):
			sid-=1
		else:
			raise Exception("{}{}{}".format("\n", "Error: ["+msqc_info+"] ID format in MSQC sites list","should be like 'XXXX.YY', whose 'Y' should be all number\n"))
		if itms[1] not in msqc_sites:
			msqc_sites[itms[1]]=dict()
		if sid not in msqc_sites[itms[1]]:
			msqc_sites[itms[1]][sid]=dict()
		msqc_sites[itms[1]][sid]={"pos":itms[2],"ref":itms[3].upper()}
	MSQCINFO.close()
	return msqc_sites


@LogExc
def judge_titv(r, a):
	to_check ="".join(map(str,sorted([r,a])))
	if to_check=='AG' or to_check=='CT':
		return 1
	else:
		return 0
@LogExc
def checkYPAR(begin,end):
	for rypar in ypar_reg:
		if begin<rypar[1] and end>rypar[0]:
			return 0
	return 1
@LogExc
def checkXPAR(begin,end):
	for rxpar in xpar_reg:
		if begin<rxpar[1] and end>rxpar[0]:
			return 0
	return 1
@LogExc
def subsame(val_ori, *string_array):
	same_opt = 1
	for item in string_array:
		if val_ori!=item:
			same_opt = 0
			break
	if 0==len(string_array) or same_opt == 1:
		return val_ori
	else:
		string_array=list(string_array)
		string_array=[val_ori]+string_array
		return "|".join(map(str,string_array))

@LogExc
def uniform_chrVarOut(varanno, sampleid):
	varout = dict()
	chrvar = varanno["var"]
	sample_h = chrvar['sample'][sampleid]
	ExAC_pop=config['ExAC_POP']
	GAD_pop = 'EAS'
	
	# must have keys in var
	varout=dict(zip('Chr Start Stop Ref VarType Call'.split(),[chrvar[key] for key in "chr pos end ref guess alt".split()]))
	
	tmpids = \
	'''
	Flank MapLoc RepeatTag
	dbsnpNote dbsnpAC dbsnpAF dbsnpMTag
	TGPpopAF TGP_AF ESP6500AC ESP6500AF
	ExAC_Filter ExAC_pop_HASC ExAC_pop_AC
	ExAC_pop_AF ExAC_HASC ExAC_AC ExAC_AF
	GAD_pop_HASC GAD_pop_AC
	GAD_pop_AF GAD_HASC GAD_AC GAD_AF
	rsID PVFD_RULE
	PVFDhomaC PVFDhomaF PVFD_AC PVFD_AF
	PanelID PanelAC PanelAF 
	phyloPpr phyloPve phyloPpm
	cosmicName cosmicSite cosmicHis cosmicStat cosmicID
	ClinSSR ClinACC ClinRevStat ClinSignificant
	'''.split()
	varout.update(dict.fromkeys(tmpids,'.'))
	
	if chrvar['guess']!="annotation-error":
		# must have if annotation OK
		varout.update(dict(zip('Flank MapLoc'.split(),[chrvar[key] for key in 'flanks cytoBand'.split()])))
	
	# for single point only
	if 'phyloPpr' in chrvar :
		varout.update(dict(zip('phyloPpr phyloPve phyloPpm'.split(),[chrvar[key] for key in 'phyloPpr phyloPve phyloPpm'.split()])))
	
	# keys for sample
	varout['SampleID']=sampleid
	if sample_list is not None:
		if sampleid in sample_info:
			varout['FamID']=sample_info[sampleid]['fam']
		else:
			arout['FamID']='.'
		if sampleid in sample_info:
			varout['Sex']=sample_info[sampleid]['sex']
		else:
			varout['Sex']='.'
	
	varout.update(dict(zip('Zygosity PhasedGID A.Index Filter'.split(),[sample_h[key] for key in 'zygosity PB AI filterTag'.split()])))
	
	# fix zygosity info for male's sex chromosome
	if 'Sex' in varout and varout['Sex']=='M' and ((re.search(r'^Y',varout['Chr'],flags=re.IGNORECASE) and checkYPAR(varout['Start'],varout['Stop'])) or (re.search(r'^X',varout['Chr'],flags=re.IGNORECASE) and checkXPAR(varout['Start'],varout['Stop']))):
		varout['Zygosity']=re.sub(r'^hom-','hem-',varout['Zygosity'])
	
	if 'NB' in sample_h:
		varout['NbGID']=sample_h['NB']
	else:
		varout['NbGID']='.'
	if 'AD' in sample_h:
		varout['A.Depth']=sample_h['AD']
	else:
		varout['A.Depth']='.'
	if 'AR' in sample_h:
		varout['A.Ratio']=sample_h['AR']
	else:
		varout['A.Ratio']='.'
	
	if 'reptag' in chrvar:
		varout['RepeatTag']=chrvar['reptag']
	
	# keys for AF info
	if 'dbsnp' in chrvar and 0<len(chrvar['dbsnp'].keys()):
		rs_ids=sorted(chrvar['dbsnp'].keys())
		dbsnpNotes, dbsnpACs, dbsnpAFs, dbsnpMTags  = list(),list(),list(),list()
		for rs in rs_ids:
			cur_rs_h = chrvar['dbsnp'][rs]
			tmp=None
			if cur_rs_h['exception']==".":
				tmp=cur_rs_h['bitfields']
			else:
				if cur_rs_h['bitfields']==".":
					tmp=cur_rs_h['exception']
				else:
					tmp=cur_rs_h['exception']+";"+cur_rs_h['bitfields']
			dbsnpNotes.append(tmp)
			dbsnpACs.append(cur_rs_h['AN'])
			dbsnpAFs.append(cur_rs_h['AF'])
			dbsnpMTags+=["U" if int(cur_rs_h['weight'])==1 else "M"]
		varout['rsID']=subsame(*rs_ids)
		varout['dbsnpNote']=subsame(*dbsnpNotes)
		varout['dbsnpAC']=subsame(*dbsnpACs)
		varout['dbsnpAF']=subsame(*dbsnpAFs)
		varout['dbsnpMTag']=subsame(*dbsnpMTags)
		
	if 'tgp' in chrvar and 'AF' in chrvar['tgp']:
		varout['TGP_AF']=chrvar['tgp']['AF']
		if 'POP_AF' in chrvar['tgp']:
			varout['TGPpopAF']=chrvar['tgp']['POP_AF']
	
	if 'esp6500' in chrvar and 'AF' in chrvar['esp6500']:
		varout['ESP6500AC']=chrvar['esp6500']['AN']
		varout['ESP6500AF']=chrvar['esp6500']['AF']
	
	if 'exac' in chrvar and 'AF' in chrvar['exac']:
		exac_hash=chrvar['exac']
		varout['ExAC_Filter']=exac_hash['Filter']
		if "HASC_" + ExAC_pop in exac_hash:
			varout['ExAC_pop_HASC']=exac_hash["HASC_" + ExAC_pop]
		if "AN_" + ExAC_pop in exac_hash:
			varout['ExAC_pop_AC']=exac_hash["AN_" + ExAC_pop]
		if "AF_" + ExAC_pop in exac_hash:
			varout['ExAC_pop_AF']=exac_hash["AF_" + ExAC_pop]
		if 'HASC' in exac_hash:
			varout['ExAC_HASC']=exac_hash['HASC']
		if 'AN' in exac_hash:
			varout['ExAC_AC']=exac_hash['AN']
		if 'AF' in exac_hash:
			varout['ExAC_AF']=exac_hash['AF']
	
	if 'gnomAD' in chrvar and 'AF' in chrvar['gnomAD']:
		gnomAD_hash=chrvar['gnomAD']
		
		
		if "HASC_" + GAD_pop in gnomAD_hash:
			varout['GAD_pop_HASC']=gnomAD_hash["HASC_"+GAD_pop]
		if "AN_" + GAD_pop in gnomAD_hash:
			varout['GAD_pop_AC']=gnomAD_hash["AN_" + GAD_pop]
		if "AF_" + GAD_pop in gnomAD_hash:
			varout['GAD_pop_AF']=gnomAD_hash["AF_" + GAD_pop]
		if 'HASC' in gnomAD_hash:
			varout['GAD_HASC']=gnomAD_hash['HASC']
		if 'AN' in gnomAD_hash:
			varout['GAD_AC']=gnomAD_hash['AN']
		if 'AF' in gnomAD_hash:
			varout['GAD_AF']=gnomAD_hash['AF']
	
	if 'PanelSQL_Fail' in chrvar:
		varout['PanelAF']='error'
	
	if 'panelctrl' in chrvar and 'AF' in chrvar['panelctrl']:
		varout['PanelID']=config['Panel_ID']
		varout['PanelAC']=chrvar['panelctrl']['AN']
		varout['PanelAF']=chrvar['panelctrl']['AF']
	
	# clinVar
	if 'clinVar_Fail' in chrvar:
		varout['ClinSignificant']='error'
	if 'clinVar' in chrvar and 'SSR' in chrvar['clinVar']:
		varout.update(dict(zip('ClinSSR ClinACC ClinRevStat ClinSignificant'.split(),[chrvar['clinVar'][key] for key in 'SSR CLNACC CLNREVSTAT CLNSIG'.split()])))
	
	# cosmic
	if 'cosmic' in chrvar and 'mutID' in chrvar['cosmic']:
		varout.update(dict(zip('cosmicName cosmicSite cosmicHis cosmicStat cosmicID'.split(),[chrvar['cosmic'][key] for key in 'mutName site histology status mutID'.split()])))
	return varout
@LogExc
def trim_pmID(pmID):
	tmp=re.findall(r'\d+',pmID)
	has_pmID = dict.fromkeys(tmp,1)
	rejoin = '|'.join(sorted(has_pmID.keys(),cmp=lambda x,y:cmp(int(x),int(y)),reverse=True))
	return rejoin


def uniform_trVarOut(varanno, varout, tr):
	#print "varout="
	#print varout
	new_var_out = copy.deepcopy(varout)
	ids= \
	'''
	mimGID mimStat mimPIDs mimInhs
	cgdCond cgdInhs cgdACond cgdManCat cgdIntCat cgdRef
	TestCode EntrezGeneID geneSym
	FuncRegion ExIn_ID pfamName pfamId Function
	Transcript Protein Strand Primary cHGVS pHGVS
	CodonChange PolarChange MutationName
	ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred
	hgmdID hgmdQMode hgmdName hgmdDis hgmdPMID hgmdPred
	'''.split()
	new_var_out.update(dict.fromkeys(ids,'.'))
	
	if tr!="":
		
		curtr=varanno["trInfo"][tr]
		new_var_out['Transcript']=tr
		
		new_var_out.update(dict(zip('EntrezGeneID geneSym TestCode FuncRegion ExIn_ID Function Protein Strand Primary cHGVS MutationName'.split(),[curtr[key] for key in 'geneId geneSym TCD r exin func prot strd major c trVarName'.split()])))
		
		geneVarNameRef=None
		custom_VarName=None
		if 'VarNameList' in config:
			if tr in VarNameList:
				geneVarNameRef=VarNameList[tr]
			elif curtr['geneSym'] in VarNameList:
				geneVarNameRef=VarNameList[curtr['geneSym']]
			elif curtr['geneId'] in VarNameList:
				geneVarNameRef=VarNameList[curtr["geneId"]]
			if geneVarNameRef is not None:
				if curtr["c"] in geneVarNameRef:
					custom_VarName=geneVarNameRef[curtr["c"]]
				else:
					tmp = curtr["c"]
					if re.search(r'\-u',curtr['c']):
						m=re.search(r'\-(\d+)\-u(\d+)',tmp)
						tmp=re.sub(r'\-(\d+)\-u(\d+)',"-"+str(int(m.group(1))+int(m.group(2))),tmp)
						tmp=re.sub(r'\-u',r'\-',tmp)
					
					if re.search(r'\+d',curtr['c']):
						m=re.search(r'\*(\d+)\+d(\d+)',tmp)
						tmp=re.sub(r'\*(\d+)\+d(\d+)',"*"+str(int(m.group(1))+int(m.group(2))),tmp)
						tmp=re.sub(r'\+d',r'\*',tmp)
					
					if re.search(r'dup',curtr['c']):
						tmp=re.sub(r'dup.+$','dup',tmp)
					
					if re.search(r'del',curtr['c']):
						tmp =re.sub(r'del.+$','del',tmp)
					
					if geneVarNameRef is not None and tmp in geneVarNameRef:
						custom_VarName=geneVarNameRef[tmp]
		
		if 'standard_cHGVS' in curtr or 'alt_cHGVS' in curtr:
			new_var_out['cHGVS']="{}{}".format(new_var_out['cHGVS'],' (')
			if 'standard_cHGVS' in curtr:
				new_var_out['cHGVS']="{}{}{}{}".format(new_var_out['cHGVS'],'std: ',curtr['standard_cHGVS'],' ')
			if 'alt_cHGVS' in curtr:
				new_var_out['cHGVS']="{}{}{}{}".format(new_var_out['cHGVS'],'alt: ',curtr['alt_cHGVS'],' ')
			new_var_out['cHGVS']="{}{}".format(new_var_out['cHGVS'],')')
			
			if geneVarNameRef is not None and custom_VarName is None:
				if curtr['standard_cHGVS'] in geneVarNameRef:
					custom_VarName=geneVarNameRef[curtr['standard_cHGVS']]
				elif curtr['alt_cHGVS'] in geneVarNameRef:
					custom_VarName=geneVarNameRef[curtr['alt_cHGVS']]
		
		if 'p' in curtr:
			new_var_out['pHGVS']=curtr['p']
			if geneVarNameRef is not None and custom_VarName is None and curtr['p'] in geneVarNameRef:
				custom_VarName=geneVarNameRef[curtr['p']]
		
		if 'standard_pHGVS' in curtr or 'alt_pHGVS' in curtr:
			new_var_out['pHGVS']="{}{}".format(new_var_out['pHGVS'],' (')
			if 'standard_pHGVS' in curtr:
				new_var_out['pHGVS']="{}{}{}{}".format(new_var_out['pHGVS'],'std: ',curtr['standard_pHGVS'],' ')
			if 'alt_pHGVS' in curtr:
				new_var_out['pHGVS']="{}{}{}{}".format(new_var_out['pHGVS'],'alt: ',curtr['alt_pHGVS'],' ')
			new_var_out['pHGVS']="{}{}".format(new_var_out['pHGVS'],')')
			
			if geneVarNameRef is not None and custom_VarName is None:
				if 'strandard_pHGVS' in curtr and curtr['standard_pHGVS'] in geneVarNameRef:
					custom_VarName=geneVarNameRef[curtr['standard_pHGVS']]
				elif 'alt_pHGVS' in curtr and curtr['alt_pHGVS'] in geneVarNameRef:
					custom_VarName=geneVarNameRef[curtr['alt_pHGVS']]
		if 'p3' in curtr and curtr['p3']!=curtr['p']:
			new_var_out['pHGVS']="{}{}{}".format(new_var_out['pHGVS'],' | ',curtr['p3'])
			if geneVarNameRef is not None and custom_VarName is None and curtr['p3'] in geneVarNameRef:
				custom_VarName=geneVarNameRef[curtr['p3']]
			if 'standard_p3' in curtr or 'alt_p3' in curtr:
				new_var_out['pHGVS']="{}{}".format(new_var_out['pHGVS'],' (')
				if 'standard_p3' in curtr:
					new_var_out['pHGVS']="{}{}{}{}".format(new_var_out['pHGVS'],'std: ',curtr['standard_p3'],' ')
				if 'alt_p3' in curtr:
					new_var_out['pHGVS']="{}{}{}{}".format(new_var_out['pHGVS'],'alt: ',curtr['alt_p3'],' ')
				new_var_out['pHGVS']="{}{}".format(new_var_out['pHGVS'],')')
				if geneVarNameRef is not None and custom_VarName is None:
					if curtr['standard_p3'] in geneVarNameRef:
						custom_VarName=geneVarNameRef[curtr['standard_p3']]
					elif curtr['alt_p3'] in geneVarNameRef:
						custom_VarName=geneVarNameRef[curtr['alt_p3']]
		
		if custom_VarName is not None:
			new_var_out['MutationName']="{}{}{}".format(new_var_out['MutationName'],' / ',custom_VarName)
		if "cc" in curtr:
			new_var_out['CodonChange']=curtr['cc']
		if "polar" in curtr:
			new_var_out['PolarChange']=curtr['polar']
		
		# pfam
		if 'pfamId' in curtr and curtr['pfamId'] is not None:
			new_var_out.update(dict(zip('pfamId pfamName'.split(),[curtr[key]  for key in 'pfamId pfamName'.split()])))
		
		# sift condel
		if 'condelPred' in curtr and curtr['condelPred'] is not None:
			new_var_out.update(dict(zip('ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred'.split(),[curtr[key] for key in 'siftScore pp2varScore pp2divScore condelScore condelPred'.split()])))
		
		if 'GeneSQL_Fail' in curtr:
			new_var_out['mimInhs']='error'
			new_var_out['cgdRef']='error'
		
		# mim
		if 'omim' in curtr and 0<len(curtr['omim'].keys()):
			genemims=sorted(curtr['omim'].keys())
			genestat, disomims, inhs=list(),list(),list()
			for gm in genemims:
				intab = "|"
				outtab = ";"
				trantab = maketrans(intab, outtab)
				curtr['omim'][gm]['genestat']=curtr['omim'][gm]['genestat'].translate(trantab)
				curtr['omim'][gm]['disomims']=curtr['omim'][gm]['disomims'].translate(trantab)
				curtr['omim'][gm]['inhs']=curtr['omim'][gm]['inhs'].translate(trantab)
				genestat.append(curtr['omim'][gm]['genestat'])
				disomims.append(curtr['omim'][gm]['disomims'])
				inhs.append(curtr['omim'][gm]['inhs'])
			new_var_out.update(dict(zip('mimGID mimStat mimPIDs mimInhs'.split(),["|".join(genemims),subsame(*genestat),subsame(*disomims),subsame(*inhs)])))
		
		# cgd
		if 'cgd' in curtr and 'geneid' in curtr['cgd']:
			new_var_out['cgdRef']=trim_pmID(curtr['cgd']['references'])
			new_var_out.update(dict(zip('cgdCond cgdInhs cgdACond cgdManCat cgdIntCat'.split(),[curtr['cgd'][key] for key in 'cond inhs allelic_cond manifest_cat intervent_cat'.split()])))
		
		if 'LocalHGMD_Fail' in curtr:
			new_var_out['hgmdPred']='error'
		
		# hgmd
		if 'hgmd' in curtr and 0< len(curtr['hgmd']):
			hgmdID, hgmdQMode, hgmdName,hgmdDis, hgmdPMID, hgmdPred=list(),list(),list(),list(),list(),list()
			all_ret=sorted(curtr['hgmd'],cmp=lambda x,y:cmp(a['id'],b['id']))
			for ret1 in all_ret:
				hgmdID+=[ret1['id']]
				if 'querylevel' in ret1:
					hgmdQMode.append(ret1['querylevel'])
				hgmdName.append(ret1['mutation_name'])
				hgmdDis.append(ret1['disease'])
				hgmdPMID.append(ret1['pmid'])
				hgmdPred.append(ret1['pred'])
			new_var_out['hgmdID']=subsame(*hgmdID)
			if len(hgmdQMode)>0:
				new_var_out['hgmdQMode']=subsame(*hgmdQMode)
			new_var_out['hgmdName']=subsame(*hgmdName)
			new_var_out['hgmdDis']=subsame(*hgmdDis)
			new_var_out['hgmdPMID']=trim_pmID(subsame(*hgmdPMID))
			new_var_out['hgmdPred']=subsame(*hgmdPred)
	else:
		if 'varName' in varanno["var"]:
			new_var_out['MutationName']=varanno['var']['varName']
	
	return new_var_out
@LogExc
def xor(a,b):
	return  (a and (not b) ) or ( (not a) and b)

# AutoInterp status comes from the following options:
# 1. whether related to a DM mutation record
# 2. whether a high quality calling
# 3. whether with a deleterious function
# 4. whether have lower allele frequency than threshold
# if all are yes: then "Disease Causing Mutation",
# if no for all without consider quality, then "Benign",
# if 1,4 or 3,4 are yes: then "Possibly Deleterious"
# if only 2,1(and/or)3, only 2,4 are yes, then "Possibly Benign",
# other case, then "Unknown"
@LogExc
def AutoInterp(rec_ent):
	global autoInterp_freq_threshold
	logger=multiprocessing.get_logger()
	Phyno_pred_opt = 0
	Func_pred_opt = 0
	freq_opt   = 1
	filter_opt = 0
	
	# phyno_pred hgmd and GaP
	if re.search(r'DM',rec_ent['hgmdPred']) or re.search(r'DFP',rec_ent['hgmdPred']) or re.search(r'pathogenic',rec_ent['ClinSignificant'],flags=re.IGNORECASE):
		Phyno_pred_opt = 1
	
	# filter
	if rec_ent['Filter']=='PASS' or rec_ent['Filter']=='.':
		filter_opt = 1
	
	# func
	if rec_ent['Function'] in deleterious_func or rec_ent['Function'] in possible_deleterious_func or rec_ent['ensCondelPred']=='deleterious':
		Func_pred_opt = 1
	
	# freq
	for freq_k in autoInterp_freq_threshold.keys():
		if freq_k in rec_ent and rec_ent[freq_k]!='.' and rec_ent[freq_k]!='error' and float(rec_ent[freq_k])>float(autoInterp_freq_threshold[freq_k]):
			freq_opt = 0
			break
	
	# if all are yes: then "Certain" (Disease Mutation),
	if Phyno_pred_opt == 1 and filter_opt == 1 and Func_pred_opt == 1 and freq_opt == 1:
		return "Certain"
	# if no for all without consider quality, then "Benign",
	elif Phyno_pred_opt == 0 and Func_pred_opt == 0 and freq_opt == 0:
		return "Benign"
	# if 1,4 are yes: then "Likely Deleterious"
	elif Phyno_pred_opt == 1 and freq_opt == 1:
		return "Likely Deleterious"
	# if 3,4 are yes: then "VUS"
	elif Func_pred_opt == 1 and freq_opt == 1:
		return "VUS"
	# if only 2,1(and/or)3, only 2,4 are yes, then "Likely Benign",
	elif filter_opt == 1 and  (xor((Phyno_pred_opt == 1 or Func_pred_opt == 1) , freq_opt )):
		return "Likely Benign"
	# other case, then "Unknown"
	else:
		return "Unknown"

# Current InExcel Strategy 
# will depend on autoInterpRules
# defined in config
# Or default rules will be implemented
@LogExc
def check_inExcel(rec_ent):
	if rec_ent['Function'] in kickoutExcel_function and (re.search(r'Benign',rec_ent['autoInterp']) or rec_ent['autoInterp']=="Unknown"):
		if rec_ent['Function']=="intron":
			if int(config['intron_edge_remain'])<=0:
				return 1
			m=re.finditer(r'\d+[\+\-](\d+)',rec_ent['cHGVS'])
			for i in m:
				if int(i.group(1))<int(config['intron_edge_remain']):
					return 1
		return 0
	
	return 1

# count while output

def print_cache(anno_cache):
	global ti
	global tv
	global var_count
	global sex_varcount
	global no_call_var_count
	global annotation_count
	global OutFp
	global out_header_keys
	def compare_rule(a,b):
		a_val, b_val = a["var"]["chr"], b["var"]["chr"]
		if a["var"]["chr"]=='X':
			a_val=23
		if b["var"]["chr"]=='X':
			b_val=23
		if a["var"]["chr"]=='Y':
			a_val=24
		if b["var"]["chr"]=='Y':
			b_val=24
		a_val, b_val=int(a_val), int(b_val)
		return cmp(a_val,b_val) or cmp(int(a["var"]["pos"]),int(b["var"]["pos"])) or cmp(int(a["var"]["end"]),int(b["var"]["end"]))
	
	# already same chromosome, only need to sort position
	sort_anno=sorted(anno_cache,cmp=compare_rule)
	
	for var_ent in sort_anno:
		sam_keys=sorted(var_ent["var"]["sample"].keys())
		for sample in sam_keys:
			curSample_h = var_ent["var"]["sample"][sample]
			invar_allele = 1
			if int(curSample_h['AI'])==0 and int(curSample_h['PLOIDY'])==2:
				invar_allele = 2
			if var_ent["var"]["guess"]=='snv':
				if judge_titv(var_ent["var"]["ref"],var_ent["var"]["alt"]):
					if sample not in ti:
						ti[sample]=0
					ti[sample]+=invar_allele
				else:
					if sample not in tv:
						tv[sample]=0
					tv[sample]+=invar_allele
			
			if re.search(r'^[xy]',var_ent["var"]["chr"],flags=re.IGNORECASE):
				# sex chromosome
				if var_ent["var"]["guess"]!='no-call':
					if sample not in sex_varcount:
						sex_varcount[sample]=dict()
					if 'total' not in sex_varcount[sample]:
						sex_varcount[sample]['total']=0
					sex_varcount[sample]['total']+=invar_allele
					if 'het' not in sex_varcount[sample]:
						sex_varcount[sample]['het']=0
					if re.search(r'het',curSample_h['zygosity']):
						sex_varcount[sample]['het']+=invar_allele
			
			if var_ent["var"]["guess"]=='no-call':
				if sample not in no_call_var_count:
					no_call_var_count[sample]=0
				no_call_var_count[sample]+=invar_allele
			else:
				if sample not in var_count:
					var_count[sample]=dict()
				if 'varType' not in var_count[sample]:
					var_count[sample]['varType']=dict()
				if 'zygosity' not in var_count[sample]:
					var_count[sample]['zygosity']=dict()
				if 'filter' not in var_count[sample]:
					var_count[sample]['filter']=dict()
				if 'total' not in var_count[sample]:
					var_count[sample]['total']=0
				if var_ent["var"]["guess"] not in var_count[sample]['varType']:
					var_count[sample]['varType'][var_ent["var"]["guess"]]=0
				if curSample_h['zygosity'] not in var_count[sample]['zygosity']:
					var_count[sample]['zygosity'][curSample_h['zygosity']]=0
				if curSample_h['filterTag'] not in var_count[sample]['filter']:
					var_count[sample]['filter'][curSample_h['filterTag']]=0
				var_count[sample]['total']+=invar_allele
				var_count[sample]['varType'][var_ent["var"]["guess"]]+=invar_allele
				var_count[sample]['zygosity'][curSample_h['zygosity']]+=invar_allele
				var_count[sample]['filter'][curSample_h['filterTag']]+=invar_allele
			
			arbitrary_remain=None
			if coreFuncRegion_h is not None:
				if coreFuncRegion_h.checkIn(var_ent["var"]["chr"],var_ent["var"]["pos"],var_ent["var"]["end"]):
					arbitrary_remain = 1
			#print "var_ent="
			#print var_ent
			cur_var_out = uniform_chrVarOut(var_ent, sample)
			if 'trInfo' in var_ent:
				for tr in sorted(var_ent["trInfo"].keys()):
					tr_var_out=uniform_trVarOut(var_ent, cur_var_out, tr)
					tr_var_out['autoInterp']=AutoInterp(tr_var_out)
					
					# must have autoInterp first
					if arbitrary_remain is not None and arbitrary_remain == 1:
						tr_var_out['InExcel']=1
					else:
						tr_var_out['InExcel']=check_inExcel(tr_var_out)
					
					# only count the major transcript for called-allele
					# to avoid multiple transcript dup
					if var_ent["var"]["guess"]!='no-call' and var_ent["var"]["varName"]==tr_var_out['MutationName']:
						if sample not in annotation_count:
							annotation_count[sample]=dict()
						if 'autoInterp' not in annotation_count[sample]:
							annotation_count[sample]['autoInterp']=dict()
						if tr_var_out['autoInterp'] not in annotation_count[sample]['autoInterp']:
							annotation_count[sample]['autoInterp'][tr_var_out['autoInterp']]=0
						if 'function' not in annotation_count[sample]:
							annotation_count[sample]['function']=dict()
						if tr_var_out['Function'] not in annotation_count[sample]['function']:
							annotation_count[sample]['function'][tr_var_out['Function']]=0
						if 'InExcel' not in annotation_count[sample]:
							annotation_count[sample]['InExcel']=dict()
						if tr_var_out['InExcel'] not in annotation_count[sample]['InExcel']:
							annotation_count[sample]['InExcel'][tr_var_out['InExcel']]=0
						if 'total' not in annotation_count[sample]:
							annotation_count[sample]['total']=0
						annotation_count[sample]['autoInterp'][tr_var_out['autoInterp']]+=invar_allele
						annotation_count[sample]['function'][tr_var_out['Function']]+=invar_allele
						annotation_count[sample]['InExcel'][tr_var_out['InExcel']]+=invar_allele
						annotation_count[sample]['total']+=invar_allele
					#global debug
					if debug :
						print >>sys.stderr,"uniformed var: "+tr_var_out
					#print "out_header_keys="
					#print out_header_keys
					outline=normalise_null([tr_var_out[key] for key in out_header_keys])
					print >>OutFp,"\t".join(outline)
			else:
				cur_var_out = uniform_trVarOut( var_ent, cur_var_out, "" )
				cur_var_out['autoInterp']=AutoInterp(cur_var_out)
				
				# must have autoInterp first
				if var_ent["var"]["guess"]!='no-call':
					if arbitrary_remain is not None and arbitrary_remain == 1:
						cur_var_out['InExcel']=1
					else:
						cur_var_out['InExcel']=check_inExcel(cur_var_out)
					if sample not in annotation_count:
						annotation_count[sample]=dict()
					if 'autoInterp' not in annotation_count[sample]:
						annotation_count[sample]['autoInterp']=dict()
					if cur_var_out['autoInterp'] not in annotation_count[sample]['autoInterp']:
						annotation_count[sample]['autoInterp'][cur_var_out['autoInterp']]=0
					if 'function' not in annotation_count[sample]:
						annotation_count[sample]['function']=dict()
					if cur_var_out['Function'] not in annotation_count[sample]['function']:
						annotation_count[sample]['function'][cur_var_out['Function']]=0
					if 'InExcel' not in annotation_count[sample]:
						annotation_count[sample]['InExcel']=dict()
					if cur_var_out['InExcel'] not in annotation_count[sample]['InExcel']:
						annotation_count[sample]['InExcel'][cur_var_out['InExcel']]=0
					if 'total' not in annotation_count[sample]:
						annotation_count[sample]['total']=0
					annotation_count[sample]['autoInterp'][cur_var_out['autoInterp']]+=invar_allele
					annotation_count[sample]['function'][cur_var_out['Function']]+=invar_allele
					annotation_count[sample]['InExcel'][cur_var_out['InExcel']]+=invar_allele
					annotation_count[sample]['total']+=invar_allele
					if debug:
						print >>sys.stderr,"uniformed var: "+str(cur_var_out)
					out_line=normalise_null([cur_var_out[key] for key in out_header_keys])
					print >>OutFp,"\t".join(out_line)

# current give filter assessment by
# judging the AD and PLtag, besides possible NB for VCF
# for tsv format(CG) using varFilter column
@LogExc
def give_filter(sample_h_in_var):
	if type=="tsv":
		# CG data
		if 'VF' in sample_h_in_var:
			if re.search(r'VQLOW',str(sample_h_in_var['VF']),flags=re.IGNORECASE):
				sample_h_in_var['filterTag']='FAIL'
			elif re.search(r'AMBIGUOUS',str(sample_h_in_var['VF']),flags=re.IGNORECASE):
				sample_h_in_var['filterTag']='DUBIOUS'
			elif re.search(r'^PASS$',str(sample_h_in_var['VF']),flags=re.IGNORECASE):
				sample_h_in_var['filterTag']='PASS'
			else:
				sample_h_in_var['filterTag']='.'
		else:
			sample_h_in_var['filterTag']='.'
	else:
		# VCF data
		if 'AD' not in sample_h_in_var or str(sample_h_in_var['AD'])=='.' or 'PLtag' not in sample_h_in_var or str(sample_h_in_var['PLtag'])=='.':
			sample_h_in_var['filterTag']='.'
		else:
			if int(sample_h_in_var['AD'])<=AD_DN_THRESHOLD or int(sample_h_in_var['PLtag'])<0:
				sample_h_in_var['filterTag']='FAIL'
			elif int(sample_h_in_var['AD'])>=AD_UP_THRESHOLD and int(sample_h_in_var['PLtag'])>0 and str(sample_h_in_var['VF'])=='.' and ('NB' not in sample_h_in_var or str(sample_h_in_var['NB'])=='.'):
				sample_h_in_var['filterTag']='PASS'
			else:
				sample_h_in_var['filterTag']='DUBIOUS'
	return sample_h_in_var
@LogExc
def tsv_var_parser(ritms, sample):
	var= {
		'locus'             : ritms[0],
		'chr'               : ritms[3],
		'begin'             : ritms[4],
		'end'               : ritms[5],
		'callerVarType'     : ritms[6],
		'referenceSequence' : ritms[7],
		'variantSequence'   : ritms[8]
	}
	
	AI, PB, VF, PLOIDY=None,None,None,None
	PB = '.' if str(ritms[12]) == '' else ritms[12]
	VF = '.' if str(ritms[11]) == '' else ritms[11]
	AI = '.' if str(ritms[2]) == '' else ritms[2]
	AI = 0  if str(AI) == 'all' else AI
	PLOIDY = ritms[1]
	
	var['sample'][sample]= {
		'PLOIDY':PLOIDY,
		'VF':VF,
		'AI':AI,
		'PB':PB
	}
	var['sample'][sample]=give_filter(var['sample'][sample])
	
	return var

# combind hom-alt mutation
@LogExc
def check_combin_group(rlg, rlsi):
	if 2 !=len(rlsi.keys()) or 2!=len(rlg):
		return [rlg, rlsi]
	if 1 not in rlsi or 2 not in rlsi:
		return [rlg, rlsi]
	v1 =rlg[0]
	v2 = rlg[1]
	samples=v1['sample'].keys()
	if v1['begin']==v2['begin'] and v1['end']==v2['end'] and v1['variantSequence']==v2['variantSequence']:
		for sam in samples:
			if v1['sample'][sam]['filterTag']!=v2['sample'][sam]['filterTag']:
				if (str(v1['sample'][sam]['filterTag'])!='PASS' and str(v1['sample'][sam]['filterTag'])!='.') or (str(v2['sample'][sam]['filterTag'])!='PASS' and str(v2['sample'][sam]['filterTag'])!='.'):
					rlg[0]['sample'][sam]['filterTag']='DUBIOUS'
			
			if v1['sample'][sam]['VF']!=v2['sample'][sam]['VF']:
				rlg[0]['sample'][sam]['VF']="{}{}{}".format(rlg[0]['sample'][sam]['VF'],"|",rlg[1]['sample'][sam]['VF'])
			
			rlg[0]['sample'][sam]['AI']=0
		del rlg[1]
		rlsi={0:[0]}
	
	return [rlg, rlsi]

@LogExc
def locus_group_parser(rlocus_group):
	# only 1 sample in tsv's var hash
	sample=rlocus_group[0]['sample'].keys()[:1]
	
	locus_strand_index = dict()
	for i in range(0,len(rlocus_group)):
		allele_index=rlocus_group[i]['sample'][sample]['AI']
		if allele_index not in locus_strand_index:
			locus_strand_index[allele_index]=list()
		locus_strand_index[allele_index].append(i)
	
	rlg, locus_strand_index=check_combin_group(rlocus_group,locus_strand_index)
	
	for ai in sorted(locus_strand_index.keys()):
		for var_idx in locus_strand_index[ai]:
			ai_var=rlg[var_idx]
			if ai_var['variantSequence']=='?':
				ai_var['sample'][sample]['zygosity']='no-call'
			else:
				cur_post ='ref' if re.search(r'ref',ai_var['callerVarType']) else 'alt'
				if ai == '0':
					pre, post=None,None
					if ai_var['sample'][sample]['PLOIDY']=='1':
						pre = "hap"
					else:
						pre = "hom"
					ai_var['sample'][sample]['zygosity']="{}{}{}".format(pre , '-' , cur_post)
				elif '3' not in  locus_strand_index:
					pair_ai='2' if ai == '1' else '1'
					if pair_ai not in locus_strand_index:
						if ai_var['sample'][sample]['PLOIDY']=='1':
							ai_var['sample'][sample]['zygosity']="{}{}".format('hap-',cur_post)
						else:
							# not complete locus?
							# set default to '.'
							print "[Critical Warning] may not complete locus "+rlg[0]['locus']+", set zygosity to '.' as unknown.\n"
							ai_var['sample'][sample]['zygosity']='.'
					else:
						# select pair var
						pair_var=None
						for pair_var_idx in locus_strand_index[pair_ai]:
							# select the first overlapped as pair_var
							if ai_var['begin']<rlg[pair_var_idx]['end'] and ai_var['end'] > rlg[pair_var_idx]['begin']:
								pair_var=rlg[pair_var_idx]
								break
						
						if pair_var is None:
							# not found overlapped pair var
							# assume overlapped with reference
							ai_var['sample'][sample]['zygosity']='hom-ref' if cur_post=='ref' else 'het-ref'
						else:
							
							if pair_var['variantSequence']=='?':
								ai_var['sample'][sample]['zygosity']="{}{}".format('half-' , cur_post)
							else:
								if re.search(r'ref',pair_var['callerVarType']):
									ai_var['sample'][sample]['zygosity']='hom-ref' if cur_post == 'ref' else 'het-ref'
								else:
									ai_var['sample'][sample]['zygosity']='het-alt'
				else:
					# multiple ploidy will always with unknown zygosity
					# except for ai_all and no-call
					ai_var['sample'][sample]['zygosity']='.'
	
	for i in range(len(rlg)-1,-1,-1):
		# delete long ref call for current protocol
		if rlg[i]['referenceSequence']=='=' and rlg[i]['variantSequence']=='=':
			del rlg[i]
	
	return rlg

# tsv file should be carefully sorted
# and formatted. 
# Here is for var-[ASM_ID].tsv.bz2
@LogExc
def read_tsv(fh,sample):
	line=None
	for line in fh:
		if re.search(r'^#',line) or re.search(r'^>',line) or re.search(r'^\s*$',line):
			continue
	if line is None:
		return [0]
	
	# headers
	# 1. locus
	# 2. ploidy
	# 3. alleleIndex (all, 1, 2)
	# 4. chr
	# 5. begin
	# 6. end
	# 7. varType
	# 8. ref
	# 9. Call
	# 10. varscoreVAF
	# 11. varscoreEAF
	# 12. varFilter
	# 13. haplink
	# 14. xRef
	# 15. AF
	# 16. altCall
	
	line = line.strip('\n')
	itms = line.split("\t",16)
	locus_group = list()
	cur_var=tsv_var_parser(itms,sample)
	cur_locus=cur_var['locus']
	if sample in prev_var:
		# not the first var for this sample
		if prev_var[sample]['locus']==itms[0]:
			# same to current reading
			locus_group.append(copy.deepcopy(prev_var[sample]))
		else:
			if prev_var[sample]['referenceSequence']==prev_var[sample]['variantSequence'] and prev_var[sample]['referenceSequence']=='=':
				# current protocol will skip this kind of
				# reference call, due to not reparse reference
				# bases from fasta.
				prev_var[sample]=cur_var
				return read_tsv(fh,sample)
			return_var=copy.deepcopy(prev_var[sample])
			if return_var['sample'][sample]['AI']=='0' or return_var['sample'][sample]['PLOIDY']==1:
				if re.search(r'^no\-',return_var['callerVarType']):
					return_var['sample'][sample]['zygosity']='no-call'
				else:
					pre, post=None,None
					if return_var['sample'][sample]['PLOIDY']==1:
						pre = "hap"
					else:
						pre = "hom"
					if re.search(r'ref',return_var['callerVarType']):
						post = 'ref'
					else:
						post = 'alt'
					
					return_var['sample'][sample]['zygosity']="{}{}{}".format(pre,'-',post)
					return [1, return_var]
			else:
				# not complete locus?
				print "[Critical Warning] not complete locus "+return_var['locus']+" for sample "+sample+". skipping ...\n"
				prev_var[sample]=cur_var
				return read_tsv( fh, sample )
	
	while cur_locus==cur_var['locus']:
		locus_group.append(cur_var)
		newline=None
		try:
			newline =next(fh)
		except StopIteration,e:
			break
		new_itms=line.split('\t',16)
		cur_var=tsv_var_parser(new_itms,sample)
	
	if cur_locus !=cur_var['locus']:
		prev_var[sample]=cur_var
	ref_locus_group=locus_group_parser(locus_group)
	return [1,ref_locus_group]

# tsv file should be carefully sorted
# and formatted. 
# Here is for masterVarBeta-[ASM_ID].tsv.bz2
@LogExc
def read_master(fh, sample):
	line=None
	for line in fh:
		if re.search(r'^#',line) or re.search(r'^>',line) or re.search(r'^\s*$',line):
			continue
	if line is None:
		return [0]
	# headers
	# 1. locus
	# 2. ploidy
	# 3. chromosome
	# 4. begin
	# 5. end
	# 6. zygosity
	# 7. varType
	# 8. reference
	# 9. allele1Seq
	# 10. allele2Seq
	# 11. allele1VarScoreVAF
	# 12. allele2VarScoreVAF
	# 13. allele1VarScoreEAF
	# 14. allele2VarScoreEAF
	# 15. allele1VarFilter
	# 16. allele2VarFilter
	# 17. allele1HapLink
	# 18. allele2HapLink
	# 19. allele1XRef
	# 20. allele2XRef
	# 21. allele1Freq
	# 22. allele2Freq
	# 23. allele1AlternativeCalls
	# 24. allele2AlternativeCalls
	# 25. evidenceIntervalId
	# 26. allele1ReadCount
	# 27. allele2ReadCount
	# 28. referenceAlleleReadCount
	# 29. totalReadCount
	# 30. allele1Gene
	# 31. allele2Gene
	# 32. pfam
	# 33. miRBaseId
	# 34. repeatMasker
	# 35. segDupOverlap
	# 36. relativeCoverageDiploid
	# 37. calledPloidy
	# 38. relativeCoverageNondiploid
	# 39. calledLevel
	# 40. bestLAFsingle
	# 41. lowLAFsingle
	# 42. highLAFsingle
	
	line = line.strip('\n')
	itms = line.split("\t",49)
	
	if itms[7]=='=' and itms[8]=='=' and re.search(r'hom|hap',itms[5]):
		# skip this kind of reference call
		# due to not reparse reference bases from fasta
		return read_master(fh,sample)
	# uniform record
	itms[2]=re.sub(r'^chr','',flags=re.IGNORECASE)
	ploidy=itms[1]
	vars = list()
	
	if itms[8] != '=':
		a1_info = dict()
		a1_info.update(dict(zip("LocID VAF EAF AltCall refAD totalAD".split(),[itms[index] for index in [0, 10, 12, 22, 27, 28]])))
		
		var1={ 
			'chr'              :itms[2],
			'begin'            :itms[3],
			'end'              :itms[4],
			'referenceSequence':itms[7],
			'variantSequence'  :itms[8]
		}
		
		ai_1=0 if itms[5] == 'no-call' or ploidy == '1' or itms[8] == itms[9] else 1
		pb_1=itms[16] if re.search(r'^\d+$',itms[16]) else '.'
		ad_1=itms[25] if re.search(r'^\d+$',itms[25]) else 0
		ar_1=None
		if re.search(r'^\d+$',a1_info['totalAD']) and int(a1_info['totalAD'])>0:
			ar_1 ='%.3f' % float(ad_1) / a1_info['totalAD']
		else:
			ar_1 ='.'
		vf_1=itms[14] if itms[14] != "" else '.'
		sample_h1={
			'PB' : pb_1,
			'AD' : ad_1,
			'AR' : ar_1,
			'AI' : ai_1,
			'VF' : vf_1,
			'PLOIDY' : ploidy
		}
		
		zy_1 = itms[5]
		
		if itms[8]=='?' or itms[8] == 'N':
			zy_1 = 'no-call'
		elif re.search(r'hap|hom|half',zy_1):
			tmp=None
			if itms[7] == itms[8]:
				tmp='-ref'
			else:
				tmp='-alt'
			zy_1="{}{}".format(zy_1,tmp)
		elif re.search(r'het-ref',zy_1) and itms[7] == itms[8]:
			zy_1 = 'het-alt'
		
		sample_h1['zygosity']=zy_1
		sample_h1=give_filter(sample_h1)
		var1['sample'][sample]=sample_h1
		vars.append(var1)
	
	if ploidy == '2' and itms[9]!='=' and not re.search(r'no\-call|hom|hap',itms[5]):
		a2_info = dict()
		a2_info.update(dict(zip("LocID VAF EAF AltCall refAD totalAD".split(),[itms[i] for i in [0, 11, 13, 23, 27, 28]])))
		
		var2 = {
			'chr'               : itms[2],
			'begin'             : itms[3],
			'end'               : itms[4],
			'referenceSequence' : itms[7],
			'variantSequence'   : itms[9]
		}
		
		ai_2 = 2
		pb_2 = itms[17] if re.search(r'^\d+$',itms[17]) else '.'
		ad_2 = itms[26] if re.search(r'^\d+$',itms[26]) else 0
		ar_2 =None
		if re.search(r'^\d+$',a2_info['totalAD']) and a2_info['totalAD']>0:
			ar_2='%.3f' % float(ad_2)/a2_info['totalAD']
		else:
			ar_2='.'
		vf_2 = itms[15] if itms[15] != "" else '.'
		
		sample_h2 = {
			'PB' : pb_2,
			'AD' : ad_2,
			'AR' : ar_2,
			'AI' : ai_2,
			'VF' : vf_2,
			'PLOIDY' : ploidy
		}
		
		zy_2=None
		if itms[9]=='?' or itms[9] == 'N':
			zy_2 = 'no-call'
		elif itms[7] == itms[8]:
			zy_2 = 'het-ref'
		else: # half call case will empty the allele2
			zy_2 = 'het-alt'
		
		sample_h2['zygosity']=zy_2
		sample_h2=give_filter(sample_h2)
		var2['sample'][sample]=sample_h2
		vars.append(var2)
	
	return [1, vars]

# generate var entry for homo ref variant
# if all samples specified in sample_ids
# have been annotated on the homoRef_var
# then return undef, otherwise generate 
# a var entry ready to push to be annotated
@LogExc
def gen_homoRef(homoRef_var, sample_ids):
	
	return_homoRef=copy.deepcopy(homoRef_var)
	if 'sample' in return_homoRef:
		del return_homoRef['sample']
	
	for sample in sample_ids:
		if 'sample' in homoRef_var and sample in homoRef_var['sample']:
			continue
		homoRef_var['sample'][sample]=1
		
		preAD=None
		if sample in pre_total_AD:
			preAD=pre_total_AD[sample]
		if not preAD:
			preAD=20
		return_homoRef['sample'][sample]={
		'zygosity'  : 'hom-ref',
		'NB'        : '.',
		'PB'        : '.',
		'AD'        : preAD,
		'AR'        : '1.00',
		'AI'        : 0,
		'PLOIDY'    : 2,
		'VF'        : '.',
		'PLtag'     : 1,
		'filterTag' : 'PASS',
		'AutoCurated' : 1,
		}
	
	if 'sample' in return_homoRef:
		return return_homoRef
	else:
		return None
@LogExc
def fix_no_call_to_ref(rReady_vars):
	ref_var_idx = -1
	no_call_var_idx = -1
	for i in range(len(rReady_vars)):
		cur_var=rReady_vars[i]
		if cur_var['variantSequence']=='?':
			no_call_var_idx = i
		if cur_var['referenceSequence']==cur_var['variantSequence']:
			ref_var_idx=i
	if no_call_var_idx < 0:
		return rReady_vars
	elif ref_var_idx < 0: # no ref-call then fix no-call var
		no_call_var=rReady_vars[no_call_var_idx]
		no_call_var['variantSequence']=no_call_var['referenceSequence']
		for smp in no_call_var['sample'].keys():
			no_call_var['sample'][smp]['zygosity']='hom-ref'
			preAD=None
			if smp in pre_total_AD:
				preAD=pre_total_AD[smp]
			if not preAD:
				preAD=20
			no_call_var['sample'][smp]['AD']=preAD
			no_call_var['sample'][smp]['AR']='1.00'
			no_call_var['sample'][smp]['PLtag']=1
			no_call_var['sample'][smp]['filterTag']='PASS'
			no_call_var['sample'][smp]['AutoCurated']=1
		return rReady_vars
	else: # with no-call ref-all, correct ref-all and delete-no-call
		no_call_var=rReady_vars[no_call_var_idx]
		ref_var=rReady_vars[ref_var_idx]
		for smp in no_call_var['sample'].keys():
			ref_var['sample'][smp]=no_call_var['sample'][smp]
			ref_var['sample'][smp]['zygosity']='hom-ref'
			preAD=None
			if smp in pre_total_AD:
				preAD=pre_total_AD[smp]
			if not preAD:
				preAD=20
			ref_var['sample'][smp]['AD']=preAD
			ref_var['sample'][smp]['AR']='1.00'
			ref_var['sample'][smp]['PLtag']=1
			ref_var['sample'][smp]['filterTag']='PASS'
			ref_var['sample'][smp]['AutoCurated']=1
		del rReady_vars[no_call_var_idx]
		return rReady_vars

# check to add homo Ref vars to variant list
@LogExc
def check_to_add_homoRef(ready_anno_vars, cur_file, cur_samples):
	test_var=ready_anno_vars[0]
	chr=test_var['chr']
	if chr in HomoRef_Var and (chr not in HRVanno_opt or cur_file not in HRVanno_opt[chr]):
		start=test_var['begin']
		if type !='tsv':
			start -= 1
		ref=test_var['referenceSequence']
		end=start+len(ref)
		
		add_all_opt = 1
		add_before = list()
		for var_in_homoR in HomoRef_Var[chr]:
			if var_in_homoR['begin']>=end:
				add_all_opt = 0
				break
			elif 'varfile' not in var_in_homoR or cur_file not in var_in_homoR['varfile']:
				if var_in_homoR['end']<=start:
					new_gen =gen_homoRef( var_in_homoR, cur_samples )
					if new_gen:
						add_before.apend(new_gen)
					var_in_homoR['varfile'][cur_file]=1
				elif var_in_homoR['begin']==start and var_in_homoR['end']==end:
					# convert no-call to ref for identical match
					ready_anno_vars=fix_no_call_to_ref(ready_anno_vars)
					var_in_homoR['varfile'][cur_file]=1
				elif var_in_homoR['begin']<end:
					# skip overlapped case
					var_in_homoR['varfile'][cur_file]=1
		
		ready_anno_vars=add_before+ready_anno_vars
		
		if add_all_opt == 1:
			HRVanno_opt[chr][cur_file]=1
	return ready_anno_vars

@LogExc
def convert_null_var(rquery_var):
	anno_error_hash=dict()
	anno_error_hash['var']=rquery_var
	
	anno_error_hash['var']["ref"]=rquery_var["referenceSequence"]
	anno_error_hash['var']["alt"]=rquery_var["variantSequence"]
	del anno_error_hash['var']["referenceSequence"]
	del anno_error_hash['var']["variantSequence"]
	
	anno_error_hash['var']["guess"]="annotation-error"
	
	if "end" in rquery_var:
		anno_error_hash['var']["pos"]=rquery_var['begin']
	else:
		anno_error_hash['var']["pos"]=rquery_var['begin']-1
		anno_error_hash['var']["end"]=anno_error_hash['var']["pos"]+len(anno_error_hash['var']["ref"])
	del anno_error_hash['var']["begin"]
	return anno_error_hash
@LogExc
def get_flanks(fa_fai,qvar):
	chr = qvar["chr"]
	if not re.search(r'^chr',chr):
		chr = "{}{}".format('chr',chr)
	if chr=='chrMT':
		chr = 'chrM_NC_012920.1'
	left_flank_region=[chr,qvar["start"]-config["FLANK_LEN"],qvar["start"]]
	right_flank_region=[chr,qvar["end"],qvar["end"]+config["FLANK_LEN"]]
	lr_flanks = list()
	for rgn in [left_flank_region, right_flank_region]:
		cut_off_seq=fa_fai.fetch(rgn[0],rgn[1],rgn[2])
		if  cut_off_seq is None or cut_off_seq == "":
			if  quiet is None:
				print "Warning: ["+rgn+"] not available."
				lr_flanks.append("[NA]")
		else:
			lr_flanks.append(cut_off_seq)
	flankseq =".".join(lr_flanks)
	return flankseq
@LogExc
def get_ref(fai_h, chr, start, end):
	if start >= end:
		return ""
	if MAX_REF_LEN > 0 and end - start > MAX_REF_LEN:
		return "="
	if not re.search(r'^chr',chr):
		chr="{}{}".format('chr',chr)
	if re.search(r'^chrM',chr,flags=re.IGNORECASE):
		chr = 'chrM_NC_012920.1'
	subseq =fai_h.fetch(chr,start,end) 
	if subseq is None or subseq=="":
		if quiet is None:
			sys.stderr.write("[Warning] Fail to get reference seq for "+str(region)+"\n")
			return "="
	return subseq.upper()

# check necessary dependencies
@LogExc
def check_dep(*args):
	for i in args:
		if  not os.access(i, os.F_OK) or not os.access(i, os.R_OK) or (os.path.isfile(i) and os.path.getsize(i)==0):
			raise Exception("Error: ["+str(i)+"] no exists or can not be read or empty.")
	return 1
# read gene symbol's test code
@LogExc
def read_TestCode(file):
	logger=multiprocessing.get_logger()
	TCD_h=None
	if re.search(r'\.gz$',file):
		try:
			TCD_h = gzip.open(file, 'rb')
		except Exception, e:
			raise Error("Error: [" + file + "] GunzipError+\n"+e)
	else:
		TCD_h = open(file, 'rb')
	geneTCD = dict()
	for line in TCD_h:
		if re.search(r'^#|^\s*$',line):
			continue
		line=re.sub(r'\s+$','',line)
		itm=re.split(r'\s+',line,1)
		if 1==len(itm):
			continue
		geneTCD[itm[0]]=itm[1]
	TCD_h.close()
	return geneTCD

# print the header line, assign out_header_keys
# and output the format of variation sheet if needed

# print the header line to the Output
@LogExc
def print_header():
	global HDR
	global OutFp
	global out_header_keys
	#out_header_keys = list() # clean output header keys
	TGP_pop    = config['TGP_POP']
	ExAC_pop   = config['ExAC_POP']
	GAD_pop = 'EAS' #current api only support EAS population
	
	# =========== Location group ========= #
	chrpos_header = [ "Chr", "Start", "Stop", "MapLoc" ]
	inExcel = "InExcel"    # firm specified.
	sample_header=list()
	if sample_list is not None:
		sample_header = [ "FamID", "SampleID", "Sex" ]
	else:
		sample_header = ["SampleID"]
	
	# ========== Genotyping Group ========= #
	geno_header = list()
	geno_header+=[
		"NbGID",     "Ref",      "VarType",   "Call",
		"Flank",     "Zygosity", "A.Depth",   "A.Ratio",
		"PhasedGID", "A.Index",  "RepeatTag", "Filter" 
	]
	out_header_keys+=chrpos_header+[inExcel]+sample_header+geno_header
	
	# ========== Gene Information ========= #
	gene_header = [
		"MIM Gene ID",
		"MIM Stat",
		"MIM Pheno IDs",
		"MIM Inheritance",
		"CGD Condition",
		"CGD Inheritance",
		"CGD Allelic Condition",
		"CGD Manifestation Categories",
		"CGD Intervention Categories",
		"CGD References",
		"TestCode",
		"EntrezGeneID",
		"Gene Symbol"
	]
	
	out_header_keys+='''mimGID mimStat mimPIDs mimInhs
					cgdCond cgdInhs cgdACond cgdManCat cgdIntCat cgdRef
					TestCode EntrezGeneID geneSym'''.split()
	
	# ======== Function Information ======== #
	func_header =[ "FuncRegion", "ExIn_ID", "pfamName", "pfamId", "Function" ]
	
	# ========== HGVS Information ========== #
	hgvs_header = ["Transcript", "Protein", "Strand",      "Primary", "cHGVS",      "pHGVS",   "CodonChange", "PolarChange","MutationName"]
	
	out_header_keys+=func_header+hgvs_header
	
	# ========= Public Frequency DB ======== #
	freq_header = list()
	
	# dbSNP
	freq_header+=["dbSNP Note","dbSNP Allele Count","dbSNP Allele Freq","dbSNP Mul Tag"]
	
	# 1000 genomes
	freq_header+=["1000G "+TGP_pop+" AF", "1000G AF"]
	
	# ESP6500
	freq_header+=["ESP6500 AC", "ESP6500 AF"]
	
	# ExAC
	freq_header+=[
		"ExAC Filter",
		"ExAC "+ExAC_pop+" HomoAlt Count",
		"ExAC "+ExAC_pop+" AC",
		"ExAC "+ExAC_pop+" AF",
		"ExAC HomoAlt Count",
		"ExAC AC",
		"ExAC AF" 
	]
	
	# GAD
	freq_header+=[
	"GAD "+GAD_pop+" HomoAlt Count",
	"GAD "+GAD_pop+" AC",
	"GAD "+GAD_pop+" AF",
	"GAD HomoAlt Count",
	"GAD AC",
	"GAD AF"
	]
	
	# rsID from dbSNP
	freq_header+=["rsID"]
	
	out_header_keys+='''dbsnpNote dbsnpAC dbsnpAF dbsnpMTag
						TGPpopAF TGP_AF ESP6500AC ESP6500AF
						ExAC_Filter ExAC_pop_HASC ExAC_pop_AC
						ExAC_pop_AF ExAC_HASC ExAC_AC ExAC_AF
						GAD_pop_HASC GAD_pop_AC
						GAD_pop_AF GAD_HASC GAD_AC GAD_AF
						rsID'''.split()
	
	# ========= Panel Frequency DB ========= #
	panel_header = [ "Panel ID", "Panel AlleleCount", "Panel AlleleFreq" ]
	
	out_header_keys+="PanelID PanelAC PanelAF".split()
	'''
	gwas_header=["GWAS Risk Allele", "GWAS pmID", "GWAS Trait"]
	phyloP_header =["PhyloP Primates", "PhyloP Vertebrates", "PhyloP Placental Mammals"]
	out_header_keys+='gwasRiskA gwasPMID gwasTrait phyloPpr phyloPve phyloPpm'.split()
	'''
	# ====== Mutation Effect Prediction ==== #
	
	# PhyloP scores
	phyloP_header=["PhyloP Primates", "PhyloP Vertebrates", "PhyloP Placental Mammals"]
	
	out_header_keys+="phyloPpr phyloPve phyloPpm".split()
	
	# Condel prediction
	condel_header =[
		"Ens SIFT Score",
		"Ens Polyphen2HumVar Score",
		"Ens Polyphen2HumDiv Score",
		"Ens Condel Score",
		"Ens Condel Pred"
	]
	
	# Cosmic prediction
	cosmic_header =[
		"Cosmic MutName",
		"Cosmic Site",
		"Cosmic Histology",
		"Cosmic MutStat",
		"Cosmic ID"
	]
	
	out_header_keys+='ensSIFT ensPP2hvar ensPP2hdiv ensCondel ensCondelPred cosmicName cosmicSite cosmicHis cosmicStat cosmicID'.split()
	
	# HGMD prediction
	hgmd_header = [ "HGMD ID" ]
	out_header_keys+=["hgmdID"]
	
	hgmd_header+=[
		"HGMD MutName",
		"HGMD Disease",
		"HGMD pmID",
		"HGMD Pred"
	]
	
	out_header_keys+="hgmdName hgmdDis hgmdPMID hgmdPred".split()
	
	# ClinVar prediction
	clinVar_header = ["ClinVar SuspectReason", "ClinVar Accession", "ClinVar RevStat", "ClinVar Significant"]
	out_header_keys+=["ClinSSR", "ClinACC", "ClinRevStat", "ClinSignificant"]
	
	# Auto Interpretation
	interpret_header = ["AutoInterpStatus"]
	out_header_keys+=["autoInterp"]
	
	all_headers=chrpos_header+ [inExcel]+ sample_header+ geno_header+   gene_header+       func_header+ hgvs_header+   freq_header+       panel_header+  phyloP_header+ condel_header+     cosmic_header+ hgmd_header+       clinVar_header+ interpret_header
	tmp="\n".join( [
	"## NCanno Version    : "+VERSION,
	"## BedAnno Version   : "+BedAnno.BedAnno.VERSION,
	"## Configure File    : " + os.path.abspath(CONFIG_FILE),
	"## Input File Format : "+ ( "vcf" if type !='tsv' else "tsv" ),
	"## Format Options    : "
	+( 'offline, ' if  offline is not None else 'online, ' )
	+ ( 'famlily' if sample_list is not None else 'single' ) + " MODE",
	"## DB Version File   : " + os.path.abspath(verList),
	"##" ])
	print >>OutFp, tmp
	print >>OutFp, "#"+("\t".join(all_headers))
	
	if headrule:
		print >>HDR, "#"+("\t".join("title width hidden level collapsed".split()))+"\n"
		for title in all_headers:
			if title=="InExcel" or title=="FamID":
				continue
			if sample_list is None and title=="SampleID":
				continue
			
			# assign default stat
			width, hidden, level, collapsed=10, 1, 1, 0
			
			if (   
			title == "#Chr" or \
			title == "Ref" or \
			title == "Call" or \
			title == "VarType" or \
			title == "A.Depth" or \
			title == "A.Ratio" or \
			title == "A.Index" or \
			title == "Strand" ):
				width = 6
			elif \
			title=="Sex" or \
			title=="Primary":
				width = 6
				hidden, level, collapsed= 0, 0, 1 
			elif \
			title == "Flank" or \
			title == "RepeatTag" or \
			title == "TestCode" or \
			title == "FuncRegion" or \
			title == "ExIn_ID" or \
			title == "pfamId" or \
			title == "EntrezGeneID":
				width = 8
			elif \
			title == "Filter" or \
			title == "Zygosity" or \
			title == "MapLoc":
				width = 8
				hidden, level, collapsed  =  0, 0, 1 
			elif \
			title == "MIM Inheritance" or \
			title == "Function" or \
			title == "rsID" or \
			title == "PVFD AF" or \
			title == "Panel AlleleFreq" or \
			title == "Panel Pred" or \
			title == "AutoInterpStatus" :
				hidden, level, collapsed  =  0, 0, 1 
			elif \
			title == "Transcript" or \
			title == "Protein" or \
			title == "pHGVS":
				width = 12
			elif title == "Gene Symbol":
				width = 12
				hidden, level, collapsed  =  0, 0, 1 
			elif re.search(r'^HGMD |^ClinVar ',title):
				width = 12
				if title == "HGMD Pred" or \
				title == "ClinVar Significant":
					hidden, level, collapsed  =  0, 0, 1 
			elif re.search(r'^dbSNP |^1000G |^ESP6500 |^ExAC |^GAD ',title) or re.search(r'^PhyloP |^Ens |^Cosmic ',title):
				width = 15
				if \
				title == "dbSNP Mul Tag" or \
				title == "1000G AF" or \
				title == "ESP6500 AF" or \
				title == "ExAC AF" or \
				title == "GAD AF" or \
				title == "PhyloP Placental Mammals" or \
				title == "Ens Condel Pred" or \
				title == "Cosmic ID":
					collapsed = 1
				else:
					level = 2
			elif title == "MutationName":
				width = 30
				hidden, level, collapsed  =  0, 0, 1 
			print >>HDR,"\t".join('"'+title+'"',width, hidden, level, collapsed)+"\n"
		print >>HDR,"\t".join('"Interpretation"', 15, 0, 0, 1)
		HDR.close()
@LogExc
def set_genome(self,genome_rz):
	if not os.access(genome_rz,os.F_OK) or not os.access(genome_rz,os.R_OK):
		raise Exception("Error: cannot read "+str(genome_rz)+".")
	setattr(self, "genome" , genome_rz)
	setattr(self, "genome_h" ,pysam.Fastafile(genome_rz))
@LogExc
def set_cytoBand(self, cytodb):
	setattr(self, "cytoBand" , cytodb)
	cytoBand_h = GetCytoBand.GetCytoBand({"db": cytodb})
	setattr(self, "cytoBand_h", cytoBand_h)
@LogExc
def set_pfam(self, pfamdb):
	setattr(self, "pfam", pfamdb)
	pfam_h = GetPfam.GetPfam({"db" : pfamdb})
	setattr(self, "pfam_h", pfam_h)
@LogExc
def set_prediction(self, predictiondb):
		logger=multiprocessing.get_logger()
		setattr(self, "prediction", predictiondb)
		common_opts = dict()
		if hasattr(self, "quiet"):
			common_opts['quiet'] = 1
		common_opts.update({"db":predictiondb})
		prediction_h = GetPrediction.GetPrediction(common_opts)
		setattr(self, "prediction_h", prediction_h)
@LogExc
def set_condel(self, condelConfig):
	setattr(self, "condel", condelConfig)
	condel_h = CondelPred.CondelPred(condelConfig)
	setattr(self, "condel_h", condel_h)
@LogExc
def set_phyloP(self, phyloPdb):
	setattr(self, "phyloP", phyloPdb)
	phyloP_h = GetPhyloP46wayScore.GetPhyloP46wayScore({"db" : phyloPdb})
	setattr(self, "phyloP_h", phyloP_h)
@LogExc
def set_cosmic(self, cosmic_db):
	setattr(self, "cosmic", cosmic_db)
	common_opts = dict()
	if hasattr(self, 'quiet'):
		common_opts['quiet'] = 1
	common_opts.update({"db":cosmic_db})
	cosmic_h = GetCOSMIC.GetCOSMIC(common_opts)
	setattr(self, "cosmic_h", cosmic_h)
@LogExc
def set_dbSNP(self, dbSNPdb):
	setattr(self, "dbSNP", dbSNPdb)
	dbSNP_h = GetDBSNP.GetDBSNP({"db" : dbSNPdb})
	setattr(self, "dbSNP_h", dbSNP_h)
@LogExc
def set_tgp(self, tgpdb):
	setattr(self, "tgp", tgpdb)
	BedAnno.load_opt_vcfaf = 1
	common_opts = dict()
	if hasattr(self, "quiet"):
		common_opts['quiet'] = 1
	common_opts.update({"db":tgpdb})
	tgp_h = GetVcfAF.GetVcfAF(common_opts)
	setattr(self, "tgp_h", tgp_h)
@LogExc
def set_esp6500(self, esp6500db):
	setattr(self, "esp6500", esp6500db)
	load_opt_vcfaf = 1
	common_opts = dict()
	if hasattr(self,'quiet'):
		common_opts['quiet'] = 1
	common_opts.update({"db":esp6500db})
	esp6500_h = GetVcfAF.GetVcfAF(common_opts)
	setattr(self, "esp6500_h", esp6500_h)
@LogExc
def set_exac(self, exac_db):
	setattr(self, "exac", exac_db)
	common_opts = dict()
	if hasattr(self, 'quiet'):
		common_opts['quiet']=1
	common_opts.update({"db":exac_db})
	exac_h = GetExAC.GetExAC(common_opts)
	setattr(self, "exac_h", exac_h)
@LogExc
def set_gnomAD(self, gnomAD_db):
	setattr(self, "gnomAD", gnomAD_db)
	common_opts = dict()
	if hasattr(self, "quiet"):
		common_opts["quiet"] = 1
	common_opts.update({"db":gnomAD_db})
	gnomAD_h = GetGAD.GetGAD(common_opts)
	setattr(self, "gnomAD_h", gnomAD_h)
@LogExc
def set_customdb(self, cusdb, dbID):
	load_opt_vcfaf = 1
	common_opts = dict()
	if hasattr(self, "quiet"):
		common_opts['quiet'] = 1
	common_opts.update({"db":cusdb})
	cusdb_h = GetVcfAF.GetVcfAF(common_opts)
	setattr(self, "cusdb_" + str(dbID) + "_h", cusdb_h)
@LogExc
def set_cg54(self, cg54db):
	setattr(self, "cg54", cg54db)
	common_opts = dict()
	if hasattr(self, "quiet"):
		common_opts['quiet'] = 1
	common_opts.update({"db":cg54db})
	cg54_h = GetCGpub.GetCGpub(common_opts)
	setattr(self, "cg54_h", cg54_h)
@LogExc
def set_wellderly(self, wellderlydb):
	setattr(self, "wellderly", wellderlydb)
	common_opts = dict()
	if hasattr(self, "quiet"):
		common_opts['quiet'] = 1
	common_opts.update({"db":wellderlydb})
	wellderly_h = GetCGpub.GetCGpub(common_opts)
	setattr(self, "wellderly_h", wellderly_h)
@LogExc
def set_rmsk(self, rmskdb):
	if rmskdb:
		setattr(self, "rmsk", rmskdb)
	rmsk_h = GetRepeatTag.GetRepeatTag({"db" : self.rmsk})
	setattr(self, "rmsk_h", rmsk_h)
@LogExc
def set_gwas(self, gwasdb):
	if gwasdb:
		setattr(self, "gwas", gwasdb)
	gwas_h = GetGWAS.GetGWAS({"db" : self.gwas})
	setattr(self, "gwas_h", gwas_h)
@LogExc
def var_NClib_anno(BedAnno_ins, var):
	if hasattr(BedAnno_ins, "cytoBand"):
		var.cytoBand = BedAnno_ins.cytoBand_h.getCB(var.chr, var.pos, var.end)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "rmsk"):
		var.reptag = BedAnno_ins.rmsk_h.getRepTag(var.chr, var.pos, var.end)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "gwas"):
		var.gwas = BedAnno_ins.gwas_h.getGWAS(var.chr, var.pos, var.end)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "phyloP"):
		if var.sm == 1:
			var.phyloPpm, var.phyloPpr, var.phyloPve = BedAnno_ins.phyloP_h.getPhyloP46wayScore(var.chr, (
			int(var.pos) + 1))  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "dbSNP"):
		
		if hasattr(var, "sep_snvs"):
			
			new_sqls = list()
			cur_start = var.sep_snvs[0]
			for i in range(len(var.sep_snvs)):
				if i == (len(var.sep_snvs) - 1) or var.sep_snvs[i + 1] - cur_start > 1:
					new_ref = var.ref[(cur_start - int(var.pos) - 1):(int(var.sep_snvs[i]) - var.pos - 1)]
					new_alt = var.alt[(cur_start - int(var.pos) - 1):(int(var.sep_snvs[i]) - var.pos - 1)]
					new_sqls.append([cur_start, var.sep_snvs[i], new_ref, new_alt])
					if i < len(var.sep_snvs) - 1:
						cur_start = var.sep_snvs[i + 1]
			for rSE in new_sqls:
				rOneSql = BedAnno_ins.dbSNP_h.getRS(var.chr, *rSE)  # 注意调用函数的实现
				for k in sorted(rOneSql.keys()):
					var.dbsnp[k] = rOneSql[k]
		else:
			var.dbsnp = BedAnno_ins.dbSNP_h.getRS(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "tgp"):
		var.tgp = BedAnno_ins.tgp_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "esp6500"):
		var.esp6500 = BedAnno_ins.esp6500_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "exac"):
		var.exac = BedAnno_ins.exac_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "gnomAD"):
		var.gnomAD = BedAnno_ins.gnomAD_h.getGAD(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	
	for dbhk in sorted(dir(BedAnno_ins)):
		m = re.match(r'cusdb_(\S+)_h', dbhk)
		if m and hasattr(BedAnno_ins, dbhk):
			dbID = m.group(1)
			setattr(var, "cusdb_" + dbID, BedAnno_ins.dbhk.getAF(var.chr, var.pos, var.end, var.ref, var.alt))
			# eval("var."+"cusdb_"+dbID)=BedAnno_ins.dbhk.getAF(var.chr,var.pos,var.end,var.ref,var.alt)#注意调用函数的实现
	if hasattr(BedAnno_ins, "cg54"):
		var.cg54 = BedAnno_ins.cg54_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "wellderly"):
		var.wellderly = BedAnno_ins.wellderly_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	if hasattr(BedAnno_ins, "cosmic"):
		var.cosmic = BedAnno_ins.cosmic_h.getCOSMIC(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现
	return var
@LogExc
def anno_threads(rquery_var):
	logger=multiprocessing.get_logger()
	common_opts = dict()
	
	global quiet
	
	if quiet is not None:
		common_opts['quiet']=1
	
	beda, panel_h, omim_h, hgmd_h, clinVar_h, vusVar_h, cgd_h, hg19_fai=[None]*8
	
	t0, t1 =None, None
	thread_timer = dict()
	
	# involke anno engine
	global hg19_fa
	global db
	global trDB
	
	hg19_fai=pysam.Fastafile(hg19_fa)
	
	arg_dict={}
	arg_dict.update({"db":db,"tr":trDB})
	global opts
	arg_dict.update(opts)
	arg_dict.update(common_opts)
	beda=BedAnno.BedAnno(arg_dict)
	
	if 'genome' in arg_dict:
		set_genome(beda,arg_dict['genome'])
	if 'cytoBand' in arg_dict:
		set_cytoBand(beda,arg_dict['cytoBand'])
	if 'pfam' in arg_dict:
		set_pfam(beda,arg_dict['pfam'])
	if 'prediction' in arg_dict:
		set_prediction(beda,arg_dict['prediction'])
	if 'condel' in arg_dict:
		set_condel(beda,arg_dict['condel'])
	if 'phyloP' in arg_dict:
		set_phyloP(beda,arg_dict['phyloP'])
	if 'cosmic' in arg_dict:
		set_cosmic(beda,arg_dict['cosmic'])
	if 'dbSNP' in arg_dict:
		set_dbSNP(beda,arg_dict['dbSNP'])
	if 'tgp' in arg_dict:
		set_tgp(beda,arg_dict['tgp'])
	if 'esp6500' in arg_dict:
		set_esp6500(beda,arg_dict['esp6500'])
	if 'exac' in arg_dict:
		set_exac(beda,arg_dict['exac'])
	if 'gnomAD' in arg_dict:
		set_gnomAD(beda,arg_dict['gnomAD'])
	for dbk in sorted(dir(beda)):
		m = re.match(r'customdb_(\S+)', dbk)
		if m:
			dbID = m.group(1)
			set_customdb(self, eval("self." + dbk), dbID)
	if 'cg54' in arg_dict:
		set_cg54(beda,arg_dict['cg54'])
	if 'wellderly' in arg_dict:
		set_wellderly(beda,arg_dict['wellderly'])
	if 'rmsk' in arg_dict:
		set_rmsk(beda,arg_dict['rmsk'])
	if 'gwas' in arg_dict:
		set_gwas(beda,arg_dict['gwas'])
	
	
	
	
	
	
	
	
	
	
	
	global config
	GetVcfAF_common_opts=copy.deepcopy(common_opts)
	
	GetVcfAF_common_opts.update({"db":config['Panel_Control']})
	if config['Panel_Control']!=".":
		panel_h = GetVcfAF.GetVcfAF(GetVcfAF_common_opts)
	
	global omim_data
	global hgmd
	global clinVar
	global vusVar
	global cgd_data
	global verList
	global DB_VERSION
	omim_h=GetOmim.GetOmim({"acdb" : omim_data})
	
	GetHGMD_common_opts=copy.deepcopy(common_opts)
	GetHGMD_common_opts.update({"db":hgmd})
	hgmd_h=GetHGMD.GetHGMD(GetHGMD_common_opts)
	clinVar_h = GetClinVar.GetClinVar( {"db" :clinVar })
	
	vusVar_h = GetClinVar.GetClinVar({"db" : vusVar})
	cgd_h = GetCGDinfo.GetCGDinfo({"db":cgd_data, "quiet" : 1 })
	
	global write_version
	if write_version is not None and 'Process-1' ==multiprocessing.current_process().name:
		try:
			VERLIST=open(verList,'w')
		except Exception,e:
			raise Exception("Error: ["+verList+"] "+e)
		for db in sorted(DB_VERSION.keys()):
			print >>VERLIST,"%20s : %s\n" % (db, DB_VERSION[db])
		VERLIST.close()
	#cur_thread_varcount = 0
	#while True:
		#rquery_var=anno_q.get()
	anno=None
	reformat_hash=dict()
	#global debug
	if debug :
		sys.stderr.write("{}{}{}{}".format("[Info from Thread ",multiprocessing.current_process().pid,"] current var: ",rquery_var))
	try:
		
		if str(rquery_var['referenceSequence'])=='=':
			rquery_var['referenceSequence']=get_ref(hg19_fai,*[rquery_var[ii] for ii in "chr begin end".split()])
		
		#anno=beda.anno(rquery_var)
		
		var_tmp = BedAnnoVar.BedAnnoVar(rquery_var)
		#if "cytoBand" in beda:
		#	anno.var.cytoBand=beda.cytoBand_h.getCB(var_tmp.chr, var_tmp.pos, var_tmp.end)
		#logger.info("1938 anno=")
		#logger.info(anno.var.__dict__)
		
		var_tmp=var_NClib_anno(beda, var_tmp)
		anno=beda.anno(var=var_tmp)
		
		
		
	except Exception,e:
		sys.stderr.write( "{}{}{}{}{}{}{}{}{}".format("[Annotation Failed] ","[Thread ",multiprocessing.current_process().pid," Var ",rquery_var['chr']+",",str(rquery_var['begin'])+",",rquery_var['referenceSequence']+",",rquery_var['variantSequence'],e))
		reformat_hash = convert_null_var(rquery_var)
	else:
		
		# remove object properties for shared_clone
		reformat_hash['var']=BedAnnoVar.BedAnnoVar.TO_JSON(anno.var)
		query_var_hash={
			"chr":reformat_hash['var']['chr'],
			"start":reformat_hash['var']['pos'],
			"end":reformat_hash['var']['end'],
			"ref":reformat_hash['var']['ref'],
			"alt":reformat_hash['var']['alt']
		}
		
		if panel_h is not None and reformat_hash['var']['guess']!='no-call':
			try:
				reformat_hash['var']['panelctrl']=panel_h.getAF(*[query_var_hash[ii] for ii in "chr start end ref alt".split()])
			except Exception,e:
				sys.stderr.write( "{}{}{}{}{}{}{}{}{}{}".format("[PanelSQL failed]","[Thread" ,multiprocessing.current_process().pid," Var ",query_var_hash['chr']+",",str(query_var_hash['start'])+",",str(query_var_hash['end'])+",",query_var_hash['ref']+",",query_var_hash['alt'],e))
				reformat_hash['var']['PanelSQL_Fail']=1
		
		try:
			clinVar_hit=clinVar_h.getCL(*[query_var_hash[ii] for ii in "chr start end ref alt".split()])
			if clinVar_hit: # hit record in clinVar db
				reformat_hash['var']['clinVar']=clinVar_hit
			else:
				
				reformat_hash['var']['clinVar']=vusVar_h.getCL(*[query_var_hash[ii] for ii in "chr start end ref alt".split()])
		except Exception,e:
			sys.stderr.write( "{}{}{}{}{}{}{}{}{}{}".format("[clinVar failed] ","[Thread" ,multiprocessing.current_process().pid," Var ",query_var_hash['chr']+",",query_var_hash['start']+",",query_var_hash['end']+",",query_var_hash['ref']+",",query_var_hash['alt']),e)
			reformat_hash["var"]["clinVar_Fail"] = 1
		
		try:
			reformat_hash['var']['flanks']=get_flanks(hg19_fai, query_var_hash)
			
		except Exception,e:
			sys.stderr.write( "{}{}{}{}{}{}{}{}{}".format("[get_flanks] ","[Thread" ,multiprocessing.current_process().pid," Var ",query_var_hash['chr']+",",query_var_hash['start']+",",query_var_hash['end']+",",query_var_hash['ref']+",",query_var_hash['alt']+" ")+e)
			reformat_hash['var']['flanks']="."
		
		if hasattr(anno,'trInfo'):
			major_var_info=anno.var.varName
			reformat_hash['trInfo']=anno.trInfo
			
			for tr in sorted(reformat_hash['trInfo'].keys()):
				curtr_ent=reformat_hash['trInfo'][tr]
				if re.search(r"^"+re.escape(tr),major_var_info):
					curtr_ent['major']="Y"
				else:
					curtr_ent['major']="N"
				
				try:
					curtr_ent['omim']=omim_h.getAnnoComb(curtr_ent['geneId'])
					
					curtr_ent['cgd']=cgd_h.getCGD(curtr_ent['geneId'])
					'''
					if curtr_ent['geneId'] not in GetCGDinfo_dict:
						curtr_ent['cgd']=dict()
					else:
						curtr_ent['cgd']=GetCGDinfo_dict[curtr_ent['geneId']]
					'''
				except Exception,e:
					sys.stderr.write("{}{}{}{}".format("[GeneSql Fail] ","[Thread ",multiprocessing.current_process().pid," for "+tr+"] ")+e)
					curtr_ent['omim']=dict()
					curtr_ent['cgd']=dict()
					curtr_ent['GeneSQL_Fail']=1
				global GeneTestCode
				if GeneTestCode is not None and curtr_ent["geneSym"] in GeneTestCode:
					curtr_ent['TCD']=GeneTestCode[curtr_ent['geneSym']]
				else:
					curtr_ent['TCD']="."
				
				# There's no need to query tr-extra database for no-call var
				if reformat_hash['var']['guess']=='no-call':
					continue
				
				trQuery=copy.deepcopy(query_var_hash)
				trQuery.update({"transcript":tr,"geneid":curtr_ent["geneId"],"genesym":curtr_ent["geneSym"],"cHGVS":curtr_ent["c"]})
				if "p" in curtr_ent :
					trQuery['pHGVS']=curtr_ent['p']
				if "p3" in curtr_ent :
					trQuery['pHGVS3']=curtr_ent['p3']
				
				query_alt_opt = 0
				altTrQ=copy.deepcopy(trQuery)
				if "alt_cHGVS" in curtr_ent or "alt_pHGVS" in curtr_ent:
					query_alt_opt = 1
					if "alt_cHGVS" in curtr_ent:
						altTrQ['cHGVS']=curtr_ent['alt_cHGVS']
					if 'alt_pHGVS' in curtr_ent:
						altTrQ['pHGVS']=curtr_ent['alt_pHGVS']
					if 'alt_p3' in curtr_ent:
						altTrQ['pHGVS3']=curtr_ent['alt_p3']
				
				query_std_opt = 0
				stdTrQ=copy.deepcopy(trQuery)
				if "standard_cHGVS" in curtr_ent or "standard_pHGVS" in curtr_ent:
					query_std_opt = 1
					if "standard_cHGVS" in curtr_ent:
						stdTrQ['cHGVS']=curtr_ent['standard_cHGVS']
					if "standard_pHGVS" in curtr_ent:
						stdTrQ['pHGVS']=curtr_ent['standard_pHGVS']
					if "standard_p3" in curtr_ent:
						stdTrQ['pHGVS3']=curtr_ent['standard_p3']
				
				if "prRef" in curtr_ent \
				and "prAlt" in curtr_ent \
				and len(curtr_ent['prRef'])==1 \
				and len(curtr_ent['prAlt'])==1:
					trQuery['aaref']=curtr_ent['prRef']
					trQuery['aaalt']=curtr_ent['prAlt']
				
				try:
					curtr_ent['hgmd']=hgmd_h.getHGMD(trQuery)
					if query_std_opt and 0 >= len(curtr_ent['hgmd']):
						curtr_ent['hgmd']=hgmd_h.getHGMD(stdTrQ)
					if query_alt_opt and 0 >= len(curtr_ent['hgmd']):
						curtr_ent['hgmd']=hgmd_h.getHGMD(altTrQ)
				except Exception,e:
					sys.stderr.write("{}{}{}{}{}{}".format("[localhgmd Fail] ","[Thread ",multiprocessing.current_process().pid," "+str(tr)+": ",trQuery['cHGVS']," ("+trQuery["pHGVS"]+")] "))
					curtr_ent['hgmd']=list()
					curtr_ent['LocalHGMD_Fail']=1
	
	#cur_thread_varcount+=1
	
	if debug :
		sys.stderr.write("{}{}{}{}".format("[Info from Thread ",multiprocessing.current_process().pid,"] var after anno:",reformat_hash))
	#print 'reformat_hash=',reformat_hash
	#anno_cache.append(reformat_hash)
	#cached_count +=1
	'''
	if cur_thread_varcount>=BUFFER_SIZE:
		if timer is not None:
			for timekey in sorted(thread_timer.keys()):
				sys.stderr.write("{}{}{}".format("[TIMER ",multiprocessing.current_process().pid,"] "+str(timekey)+" : "+str(thread_timer[timekey])+"\n"))
				thread_timer[timekey]=0
		cur_thread_varcount = 0
	'''
	hg19_fai.close()
	return reformat_hash

# check PL score
# return PLtag:
# -1:   FAIL
#  0:   DUBIOUS
#  1:   PASS
@LogExc
def check_PL(PL, a1, a2):
	PLs = PL
	PLtag = 1
	tmp=sorted([a1, a2])
	PLindex=pIndex_cal(*tmp)
	gtype_PL=PLs.pop(int(PLindex))
	if str(gtype_PL)==".":
		PLtag = 1
	elif int(gtype_PL)>int(PL_THRESHOLD):
		PLtag = -1
	elif int(gtype_PL) > 0:
		PLtag = 0
	else:
		for otherPL in PLs:
			if str(otherPL)==".":
				continue
			elif int(otherPL) == 0:
				PLtag = -1
				break
			elif int(otherPL)<int(PL_UP_THRESHOLD):
				PLtag = 0
	return PLtag


# calculate the index in GL/PL string, when geno pair is from $p1 and $p2 alleles
@LogExc
def pIndex_cal(p1, p2):
	p1=int(p1)
	p2=int(p2)
	return (p2 * (p2 + 1) / 2 + p1)

# vcf reader, input a vcf handler, return a var array.
# IMPRECISE variant will be ignored.
@LogExc
def read_vcf(vcf_h):
	rec=None
	try:
		rec = next(vcf_h)
	except StopIteration:
		return [0,None]
	if "IMPRECISE" in rec.INFO:
		return read_vcf(vcf_h)
	if len(rec.FILTER)==0:
		rec.FILTER=['PASS']
	filter_tag =";".join(rec.FILTER)
	filter_tag =re.sub(r'PASS',".",filter_tag,flags=re.IGNORECASE)
	
	chr, vcf_pos, vcf_ref, vcf_alts = rec.CHROM,rec.POS,rec.REF,rec.ALT
	
	chr =re.sub(r'^chr',"",chr,flags=re.IGNORECASE)
	
	# unshift vcf ref to all case
	all_case=[vcf_ref]+vcf_alts
	nocall_var={
		'chr':chr,
		'begin':int(vcf_pos) - 1,
		'end':(int(vcf_pos)+len(vcf_ref)-1),
		'referenceSequence':vcf_ref,
		'variantSequence':'?'
	}
	
	vars = list()
	for i in range(len(all_case)):
		# ignore the non precise variation for vcf format
		# ignore large sv or breakend for vcf format
		if re.search(r'[<>\[\]]',str(all_case[i])) or (str(all_case[i])!="." and re.search(r'\.',str(all_case[i]))):
			return read_vcf(vcf_h)
		if str(all_case[i])=='.':
			all_case[i]=all_case[0]
		
		# this is just indicate it's from vcf
		# DO NOT add "end" key or you will have to
		# change begin to 0-based position
		var = {
			'chr':chr,
			'begin':vcf_pos,
			'referenceSequence':vcf_ref,
			'variantSequence':all_case[i]
		}
		vars.append(var)
	
	# assign sample info to var
	# "sample" => {
	#   $sample => {
	#       zygosity => $zygosity,
	#       AD => $allelicDepth,
	#       AR => $allelicRatio,
	#       VF => $CallerFilterTag,
	#       AI => $AllelicIndex
	#       PLOIDY => $ploidy,
	#       PLTAG => $judge_reliablity_from_PL,
	#    
	#       # Here is customized keys
	#       PB => $Phased_Group_ID,
	#       NB => $Neighbour_Group_ID,
	#   },
	#   ...
	# }
	
	involve_nocall = 0
	for sample in rec.samples:
		sample_hash=dict(sample.data._asdict())
		sample=sample.sample
		
		#print 'sample=',sample.sample
		#print 'sample_hash=',dict(sample_hash)
		#print sample.data
		#print sample_hash
		#sys.exit(1)
		m=re.search(r'^(\.|\d+)([\\|/]?)(\.?|\d*)$',sample_hash["GT"])
		if not m:
			raise "Could not parse gtype string ["+sample_hash["GT"]+"] ["+chr+":"+vcf_pos+"]"
		a1, sep, a2 = m.group(1),m.group(2),m.group(3)
		Zyg, AD, AR, PLtag, PB, NB=['.'] * 6
		if 'PB' in sample_hash:
			PB = sample_hash['PB']
		if 'NB' in sample_hash:
			NB = sample_hash['NB']
		AI = 0 # default homo/hap
		PLOIDY = None
		if a2 is not None and a2 !='':
			PLOIDY = 2
		else:
			PLOIDY = 1
		# no calling case will non-exists the tag.
		ADs = list()
		total_Dp = 0
		if 'AD' in sample_hash and sample_hash['AD'] is not None and sample_hash['AD'] !='.':
			ADs = sample_hash['AD']
			total_Dp=sum(ADs)
			pre_total_AD[sample]=total_Dp
		
		if 'PL' in sample_hash and a2 is not None and a2!='':
			PLtag = check_PL(sample_hash['PL'], a1, a2)
		if a2 is None or a2=='' or a1 == a2:
			AI = 0
			cur_sample_h=None
			if a1 == '.':
				if 'sample' not in nocall_var:
					nocall_var['sample']=dict()
				if sample not in nocall_var['sample']:
					nocall_var['sample'][sample]=dict()
				cur_sample_h = nocall_var['sample'][sample]
				involve_nocall = 1
				Zyg = 'no-call'
			else:
				if 'sample' not in vars[int(a1)]:
					vars[int(a1)]['sample']=dict()
				vars[int(a1)]['sample'][sample]=dict()
				cur_sample_h = vars[int(a1)]['sample'][sample]
				if a2 is not None and a2 == '0':
					Zyg = 'hom-ref'
				elif a2 is not None and a2!='':
					Zyg = 'hom-alt'
				elif a1 == '0':
					Zyg = 'hap-ref'
				else:
					Zyg = 'hap-alt'
				
				if total_Dp > 0:
					if ADs[int(a1)] is None:
						print "{}{}{}{}{}".format("AD index error! [",chr, vcf_pos, sample,"]")
					else:
						AD = ADs[int(a1)]
						AR = '%.2f' % (AD/total_Dp)
			cur_sample_h['zygosity']= Zyg
			cur_sample_h['NB']      = NB
			cur_sample_h['PB']      = PB
			cur_sample_h['AD']      = AD
			cur_sample_h['AR']      = AR
			cur_sample_h['VF']      = filter_tag
			cur_sample_h['AI']      = AI
			cur_sample_h['PLOIDY']  = PLOIDY
			cur_sample_h['PLtag']   = PLtag
			cur_sample_h=give_filter(cur_sample_h)
		else:
			two_sample_h = list()
			aids = [a1, a2]
			for i in ( 0, 1 ):
				Zyg, AD, AR=['.']*3
				cur_aid  = aids[i]
				pair_aid = aids[ i - 1 ]
				AI = i + 1
				if cur_aid == '.':
					nocall_var['sample'][sample]=dict()
					two_sample_h[i]=nocall_var['sample'][sample]
					involve_nocall = 1
					Zyg = 'no-call'
				else:
					if 'sample' not in vars[int(cur_aid)]:
						vars[int(cur_aid)]['sample']=dict()
					vars[int(cur_aid)]['sample'][sample]=dict()
					two_sample_h.insert(i,vars[int(cur_aid)]['sample'][sample])
					if pair_aid == '.':
						Zyg='half-ref' if cur_aid == '0' else 'half-alt'
					elif pair_aid != '0': # pair_aid eq '0'
						Zyg = 'het-alt'
					else:
						Zyg = 'het-ref'
					if total_Dp > 0:
						if int(cur_aid) >= len(ADs):
							print "{}{}{}{}{}".format("AD index error! [",chr, vcf_pos, sample,"]")
						else:
							AD = ADs[int(cur_aid)]
							AR = '%.2f' % (AD/total_Dp)
				two_sample_h[i]['zygosity']= Zyg
				two_sample_h[i]['NB']      = NB
				two_sample_h[i]['PB']      = PB
				two_sample_h[i]['AD']      = AD
				two_sample_h[i]['AR']      = AR
				two_sample_h[i]['AI']      = AI
				two_sample_h[i]['PLOIDY']  = PLOIDY
				two_sample_h[i]['VF']      = filter_tag
				two_sample_h[i]['PLtag']   = PLtag
				two_sample_h[i]=give_filter(two_sample_h[i])
	
	if involve_nocall:
		vars.append(nocall_var)
	
	available_vars = list()
	refVar = 1
	for v in vars:
		zyg_tmp=None
		if refVar:
			zyg_tmp = 'ref'
			refVar = 0
		else:
			zyg_tmp = 'alt'
		if 'sample' not in v: # those dabase like vcf don't have sample infos
			v['sample']={
				"nullSample":{
					'zygosity':'{}{}'.format('homo-',zyg_tmp),
					'AD' : '.',
					'AR' : '.',
					'VF' : '.',
					'AI' : 0,
					'PLOIDY' : 2,
					'PLtag' : '.',
					'PB' : '.',
					'NB' : '.',
					'filterTag' : '.'
				}
			}
	return [1, vars]

NOAHCARE=None
NOAHCARE=os.path.dirname(os.path.abspath(sys.argv[0]))

VERSION = '1.10'

beda=None


# database file default name configuration




hg19_fa   = NOAHCARE+'/hg19_chM_db/hg19_chM.fa'
db        = NOAHCARE + '/annodb/annodb/ncbi_anno_rel104_db.bed.gz'
trDB      = NOAHCARE + '/annodb/annodb/ncbi_anno_rel104_trSeq.fas'
cytoband  = NOAHCARE + '/annodb/cytoBand/cytoBand_hg19_grch37.txt.gz'
pfam      = NOAHCARE + '/annodb/pfam/Pfam-A-ncbi_2012-12-21.bed.gz'
preds     = NOAHCARE + '/annodb/predictions/predictDB_for_anno104.tab.gz'
condel    = NOAHCARE + '/config/Condel/'
phyloP    = NOAHCARE + '/annodb/phyloP/phyloP_all3class_combin_2013-09-25.bed.gz'
dbsnp     = NOAHCARE + '/annodb/dbsnp/snp147.bed.gz'
esp6500   = NOAHCARE + '/annodb/NHLBI/ESP6500SI-V2-SSA137.NHLBI.bed.rmanchor.uniq.gz'
tgp       = NOAHCARE + '/annodb/tgp/TGP.phase3.bed.gz'
repeat_db = NOAHCARE + '/annodb/rep_dup_db/dup_trf_rmsk.bed.gz'
omim_data = NOAHCARE + '/annodb/omim/anno_combin.tsv.gz'
cgd_data  = NOAHCARE + '/annodb/CGD/CGD.txt.gz'
cosmic    = NOAHCARE + '/annodb/cosmic/COSMICv80_GRCh37_20170213.bed.gz'
exac_data = NOAHCARE + '/annodb/exac/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz'
gnomAD_db = NOAHCARE + '/annodb/gnomAD/gnomad.exomes.r2.0.1.sites.vcf.gz'
hgmd      = NOAHCARE + '/annodb/hgmd/localhgmd.db.gz'
verList   = NOAHCARE + '/config/db_version.list'
msqc_info = NOAHCARE + '/analyze/snv/MSQC_info.list'
clinVar   = NOAHCARE + '/annodb/clinVarDB/clinVarDB.vcf.gz'
vusVar    = NOAHCARE + '/annodb/clinVarDB/vusVarDB.vcf.gz'



DB_VERSION = {
"REFERENCE_BUILD"  : 'GRCh37',
"ANNOTATION_BUILD" : 'NCBI Annotation Release 104',
"ENSDB_BUILD"      : '73',
"DBSNP_BUILD"      : '147',
"COSMIC_BUILD"  : '80',
"TGP_BUILD" : 'phase3_v1',
"CGD_Date" : '2013-12-30',
"OMIM_Date" : '2015-01-19',
"PFAM_Date" : '2012-12-21',
"LOCALHGMD_YEAR" : '2012',
"ESP6500_VERSION" : '2',
"ExAC_VERSION" :'0.3.1',
"GAD_VERSION" : '2.0.1',
"UCSC_HG19_Date" : '2013-12-06',
"CLIVAR_Date" : '2016-11-28'
}

finished=None  # don't define it any where until finished.



check_dep( db, trDB, hg19_fa ) 

help,         outfile,       sample_list, quiet          = None, None, None, None
headrule,     type,          offline,     number_threads = None, None, None, None
buffer_cache, write_version, debug,       timer          = None, None, None, None
in_msqc,      msqc_rst ,    logfile                      = None, None, None

# running mode
type        = "vcf"
NUM_THR     = 2
#BUFFER_SIZE = 5000
MAX_REF_LEN = 200

# filter rule thresholds
AD_UP_THRESHOLD = 8  # used to decide PASS
AD_DN_THRESHOLD = 2  # used to decide FAIL
PL_THRESHOLD    = 3  # used to decide PLTAG (FAIL)
PL_UP_THRESHOLD = 10 # used to decied PLTAG (DUBIOUS)

deleterious_func=dict()
deleterious_func=dict.fromkeys("cds-loss init-loss stop-gain stop-loss frameshift knockout nonsense".split(),1)

possible_deleterious_func=dict()
possible_deleterious_func=dict.fromkeys("missense stop-retained splice splice-3 splice-5 span abnormal-intron".split(),1)

autoInterp_freq_threshold = {
'TGP_AF'    : 0.01,
'ESP6500AF' : 0.01,
'PVFD_AF'   : 0.01,
'ExAC_AF'   : 0.01,
'GAD_AF'    : 0.01,
'PanelAF'   : 0.05
}

kickoutExcel_function=None
kickoutExcel_function=dict.fromkeys("utr-3 utr-5 promoter intron annotation-fail unknown-no-call unknown no-change .".split(),1)

AutoInterpWords = {
	"Certain" : 1,
	"Likely Deleterious" : 2,
	"VUS" : 3,
	"Likely Benign" : 4,
	"Benign" : 5,
	"Unknown" : 6,
}

xpar_reg = [
	[ 60000, 2699520 ],
	[ 154931043, 155260560 ]
]

ypar_reg = [
	[ 10000, 2649520 ],
	[ 59034049, 59363566 ]
]

VarNameList=dict()

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('config_file')
parser.add_argument('varfile',nargs="+")
parser.add_argument('-t','--type', help='type of input files vcf/tsv, default vcf,\nall input files should be in the same\ntype, mixed cases are not supported.',required=True,type=str)
parser.add_argument('-n','--numthreads', help='Number of threads to be used for annotation',required=True,type=int)
#parser.add_argument('-b','--buffersize', help='Number of record to be cached, ready for\nsorting and output.',required=True,type=int)
parser.add_argument('-l','--logfile', help='multiprocessing log file.',required=True,type=str)
parser.add_argument('-s','--samplist', help='sample list used in pipeline, if sample\nlist is specified, extra field "FAM_ID", "SAMPLE_ID", "Sex" will be added into output, default only "SAMPLE_ID" got from vcf col header or from CG\'s input file name is added.',type=str)
parser.add_argument('-o','--outfile', help='Output file, default STDOUT.',required=True,type=str)
parser.add_argument('-r','--headrule', help='Output header format file, default \nno output.',type=str)
parser.add_argument('-i','--inmsqc', help='Result of MS Genotyping result for QC sites',type=str)
parser.add_argument('-c','--qcresult', help='Output MSQC result to FILE, default STDERR',type=str)
parser.add_argument('-f','--offline', help='Using an offline version of database \nfor emergency, network error or any other \nspecial situation.',action="store_true")
parser.add_argument('-q','--quiet', help='Suppress the warning information, default \nnot suppress.',action="store_true")
parser.add_argument('-d','--debug', help='Show debug information',action="store_true")
parser.add_argument('-m','--timer', help='Show timer information',action="store_true")
parser.add_argument('-w','--writever', help='write DB version information to verList, \ndefault no write.',action="store_true")

args = parser.parse_args()

type           = args.type
number_threads = args.numthreads
#buffer_cache   = args.buffersize
sample_list    = args.samplist
outfile        = args.outfile
headrule       = args.headrule
in_msqc        = args.inmsqc
msqc_rst       = args.qcresult
offline        = args.offline
quiet          = args.quiet
debug          = args.debug
timer          = args.timer
write_version  = args.writever
logfile        = args.logfile

#log_file="/home/zhouzhq/workdir/project/NCanno/python_BNC/BNC/python_multiprocessing_log"
f=logging.Formatter('[%(levelname)s %(processName)s %(asctime)s %(funcName)s] %(message)s')
h=logging.FileHandler(logfile, 'w')
h.setFormatter(f)
logger=multiprocessing.get_logger()
logger.addHandler(h)
logger.setLevel(logging.INFO)
logger.info('program NCanno started in %s with command: %s', os.getcwd(), ' '.join(sys.argv))
logger.info('\n[ the number of process is  %s]',number_threads)
logger.info('\n[ start time: %s]',time.asctime( time.localtime(time.time()) ))
opts=dict() # BedAnno options
# some db is not available online currently, 
# these local version database APIs have been 
# integrated into BedAnno
check_dep(cytoband,  pfam,   preds,  condel, repeat_db, esp6500,omim_data, phyloP, cosmic, dbsnp,  exac_data, gnomAD_db, tgp)


opts['cytoBand']   = cytoband
opts['pfam']       = pfam
opts['prediction'] = preds
opts['condel']     = condel
opts['esp6500']    = esp6500
opts['phyloP']     = phyloP
opts['cosmic']     = cosmic
opts['rmsk']       = repeat_db
opts['dbSNP']      = dbsnp
opts['genome']     = hg19_fa
opts['exac']       = exac_data
opts['gnomAD']     = gnomAD_db
opts['tgp']        = tgp

# Omim and CGD is outof the scope of BedAnno
check_dep(cgd_data, hgmd, clinVar, vusVar)

OutFp=sys.stdout
OUT, QCSTAT, HDR = None, None, None

if outfile is not None:
	OUT=open(outfile,"w")
	OutFp=OUT

if msqc_rst is not None:
	QCSTAT=open(msqc_rst,"w")

if headrule is not None:
	HDR=open(headrule,"w")

# read the config items
#config=dict()

sys.path.append(os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(args.config_file)))))
(filepath,tempfilename)=os.path.split(args.config_file)
(shortname,extension)=os.path.splitext(tempfilename)
config_module=__import__(shortname)
config= config_module.config

# Annotation rules
# available pop: ASN,AMR,AFR,EUR
# available pop: EAS,SAS,AMR,AFR,EUR,FIN,NFE,OTH
for k,v in zip(['TGP_POP','ExAC_POP','PVFD_rule','GaP_rule','Panel_Control','Panel_ID','FLANK_LEN','intron_edge_remain'],['ASN','EAS','C1125','.','.','.',2,0]):
	if k not in config:
		config[k]=v


if "DBVER_LIST" in config:
	verList =config['DBVER_LIST']
if "MSQC_INFO" in config:
	msqc_info =config['MSQC_INFO']
if "ANNO_THR_NUM" in config:
	NUM_THR =config['ANNO_THR_NUM']
#if "ANNO_BUFFER_SIZE" in config:
#	BUFFER_SIZE =config['ANNO_BUFFER_SIZE']
if number_threads is not None:
	NUM_THR = number_threads
#if buffer_cache is not None:
#	BUFFER_SIZE = buffer_cache

# command line option have the highest priority
if NUM_THR < 2:
	NUM_THR = 2
#if BUFFER_SIZE < 100:
#	BUFFER_SIZE = 100
#SINGLE_TIMEOUT = 60 * (NUM_THR - 1)

if 'autoInterp_freq_threshold' in config:
	autoInterp_freq_threshold=copy.deepcopy(config['autoInterp_freq_threshold'])

if 'kickoutExcel_function' in config:
	kickoutExcel_function=copy.deepcopy(config['kickoutExcel_function'])

possible_deleterious=dict()
if 'possible_deleterious' in config:
	possible_deleterious_func.update(dict.fromkeys(config['possible_deleterious'],1))

if 'VarNameList' in config:
	VNL=None
	try:
		VNL=open(config['VarNameList'],"r")
	except e:
		raise BaseException("Error: ["+config['VarNameList']+"]"+e)
	for line in VNL:
		if re.search(r'^#',line) or re.search(r'^\s*$',line):
			continue
		line=re.sub(r'\s+$','',line)
		g, v, n = line.split('\t')[0:3]
		if g not in VarNameList:
			VarNameList[g]=dict()
		VarNameList[g][v]=n
	VNL.close()

if 'CoreFuncRegion' in config:
	check_dep(config['CoreFuncRegion'])

GeneTestCode=None
if 'GeneTestCode' in config:
	GeneTestCode=read_TestCode(config['GeneTestCode'])

# init Homo reference var list
HomoRef_Var, HRVanno_opt=None,None
if 'HomoRef_Var' in config:
	HomoRef_Var = read_HomoRef_var(config['HomoRef_Var'])

# read and init MS QC info
QCsites, msQCrst, ngsQCrst, nQCs=None,None,None,None
if in_msqc is not None:
	QCsites = read_MSQC_info(msqc_info)
	msQCrst = read_inmsqc(in_msqc)
	msqc_samples = msQCrst.keys()
	sids = dict()
	for chr in QCsites.keys():
		sids+=sorted(QCsites[chr].keys())
	qcsids=dict.fromkeys(sids,1)
	nQCs=len(qcsids.keys())
	if not (0 < len(msqc_samples) and nQCs==len(msQCrst[msqc_samples[0]])):
		raise Exception("\n".join(["Error: inconsistent - QC sites info file and the MS QC result","       QC sites info file : "+msqc_info,"       MS QC result file  : "+in_msqc]))

sample_info = dict()
if sample_list is not None:
	sample_info = read_sample_list(sample_list)

if "genes" in config:
	check_dep(config['genes'])
	opts['genes']=config['genes']

if "trans" in config:
	check_dep(config["trans"])
	opts['trans']=config['trans']


sys.stderr.write("{}{}{}".format("== Annotation Start at ", time.asctime( time.localtime(time.time()) )," ==\n"))

# shift the first args of config file
CONFIG_FILE = args.config_file

var_count = dict()
no_call_var_count = dict()
sex_varcount = dict()
ti = dict() # A <=> G, C <=> T
tv = dict() # others
annotation_count = dict()
pre_total_AD = dict()

out_header_keys=list()



print_header()

#manager = multiprocessing.Manager()
anno_cache =list()

# The input variation file must be carefully sorted
#cached_count=multiprocessing.Value('i', 0)
#cached_count_lock=multiprocessing.Lock()

# Annotation Queue
anno_q=list()

#for i in range(1,NUM_THR ):
'''
anno_thrs.apply_async(anno_threads,anno_cache,anno_q)
anno_thrs.close()
anno_thrs.join()
'''
prev_var=dict()
#enqueue_count = 0
read_record_time = 0
#main_record_count = 0
total_record_count = 0
tsv_type=None
coreFuncRegion_h=None
if "CoreFuncRegion" in config:
	coreFuncRegion_h=CheckInOut.CheckInOut({"db":config['CoreFuncRegion']})

for var_file in args.varfile:
	varf_h=None
	tsv_sample_name=None
	if type!="tsv":
		varf_h =vcf.Reader(open(var_file, 'r'))
	else:
		base_filename = os.path.basename(var_file)
		m=re.search(r'^var\-([^\\\/]+)\.tsv(\.bz2)?$',base_filename)
		m2=re.search(r'^masterVarBeta\-([^\\\/]+)\.tsv(\.bz2)?$',base_filename)
		if m:
			tsv_sample_name = m.group(1)
			if m.group(2) is None:
				try:
					varf_h=open(var_file)
				except Exception,e:
					raise Exception("Error: ["+var_file+"]"+e)
			else:
				try:
					varf_h=bz2.BZ2File(var_file)
				except Exception,e:
					raise Exception("Error: ["+var_file+"]"+e)
			tsv_type = "var"
		elif m2:
			tsv_sample_name = m2.group(1)
			if m2.group(2) is None:
				try:
					varf_h=open(var_file)
				except Exception,e:
					raise Exception("Error: ["+var_file+"]"+e)
			else:
				try:
					varf_h=bz2.BZ2File(var_file)
				except Exception,e:
					raise Exception("Error: ["+var_file+"]")
			tsv_type = "master"
		else:
			raise Exception("\n".join(["Error: cannot recognized the file name of tsv result.\n","      ["+var_file+"] should be in format of var-[ASM_ID].tsv.bz2,\n","      of in format of masterVarBeta-[ASM_ID].tsv.bz2.\n","      ASM_ID is the sample ID in CG's system. which defined in\n","      Standard Sequencing Service Data Format v2.4"]))
	
	read_stat = 1
	pre_chr=None
	sample_ids=None
	while(read_stat):
		tobe_anno_vars = list()
		t0, t1=None,None
		if timer is not None:
			t0 = time.time()
		if type !="tsv":
			read_stat, tobe_anno_vars = read_vcf(varf_h)
		else:
			if tsv_type == "var":
				read_stat, tobe_anno_vars=read_tsv(varf_h, tsv_sample_name)
			elif tsv_type == "master":
				read_stat, tobe_anno_vars=read_master(varf_h, tsv_sample_name)
		if timer is not None:
			t1 = time.time()
			read_record_time+=t1-t0
			'''
			if read_stat:
				#main_record_count+=1
				pass
			if main_record_count>BUFFER_SIZE:
				sys.stderr.write("{}{}{}{}{}".format("[TIMER RECORD_READ ",BUFFER_SIZE,"] : ",read_record_time,"\n"))
				read_record_time = 0
				main_record_count = 0
			'''
		
		if not read_stat or (pre_chr is not None and pre_chr!=tobe_anno_vars[0]['chr']):
			# check homo ref vars to be add before chr changing or end of file
			if HomoRef_Var is not None \
			and sample_ids is not None \
			and pre_chr is not None \
			and pre_chr in HRVanno_opt \
			and ( \
			pre_chr  not in HRVanno_opt \
			or var_file not in HRVanno_opt[pre_chr] \
			):
				chrom_final_add = list()
				for var_in_homoR in HomoRef_Var[pre_chr]:
					if 'varfile' not in var_in_homoR or var_file not in var_in_homoR['varfile']:
						new_gen =gen_homoRef( var_in_homoR, sample_ids )
						if new_gen:
							chrom_final_add.append(new_gen)
						var_in_homoR['varfile'][var_file]=1
				#enqueue_count+=len(chrom_final_add)
				anno_q.append(chrom_final_add)
				HRVanno_opt[pre_chr][var_file]=1
		'''
		if enqueue_count >= BUFFER_SIZE   \
		or not read_stat \
		or (pre_chr is not None and pre_chr !=tobe_anno_vars[0]['chr']):
			print_cache()
			
			enqueue_count =enqueue_count- cached_count
			anno_cache   = list()
			cached_count = 0
		'''
		if read_stat:
			
			total_record_count+=1
			# intilize sample ids in current vcf file
			if sample_ids is None:
				sample_ids=sorted(tobe_anno_vars[0]['sample'].keys())
			
			if HomoRef_Var is not None:
				# add reference allele to vcf no-called region due to
				# by default only mutant site will be output in VCF
				check_to_add_homoRef(tobe_anno_vars,var_file,sample_ids)
			#enqueue_count += len(tobe_anno_vars)
			pre_chr = tobe_anno_vars[0]['chr']
			map(anno_q.append,tobe_anno_vars)
			
			if not re.search(r'^chr',pre_chr):
				qcusing_chr ="{}{}".format('chr',pre_chr)
			if in_msqc is not None and msqc_rst is not None and qcusing_chr in QCsites:
				# due to the co-presence of variants for one locus
				# we check ms-qc sites here
				for qcidx in QCsites[qcusing_chr].keys():
					cur_qc_pos = QCsites[qcusing_chr][qcidx]["pos"]
					cur_qc_ref = QCsites[qcusing_chr][qcidx]["ref"]
					for tba in tobe_anno_vars:
						if 'end' in tba:    # CG var format
							if tba['begin']< cur_qc_pos and tba['end']>=cur_qc_pos:
								if tba['variantSequence']==tba['referenceSequence']:    # reference called allele
									for samp in tba['sample'].keys():
										if samp not in ngsQCrst or qcidx not in ngsQCrst[samp] or ngsQCrst[samp][qcidx]!='.':
											ngsQCrst[samp][qcidx]="{}{}".format(ngsQCrst[samp][qcidx],cur_qc_ref)
								elif tba['end']-tba['begin']==1 and re.search(r'^[ACGT]$',tba['variantSequence'],flags=re.IGNORECASE):
									for samp in tba['sample'].keys():
										if samp not in ngsQCrst or qcidx not in ngsQCrst[samp] or ngsQCrst[samp][qcidx]!='.':
											ngsQCrst[samp][qcidx]="{}{}".format(ngsQCrst[samp][qcidx],tba['variantSequence'])
								else:
									# other cases all equal to not available
									# This will overide any called allele
									for samp in tba['sample'].keys():
										ngsQCrst[samp][qcidx]='.'
						else:                         # vcf var format
							if tba['begin']<=cur_qc_pos and (tba['begin']+len(tba['referenceSequence'])-1)>= cur_qc_pos:
								if tba['variantSequence'] == tba['referenceSequence']:    # reference called allele
									for samp in tba['sample'].keys():
										if samp not in ngsQCrst or qcidx not in ngsQCrst[samp] or ngsQCrst[samp][qcidx]!='.':
											ngsQCrst[samp][qcidx]="{}{}".format(ngsQCrst[samp][qcidx],cur_qc_ref)
								elif tba['begin']==cur_qc_pos and re.search(r'^[ACGTN]+$',tba['variantSequence']):
									# assume anchor sites to be available
									# for qc using
									called_allele=tba['variantSequence'][0:1]
									for samp in tba['sample'].keys():
										if called_allele!='N':
											if samp not in ngsQCrst or qcidx not in ngsQCrst[samp] or ngsQCrst[samp][qcidx]!='.':
												ngsQCrst[samp][qcidx]="{}{}".format(ngsQCrst[samp][qcidx],called_allele)
										else:
											ngsQCrst[samp][qcidx]='.'
								else:    # any ranged or complex variants
									for samp in tba['sample'].keys():
										ngsQCrst[samp][qcidx]='.'
		elif in_msqc is not None and msqc_rst is not None: # finished reading
			# complete qc-sites info for each sample
			for samp in sample_ids:
				for i in range(nQCs):
					for ch in sorted(QCsites.keys()):
						if i not in QCsites[ch]:
							continue
						if samp not in ngsQCrst or i not in ngsQCrst[samp]:
							ngsQCrst[samp][i]=QCsites[ch][i]['ref']

anno_thrs= Pool(processes=NUM_THR)
'''
GetCGDinfo_obj=GetCGDinfo.GetCGDinfo({"db" : cgd_data, "quiet" : 1})
GetCGDinfo_dict=Manager().dict(GetCGDinfo_obj.cgd_h)
'''
result = list()
for line in anno_q:
	result.append(anno_thrs.apply_async(anno_threads,(line,)))
anno_thrs.close()
anno_thrs.join()
for res in result:
	anno_cache.append(res.get())
print_cache(anno_cache)
OutFp.close()
# Output all samples qc comparison into a single file
if msqc_rst is not None and in_msqc is not None:
	sys.stderr.write("\t".join(['#SampleID','Class']+[str(i) for i in range(1,22)])+"\n")
	for samp in sorted(msQCrst.keys()):
		
		if samp not in ngsQCrst:
			sys.stderr.write("{}{}{}".format("[Warning] the MS QC result contains some samples not in input. (",samp,") \n"))
			continue
		
		qc_cmp_stat = list()
		for i in range(nQCs):
			if msQCrst[samp][i]=='.' or \
			ngsQCrst[samp][i]=='.' or \
			re.search(r'[^ACGT]',msQCrst[samp][i]) or \
			re.search(r'[^ACGT]',ngsQCrst[samp][i]):
				qc_cmp_stat[i] = '?'
			else:
				uniform_ms ="".join(sorted(list(msQCrst[samp][i]))).upper()
				uniform_ngs ="".join(sorted(list(ngsQCrst[samp][i]))).upper()
				if uniform_ms==uniform_ngs:
					qc_cmp_stat[i]='ok'
				else:
					qc_cmp_stat[i]='unmatch'
		
		print>>QCSTAT,"\t".join([samp,"(MS)"]+msQCrst[samp])+"\n"
		print>>QCSTAT,"\t".join([samp,"(NGS)"]+ngsQCrst[samp])+"\n"
		print>>QCSTAT,"\t".join([samp,"(RST)"]+qc_cmp_stat)+"\n"
	QCSTAT.close()
'''
tmp_manager = multiprocessing.Manager()
temp_queue=tmp_manager.Queue()
temp_queue.put([None]*(NUM_THR - 1))
while True:
	if not anno_q.empty():
		temp_queue.put(anno_q.get())
	else:
		break
anno_q.put([None]*(NUM_THR - 1))
while 1:
	anno_q.put(temp_queue.get())
#anno_q=copy.deepcopy(temp_queue)
time.sleep(1)
'''
# Output Statistic Information
sys.stderr.write("\n\n============= Variation Statistics ============\n")
for sampleID in sorted(var_count.keys()):
	sample_total=var_count[sampleID]['total']
	no_call_total=no_call_var_count[sampleID] if sampleID in no_call_var_count else 0
	real_total=sample_total+no_call_total
	sys.stderr.write("{}{}{}{}".format("[INFO] "+str(sampleID)+" - called ratio : ",str(sample_total)+" / "+str(real_total)+" (","%.3f" %(sample_total/real_total * 100),"%)\n"))
	del var_count[sampleID]['total']
	for classID in sorted(var_count[sampleID].keys()):
		for sub_classID in sorted(var_count[sampleID][classID].keys()):
			sys.stderr.write("{}{}{}{}".format("[INFO] "+str(sampleID)+" - "+str(classID)+" - "+str(sub_classID)+" : ",str(var_count[sampleID][classID][sub_classID])+" / "+str(sample_total)+" (","%.3f" %(var_count[sampleID][classID][sub_classID]/sample_total*100),"%)\n"))
	
	if sampleID not in tv:
		tv[sampleID]=0
	if sampleID not in ti:
		ti[sampleID]=0
	tmp=None
	if tv[sampleID]==0:
		tmp="inf"
	else:
		tmp="%.3f" %(ti[sampleID]/tv[sampleID])
	sys.stderr.write("{}{}{}{}".format("[INFO] "+str(sampleID)+" - ti/tv : ",str(ti[sampleID])+" /  "+str(tv[sampleID])+" = ",tmp,"\n"))
	if sampleID in sex_varcount:
		if "het" not in sex_varcount[sampleID]:
			sex_varcount[sampleID]['het']=0
		sys.stderr.write("{}{}{}{}".format("[INFO] "+str(sampleID)+" - het in sex : ",str(sex_varcount[sampleID]['het'])+" / "+str(sex_varcount[sampleID]['total'])+" (","%.3f" %(sex_varcount[sampleID]['het']/sex_varcount[sampleID]['total']*100),"%)\n"))
	
	anno_total=annotation_count[sampleID]['total']
	del annotation_count[sampleID]['total']
	for classID in sorted(annotation_count[sampleID].keys()):
		for sub_classID in sorted(annotation_count[sampleID][classID].keys()):
			sys.stderr.write("{}{}{}{}".format("[INFO] "+str(sampleID)+" - "+str(classID)+" - "+str(sub_classID)+" : ",str(annotation_count[sampleID][classID][sub_classID])+" / "+str(anno_total)+" (","%.3f" %(annotation_count[sampleID][classID][sub_classID]/anno_total*100),"%)\n"))
sys.stderr.write("{}{}{}".format("== Annotation Done at ",time.asctime( time.localtime(time.time()) )," ==\n"))

finished = 1
time.sleep(1)
logger.info('\n[ end time: %s]',time.asctime( time.localtime(time.time()) ))




