import os, pysam
import multiprocessing
import functools
def LogExc(f):
    @functools.wraps(f)
    def wrapper(*largs, **kwargs):
        try:
            res=f(*largs, **kwargs)
        except Exception, e:
            multiprocessing.get_logger().exception(e)
            raise
        return res
    return wrapper
class GetClinVar(object):
	SSRINF = {
		"0"    : "unspecified",
		"1"    : "Paralog",
		"2"    : "byEST",
		"4"    : "oldAlign",
		"8"    : "Para_EST",
		"16"   : "1kg_failed",
		"1024" : "other",
		"2048"   : "no_record"
	}
	SSRSET = [ 1, 2, 4, 8, 16 ]
	CLNSIGINF = {
		"0"   : "Uncertain significance",
		"1"   : "not provided",
		"2"   : "Benign",
		"3"   : "Likely benign",
		"4"   : "Likely pathogenic",
		"5"   : "Pathogenic",
		"6"   : "drug response",
		"7"   : "histocompatibility",
		"255" : "other",
		"254"   : "unknow_impact"
	}
	@LogExc
	def __init__(self, opts):
		if "db" not in opts or not os.path.exists(opts["db"]) or not os.access(opts["db"],os.R_OK):
			raise "db option error!"
		if not os.path.exists(opts["db"]+".tbi") or not os.access(opts["db"]+".tbi",os.R_OK):
			raise "Error: ["+str(opts["db"])+"] index (.tbi) file not found, please build it."
		opts['tidb'] = pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getCL(self, chrom, start, end, ref, alt):
		logger=multiprocessing.get_logger()
		qstart, qend  = start, end
		if ref == "":
			ref = "."
		if alt == "":
			alt = "."
		if int(start) == int(end):
			qstart-=1
			qend+=1
		q = None
		try:
			q = self.tidb.fetch(chrom, qstart, qend)
		except Exception,e:
			return dict()
		tmp=None
		try:
			tmp=next(q)
		except Exception, e:
			return dict()
		@LogExc
		def process_rec(rec):
			logger=multiprocessing.get_logger()
			rec=rec.strip("\r\n")
			linesInfo=rec.split('\t')
			vcfref    = linesInfo[3]
			vcfalt    = linesInfo[4]
			altcount  = 0
			altsList  = vcfalt.split(',')
			for inum in range(len(altsList)):
				talt = altsList[inum]
				tstart, tend = linesInfo[1], linesInfo[1]
				tref = vcfref
				if len(talt) == len(tref):
					refbaseList = list(vcfref)
					altbaseList = list(talt)
					if refbaseList[0] != altbaseList[0]:
						talt = altbaseList[0]
						tref = refbaseList[0]
						tstart=int(tstart)-1
				elif len(talt) > len(tref) and len(tref) == 1:
					# ins
					talt = talt[1:]
					tref = "."
				elif len(talt) < len(tref) and len(talt) ==1:
					# del
					tref = tref[1:]
					talt = "."
					tend = int(tstart) + len(tref)
				else:
					#delins
					tref = tref[1:]
					talt = talt[1:]
					tend = int(tstart) + len(tref)
				if int(tstart) == int(start) \
					and int(tend) == int(end) \
					and str(tref) == str(ref) \
					and str(talt) == str(alt) :
					outinf = GetClinVar.fetch_inf(linesInfo[7], inum)
					return outinf
			return None
		ret_val=process_rec(tmp)
		if ret_val is None:
			for rec in q:
				ret_val=process_rec(rec)
				if ret_val is not None:
					return ret_val
		else:
			return ret_val
		return dict()
	@staticmethod
	@LogExc
	def fetch_inf(inf, inum):
		logger=multiprocessing.get_logger()
		ainf = inf.split(';')
		output={"SSR":"2048","CLNSIG":"254","CLNREVSTAT":"no_record","CLNACC":"no_record"}
		for tags in ainf:
			arr = tags.split('=')
			header, svalue=None, None
			if len(arr)==1:
				header, svalue=arr[0],""
			elif len(arr)==2:
				header, svalue=arr[0],arr[1]
			
			if header == "SSR":
				output["SSR"] = GetClinVar.split_inf(svalue, inum)
			elif header == "CLNSIG":
				output["CLNSIG"] = GetClinVar.split_inf(svalue, inum)
			elif header == "CLNREVSTAT":
				output["CLNREVSTAT"] = GetClinVar.split_inf(svalue, inum)
			elif header == "CLNACC":
				output["CLNACC"] = GetClinVar.split_inf(svalue, inum)
		clnr = output["CLNSIG"].split('|')
		outcln=list()
		for stag in clnr:
			if stag not in GetClinVar.CLNSIGINF:
				outcln.append('')
			else:
				outcln.append(GetClinVar.CLNSIGINF[stag])
		output["CLNSIG"]="|".join(outcln)
		output["SSR"]=GetClinVar.analyse_ssr(output["SSR"])
		return output
	@staticmethod
	@LogExc
	def analyse_ssr(inf):
		if inf in GetClinVar.SSRINF:
			return GetClinVar.SSRINF[inf]
		else:
			output = list()
			for tag in GetClinVar.SSRSET:
				if tag & inf ==tag:
					output.append(GetClinVar.SSRINF[tag])
			if len(output)==0:
				return GetClinVar.SSRINF["2048"]
			else:
				return "|".join(output)
	@staticmethod
	@LogExc
	def split_inf(inf, inum):
		tagInf = inf.split(',')
		if len(tagInf)>inum:
			return tagInf[inum]
		else:
			return inf