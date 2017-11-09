import pysam
import os
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
class GetGAD(object):
	@LogExc
	def __init__(self, opts):
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getGAD(self, chrom, start, end, ref, alt):
		logger=multiprocessing.get_logger()
		
		qstart, qend  = start, end
		if ref == "":
			ref="."
		if alt == "":
			alt="."
		if start == end:
			qstart -= 1
			qend += 1
		q=None
		try:
			q = self.tidb.fetch(chrom, qstart, qend)
		except Exception, e:
			return dict()
		n=None
		try:
			n=next(q)
		except Exception, e:
			return dict()
		def process_rec(rec):
			logger=multiprocessing.get_logger()
			rec =rec.strip("\r\n")
			linesInfo = rec.split('\t')
			vcfref = linesInfo[3]
			vcfalt = linesInfo[4]
			altcount = 0
			altsList = vcfalt.split(',')
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
						tstart =int(tstart) - 1
				elif len(talt) > len(tref) and len(tref)==1:
					talt =talt[1:]
					tref = "."
				elif len(talt) < len(tref) and len(talt)==1:
					tref = tref[1:]
					talt = "."
					tend = int(tstart) + len(tref)
				else:
					tref = tref[1:]
					talt = talt[1:]
					tend = int(tstart) +len(tref)
				if int(tstart) == int(start) and int(tend) == int(end) and str(tref) == str(ref) and str(talt) == str(alt):
					outinf =GetGAD.fetch_inf(linesInfo[7], inum)
					return outinf
			return None
		ret_val=process_rec(n)
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
		output = {"AN" : 0, "AF" : 0, "AN_EAS" : 0, "AF_EAS" : 0, "Hom" : 0, "Hom_EAS" : 0}
		for tags in ainf:
			tmp=tags.split('=')
			header, svalue="",""
			if len(tmp) ==2:
				header, svalue = tmp[0], tmp[1]
			elif len(tmp) ==1:
				header, svalue = tmp[0], ""
			if header == "AC":
				output['AN'] = GetGAD.split_inf(svalue, inum)
			elif header == "AF":
				output["AF"] = GetGAD.split_inf(svalue, inum)
			elif header == "AC_EAS":
				output["AN_EAS"] = GetGAD.split_inf(svalue, inum)
			elif header == "AF_EAS":
				output["AF_EAS"] = GetGAD.split_inf(svalue, inum)
			elif header == "Hom":
				output["HASC"] = GetGAD.split_inf(svalue, inum)
			elif header == "Hom_EAS":
				output["HASC_EAS"] = GetGAD.split_inf(svalue, inum)
			else:
				continue
		return output
	@staticmethod
	@LogExc
	def split_inf(inf, inum):
		tagInf = inf.split(',')
		if len(tagInf) > inum:
			return tagInf[inum]
		else:
			return inf