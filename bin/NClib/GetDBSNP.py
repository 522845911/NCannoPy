import pysam
import re
import os
from string import maketrans
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
class GetDBSNP(object):
	@LogExc
	def __init__(self, opts):
		logger=multiprocessing.get_logger()
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getRS(self, chr, start, end, ref, alt):
		logger=multiprocessing.get_logger()
		if not re.search(r'^chr', chr):
			chr = 'chr' + str(chr)
		qstart, qend = start, end
		if start == end:
			qstart -= 1
			qend   += 1
		q=None
		try:
			q = self.tidb.fetch(chr, qstart , qend)
		except Exception, e:
			return dict()
		n=None
		try:
			n=next(q)
		except Exception, e:
			return dict()
		
		rdbsnp_rsIDs = dict()
		@LogExc
		def process_rec(rec):
			if re.search(r'^\s*$', rec):
				return
			rec = rec.strip("\r\n")
			itms = rec.split('\t')
			if len(itms) < 25:
				raise Exception("Error: ["+str(self.db)+"] not enough columns ["+str(len(itms))+"]. ["+str(chr)+","+str(start)+","+str(end)+"].")
			if int(itms[1]) != int(start) or (int(itms[2]) != int(end) and not re.search(r'microsat', itms[10])):
				return
			alleles = GetDBSNP.group_revcom(itms[21]) if itms[5] == '-' else itms[21]
			total_AC = 0
			altAF = 0
			if alleles != '':
				alleles  = alleles.split(',')
				alleleNs = itms[22].split(',')
				alleleFs = itms[23].split(',')
				alleleFs=[alleleFs[i] for i in range(len(alleleFs)) if alleleFs[i]]
				alleleNs=[alleleNs[i] for i in range(len(alleleNs)) if alleleNs[i]]
				alleles=[alleles[i] for i in range(len(alleles)) if alleles[i]]
				if re.search(r'microsat', itms[10]):
					for i in range(len(alleles)):
						m=re.search(r'^\(([ACGT]+)\)(\d+)$', alleles[i])
						if m:
							alleles[i] = m.group(1)*int(m.group(2))
				hit_altAF = 0
				for k in range(len(alleles)):
					
					total_AC += float(alleleNs[k])
					alleles[k] = "" if alleles[k] == '-' else alleles[k]
					if alleles[k] == alt:
						altAF = alleleFs[k]
						hit_altAF = 1
				if hit_altAF == 0:
					return
			if float(altAF) > 0:
				altAF = "%.6f" %(float(altAF))
			if total_AC-int(total_AC)==0:
				total_AC=int(total_AC)
			rdbsnp_rsIDs[itms[3]] = {
				"AN"        : total_AC,
				"AF"        : altAF,
				"Alleles"   : itms[21],
				"class"     : itms[10],
				"func"      : itms[14],
				"weight"    : itms[16],
				"exception" : itms[17],
				"bitfields" : itms[24]
			}
		process_rec(n)
		for rec in q:
			process_rec(rec)
			
		return rdbsnp_rsIDs
	
	@staticmethod
	@LogExc
	def group_revcom(s):
		tmps = s
		m = re.finditer(r'([ATCGatcgRYKMSWBDHVrykmswbdhv]+)', s)
		for i in m:
			al = GetDBSNP.revcom(i.group(1))
			tmps=tmps[:i.start()]+al+tmps[i.end():]
		return tmps
	
	@staticmethod
	@LogExc
	def revcom(Seq):
		Seq=Seq[::-1]
		intab = "ATCGatcgRYKMSWBDHVrykmswbdhv"
		outtab = "TAGCtagcYRMKSWVHDByrmkswvhdb"
		trantab = maketrans(intab, outtab)
		Seq=Seq.translate(trantab)
		return Seq
