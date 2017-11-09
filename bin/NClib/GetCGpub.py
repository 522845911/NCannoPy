from __future__ import division
import pysam
import os
import re
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
class GetCGpub(object):
	@LogExc
	def __init__(self, opts):
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error! Please ensure that path is right and the absolute address!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getAF(self, *args):
		if len(args) <5:
			raise Exception("Missing input parameters. need [chr, sta, sto, ref, alt]")
		chr, sta, sto, ref, alt = args[:5]
		if alt == '?':
			return dict()
		if not re.search(r'^chr', chr, flags=re.IGNORECASE):
			chr='chr{}'.format(chr)
		qsta, qsto = sta, sto
		if sta == sto:
			qsta -= 1
			qsto += 1
		query=None
		try:
			query = self.tidb.fetch(chr, qsta, qsto)
		except Exception,e:
			return dict()
		n=None
		try:
			n=next(query)
		except Exception, e:
			return dict()
		total_AC, ret_AF = [0], [0]
		def process_rec(rec):
			ti_read = ti_read.strip('\r\n')
			ar = ti_read.split('\t',15)
			if ar[2] != sta or ar[3] != sto:
				return "continue"
			alleles = [ar[i] for i in [5, 6]]
			AC      = [ar[i] for i in [11, 12]]
			hit_ref = 0
			hit_alt = 0
			for i in range(len(alleles)):
				total_AC[0] += float(AC[i])
				if alleles[i] == alt:
					ret_AF[0] = AC[i]
					hit_alt = 1
				if alleles[i] == ref:
					hit_ref = 1
			if hit_ref == 0 and ref != '=' and not hasattr(self, 'quiet'):
				print "Warning: ref unmatch [{}, {}, {}, {}]".format(chr, sta, sto, ref)
			return "break"
		ret_val=process_rec(n)
		if ret_val=="continue":
			for ti_read in query:
				ret_val=process_rec(ti_read)
				if ret_val=="break":
					break
		total_AC=total_AC[0]
		ret_AF  =ret_AF[0]
		if total_AC-int(total_AC)==0:
			total_AC=int(total_AC)
		if total_AC > 0:
			ret_AF = "%.6f" %(ret_AF/total_AC)
			return {"AN":total_AC,"AF":ret_AF}
		else:
			return dict()