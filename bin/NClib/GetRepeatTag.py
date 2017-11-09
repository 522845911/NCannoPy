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
class GetRepeatTag():
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
	def getRepTag(self, chr, start, end):
		if not re.search(r'^chr', chr, flags=re.IGNORECASE):
			chr='chr{}'.format(chr)
		qstart, qend = start, end
		if start == end:
			qstart -= 1
			qend += 1
		q=None
		try:
			q = self.tidb.fetch(chr, qstart, qend)
		except Exception, e:
			return "."
		n=None
		try:
			n=next(q)
		except Exception, e:
			return "."
		all_rmsks = list()
		rec=n
		@LogExc
		def process_rec(rec):
			if not re.search(r'^\s*$', rec):
				rec = rec.strip("\r\n")
				itms = rec.split('\t', 4)
				if len(itms) < 4:
					raise Exception("Error: [{}] not enough columns. [{},{},{}].".format(self.db, chr, start, end))
				all_rmsks.append(itms[3])
		process_rec(n)
		for rec in q:
			process_rec(rec)
		if 0 == len(all_rmsks):
			return "."
		return ";".join(all_rmsks)