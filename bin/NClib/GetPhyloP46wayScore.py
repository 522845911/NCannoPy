import pysam
import re
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
class GetPhyloP46wayScore(object):
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
	def getPhyloP46wayScore(self, chr, pos):
		chr=re.sub(r'^chr','',chr,flags=re.IGNORECASE)
		if not re.search(r'^\d+$',str(pos)):
			raise Exception("Error: pos_string ["+str(pos)+"] should be a position number!")
		q=None
		try:
			q = self.tidb.fetch(chr, pos - 1, pos)
		except Exception, e:
			return [ ".", ".", "." ]
		n=None
		try:
			n=next(q)
		except Exception, e:
			return [ ".", ".", "." ]
		ptmScore, prmScore, vtbScore=[None], [None], [None]
		@LogExc
		def process_rec(rec):
			rec=re.sub(r'\s+$','',rec)
			rec=re.sub(r'\r\n$','',rec)
			column = rec.split("\t",5)
			if len(column) < 5:
				raise Exception("Error: ["+str(self.db)+"] not enough columns.")
			ptmScore[0] = column[2] if not re.search(r'^\s*$',column[2]) else "."
			prmScore[0] = column[3] if not re.search(r'^\s*$',column[3]) else "."
			vtbScore[0] = column[4] if not re.search(r'^\s*$',column[4]) else "."
			return "last"
		ret_val=process_rec(n)
		if ret_val !="last":
			for rec in q:
				ret_val=process_rec(rec)
				if ret_val =="last":
					break
		if ptmScore[0] is None:
			return [ ".", ".", "." ]
		else:
			return [ ptmScore[0], prmScore[0], vtbScore[0] ]
	