import pysam, re, os
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
class GetPfam(object):
	@LogExc
	def __init__(self,opts):
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getPfam(self, np, start, end):
		q=None
		try:
			q = self.tidb.fetch(str(np), int(start) - 1, int(end))
		except ValueError, e:
			return [".","."]
		n=None
		try:
			n=next(q)
		except Exception, e:
			return [".","."]
		pfam_hitted = dict()
		@LogExc
		def process_rec(rec):
			if re.search(r'^\s*$',rec):
				return
			rec=rec.strip("\r\n")
			itms = rec.split('\t',5)
			if len(itms)<5:
				raise Exception("Error: ["+str(self.db)+"] not enough columns. ["+str(chr)+","+str(start)+","+str(end)+"].")
			pfam_hitted[itms[3]]=itms[4]
		process_rec(n)
		for rec in q:
			process_rec(rec)
		pfam_ID=sorted(pfam_hitted.keys())
		if 0 == len(pfam_ID):
			return [ ".", "." ]
		else:
			return [";".join(pfam_ID),";".join([pfam_hitted[i] for i in pfam_ID])]