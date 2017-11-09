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
class GetCOSMIC(object):
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
	def getCOSMIC(self, chr, start, end, ref, alt):
		logger=multiprocessing.get_logger()
		if not re.search(r'^chr',chr, flags=re.IGNORECASE):
			chr = 'chr' + str(chr)
		if ref != "=":
			ref = ref.upper()
		alt=alt.upper()
		alt_list=[alt]
		ref_list=[ref]
		q=None
		try:
			if start == end:
				q=self.tidb.fetch(chr, start - 1, end + 1)
			else:
				q=self.tidb.fetch(chr, start , end)
		except Exception, e:
			return dict()
		n=None
		try:
			n=next(q)
		except Exception, e:
			return dict()
		results=list()
		@LogExc
		def process_rec(rec):
			rec=rec.strip("\r\n")
			ar = rec.split('\t', 11)
			if alt_list[0] =="" or alt_list[0]==".":
				alt_list[0] = '-'
			if int(start) != int(ar[1]) or int(end) != int(ar[2]) or str(alt_list[0]) != str(ar[4]):
				return "next"
			if ref_list[0] is not None and not hasattr(self,'quiet'):
				if ref_list[0]=="" or ref_list[0] == ".":
					ref_list[0] = '-'
				if ref_list[0] != '=' and ref_list[0] != ar[3]:
					print "Warning: reference unmatch."
			ret = dict()
			ret = dict(zip("mutID gene mutName site histology status".split() ,[ar[i] for i in range(5,11)]))
			return ret
		ret_val=process_rec(n)
		if ret_val=="next":
			for rec in q:
				ret_val=process_rec(rec)
				if ret_val!="next":
					return ret_val
		else:
			return ret_val
		return dict()