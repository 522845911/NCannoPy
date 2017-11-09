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
class GetGWAS(object):
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
	def getGWAS(self, chr, start, end):
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
			return dict()
		n=None
		try:
			n=next(q)
		except Exception, e:
			return dict()
		rgwas = dict()
		@LogExc
		def process_rec(rec):
			if not re.search(r'^\s*$', rec):
				rec = rec.strip('\r\n')
				itms = rec.split('\t', 6)
				if len(itms) < 6:
					raise Exception("Error: [{}] not enough columns [{}]. [{},{},{}].".format(self.db, chr, start, end))
				if not (itms[1] != start or itms[2] != end):
					if 'risk_allele' not in rgwas:
						rgwas['chr'] = chr
						rgwas['start'] = start
						rgwas['end'] = end
						rgwas['risk_alleles'] = list()
						rgwas['pmids'] = list()
						rgwas['traits'] = list()
					rgwas['risk_alleles'].append(itms[3])
					rgwas['pmids'].append(itms[4])
					rgwas['traits'].append(itms[5])
		process_rec(n)
		for rec in q:
			process_rec(rec)
		return rgwas