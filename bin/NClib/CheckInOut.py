import os, pysam, multiprocessing, re
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
class CheckInOut(object):
	@LogExc
	def __init__(self,opt):
		logger=multiprocessing.get_logger()
		if "db" not in opt or not os.access(opt["db"],os.F_OK) or not os.access(opt["db"],os.R_OK):
			raise Exception("db option error! Please ensure that path is right and the absolute address!")
		if not os.access(opt['db']+".tbi",os.F_OK) or not os.access(opt['db']+".tbi",os.R_OK):
			raise Exception("Error: ["+opt['db']+"] index (.tbi) file not found, please build it.")
		opt['tidb'] = pysam.TabixFile(opt['db'])
		for k in opt:
			setattr(self,k,opt[k])
	@LogExc
	def checkIn(self,*args):
		if 3>len(args):
			raise Exception("Missing input parameters. need [chr, sta, sto]")
		chr, sta, sto=args[:3]
		if not re.search(r'^chr',chr):
			chr="chr"+chr
		
		qsta, qsto=sta, sto
		if sta == sto:
			qsta = sta - 1
			qsto = sto + 1
		
		q=None
		try:
			q =self.tidb.fetch(chr, qsta, qsto)
		except Exception, e:
			return 0
		else:
			no_hit = 0
			for ret in q:
				if not re.search(r'^\s*$',ret):
					no_hit = 1
					break
			return no_hit

