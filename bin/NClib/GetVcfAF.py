import os,pysam,sys
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
class GetVcfAF(object):
	@LogExc
	def __init__(self,opt):
		if "db" not in opt \
		or not os.access(opt['db'],os.F_OK) \
		or not os.access(opt['db'],os.R_OK):
			raise Exception("db option error! Please ensure that path is right and the absolute address!")
		if not os.access(opt['db']+".tbi",os.F_OK) \
		or not os.access(opt['db']+".tbi",os.R_OK):
			raise Exception("Error: ["+opt['db']+"] index (.tbi) file not found, please build it.")
		opt['tidb']=pysam.TabixFile(opt['db'])
		for k in opt:
			setattr(self,k,opt[k])
	@LogExc
	def getAF(self,*args):
		logger=multiprocessing.get_logger()
		if 5>len(args):
			raise Exception("Missing input parameters. need [chr, sta, sto, ref, alt]")
		chr, sta, sto, ref, alt=args
		
		if alt == '?':
			return dict()
		
		qsta, qsto  = sta, sto
		if sta == sto:
			qsta = sta - 1
			qsto = sto + 1
		q=None
		try:
			q = self.tidb.fetch(chr,qsta,qsto)
		except ValueError, e:
			return dict()
		item=None
		try:
			item=next(q)
		except Exception,e:
			return dict()
		
		results=list()
		
		total_AC, ret_AF = [0], [0]
		@LogExc
		def process_rec(rec):
			rec=rec.strip("\r\n")
			ar = rec.split('\t',6)
			if int(ar[1]) != int(sta) or int(ar[2]) != int(sto):
				return "next"
			#1  69427   69428   G,T	327,10343,  0.030647,0.969353,
			#1  1267324 1267326 GCC,G,GC        2025,22,10197,  0.165387,0.001797,0.832816,
			alleles = ar[3].split(',')
			AC      = ar[4].split(',')
			AF      = ar[5].split(',')
			alleles=[alleles[i] for i in range(len(alleles)) if alleles[i]]
			AC=[AC[i] for i in range(len(AC)) if AC[i]]
			AF=[AF[i] for i in range(len(AF)) if AF[i]]
			hit_ref = 0
			hit_alt = 0
			for i in range(len(alleles)):
				total_AC[0]+=float(AC[i])
				if str(alleles[i]) == '-' and str(alt) == '':
					ret_AF[0] = AF[i]
					hit_alt = 1
				if str(alleles[i]) == str(alt):
					ret_AF[0] = AF[i]
					hit_alt = 1
				if str(alleles[i]) == '-' and str(ref) == '':
					hit_ref = 1
				if str(alleles[i]) == str(ref):
					hit_ref = 1
			if hit_ref == 0 and str(ref) != '=' and not hasattr(self,'quiet'):
				print "Warning: ref unmatch ["+", ".join(chr, sta, sto, ref)+"]"+" at "+str(__file__)+" line "+str(sys._getframe().f_lineno)
			if hit_alt:
				return "last"
			return None
		ret_val=process_rec(item)
		if ret_val == "next" or ret_val is None:
			for rec in q:
				ret_val=process_rec(rec)
				if ret_val == "last":
					break
		total_AC=total_AC[0]
		ret_AF=ret_AF[0]
		if total_AC-int(total_AC)==0:
			total_AC=int(total_AC)
		if float(ret_AF) > 0:
			ret_AF="%.6f" %(float(ret_AF))
		if total_AC > 0:
			return {'AN':total_AC,'AF':ret_AF}
		else:
			return dict()
		