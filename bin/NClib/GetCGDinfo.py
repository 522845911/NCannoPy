import re, gzip
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
class GetCGDinfo(object):
	@LogExc
	def __init__(self, opts):
		logger=multiprocessing.get_logger()
		for k in opts:
			setattr(self,k,opts[k])
		if "db" not in opts:
			raise Exception("Error: need db path.")
		CGD_h = None
		if re.search(r'\.gz$',opts['db']):
			try:
				CGD_h = gzip.open(opts['db'],'rb')
			except Exception, e:
				raise Exception("Error: cannot open "+str(opts["db"])+", GunzipError"+e)
		else:
			try:
				CGD_h = open(opts['db'])
			except Exception, e:
				raise Exception("Error: cannot open "+str(opts["db"])+e)
		
		# headers:
		# 1.  gene                     
		# 2.  entrez_gene_id           
		# 3.  condition                
		# 4.  inheritance              
		# 5.  age_group                
		# 6.  allelic_conditions       
		# 7.  manifestation_categories 
		# 8.  intervention_categories  
		# 9.  comments                 
		# 10. intervention             
		# 11. references        
		
		cgd = dict()
		for line in CGD_h:
			if re.search(r'^#',line):
				continue
			line=line.strip("\r\n")
			itms = line.split('\t', 11)
			if len(itms) < 11:
				raise Exception("Error: acdb format error. ["+str(opts["db"])+"]")
			infos = dict(zip("genesym geneid cond inhs age_grp allelic_cond manifest_cat intervent_cat comments intervent references".split(),[itms[i] for i in range(11)]))
			cgd[itms[0]] = infos
			cgd[itms[1]] = infos
		CGD_h.close()
		setattr(self, "cgd_h", cgd)
	@LogExc
	def getCGD(self, gene_tag):
		logger=multiprocessing.get_logger()
		if gene_tag not in self.cgd_h:
			return dict()
		else:
			return self.cgd_h[gene_tag]