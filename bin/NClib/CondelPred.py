from __future__ import division
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

class CondelPred(object):
	@LogExc
	def __init__(self, args):
		logger=multiprocessing.get_logger()
		self.config = dict()
		self.sift = dict()
		self.polyphen = dict()
		self.config['condel_dir'] = args
		config_file = self.config['condel_dir']+"config/condel_SP.conf"
		CONF=None
		try:
			CONF = open(config_file)
		except Exception, e:
			raise Exception("Could not open"+config_file)
		conf = CONF.readlines()
		CONF.close()
		
		for i in range(len(conf)):
			m = re.search(r"(cutoff\.HumVar\.\w+)=\'(\S+)\'",conf[i])
			m2=re.search(r"(max\.HumVar\.\w+)=\'(\S+)\'",conf[i])
			if m:
				self.config[m.group(1)] = m.group(2)
			if m2:
				self.config[m2.group(1)] = m2.group(2)
		if len(self.config.keys()) < 4:
			raise Exception("bad config file! please check")
		self.sift = {"tp":dict(),"tn":dict()}
		self.polyphen = {"tp":dict(),"tn":dict()}
		SIFT=None
		try:
			SIFT = open(self.config['condel_dir']+"/methdist/sift.data")
		except Exception, e:
			raise Exception("no sift.data in "+self.config['condel_dir']+"/methdist/")
		sift = SIFT.readlines()
		SIFT.close()
		for i in range(len(sift)):
			m=re.search(r'(\S+)\s+(\S+)\s+(\S+)',sift[i])
			if m:
				self.sift['tp'][m.group(1)] = m.group(2)
				self.sift['tn'][m.group(1)] = m.group(3)
		POLYPHEN=None
		try:
			POLYPHEN=open(self.config['condel_dir']+"/methdist/polyphen.data")
		except Exception, e:
			raise Exception("no polyphen.data in "+self.config['condel_dir']+"/methdist/")
		polyphen=POLYPHEN.readlines()
		POLYPHEN.close()
		for i in range(len(polyphen)):
			m=re.search(r'(\S+)\s+(\S+)\s+(\S+)',polyphen[i]) 
			if m:
				self.polyphen['tp'][m.group(1)]=m.group(2)
				self.polyphen['tn'][m.group(1)]=m.group(3)
	@LogExc
	def pred(self, tva):
		logger=multiprocessing.get_logger()
		pph_score, pph_pred, sift_score=[None] * 3
		if 'pp2varScore' in tva:
			pph_score   = tva['pp2varScore']
		if 'pp2varPred' in tva:
			pph_pred   = tva['pp2varPred']
		if 'siftScore' in tva:
			sift_score   = tva['siftScore']
		
		if pph_score is not None and sift_score is not None and str(pph_pred) != 'unknown':
			condel_pred, condel_score = self.compute_condel(sift_score, pph_score)
			condel_score = "%.3f" %(condel_score)
			return {"condelPred":condel_pred,"condelScore":condel_score}
		else:
			return dict()
	@LogExc
	def compute_condel(self, sift_score, polyphen_score):
		USE_V2 = 1
		base = 0
		int_score = 0
		sift_score     = "%.3f" %(float(sift_score))
		polyphen_score = "%.3f" %(float(polyphen_score))
		if float(sift_score) <= float(self.config['cutoff.HumVar.sift']):
			int_score+= float("%.3f" %((1-float(sift_score)/float(self.config['max.HumVar.sift']))*(1-float(self.sift['tn'][sift_score]))))
			tmp = None
			if USE_V2:
				tmp = 1
			else:
				tmp = 1-float(self.sift['tn'][sift_score])
			base += tmp
		else:
			int_score += float("%.3f" %((1-float(sift_score)/float(self.config['max.HumVar.sift']))*(1-float(self.sift['tp'][sift_score]))))
			tmp = None
			if USE_V2:
				tmp = 1
			else:
				tmp = 1-float(self.sift['tp'][sift_score])
			base += tmp
		if float(polyphen_score) >= float(self.config['cutoff.HumVar.polyphen']):
			int_score += float("%.3f" %(float(polyphen_score)/float(self.config['max.HumVar.polyphen'])*( 1 - float(self.polyphen['tn'][polyphen_score]))))
			tmp = None
			if USE_V2:
				tmp = 1
			else:
				tmp = 1 - float(self.polyphen['tn'][polyphen_score])
			base += tmp
		else:
			int_score += float("%.3f" %(float(polyphen_score) / float(self.config['max.HumVar.polyphen'])*(1-float(self.polyphen['tp'][polyphen_score]))))
			tmp = None
			if USE_V2:
				tmp = 1
			else:
				tmp = 1- float(self.polyphen['tp'][polyphen_score])
			base += tmp
		pred = None
		if base == 0:
			int_score = -1
			pred = 'not_computable_was'
		else:
			int_score = "%.3f" %(int_score/base)
			int_score=float(int_score)
		if int_score >= 0.469:
			pred = 'deleterious'
		elif int_score >= 0  and int_score < 0.469:
			pred = 'neutral'
		return [pred, int_score]