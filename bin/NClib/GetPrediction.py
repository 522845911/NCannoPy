import os, pysam, re
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
class GetPrediction(object):
	ALL_AAS   = "A C D E F G H I K L M N P Q R S T V W Y".split()
	NUM_AAS   = len(ALL_AAS)
	AA_LOOKUP = dict([(ALL_AAS[i],i) for i in range(len(ALL_AAS))])
	PRED_GROUP = "sift polyphen2_humvar polyphen2_humdiv".split()
	NUM_PREDS = len(PRED_GROUP)
	
	VAL_TO_PREDICTION={
		"polyphen2_humdiv":{
			0:'Probably damaging',
			1:'Possibly damaging',
			2:'Benign',
			3:'unknown'
		},
		"polyphen2_humvar":{
			0:'Probably damaging',
			1:'Possibly damaging',
			2:'Benign',
			3:'unknown'
		},
		"sift":{
			0:'Tolerated',
			1:'Damaging'
		}
	}
	@LogExc
	def __init__(self,opts):
		logger=multiprocessing.get_logger()
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getPredScore(self, Pid, Pos, Alt):
		logger=multiprocessing.get_logger()
		Alt = Alt.upper()
		if Alt not in GetPrediction.AA_LOOKUP:
			return dict()
		q=None
		try:
			q = self.tidb.fetch(str(Pid), int(Pos)-1, int(Pos))
		except ValueError, e:
			return dict()
		n=None
		try:
			n=next(q)
		except Exception, e:
			return dict()
		@LogExc
		def process_rec(rec):
			logger=multiprocessing.get_logger()
			rec=rec.strip("\r\n")
			preds = rec.split('\t',3)[2]
			preds=re.sub(r'\s+','',preds)
			aa_preds = preds.split(',',GetPrediction.NUM_AAS)[GetPrediction.AA_LOOKUP[Alt]]
			if aa_preds == '.':
				return dict()
			aa_pred_grp = aa_preds.split(';',GetPrediction.NUM_PREDS)
			if len(aa_pred_grp) != GetPrediction.NUM_PREDS:
				if len(aa_pred_grp)<GetPrediction.NUM_PREDS:
					raise Exception("Error: ["+str(self.db)+"] format error, not enough prediction group.\n"+"       [pid:"+str(Pid)+", pos:"+str(Pos)+", aa:"+str(Alt)+", preds:"+str(aa_preds)+"]\n")
				elif not hasattr(self,"quiet"):
					print "Warnings: ["+self.db+"] format not match, database have more prediciton groups.\n"+"       [pid:"+str(Pid)+", pos:"+str(Pos)+", aa:"+str(Alt)+", preds:"+str(aa_preds)+"]\n"
			tmpPreds = dict([(GetPrediction.PRED_GROUP[i],aa_pred_grp[i]) for i in range(len(GetPrediction.PRED_GROUP))])
			ret_preds = dict()
			for pred in tmpPreds.keys():
				if tmpPreds[pred] != '.':
					predval, score = tmpPreds[pred].split('|')[:2]
					if int(predval) not in GetPrediction.VAL_TO_PREDICTION[pred]:
						raise Exception( "Error: ["+self.db+"] format error, can not recognize the prediction value.\n"+"       [pid:"+str(Pid)+", pos:"+str(Pos)+", aa:"+str(Alt)+", preds:"+str(pred)+", val:"+str(predval)+"]\n")
					ret_preds[pred]=[GetPrediction.VAL_TO_PREDICTION[pred][int(predval)],score]
			return ret_preds
		ret_val=process_rec(n)
		return ret_val
		for rec in q:
			process_rec(rec)
		return dict()
