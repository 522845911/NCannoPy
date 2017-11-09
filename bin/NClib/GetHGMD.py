import os, pysam, sys
import multiprocessing
import functools
import re
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
class GetHGMD(object):
	@LogExc
	def __init__(self,opts):
		if "db" not in opts or not os.access(opts['db'],os.F_OK) or not os.access(opts['db'],os.R_OK):
			raise Exception("db option error!")
		if not os.access(opts['db']+".tbi",os.F_OK ) or not os.access(opts['db']+".tbi",os.R_OK ):
			raise Exception("Error: ["+str(opts['db'])+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self,k,opts[k])
	@LogExc
	def getHGMD(self, rquery):
		logger=multiprocessing.get_logger()
		# rquery is the hash ref
		if not type(rquery) is dict:
			if not hasattr(self,'quiet'):
				sys.stderr.write("check your input in 'getHGMD'!\n")
				return list()
		chr     = rquery['chr']
		start   = rquery['start']
		end     = rquery['end']
		ref,alt=None,None
		if 'ref' in rquery:
			ref=rquery['ref']
		if 'alt' in rquery:
			alt=rquery['alt']
		inchgvs=None
		if 'cHGVS' in rquery and not re.search(r'^\s*$',rquery['cHGVS']):
			inchgvs=rquery['cHGVS']
		inphgvs=None
		if 'pHGVS' in rquery and  not re.search(r'^\s*$',rquery['pHGVS']):
			inphgvs=rquery['pHGVS']
		if not (chr is not None and start is not None and end is not None):
			sys.stderr.write("check your chr postion in 'getHGMD'\n")
			return list()
		if not re.search(r'^chr',chr):
			chr = "{}{}".format("chr" , str(chr))
		qstart, qend  =  start, end
		if int(start)==int(end):
			qstart-=1
			qend+=1
		q=None
		try:
			q =self.tidb.fetch(chr, qstart, qend)
		except Exception,e:
			return list()
		item=None
		try:
			item=next(q)
		except Exception,e:
			return list()
		returns = list()
		@LogExc
		def process_rec(rec):
			if re.search(r'^\s*$',rec):
				return
			rec=rec.strip("\r\n")
			keya ="id mutation_name pmid disease pred".split()
			cols = rec.split('\t')
			qresult=dict()
			mutid,  geneid,  gene, nm, chgvs, \
			phgvs, disease,  pmid, pred = cols[3:12]
			rchr, rstart, rend=cols[:3]
			if nm is not None and not re.search(r'^\s*$',nm) and nm !=".":
				nm=nm
			else:
				nm=""
			if gene is not None and not re.search(r'^\s*$',gene) and gene !=".":
				gene="("+str(gene)+")"
			else:
				gene=""
			if chgvs is not None and not re.search(r'^\s*$',chgvs) and chgvs !=".":
				chgvs=chgvs
			else:
				chgvs=""
			if phgvs is not None and not re.search(r'^\s*$',phgvs) and phgvs !=".":
				phgvs=" ("+str(phgvs)+")"
			else:
				phgvs=""
			mutname = "{}{}{}{}{}".format(nm,gene,": ",chgvs,phgvs)
			qresult=dict(zip(keya,[mutid, mutname, pmid, disease, pred]))
			if inchgvs is not None and inchgvs== chgvs:
				returns.append(qresult)
			elif inchgvs is not None \
			and re.search(r'dup',inchgvs) \
			and not re.serach(r'del',inchgvs) \
			and not re.serach(r'^\s*$',chgvs) \
			and re.serach(r'ins',chgvs) \
			and not re.serach(r'del',chgvs):
				tempc = inchgvs
				ins1, ins2 =  "", ""
				dup2ins=self.dup2ins(inchgvs)
				if dup2ins is not None and len(dup2ins)==2:
					ins1, ins2  = dup2ins[:2]
				if (ins1 is not None and not re.search(r'^\s*$',ins1) and ins1==chgvs) \
					or (not re.search(r'^\s*$',ins2) and ins2 == chgvs):
					returns.append(qresult)
			elif not(inchgvs is not None or inphgvs is not None):
				returns.append(qresult)
			
		process_rec(item)
		for rec in q:
			process_rec(rec)
		if len(returns)==0:
				return list()
		return returns





	@LogExc
	def dup2ins( self, chgvs ):
		if not(chgvs is not None and re.search(r'dup',chgvs)):
			sys.stderr.write("cHGVS is not 'dup'\n")
			return
		dup1, dup2 =  "", ""
		m=re.serach(r'c\.([\*\-]?)(\d+)([\+\-]\d+)?_([\*\-]?)(\d+)([\+\-]\d+)?dup(.*)',chgvs)
		m2=re.search(r'c\.([\*\-]?)(\d+)([\+\-]?)(\d+?)dup(.*)',chgvs)
		if m:
			c1 = m.group(1) if m.group(1) is not None else ""
			c2 = m.group(2) if m.group(2) is not None else ""
			c3,c4 = "",""
			tf =None
			if m.group(3) is not None:
				tf=m.group(3)
			c5 = m.group(4) if m.group(4) is not None else ""
			c6 = m.group(5) if m.group(5) is not None else ""
			c7,c8="",""
			tb =None
			if m.group(6) is not None:
				tb =m.group(6)
			c9 = m.group(7) if m.group(7) is not None else ""
			m=re.search(r'([\+\-])(\d+)',tf)
			if tf is not None and not re.search(r'^\s*$',tf) and m:
				c3 = m.group(1) if m.group(1) is not None  else ""
				c4 = m.group(2) if m.group(2) is not None  else ""
			m=re.search(r'([\+\-])(\d+)',tb)
			if tb is not None or not re.search(r'^\s*$',tb) and m:
				c7 = m.group(1) if m.group(1) is not None  else ""
				c8 = m.group(2) if m.group(2) is not None  else ""
			if not re.search(r'^\s*$',c4):
				t = "{}{}".format(c3 ,int(c4)-1)
				if int(c4)-1==0:
					t = ""
				dup1 = "{}{}{}{}{}{}{}{}{}{}{}".format("c\.",c1,c2,t,"\_",c1,c2,c3,c4,"ins",c9)
			else:
				t = int(c2) - 1
				dup1 = "{}{}{}{}{}{}{}{}".format("c\." , c1 , t , "\_" , c1 , c2 , "ins" , c9)
			if not re.search(r'^\s*$',c8):
				dup2 ="{}{}{}{}{}{}{}{}{}{}{}{}".format("c\.",c5,c6,c7,c8 , "\_",c5,c6,c7,int(c8) + 1 ,"ins",c9)
			else:
				dup2 = "{}{}{}{}{}{}{}{}".format("c\." , c5 , c6 , "\_" , c5 , int(c6) + 1  , "ins" , c9)
		elif m2:
			c1 = m2.group(1) if m2.group(1) is not None else ""
			c2 = m2.group(2) if m2.group(2) is not None else ""
			c3 = m2.group(3) if m2.group(3) is not None else ""
			c4 = m2.group(4) if m2.group(4) is not None else ""
			c5 = m2.group(5) if m2.group(5) is not None else ""
			if not re.search(r'^\s*$',c4):
				t = "{}{}".format(c3,int(c4)-1)
				if int(c4)-1==0:
					t = ""
				dup1 ="{}{}{}{}{}{}{}{}{}{}{}".fortmat("c\.",c1,c2,t,"\_",c1,c2,c3,c4,"ins",c5)
				dup2 ="{}{}{}{}{}{}{}{}{}{}{}{}".format("c\.",c1,c2,c3,c4,"\_",c1,c2,c3,int(c4) + 1,"ins",c5)
			else:
				t = int(c2) - 1
				if t == 0:
					t    = ""
				dup1 = "{}{}{}{}{}{}{}{}".format("c\." ,c1 ,t , "\_" ,c1 ,c2 , "ins" ,c5)
				dup2 = "{}{}{}{}{}{}{}{}".format("c\." ,c1 ,c2 , "\_" ,c1, int(c2) + 1  , "ins" ,c5)
		return [ dup1, dup2 ]