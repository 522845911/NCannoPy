from __future__ import division
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
class GetExAC(object):
	CSQ_headers = '''
	Allele Gene Feature Feature_type Consequence cDNA_position CDS_position 
  Protein_position Amino_acids Codons Existing_variation AA_MAF EA_MAF 
  ALLELE_NUM EXON INTRON MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE 
  DISTANCE STRAND CLIN_SIG CANONICAL SYMBOL SYMBOL_SOURCE SIFT PolyPhen GMAF 
  BIOTYPE ENSP DOMAINS CCDS HGVSc HGVSp AFR_MAF AMR_MAF ASN_MAF EUR_MAF PUBMED 
  LoF_info LoF_flags LoF_filter LoF
	'''.split()
	Pops = "AFR AMR EAS EUR FIN NFE OTH SAS".split()
	@LogExc
	def __init__(self,opts):
		if "db" not in opts or not os.access(opts["db"],os.F_OK) or not os.access(opts["db"],os.R_OK):
			raise Exception("db option error! Please ensure that path is right and the absolute address!")
		if not os.access(opts["db"]+".tbi",os.F_OK) or not os.access(opts["db"]+".tbi",os.R_OK):
			raise Exception("Error: ["+opts["db"]+"] index (.tbi) file not found, please build it.")
		opts['tidb']=pysam.TabixFile(opts['db'])
		for k in opts:
			setattr(self, k, opts[k])
	@LogExc
	def getAF(self, *args):
		logger=multiprocessing.get_logger()
		if len(args) < 5:
			raise Exception("Missing input parameters. need [chr, sta, sto, ref, alt]")
		chr, sta, sto, ref, alt = args[:5]
		if str(alt) == '?':
			return dict()
		chr = re.sub(r'^chr','', chr, flags=re.IGNORECASE)
		qsta, qsto  = sta, sto
		if sta == sto:
			qsta = sta - 1
			qsto = sto + 1
		q=None
		try:
			q = self.tidb.fetch(chr, qsta, qsto)
		except Exception,e:
			return dict()
		item=None
		try:
			item = next(q)
		except Exception,e:
			return dict()
		
		results = []
		@LogExc
		def process_rec(rec):
			logger=multiprocessing.get_logger()
			rec = rec.strip("\r\n")
			ar = rec.split('\t', 8)
			alts = ar[4].split(',')
			cur_rec_hit = GetExAC.vcfhit(sta, sto, ref, alt, ar[1], ar[3], alts)
			if int(cur_rec_hit) < 0:
				return 
			INFO_hash = GetExAC.parse_vcfINFO(ar[7])
			infos =dict()
			infos['INFO'] = INFO_hash
			infos['dbSNPid'] = ar[2]
			infos['Qual'] = ar[5]
			infos['Filter'] = ar[6]
			
			AN, rAC = None, None
			if 'AN_Adj' in infos['INFO'] and 'AC_Adj' in infos['INFO']:
				AN  = infos['INFO']['AN_Adj']
				rAC = infos['INFO']['AC_Adj']
			elif 'AN' in infos['INFO'] and 'AC' in infos['INFO']:
				AN  = infos['INFO']['AN']
				rAC = infos['INFO']['AC']
			if AN is None:
				if not hasattr(self, 'quiet'):
					print "No AN information available for "+"{}:{}-{}({}>{})".format(chr, sta, sto, ref, alt)
				return 
			else:
				if int(cur_rec_hit) == 0:
					infos['AN'] = AN
					other_allele = 0
					for othac in rAC:
						other_allele += int(othac)
					ref_AC = int(AN) - other_allele
					if int(AN) == 0:
						infos['AF'] = 0
					else:
						infos['AF'] = "%.6f" %((int(AN) - other_allele) / int(AN))
					if 'AC_Het' in infos['INFO']:
						other_allele = 0
						for i in range(len(rAC)):
							other_allele += int(infos['INFO']['AC_Het'][i])
						ref_HASC = ref_AC - other_allele
						if int(ref_HASC) % 2 == 0 and not re.search(r'[XY]', chr):
							infos['HASC'] = ref_HASC / 2
						elif re.search(r'[XY]', chr):
							infos['HASC'] = ref_HASC
						else:
							if not hasattr(self, 'quiet'):
								print "Warning: non-even number "+"[Hom Allele Count] for refcall of "+"{}:{}-{}".format(chr, sta, sto)
					for pop in GetExAC.Pops:
						if "AN_" + pop in infos['INFO'] and "AC_" + pop in infos['INFO']:
							infos["AN_" + pop] = infos['INFO']["AN_" + pop]
							other_allele = 0
							for othac in infos['INFO']["AC_" + pop]:
								other_allele += int(othac)
							if int(infos["AN_" + pop]) == 0:
								infos["AF_" + pop] = 0
							else:
								infos["AF_" + pop] = "%.6f" %((int(infos["AN_" + pop]) - other_allele) / int(infos["AN_" + pop]))
							if "Het_"+pop in infos['INFO']:
								other_allele = 0
								for i in range(len(rAC)):
									other_allele += int(infos['INFO']["Het_"+pop][i])
								if "AC_"+pop not in infos:
									pop_ref_HASC = 0 - other_allele
								else:
									pop_ref_HASC = int(infos["AC_"+pop]) - other_allele
								if int(pop_ref_HASC) % 2 == 0 and not re.search('[XY]', chr):
									infos["HASC_"+pop] = int(pop_ref_HASC) / 2
								elif re.search(r'[XY]', chr):
									infos["HASC_"+pop] = pop_ref_HASC
								else:
									if not hasattr(self,"quiet"):
										print "Warning: non-even number" +" ["+str(pop)+" Hom Allele Count] for refcall of "+"{}:{}-{}".format(chr, sta, sto)
				else:
					hitted = cur_rec_hit - 1
					infos['AN'] = AN
					if int(AN) == 0:
						infos['AF'] = 0
					else:
						infos['AF'] = "%.6f" %(int(rAC[hitted])/int(AN))
					if 'AC_Hom' in infos['INFO']:
						infos['HASC'] = infos['INFO']['AC_Hom'][hitted]
					for pop in GetExAC.Pops:
						if "AN_" + pop in infos['INFO'] and "AC_" + pop in infos['INFO']:
							infos["AN_" + pop] = infos['INFO']["AN_" + pop]
							if int(infos["AN_" + pop]) == 0:
								infos["AF_" + pop] = 0
							else:
								infos["AF_" + pop] = "%.6f" %(int(infos['INFO']["AC_" + pop][hitted])/int(infos["AN_" + pop]))
							if "Hom_" + pop in infos['INFO']:
								infos["HASC_" + pop] = infos['INFO']["Hom_" + pop][hitted]
			
			results.append(infos)
		process_rec(item)
		for rec in q:
			process_rec(rec)
		if 0 == len(results):
			return dict()
		elif 1 == len(results):
			return results[0]
		else:
			if not hasattr(self, "quiet"):
				print "Multiple record for "+"{}:{}-{}({}>{})".format(chr, sta, sto, ref, alt)+"only return the last one."
			return results[-1]
	@staticmethod
	@LogExc
	def parse_vcfINFO(arg):
		logger=multiprocessing.get_logger()
		all_infos = arg.split(';')
		info_hash = dict()
		for inf in all_infos:
			m = re.search(r'^(\S+)=(.+)$', inf)
			if m:
				tag = m.group(1)
				content = m.group(2)
				if re.search(r'^AC|^AF|^Het_|^Hom_|^CSQ$|_HIST$', tag):
					info_hash[tag] = content.split(',')
				else:
					info_hash[tag] = content
				if tag == "CSQ":
					for i in range(len(info_hash[tag])):
						sub_content = info_hash[tag][i]
						csq_hash =dict(zip(GetExAC.CSQ_headers, sub_content.split('|',len(GetExAC.CSQ_headers))))
						info_hash[tag][i]=csq_hash
			else:
				info_hash[inf]=1
		
		if "AN_FIN" in info_hash and "AN_NFE" in info_hash:
			info_hash['AN_EUR']=int(info_hash['AN_FIN'])+int(info_hash['AN_NFE'])
			if "AC_FIN" in info_hash and "AC_NFE" in info_hash:
				alta_count = len(info_hash['AC_FIN'])
				info_hash['AC_EUR'] = [0]*alta_count
				for k in range(alta_count):
					info_hash['AC_EUR'][k]=int(info_hash['AC_EUR'][k])+int(info_hash['AC_FIN'][k])+int(info_hash['AC_NFE'][k])
			if 'Hom_FIN' in info_hash and "Hom_NFE" in info_hash:
				homa_count = len(info_hash['Hom_FIN'])
				info_hash['Hom_EUR'] = [0]*homa_count
				for k in range(homa_count):
					info_hash['Hom_EUR'][k] =int(info_hash['Hom_EUR'][k]) +int(info_hash['Hom_FIN'][k])+int(info_hash['Hom_NFE'][k])
			if 'Het_FIN' in info_hash and 'Het_NFE' in info_hash:
				heta_count = len(info_hash['Het_FIN'])
				info_hash['Het_EUR']=[0]*heta_count
				for k in range(heta_count):
					info_hash['Het_EUR'][k]=int(info_hash['Het_EUR'][k])+int(info_hash['Het_FIN'][k])+int(info_hash['Het_NFE'][k])
		return info_hash
	
	@staticmethod
	@LogExc
	def vcfhit(sta, sto, ref, alt, pos_in_vcf, ref_in_vcf, ralts_in_vcf):
		logger=multiprocessing.get_logger()
		pos_in_vcf=int(pos_in_vcf)
		len_vcf_ref = len(ref_in_vcf)
		end_vcf_ref = pos_in_vcf + len_vcf_ref - 1
		if sta < pos_in_vcf - 1 and sto >= pos_in_vcf and ref[:pos_in_vcf - 1 - sta] == alt[:pos_in_vcf - 1 - sta]:
			sta = pos_in_vcf - 1
			ref = ref[sta - pos_in_vcf -1:]
			alt = alt[sta - pos_in_vcf -1:]
		if sto > end_vcf_ref and sta <= end_vcf_ref and ref[end_vcf_ref - sto:] == alt[end_vcf_ref - sto:]:
			sto = end_vcf_ref
			ref = ref[:end_vcf_ref - sta]
			alt = alt[:len(alt)-(sto - end_vcf_ref)]
		if sta > end_vcf_ref or (sta < pos_in_vcf - 1) or sto > end_vcf_ref or sto < pos_in_vcf:
			return -1
		
		check_ref, check_alt = ref_in_vcf, ref_in_vcf
		check_alt =check_alt[:sta - pos_in_vcf + 1]+alt+check_alt[sto - pos_in_vcf+1:]
		if check_ref == ref and ref == alt:
			return 0
		for i in range(len(ralts_in_vcf)):
			if ralts_in_vcf[i] == check_alt:
				return i + 1
		return -1