import gzip
import functools
import multiprocessing
import re
import copy
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
'''
=head1 SYNOPSIS

    # gene_tag can be entrez gene id or gene symbol
    # synopsis case 1: annotation
    use GetOmim;
    my $omim = GetOmim->new( "acdb" => $anno_combin_db );
    my $rAnnoComb = $omim->getAnnoComb( $gene_tag );

    # synopsis case 2: display
    use GetOmim;
    my $omim = GetOmim->new(
        "gndb" => $gene_mim_db,
        "mbdb" => $morbid_mim_db,
        "gvdb" => $gene_vars_db
    );
    my $rGeneMimInfo = $omim->getGeneMim( $gene_mim );
    my $rMorbMimInfo = $omim->getMorbMim( $morb_mim, $gene_mim );
    my $rGeneVarInfo = $omim->getGeneVar( $gene_mim );
    my $rSiglVarInfo = $omim->getGeneVar( $gene_mim, $varID );

=cut
'''
class GetOmim(object):
	@LogExc
	def __init__(self,opts):
		for k in opts:
			setattr(self,k,opts[k])
		if "acdb" in opts:
			self.set_acdb(opts["acdb"])
		if "gndb" in opts:
			self.set_gndb(opts["gndb"])
		if "mbdb" in opts:
			self.set_mbdb(opts["mbdb"])
		if "gvdb" in opts:
			self.set_gvdb(opts["gvdb"])
	@LogExc
	def set_acdb(self,dbf):
		setattr(self,'acdb',dbf)
		ACDB_h=None
		if re.search(r'\.gz$',dbf):
			try:
				ACDB_h = gzip.open(dbf,'rb')
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+", GunzipError"+e
		else:
			try:
				ACDB_h=open(dbf)
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+e
		
		# headers:
		# 1. gene symbol in annotation database
		# 2. Entrez Gene id
		# 3. Gene MIM id
		# 4. Gene Status
		# 5. Possible Inheritance
		# 6. Related disorder MIM id
		acdb = dict()
		for line in ACDB_h:
			line=line.strip("\r\n")
			itms = line.split('\t',6)
			if len(itms)<6:
				raise Exception("Error: acdb format error. ["+str(dbf)+"]")
			infos = dict()
			infos=dict(zip("genesym geneid genemim genestat inhs disomims".split(),itms))
			logger=multiprocessing.get_logger()
			if itms[1] not in acdb:
				acdb[itms[1]]=dict()
			if itms[0] not in acdb:
				acdb[itms[0]]=dict()
			acdb[itms[1]][itms[2]]=infos
			acdb[itms[0]][itms[2]]=infos
		ACDB_h.close()
		setattr(self,'acdb_h',copy.deepcopy(acdb))
	@LogExc
	def set_gndb(self,dbf):
		setattr(self,'gndb',dbf)
		GNDB_h=None
		if re.search(r'\.gz$',dbf):
			try:
				GNDB_h = gzip.open(dbf,'rb')
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+", GunzipError"+e
		else:
			try:
				GNDB_h=open(dbf)
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+e
		
		# headers:
		# 1.  gene mim id
		# 2.  entrez gene id
		# 3.  appoved symbol
		# 4.  gene symbols
		# 5.  gene status
		# 6.  all possible inheritance
		# 7.  all related disorders mims
		# 8.  description in genemap
		# 9.  main title in omim
		# 10. alternative titles in omim
		gndb = dict()
		for line in ACDB_h:
			line=line.strip("\r\n")
			itms = line.split('\t',maxsplit=10)
			if len(itms)<10:
				raise "Error: gndb format error. ["+str(dbf)+"]"
			infos = dict()
			infos=dict(zip("geneid appsym genesyms genestat inhs disomims desc main_title alt_titles".split(),itms[1:10]))
			gndb[itms[0]]=infos
		GNDB_h.close()
		setattr(self,'gndb_h',copy.deepcopy(gndb))
	@LogExc
	def set_mbdb(self,dbf):
		setattr(self,'mbdb',dbf)
		MBDB_h=None
		if re.search(r'\.gz$',dbf):
			try:
				MBDB_h = gzip.open(dbf,'rb')
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+", GunzipError"+e
		else:
			try:
				MBDB_h=open(dbf)
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+e
		
		# headers:
		# 1. disorder mim id
		# 2. gene mim id
		# 3. disorder class
		# 4. disorder map status
		# 5. all possible inheritances
		# 6. disorder description in morbidmap
		# 7. disorder main title
		# 8. disorder alternative titles
		mbdb = dict()
		for line in MBDB_h:
			line=line.strip("\r\n")
			itms = line.split('\t',maxsplit=9)
			if len(itms)<9:
				raise "Error: mbdb format error. ["+str(dbf)+"]"
			infos = dict()
			infos=dict(zip("disoclass disomapstat maploc inhs desc main_title alt_titles".split(),itms[2:9]))
			mbdb[itms[0]][itms[1]]=infos
		MBDB_h.close()
		setattr(self,'mbdb_h',copy.deepcopy(mbdb))
	@LogExc
	def set_gvdb(self,dbf):
		setattr(self,'gvdb',dbf)
		GVDB_h=None
		if re.search(r'\.gz$',dbf):
			try:
				GVDB_h = gzip.open(dbf,'rb')
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+", GunzipError"+e
		else:
			try:
				GVDB_h=open(dbf)
			except Exception,e:
				raise "Error: cannot open "+str(dbf)+e
		
		# headers:
		# 1. var_mim
		# 2. genesym
		# 3. pHGVS
		# 4. var_desc
		# 5. ext_desc (usually dbsnp id)
		# 6. var_titles
		gvdb = dict()
		for line in GVDB_h:
			line=line.strip("\r\n")
			itms = line.split('\t',maxsplit=6)
			if len(itms)<6:
				raise "Error: mbdb format error. ["+str(dbf)+"]"
			infos = dict()
			infos=dict(zip("varmim gsym pHGVS desc ext_desc titles".split(),itms))
			two_key = itms[0].split('\.')
			gvdb[two_key[0]][two_key[1]]=infos
		GVDB_h.close()
		setattr(self,'gvdb_h',copy.deepcopy(gvdb))
	@LogExc
	def getAnnoComb(self,gene_tag):
		if gene_tag not in self.acdb_h:
			return dict()
		else:
			return self.acdb_h[gene_tag]
	@LogExc
	def getGeneMim(self,genemim):
		if genemim is None or genemim not in self.gndb_h:
			return dict()
		else:
			return self.gndb_h[genemim]
	@LogExc
	def getMorbMim(self,morbmim, genemim):
		if morbmim is None \
		or genemim is None \
		or morbmim not in self.mbdb_h \
		or genemim not in self.mbdb_h[morbmim]:
			return dict()
		else:
			return self.mbdb_h[morbmim][genemim]
	@LogExc
	def getGeneVar(self,genemim,varid):
		if genemim is None or genemim not in self.gvdb_h:
			return dict()
		if varid is not None:
			if varid not in self.gvdb_h[genemim]:
				return dict()
			else:
				return self.gvdb_h[genemim][varid]
		return self.gvdb_h[genemim]