args <- commandArgs(trailingOnly = TRUE)
f = as.numeric(args[1])
chrnum = as.numeric(args[2])
nsites = as.numeric(args[3])

# Extract candidate recombination locations in haplo data
maxw = as.numeric(args[4])

chr=paste0("chr",chrnum)
# Number of replicates of haplotype sampling
nr=10
# Window width in number of variants
ww = 80
# Lag between windows
decal = 8

path = "/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/recombinations/BG/"
pathIBD = "/dcl01/mathias1/data/TOPMED_WGS/TOPMed_Genomes_Freeze8/temp/OmniExpress/freeze.8."

# Function to recode haplotypes into consecutive numbers in order of appearance of the haplotypes
recodeHap = function(vec)
{
	lev = unique(vec)
	fac = factor(vec,levels=lev)
	unclass(fac)
}

# Function to turn homozygous subjects into heterozygous subjects by assigning them a new allele nlev + 1
het_subject = function(vec)
{
	nlev = length(unique(vec))
	vec[seq(2,length(vec),by=2)] = ifelse(vec[seq(2,length(vec),by=2)]==vec[seq(1,length(vec)-1,by=2)],nlev+1,vec[seq(2,length(vec),by=2)])
	vec
}
common.allele = function(all1, all2)
{
    p1 = which(!is.na(all1))
	if (length(p1)>0) c1 = sum(all1[p1[1]] == all1[-p1[1]] | all1[p1[1]] == all2[-p1[1]],na.rm=T)
	else c1 = 0
    p2 = which(!is.na(all2))
	if (length(p2)>0) c2 = sum(all2[p2[1]] == all1[-p2[1]] | all2[p2[1]] == all2[-p2[1]],na.rm=T)
	else c2 = 0
	if (c1>0 | c2>0)
	{
		if (c1>=c2) return (all1[p1[1]])
		else return (all2[p2[1]])
	}
	else return (NA)
}

# name of sample file
name = ".pass.BG.sim"
sample.file = paste0(path,chr,name,"1.sample")
sample.dat = read.table(sample.file,header=T)
n = nrow(sample.dat)-1

# Identifying the family processed
uf.vec = unique(as.character(sample.dat$ID_1[-1]))
fam = uf.vec[f]

# name of IBD cluster file
#IBDclust.name = ".pass_and_fail.gtonly.minDP10.dbGaPID.parents.Mendelcompat.geno05.OmniExpress.dash_cc.hcl"
IBDclust.name = ".dash_cc.IncluPairs.clst"
IBDclust.file = paste0(pathIBD,chr,IBDclust.name)
IBDclust.dat = read.table(IBDclust.file,header=F,col.names = paste0("V",seq_len(maxw)),fill=T)
kc = apply(IBDclust.dat,2,function(vec) !any(is.na(vec)))
IBDclust.dat = as.matrix(IBDclust.dat[,kc])
# Remove second set of beginning and end of interval for .clst files
IBDclust.dat = IBDclust.dat[,-c(4,5)]
ncol(IBDclust.dat)

# Keeping only cluster involving current family
clustwithf = apply(IBDclust.dat,1,function(vec,fam) any(vec[!is.na(vec)]==fam),fam=fam)
IBDclust.dat = IBDclust.dat[clustwithf,]
clustbegin = as.numeric(IBDclust.dat[,2])
clustend = as.numeric(IBDclust.dat[,3])

# Vector of columns in the data corresponding to children and vector for the family analyzed
children.bool = sample.dat$father[-1] != "0" & sample.dat$ID_1[-1] == fam
col.to.keep = c(rep(F,5),rep(children.bool,rep(2,length(children.bool))))
ID.vec = as.character(sample.dat$ID_2[-1][children.bool])
IDrep.vec = rep(ID.vec,rep(2,length(ID.vec)))
#nucfams = rep(paste0(sample.dat$father[-1],sample.dat$mother[-1]),rep(2,nrow(sample.dat)-1))
nucfams = rep(paste0(sample.dat$father[-1],sample.dat$mother[-1])[children.bool],rep(2,sum(children.bool)))

# Initializing the list of list of values of the index where haplotype configuration of the children changes
change.list = list()
for (r in 1:nr) change.list[[r]] = numeric(0)

# Loop over the replicates
for (r in 1:nr)
{
cat(r,"\n")
	# Reading the haplotype file
haplo.file = paste(path,chr,name,r,".haps",sep="")
haplo.dat = matrix(scan(haplo.file,what=character()),nsites,5+2*n,byrow=TRUE)

# Initialization of the haplotype configuration of the children over the first window
haplotypes = haplo.dat[1:ww,col.to.keep]
vpos = as.numeric(haplo.dat[c(1,ww),3])
# Recoding haplotypes
oldhaplo_string=apply(haplotypes,2,function(vec) paste(vec,collapse=""))
oldgroup = recodeHap(oldhaplo_string)

# Distinguishing alleles of the same subject
oldgrouphet = het_subject(oldgroup)
		# Tabulating count of each allele per nuclear family
		tab = table(oldgrouphet,nucfams)
		nfam.allele = apply(tab,1,function(vec) sum(vec>0))
		# If there is more than one allele seen in multiple nuclear families
		if (sum(nfam.allele>1) > 1)
		{
			# Identify allele to keep
			ak = which.max(apply(tab,1,sum))[1]
			# Other alleles present in more than one nuclear family
			for (l in which(names(nfam.allele)!=ak & nfam.allele>1))
			{
				# Identifying nuclear families where allele l is present
				nf = colnames(tab)[tab[names(nfam.allele)==l]>0]
				# Assigning a new allele label to nuclear families other than the first one
				oldgrouphet[oldgrouphet==l][nucfams[oldgrouphet==l]%in%nf[-1]] = max(oldgrouphet) + rep(1:(length(nf)-1),rep(tab[l,][tab[l,]>0][-1]))
			}
		}
	# Relabel the alleles
	oldgrouphet = recodeHap(oldgrouphet)	

# Removing alleles associated to clusters beginning before the end of the current window
cliv = which (clustbegin<vpos[2])
if (length(cliv)>0)
{
	for (cli in cliv)
	{
		tmp = IBDclust.dat[cli,IBDclust.dat[cli,]!=""]
		faml = tmp[seq(4,length(tmp)-1,by=2)]
		tmp2 = strsplit(tmp[seq(5,length(tmp),by=2)],".",fixed=T)
		IDl = sapply(tmp2,function(vec) vec[1])
		IDbyfam = split(IDl,faml)
			ll = length(oldgrouphet)
			all1 = oldgrouphet[seq(1,ll-1,by=2)]
			all2 = oldgrouphet[seq(2,ll,by=2)]
			# Set to missing the allele shared by the most subjects in the cluster if not already missing
			ca = common.allele(all1[ID.vec%in%IDbyfam[[fam]]],all2[ID.vec%in%IDbyfam[[fam]]])
			if (!is.na(ca)) oldgrouphet[IDrep.vec%in%IDbyfam[[fam]]][oldgrouphet[IDrep.vec%in%IDbyfam[[fam]]]==ca] = NA
	}
}

# Loop over the variant positions
#for (i in 2:(nrow(haplo.dat)-ww))
#for (i in 2:100)
for (i in seq(decal+1,(nrow(haplo.dat)-ww),by=decal))
{
#	cat(i,"\n")
	haplotypes = haplo.dat[i:(i-1+ww),col.to.keep]
	vpos = as.numeric(haplo.dat[c(i,i-1+ww),3])
	# Recoding haplotypes
	newhaplo_string=apply(haplotypes,2,function(vec) paste(vec,collapse=""))
	if (!all(newhaplo_string==oldhaplo_string))
	{
	newgroup = recodeHap(newhaplo_string)
	newgrouphet = het_subject(newgroup)

		# Tabulating count of each allele per nuclear family
		tab = table(newgrouphet,nucfams)
		nfam.allele = apply(tab,1,function(vec) sum(vec>0))
		# If there is more than one allele seen in multiple nuclear families
		if (sum(nfam.allele>1) > 1)
		{
			# Identify allele to keep
			ak = which.max(apply(tab,1,sum))[1]
			# Other alleles present in more than one nuclear family
			for (l in which(names(nfam.allele)!=ak & nfam.allele>1))
			{
				# Identifying nuclear families where allele l is present
				nf = colnames(tab)[tab[names(nfam.allele)==l]>0]
				# Assigning a new allele label to nuclear families other than the first one
				newgrouphet[newgrouphet==l][nucfams[newgrouphet==l]%in%nf[-1]] = max(newgrouphet) + rep(1:(length(nf)-1),rep(tab[l,][tab[l,]>0][-1]))
			}
		}
		# Relabel the alleles
		newgrouphet = recodeHap(newgrouphet)	
    
	# Removing alleles associated to clusters beginning before the end of the current window 
	# and ending after the beginning of the current window
	cliv = which (clustbegin<vpos[2] & clustend>vpos[1])
	if (length(cliv)>0)
	{
		for (cli in cliv)
		{
			tmp = IBDclust.dat[cli,IBDclust.dat[cli,]!=""]
			faml = tmp[seq(4,length(tmp)-1,by=2)]
			tmp2 = strsplit(tmp[seq(5,length(tmp),by=2)],".",fixed=T)
			IDl = sapply(tmp2,function(vec) vec[1])
			IDbyfam = split(IDl,faml)
				ll = length(newgrouphet)
				all1 = newgrouphet[seq(1,ll-1,by=2)]
				all2 = newgrouphet[seq(2,ll,by=2)]
				# Set to missing the allele shared by the most subjects in the cluster if not already missing
				ca = common.allele(all1[ID.vec%in%IDbyfam[[fam]]],all2[ID.vec%in%IDbyfam[[fam]]])
				if (!is.na(ca)) newgrouphet[IDrep.vec%in%IDbyfam[[fam]]][newgrouphet[IDrep.vec%in%IDbyfam[[fam]]]==ca] = NA
		}
	}

		# Recording group changes
		if (sum(abs(newgrouphet-oldgrouphet),na.rm=T)!=0) change.list[[r]][length(change.list[[r]])+1] = i
	oldgrouphet = newgrouphet
	oldhaplo_string = newhaplo_string

}
}
}
outfile = paste0("change_",chr,"_",fam,"_d",decal,".RData")
save(change.list,file=outfile)


