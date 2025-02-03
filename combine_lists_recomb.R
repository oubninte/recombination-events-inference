# Function to shift positions of candidate recombination sites by the window width
# and remove duplicated candidates (candidates detected at start and end of the window)
intervalf = function(pos,ww)
{
	posd = pos + ww
	# positions in 1st window + positions shifted by window length, if not in 
	#c(pos[pos<ww],posd[!(pos%in%posd | posd%in%pos)])
	c(pos[pos<ww & !(posd%in%pos)],posd[!(pos%in%posd)])
}

# Test of the function
intervalf(change22.list40$GSF5519,40)

len.init = sapply(change22.list40,function(l) sapply(l,length))

# Applying intervalf to list of list of values of the index where haplotype configuration changes
change22f.list40 = lapply(change22.list40,function(l) sapply(l,intervalf,ww=40))

len.filtered = sapply(change22f.list40,function(l) sapply(l,length))

# Tallying the number of replicates where any haplotype configuration change was observed
tab22f.list40 = sapply(change22f.list40,function(l) table(unlist(l)))
table(tab22f.list40[[1]])
sapply(tab22f.list40,length)

# Creating the candidate positions file by looping over the families and concatenating their positions
candidatepos22 = cbind(names(tab22f.list40)[1],"22",haplo.dat[as.numeric(names(tab22f.list40[[1]]))-1,3],haplo.dat[as.numeric(names(tab22f.list40[[1]])),3],tab22f.list40[[1]]/10)
for (i in 2:length(tab22f.list40))
{
	candidatepos22 = rbind(candidatepos22,cbind(names(tab22f.list40)[i],"22",haplo.dat[as.numeric(names(tab22f.list40[[i]]))-1,3],haplo.dat[as.numeric(names(tab22f.list40[[i]])),3],tab22f.list40[[i]]/10))
}
dimnames(candidatepos22)[[2]] = c("FAMILY","CHR","START","END","PROB")

# Writing candidate positions file
write.table(candidatepos22,"candidatepos22.txt",quote=F,row.names=F)
