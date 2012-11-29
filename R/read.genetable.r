
#must provide a filename and if pop, ind and in 2 or two columns per allele
read.genetable <- function(filename, pop=NULL, ind=NULL,lat=NULL, long=NULL, other.min=NULL, other.max=NULL,oneColPerAll, missing=NA,  sep=NULL, ncode=NULL  )
{
gfile <-read.csv(filename)

popnr=NULL
pops=NULL
indnr=NULL
inds=NULL

latlong=NULL

if ("pop" %in% colnames(gfile))            
  {
  pops <- as.factor(gfile$pop)
  popnr <- which(colnames(gfile)=="pop")
  }
if ("ind" %in% colnames(gfile)) 
  {
  inds <- as.factor(gfile$ind)
  indnr <- which(colnames(gfile)=="ind")
  }





genes <- gfile
rem=NULL

if (!is.null(other.min) & !is.null(other.max)) rem <- other.min:other.max         
if (!is.null(popnr)) rem <- c(rem, popnr)
if (!is.null(indnr)) rem <- c(rem, indnr)
if (!is.null(lat))  rem <- c(rem, lat)
if (!is.null(long)) rem <- c(rem, long)

if (!is.null(rem)) genes <- gfile[, -c(rem)]   else genes <- gfile


if (oneColPerAll==F)
{
res <-data.frame(allele= rep(NA,dim(genes)[1]) )

for (i in seq(1,dim(genes)[2]-1,2))
  {
  res[,ceiling(i/2)] <- paste(genes[,i], genes[,i+1],sep="/")
  
  colnames(res)[ceiling(i/2)] <- paste( colnames(genes)[i], colnames(genes)[i+1],sep="/" )
  }
  
sep="/"

genes <- res
} 

  
  
  df <- df2genind(genes, pop=pops, ind.names=inds, missing=missing,     ncode=ncode,loc.names=colnames(genes), sep=sep)
  
  if (!is.null(lat) & !is.null(long))
  df@other$latlong <- data.frame(lat=gfile[,lat], long=gfile[,long])
  if (!is.null(other.min) & !is.null(other.max)) 
  {
  df@other$data  <- data.frame(gfile[,c(other.min:other.max)])
  }

return(df)



}
