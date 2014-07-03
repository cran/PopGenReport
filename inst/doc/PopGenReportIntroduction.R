## ----c1,echo=FALSE,results='hide', message=FALSE, warning=FALSE----------
ip <- installed.packages()
ver <- ip[which(ip[,1]=="PopGenReport"),"Version"]
require(xtable)
require(adegenet)

## ----c1b,echo=FALSE, results='asis'--------------------------------------
cat(paste("\\subtitle {using PopGenReport Ver.",ver,"}\n"))

## ----c2,echo=TRUE, results='markup'--------------------------------------
Sys.setenv(PATH = paste(Sys.getenv("PATH"),
"E:\\PopGenPack\\Miktex\\miktex\\bin\\", sep=.Platform$path.sep))

## ----c3,echo=TRUE, results='markup', message=FALSE, warning=FALSE--------
# Uncomment by deleting the hash from the install.packages command
# Needs to be done only once then it is permanently on your computer.
#install.packages("PopGenReport", repos="http://cran.rstudio.com/")
require(PopGenReport)

## ----c4,echo=TRUE, results='markup', message=FALSE, warning=FALSE--------
#a short test
data(bilby)
summary(bilby)

## ----c5,echo=TRUE, results='hide'----------------------------------------
paste(.libPaths()[1],"/PopGenReport/extdata/platypus1c.csv", sep="" )

## ----c6,echo=TRUE, results='markup'--------------------------------------
platy <- read.csv(paste(.libPaths()[1],"/PopGenReport/extdata/platypus1c.csv", sep="" ))
head(platy)

## ----c6b, echo=FALSE, results='asis'-------------------------------------
platy.csv <- read.csv( paste(.libPaths()[1], "/PopGenReport/extdata/platypus1c.csv",sep=""))[1:4,]
xt <- xtable(platy.csv, caption="Example of a correctly formatted data set in a spreadsheet calculator (e.g. Excel). Please make sure to use the exactly same headings, such as 'ind', 'pop', ['lat', 'long' optional] for the first columns and use unique identifiers for your allel or loci headings. Save your data set as csv file. Please be aware to use a ', as separator. In some countries Excel uses ';' as default (e.g. Germany).")
print(xt, include.rownames=FALSE, size="\\small" )


## ----c7,echo=TRUE, results='markup'--------------------------------------
platy.gen <- read.genetable( paste(.libPaths()[1], "/PopGenReport/extdata/platypus1c.csv"
, sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE, 
sep="/", ploidy=2)

platy.gen  #to check the data. 

## ----c8,echo=TRUE, results='markup'--------------------------------------
platy.gen@ind.names  #check the number and names of individuals
pop(platy.gen) #similar to platy.gen@pop
platy.gen@loc.names #names of all loci
platy.gen@loc.fac  #the number of allels per loci. e.g. the first locus has 4 alleles
table(platy.gen@loc.fac) # a better way to check the number of alleles per loci
platy.gen@all.names #allele names for each loci

## ----c9,echo=TRUE, results='markup'--------------------------------------
platy.gen@tab

## ----c10,echo=TRUE, results='markup'-------------------------------------
str(platy.gen@other)  #lists the content of the slot

## ----xyplot, echo=TRUE, results='as.is', fig.width=4,fig.height=4--------
plot(platy.gen@other$latlong, pch=as.numeric(platy.gen@other$data$group)+15, col=platy.gen@other$data$age)

## ----c11,echo=TRUE,  results='markup', tidy=FALSE------------------------
### for technical reasons we set here mk.pdf=FALSE 
### please set mk.pdf=TRUE if you want to have a report !!!!!!
platy.out1 <- popgenreport(platy.gen, mk.counts=TRUE, mk.pdf=FALSE, 
                           foldername="platy")

## ----c12,echo=TRUE, results='markup', tidy=FALSE-------------------------
platy.out1 <- popgenreport(platy.gen, mk.counts=TRUE, mk.pdf=FALSE, 
                           foldername="platy")

## ----c13,echo=TRUE, results='markup'-------------------------------------
summary(platy.out1)
platy.out1

## ----allperpop, echo=TRUE, fig.width=6, fig.height=4---------------------
barplot(platy.out1$counts$nallelesbypop)

## ----c14,echo=TRUE, results='markup'-------------------------------------
dir(paste(tempdir(),"platy",sep=.Platform$file.sep))

## ----c15,echo=TRUE,  results='markup', tidy=FALSE------------------------
data(bilby)
### for technical reasons we set here mk.pdf=FALSE 
### please set mk.pdf=TRUE if you want to have a report !!!!!!
bilby.complete <- popgenreport(bilby, mk.complete=TRUE, mk.Rcode=TRUE,
                               mk.pdf=FALSE)

## ----c16,echo=TRUE,  results='markup'------------------------------------
### for technical reasons we set here mk.pdf=FALSE 
### please set mk.pdf=TRUE if you want to have a report !!!!!!
bilby.fem <- popgenreport(bilby[bilby@other$sex=="Female"], 
path.pgr=tempdir(), fname="Female", mk.fst=TRUE, mk.counts=FALSE,mk.pdf=FALSE)

### for technical reasons we set here mk.pdf=FALSE 
### please set mk.pdf=TRUE if you want to have a report !!!!!!
bilby.mal <- popgenreport(bilby[bilby@other$sex=="Male"], 
path.pgr=tempdir(), fname="Male", mk.fst=TRUE, mk.counts=FALSE,mk.pdf=FALSE)

## ----c17,echo=TRUE, results='markup'-------------------------------------
bilby.mal
bilby.fem
bilby.mal$fst$FSTpairwise> bilby.fem$fst$FSTpairwise

## ----c18,echo=TRUE, results='markup',message=FALSE, tidy=FALSE-----------

#### convert from utm 55S to latlong WGS84

tiger.gen <- read.genetable( paste(.libPaths()[1],"/PopGenReport/extdata/tiger.csv",
sep="" ), ind=1, pop=2, other.min=3, other.max=6, oneColPerAll=TRUE)


require(rgdal) #load package
xy <- as.matrix(tiger.gen@other$data[,2:1])  #y first then x
#projection from utm 55S to latlong WGS84
latslongs <-project(xy, 
"+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs",inv=T) 

#add it to tiger.gen at the right place (again lats)
tiger.gen@other$latlong <- latslongs[,2:1] #again lat and then long

## ----c19,echo=TRUE,  results='markup'------------------------------------
#customised map
### for technical reasons we set here mk.pdf=FALSE 
### please set mk.pdf=TRUE if you want to have a report !!!!!!
popgenreport(tiger.gen, mk.map=TRUE,mk.counts=FALSE, path.pgr=tempdir(),
mapdotcolor="orange", mapdotsize=as.numeric(tiger.gen@other$data$sex), 
maptype="roadmap", mapdotalpha=0.9,mk.pdf=FALSE)

## ----c20,echo=TRUE, results='hide'---------------------------------------
#custom.snw file
snw <- readLines(paste(.libPaths()[1],"/PopGenReport/swchunks/custom.snw", sep="" ))

## ----c21,echo=FALSE, results='asis'--------------------------------------
snwtab <- xtable(data.frame(snw))
print(snwtab )

