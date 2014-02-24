#leatcost function
genleastcost <- function(cats, fr.raster, gen.dist)
{
#create friction matrix
 fric.mat <- transition(fr.raster,function(x) 1/x[2],4)
# fric.mat <- transition(fr.raster,function(x) 1/(abs(x[1]-x[2])),8)
#set distances to meters
fric.mat@crs@projargs<- "+proj=merc +units=m"
fric.mat.cor <- geoCorrection(fric.mat)

dist.type<-NA
if (gen.dist=="D" || gen.dist=="Gst.Hedrick" || gen.dist=="Gst.Nei") dist.type<- "pop" else dist.type<-"ind" 


if (dist.type=="pop")
{
#calculate the centers if population meassurment is wanted
c.x <- tapply(cats@other$xy[,1],cats@pop, mean)
c.y <- tapply(cats@other$xy[,2],cats@pop, mean)
cp<-cbind(c.x, c.y)
} else 
{
cp <- cbind(cats@other$xy[,1], cats@other$xy[,2])
}

plot(fr.raster)
 #image(fr.raster, col=fr.raster@legend@colortable, asp=1)

 mapcolor <-   col2rgb(cats@other$mapdotcolor)/255
 colX = rgb(mapcolor[1,],mapcolor[2,], mapcolor[3,],cats@other$mapdotalpha)
 points(cats@other$xy,cex=cats@other$mapdotsize, pch= cats@other$mapdottype, col=colX)
if (dist.type=="pop")  points(cp,cex=cats@other$mapdotsize*1.5 , pch= 16, col=rgb(1,0.7,0,0.7))

eucl.mat <- round(as.matrix(dist(cp)),3)
cd.mat <-costDistance(fric.mat.cor, cp, cp)

if (dist.type=="pop") 
{
  dimnames(cd.mat) <- list(cats@pop.names, cats@pop.names) 
  dimnames(eucl.mat) <- list(cats@pop.names, cats@pop.names) 
  npop <- length(levels(cats@pop))
} else
  
{
  dimnames(cd.mat) <- list(cats@ind.names, cats@ind.names) 
  dimnames(eucl.mat) <- list(cats@ind.names, cats@ind.names) 
  npop <- length(cats@ind.names)
}

comb <- t(combn(1:npop,2))

#pathlength matrix
pathlength.mat <- cd.mat
pathlength.mat[,] <- 0


paths <- list()
cols <- rainbow(dim(comb)[1], alpha=0.5)
for (i in 1:dim(comb)[1])
{
sPath <- shortestPath(fric.mat.cor, cp[comb[i,1],], cp[comb[i,2],], output="SpatialLines")
lines(sPath, lwd=1.5, col=cols[i])
paths[[i]] <- sPath
ll <-  round(SpatialLinesLengths(sPath),3)
pathlength.mat[comb[i,1],comb[i,2]] <- ll
pathlength.mat[comb[i,2],comb[i,1]] <- ll
}



#put other calculations here....
# Calculate genetic distances across subpopulations

if (gen.dist=="Gst.Nei")
{
gendist.mat<-as.matrix(pairwise_Gst_Nei(cats))
}
if (gen.dist=="Gst.Hedrick")
{
gendist.mat<-as.matrix(pairwise_Gst_Hedrick(cats))
}
if (gen.dist=="D")
{
gendist.mat<-round(as.matrix(pairwise_D(cats)),4)
}

if (gen.dist=="Smouse")
{
gendist.mat <- as.matrix(gd.smouse(cats,verbose=FALSE))
}
if (gen.dist=="Kosman")
{
gendist.mat <-as.matrix(as.dist(gd.kosman(cats)$geneticdist))
}

if (dist.type=="pop") 
{
  colnames(gendist.mat)<-cats@pop.names
  rownames(gendist.mat)<-cats@pop.names
} else 
{
  colnames(gendist.mat)<-cats@ind.names
  rownames(gendist.mat)<-cats@ind.names
}


return(list(eucl.mat=eucl.mat, cost.mat=cd.mat,pathlength.mat= pathlength.mat, gen.mat=gendist.mat,  paths=paths))



}