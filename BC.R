
require(vegan)
require(cluster)
bc=vegdist(community_nona,"bray")

library(sp)
dist<-spDists(as.matrix(env_var[,10:11]),longlat=TRUE) # distance g?ographique
pluie<-vegdist(env_var$PREC,"euclidean")

pluie<-vegdist(env_var$PREC,"euclidean")

twi<- vegdist(env_var$TWI,"euclidean")

plot(bc~pluie,pch=20)
plot(as.matrix(bc)~dist,ylab="Indice de dissimilarit? de Bray-Curtis",xlab="Distance g?ographique (km)",pch=20)

legend("bottomright","Mantel r = 0.33***",cex=1.2,bty="n")



par(mfrow=c(3,1), mar=c(4,4,2,2))

plot(bc~pluie,pch=20,xlab="delta rainfall",ylab="Bray-Curtis dissimilarity index", main="All data")
legend("bottomright","Mantel r = 0.65***",cex=1.2,bty="n")



mantel(xdis = as.matrix(bc), ydis = dist, method="pearson",permutations = 999)
mantel(xdis = bc, ydis = twi, method="pearson",permutations = 999)
mantel(xdis = bc, ydis = pluie, method="pearson",permutations = 999)



bc_west=vegdist(community_west,"bray")
bc_east=vegdist(community_east,"bray")

twi_west <-vegdist(env_var$TWI[env_var$SIDE=="West"],"euclidean")
twi_east <-vegdist(env_var$TWI[env_var$SIDE=="East"],"euclidean")



plot(bc_west~twi_west,pch=20,xlab="delta TWI",ylab="Bray-Curtis dissimilarity index", main="West")

mantel(xdis = bc_west, ydis = twi_west, method="pearson",permutations = 999)

legend("bottomright","Mantel r = NS",cex=1.2,bty="n")


plot(bc_east~twi_east,pch=20,xlab="delta TWI",ylab="Bray-Curtis dissimilarity index", main="East")

mantel(xdis = bc_east, ydis = twi_east, method="pearson",permutations = 999)

legend("bottomright","Mantel r = 0.39**",cex=1.2,bty="n")



par(mfrow=c(3,1), mar=c(4,4,2,2))


dist_trait <- vegdist(scale(FD_test$CWM),"euclidean")

dist_trait <- vegdist(stat_obs_all$tab_stats_obs[,c("CWM.AxcQ1","CWM.AxcQ2")],"euclidean")

plot(dist_trait~bc,pch=20,ylab="Functional distance",xlab="Bray-Curtis dissimilarity index", main="All data")

mantel(xdis = bc, ydis = dist_trait, method="pearson",permutations = 999)



  legend("topleft","Mantel r = 0.54**",cex=1.2,bty="n")

  dist_trait_west <- vegdist(stat_obs_all$tab_stats_obs[W,c("CWM.AxcQ1","CWM.AxcQ2")],"euclidean")
  

plot(dist_trait_west~bc_west,pch=20,ylab="Functional distance",xlab="Bray-Curtis dissimilarity index", main="West")

mantel(xdis = bc_west, ydis = dist_trait_west, method="pearson",permutations = 999)


  legend("topleft","Mantel r = 0.55**",cex=1.2,bty="n")

  dist_trait_east <- vegdist(stat_obs_all$tab_stats_obs[E,c("CWM.AxcQ1","CWM.AxcQ2")],"euclidean")
  
  

plot(dist_trait_east~bc_east,pch=20,ylab="Functional distance",xlab="Bray-Curtis dissimilarity index", main="East")

mantel(xdis = bc_east, ydis = dist_trait_east, method="pearson",permutations = 999)


  legend("topleft","Mantel r = 0.42**",cex=1.2,bty="n")

