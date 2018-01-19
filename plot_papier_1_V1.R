########################################################################################################################

########################################################################################################################




#####  Version 1 : on ne prend en compte que l'abondance, pas de correction de Pvalues sur le 4eme coin #####


##### init data #####




##### initialise data #####


library(vegan)
library(ade4)




# modif function from ade4 : color arrows and labels
scatterutil.eti.circ2 <- 
  function (x, y, label, clabel, origin = c(0, 0), boxes = TRUE, colboxes="black", collab = "black", labfont = 1) 
  {# modif function from ade4 : color arrows and labels
    
    if (is.null(label)) 
      return(invisible())
    if (any(is.na(label))) 
      return(invisible())
    if (any(label == "")) 
      return(invisible())
    xref <- x - origin[1]
    yref <- y - origin[2]
    for (i in 1:(length(x))) {
      cha <- as.character(label[i])
      cha <- paste(" ", cha, " ", sep = "")
      cex0 <- par("cex") * clabel
      xh <- strwidth(cha, cex = cex0)
      yh <- strheight(cha, cex = cex0) * 5/6
      if ((xref[i] > yref[i]) & (xref[i] > -yref[i])) {
        x1 <- x[i] + xh/2
        y1 <- y[i]
      }
      else if ((xref[i] > yref[i]) & (xref[i] <= (-yref[i]))) {
        x1 <- x[i]
        y1 <- y[i] - yh
      }
      else if ((xref[i] <= yref[i]) & (xref[i] <= (-yref[i]))) {
        x1 <- x[i] - xh/2
        y1 <- y[i]
      }
      else if ((xref[i] <= yref[i]) & (xref[i] > (-yref[i]))) {
        x1 <- x[i]
        y1 <- y[i] + yh
      }
      if (boxes) {
        rect(x1 - xh/2, y1 - yh, x1 + xh/2, y1 + yh, border = colboxes )
      }
      text(x1, y1, cha, cex = cex0, col = collab, font = labfont)
    }
  }

s.arrow2 <- 
  function (dfxy, xax = 1, yax = 2, label = row.names(dfxy), clabel = 1, 
            pch = 20, cpoint = 0, boxes = TRUE, edge = TRUE, origin = c(0, 
                                                                        0), xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
            cgrid = 1, sub = "", csub = 1.25, possub = "bottomleft", 
            pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE, colboxes="black", collab = "black" , arrowcol = "black", 
            labfont = 1,  box = T) 
  {# modif function from ade4 : color arrows and labels
    
    arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1, 
                       edge, arrowcol = 'black') {
      d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
      if (d0 < 1e-07) 
        return(invisible())
      segments(x0, y0, x1, y1, lty = lty, col = arrowcol)
      h <- strheight("A", cex = par("cex"))
      if (d0 > 2 * h) {
        x0 <- x1 - h * (x1 - x0)/d0
        y0 <- y1 - h * (y1 - y0)/d0
        if (edge) 
          arrows(x0, y0, x1, y1, angle = ang, length = len, 
                 lty = 1, col = arrowcol)
      }
    }
    dfxy <- data.frame(dfxy)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
                            xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
                            cgrid = cgrid, include.origin = TRUE, origin = origin, 
                            sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
                            contour = contour, area = area, add.plot = add.plot)
    if (grid & !add.plot) 
      scatterutil.grid(cgrid)
    if (addaxes & !add.plot) 
      abline(h = 0, v = 0, lty = 1)
    if (cpoint > 0) 
      points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
    for (i in 1:(length(coo$x))) arrow1(origin[1], origin[2], 
                                        coo$x[i], coo$y[i], edge = edge, arrowcol = arrowcol)
    if (clabel > 0) 
      scatterutil.eti.circ2(coo$x, coo$y, label, clabel, origin, 
                            boxes,colboxes=colboxes, collab = collab, labfont = labfont)
    if (csub > 0) 
      scatterutil.sub(sub, csub, possub)
    if (box == T)
      box()
    invisible(match.call())
  }
arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1, 
                   edge, arrowcol = 'black') {
  d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
  if (d0 < 1e-07) 
    return(invisible())
  segments(x0, y0, x1, y1, lty = lty, col = arrowcol)
  h <- strheight("A", cex = par("cex"))
  if (d0 > 2 * h) {
    x0 <- x1 - h * (x1 - x0)/d0
    y0 <- y1 - h * (y1 - y0)/d0
    if (edge) 
      arrows(x0, y0, x1, y1, angle = ang, length = len, 
             lty = 1, col = arrowcol)
  }
}

#####  load data ##### 

# occurences

FBA <- read.csv("/home/thesardfou/Documents/data/tabs/Fiche terrain FBA_saisie_10_10_17.csv", sep = ",")

# env var


tab_extract3 <- read.csv('/home/thesardfou/Documents/data/tab_extract3-all_fin2017.csv')


# TF

FT_all_kone = read.csv( '/home/thesardfou/Documents/data/FT_all_kone_111017.csv')


##### make clean ####


#enlève les lignes avec les morts
FBA<-subset(FBA,!FBA$strate==0)
FBA<-subset(FBA,!FBA$Taxon=="")


# tout les points  
env_var <- tab_extract3

names(env_var) <- c("X","id","ELEV", "INS", "PREC", "TWI", "CURV", "ASP","SLOP","LONG","LAT","EDGE","TWIxPREC","AREA")


# omit sites 2* : pour un nombre identique de pts par cote

# supprimer les pts 1 et 2 car en zone plate !!! 

###### supprimer tiawe tribe et flats? #####

omit_sites = c(1,2,11,27,28,29)

# supprimer tiawe tribe

#  omit_sites = c(11,27,28,29)


# moins le pts 12 ?


omit_rows <- which(FBA$id_pts %in% omit_sites == T)

FBA <- FBA[-omit_rows,]



env_var <- tab_extract3[-omit_sites,] # enleve les parcelles tiawe tribe & flat




names(env_var) <- c("X","id","ELEV", "INS", "PREC", "TWI", "CURV", "ASP","SLOP","LONG","LAT","EDGE","TWIxPREC","AREA", "PMA")


##### pour récupérer les points  avec les valeurs de TWI les plus extrmes à kone (après omit sites !!!!!!!!! ) * #####

# autant de pts dans chaque sites à kone !!!! 


order_twi_kone <- order(env_var$TWI[env_var$id %in% 1:32]) # plaine de kone sans tiawe tribe



##  dans chaque paysages : 10 low et 10 high  # en gardant 26 points
kone_pts_low_twi <- order_twi_kone[1:10]
kone_pts_high_twi <- order_twi_kone[(length(order_twi_kone)-9):length(order_twi_kone)]
compare_EW <- sort(c(kone_pts_low_twi,kone_pts_high_twi, (length(order_twi_kone)+1):nrow(env_var)))


id_sites_compare_EW <- env_var$id[compare_EW]

FBA <- FBA[FBA$id_pts %in% id_sites_compare_EW,]

env_var <- env_var[compare_EW,][1:length(compare_EW),]

dim(env_var)
length(unique(compare_EW))

length(unique(FBA$id_pts))

## *****



# corrige NA sapect = transform en southness 


env_var$ASP[is.na(env_var$ASP)] = pi/2 # NA = plat = pi/2 from south


env_var$ASP = abs(env_var$ASP-pi)



# class topo


env_var$TOPO <- ifelse(env_var$TWI < 9,"ridge", "talweg")
env_var$SIDE <- ifelse(env_var$PREC < 2000,"West", "East")

# autre class topo : 

# env_var$TOPO2 <- ifelse(env_var$TWI < 9 & env_var$CURV>0, "ridge", ifelse(env_var$TWI > 9 & env_var$CURV<0, "talweg", "NA") )


env_var$TOPO= as.factor(env_var$TOPO )

env_var$SIDE = as.factor(env_var$SIDE)

ridge_west <- which(env_var$TOPO == "ridge" & env_var$SIDE == "West")
ridge_east <- which(env_var$TOPO == "ridge" & env_var$SIDE == "East")
talweg_west <- which(env_var$TOPO == "talweg" & env_var$SIDE == "West")
talweg_east <- which(env_var$TOPO == "talweg" & env_var$SIDE == "East")
prec_topo_class <- c()
prec_topo_class[ridge_west] <- rep("ridge_west",length(ridge_west))
prec_topo_class[talweg_west] <- rep("talweg_west",length(talweg_west))
prec_topo_class[ridge_east] <- rep("ridge_east",length(ridge_east))
prec_topo_class[talweg_east] <- rep("talweg_east",length(talweg_east))
prec_topo_class <- factor(prec_topo_class)
prec_topo_class <- factor(prec_topo_class, levels = c("ridge_west", "talweg_west", "ridge_east","talweg_east"))


env_var$SIDE = as.factor(env_var$SIDE)

env_var$CLASS <- prec_topo_class



RW = env_var$CLASS == levels(env_var$CLASS)[1]
TW = env_var$CLASS == levels(env_var$CLASS)[2]
RE = env_var$CLASS == levels(env_var$CLASS)[3]
TE = env_var$CLASS == levels(env_var$CLASS)[4]

REW = env_var$CLASS == levels(env_var$CLASS)[1] | env_var$CLASS == levels(env_var$CLASS)[3]
TEW = env_var$CLASS == levels(env_var$CLASS)[2] | env_var$CLASS == levels(env_var$CLASS)[4]


W = env_var$SIDE == levels(env_var$SIDE)[2]
E = env_var$SIDE == levels(env_var$SIDE)[1]



# id_points = rownames & supprime la premiere colone en trop

rownames(env_var) = env_var$id

env_var = env_var[,-c(1)]

dim(env_var)



###########################################################################################
############################ matrice d'abondance espèce par site ##########################
###########################################################################################

community <- aggregate(FBA$Taxon, list(Site=FBA$id_pts, Taxon = FBA$Taxon), length)

# unique(community$Taxon)

library(reshape)

community <- cast(community,Site~Taxon) #Les variables site et taxon deviennent lignes et colonnes

community[is.na(community)] <- 0 #remplace les NA en 0

ordre_site <- order(as.numeric(as.character(community$Site)))

community <- community[ordre_site,]

community <- community[,2:ncol(community)]

dim(community)

###########################################################################################################
############################ Matrice espèce*traits fonctionnels ###########################################
###########################################################################################################


Taxon <- colnames(community)

###################################################################################################
######################### TRAITS FONCTIONNELS ~ COMMUNAUTÉ ############################
###################################################################################################

taxons_TF<-FT_all_kone$taxon       
TF_kone<-FT_all_kone[,2:ncol(FT_all_kone)]

community<-as.matrix(community)

colnames(community)<-Taxon

rownames(TF_kone)<-taxons_TF

# colnames(community)


TF_kone_temp <- data.frame(SLA = TF_kone$mean_SLA..m..kg..1., LA = TF_kone$LA_tax..mean_LA..cm..., LDMC = TF_kone$LDMC_tax..mean_LDMC..mg.g.1..,
                           WD =TF_kone$WD.mean  )

rownames(TF_kone_temp) <- taxons_TF

# only species with fonctionnal traits


TF_kone_temp_nona <- TF_kone_temp[complete.cases(TF_kone_temp),]


# only species with fonctionnal traits in community
Taxon[!Taxon %in% row.names(TF_kone_temp_nona)]

community_nona <-  community[,colnames(community) %in% rownames(TF_kone_temp_nona) ]
community_all = community_nona

# traits only for species with fonctionnal traits in community

TF_kone_com_nona <- TF_kone_temp_nona[rownames(TF_kone_temp_nona) %in% colnames(community_nona),]



# FT_all_kone$taxon

FBA_nona <-  FBA[FBA$Taxon %in% rownames(TF_kone_com_nona),]

##### abundances with traits #####


tf_all <- data.frame(Taxon = rownames(TF_kone_com_nona), TF_kone_com_nona)
abund_all <- data.frame(Taxon = FBA$Taxon)

abund_all_tf <-  merge(abund_all,tf_all , by = "Taxon")


abund_all_site <- data.frame(Taxon = FBA$Taxon, site = FBA$id_pts)


abund_all_tf_site <-  merge( abund_all_site, tf_all  , by = "Taxon")

####### transform LA SLA & ajoute LT #######



# TF_kone_com_nona$LA =log(TF_kone_com_nona$LA)
# TF_kone_com_nona$SLA = log(TF_kone_com_nona$SLA)
# TF_kone_com_nona$LT = 1/(TF_kone_com_nona$SLA*TF_kone_com_nona$LDMC)


###### data est west ######

community_west <- community_nona[env_var$SIDE == "West",][,!apply(community_nona[env_var$SIDE == "West",],2,sum)==0]
community_east <- community_nona[env_var$SIDE == "East",][,!apply(community_nona[env_var$SIDE == "East",],2,sum)==0]


TF_west <- TF_kone_com_nona[rownames(TF_kone_com_nona) %in% colnames(community_west),]

TF_east <- TF_kone_com_nona[rownames(TF_kone_com_nona) %in% colnames(community_east),]

env_var_west = env_var [env_var$SIDE == "West",]
env_var_east = env_var [env_var$SIDE == "East",]


#######################################################################################################################
############################################################################################################################

##### NULL MODEL version 1 = coordordinates on RLQ = ps traits ==> compute stats ==> swap community matrix (null model) ==> compare (SES)  #####

#######################################################################################################################
#######################################################################################################################


###############################################################################################################################################################

###### stats obs TF === CWV, CR..... #####

#################################################################################################################################################


# obs datas


library(ade4)

compall = list()

# with abundances

compall$spe <- data.frame(community_all, check.names = F)




compall$traits = TF_kone_com_nona

# log transform LA & SLA

# compall$traits$SLA = log(compall$traits$SLA)

# compall$traits$LA =log(compall$traits$LA)


# compall$traits$LT = 1/(compall$traits$SLA*compall$traits$LDMC)



env_var_select = env_var[,c("TWI","INS", "CURV", "SIDE")]
# env_var_select = env_var[,c( "TWI","INS", "CURV", "SIDE", "CLASS")]

#  env_var_select = env_var[,c( "TWI", "SIDE")]

# env_var_select = env_var[,c("SIDE", "CLASS")]


compall$env = env_var_select


compall_presabs = compall


compall_presabs$spe = data.frame(ifelse(compall$spe == 0, 0, 1), check.names = F) # transform direct la matrix abundance en pers / abs


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.all <- dudi.coa(compall$spe, scannf = FALSE)

acpR.all <- dudi.hillsmith(compall$env, row.w = afcL.all$lw,
                           scannf = FALSE)


acpQ.all <- dudi.pca(compall$traits, row.w = afcL.all$cw,
                     scannf = FALSE)

s.arrow(acpQ.all$co)



afcL.all_presabs <- dudi.coa(compall_presabs$spe, scannf = FALSE)

acpR.all_presabs <- dudi.hillsmith(compall_presabs$env, row.w = afcL.all_presabs$lw,
                           scannf = FALSE)


acpQ.all_presabs <- dudi.pca(compall_presabs$traits, row.w = afcL.all_presabs$cw,
                     scannf = FALSE)
#### RLQ ####

rlq.all <- rlq(acpR.all, afcL.all, acpQ.all,
               scannf = FALSE)


rlq.all_presabs <- rlq(acpR.all_presabs, afcL.all_presabs, acpQ.all_presabs,
               scannf = FALSE)

##### fourthcorner #####

nrepet <- 9999

four.comb.all <- fourthcorner(compall$env, compall$spe,
                              compall$traits, modeltype = 6, p.adjust.method.G = "none", 
                              p.adjust.method.D = "none", nrepet = nrepet) 
four.comb.all
four.comb.all_presabs <- fourthcorner(compall_presabs$env, compall_presabs$spe,
                                      compall_presabs$traits, modeltype = 6, p.adjust.method.G = "none", 
                                         p.adjust.method.D = "none", nrepet = nrepet)


Srlq <- fourthcorner2(compall$env, compall$spe, compall$traits,
                      modeltype = 6, p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ



four.comb.all_presabs

plot(four.comb.all)

plot(four.comb.all_presabs)
testrlq.all <- randtest(rlq.all, modeltype = 6, nrepet = nrepet)
testrlq.all

plot(testrlq.all)


testrlq.all_presabs <- randtest(rlq.all_presabs, modeltype = 6, nrepet = nrepet)
testrlq.all_presabs

plot(testrlq.all_presabs)

testRaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testRaxes.comb.all, stat = "D")

testQaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testQaxes.comb.all, stat = "D")

par(mfrow = c(1, 2))
plot(testQaxes.comb.all, alpha = 0.05, type = "biplot",
     stat = "D2")
plot(testRaxes.comb.all, alpha = 0.05, type = "biplot",
     stat = "D2")
par(mfrow = c(1, 1))


testRaxes.comb.all_presabs <- fourthcorner.rlq(rlq.all_presabs, modeltype = 6,
                                       typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testRaxes.comb.all_presabs, stat = "D")

testQaxes.comb.all_presabs <- fourthcorner.rlq(rlq.all_presabs, modeltype = 6,
                                       typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testQaxes.comb.all_presabs, stat = "D")

par(mfrow = c(1, 2))
plot(testRaxes.comb.all_presabs, alpha = 0.05, type = "table",
     stat = "D2")
plot(testQaxes.comb.all_presabs, alpha = 0.05, type = "table",
     stat = "D2")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(testQaxes.comb.all_presabs, alpha = 0.05, type = "biplot",
     stat = "D2")
plot(testRaxes.comb.all_presabs, alpha = 0.05, type = "biplot",
     stat = "D2")
par(mfrow = c(1, 1))


plot(four.comb.all, x.rlq = rlq.all, alpha = 0.05,
     stat = "D2", type = "biplot")

plot(four.comb.all_presabs, x.rlq = rlq.all_presabs, alpha = 0.05,
     stat = "D2", type = "biplot")



###### plot RLQ all abund ( no class) #####

percent_axis = signif(100*rlq.all$eig/sum(rlq.all$eig), digits = 4)

par(mar = c(4,4,1,1))
plot(0, 0, xlim = c(-.53, .82), ylim = c(-.22,.22), type = "n", 
     xlab = paste("RLQ axis 1 (",percent_axis[1],"%)",sep = ""), 
     ylab = paste("RLQ axis 2 (",percent_axis[2],"%)",sep = ""), xaxt='n', yaxt='n')

axis(1, seq(-.5,.75,.25))
axis(2, seq(-.4,.4,.1))

abline(h=seq(-.4,.4,.1), v=seq(-.5,.75,.25), col="gray", lty=3)
abline(h=seq(0), v=seq(0), col="black", lty=1)

# s.arrow(rlq.all_presabs$co,add.plot = T)
s.arrow2(rlq.all$co[c(1,2,4),], add.plot = T, labfont = 2, box = F)
s.arrow2(rlq.all$co[3,], colboxes = 'grey60',  collab = 'grey60', arrowcol = 'grey60', add.plot = T, labfont = 2, box = F)

r_var = rownames(rlq.all$li)
r_var[4] = "EAST"
r_var[5] = "WEST"

text(rlq.all$li$Axis1, rlq.all$li$Axis2, r_var, col = c("black", 'grey60','black','black', 'black'), pos = c(3,4,1,3,3), font = 2)

arrow1(0,0,rlq.all$li$Axis1[1], rlq.all$li$Axis2[1], edge = T, arrowcol = 'black', lty = 2)
arrow1(0,0,rlq.all$li$Axis1[2], rlq.all$li$Axis2[2], edge = T, arrowcol = 'grey60', lty = 2)
arrow1(0,0,rlq.all$li$Axis1[3], rlq.all$li$Axis2[3], edge = T, arrowcol = 'black', lty = 2)

points(rlq.all$li$Axis1[4:5], rlq.all$li$Axis2[4:5], col = c('black', 'black'), pch = c(0,0))



###### plot RLQ all pres abs ( no class) #####

percent_axis = signif(100*rlq.all_presabs$eig/sum(rlq.all_presabs$eig), digits = 4)

par(mar = c(4,4,1,1))
plot(0, 0, xlim = c(-.53, .77), ylim = c(-.22,.22), type = "n", 
     xlab = paste("RLQ axis 1 (",percent_axis[1],"%)",sep = ""), 
     ylab = paste("RLQ axis 2 (",percent_axis[2],"%)",sep = ""), xaxt='n', yaxt='n')

axis(1, seq(-.5,.75,.25))
axis(2, seq(-.4,.4,.1))

abline(h=seq(-.4,.4,.1), v=seq(-.5,.75,.25), col="gray", lty=3)
abline(h=seq(0), v=seq(0), col="black", lty=1)

# s.arrow(rlq.all_presabs$co,add.plot = T)
s.arrow2(rlq.all_presabs$co[c(1,2,4),], add.plot = T, labfont = 2, box = F)
 s.arrow2(rlq.all_presabs$co[3,], colboxes = 'grey60',  collab = 'grey60', arrowcol = 'grey60', add.plot = T, labfont = 2, box = F)

r_var = rownames(rlq.all_presabs$li)
r_var[4] = "EAST"
r_var[5] = "WEST"

text(rlq.all_presabs$li$Axis1, rlq.all_presabs$li$Axis2, r_var, col = c("black", 'grey60','black','black', 'black'), pos = c(3,2,1,3,3), font = 2)

arrow1(0,0,rlq.all_presabs$li$Axis1[1], rlq.all_presabs$li$Axis2[1], edge = T, arrowcol = 'black', lty = 2)
arrow1(0,0,rlq.all_presabs$li$Axis1[2], rlq.all_presabs$li$Axis2[2], edge = T, arrowcol = 'grey60', lty = 2)
arrow1(0,0,rlq.all_presabs$li$Axis1[3], rlq.all_presabs$li$Axis2[3], edge = T, arrowcol = 'black', lty = 2)

 points(rlq.all_presabs$li$Axis1[4:5], rlq.all_presabs$li$Axis2[4:5], col = c('black', 'black'), pch = c(0,0))


 ###### plot RLQ all pres abs ( with topo class) #####
 
 percent_axis = signif(100*rlq.all_presabs$eig/sum(rlq.all_presabs$eig), digits = 4)
 
 par(mar = c(4,4,1,1))
 plot(0, 0, xlim = c(-.65, 1), ylim = c(-.22,.22), type = "n", 
      xlab = paste("RLQ axis 1 (",percent_axis[1],"%)",sep = ""), 
      ylab = paste("RLQ axis 2 (",percent_axis[2],"%)",sep = ""), xaxt='n', yaxt='n')
 
 axis(1, seq(-.5,1,.25))
 axis(2, seq(-.4,.4,.1))
 
 abline(h=seq(-.4,.4,.1), v=seq(-.5,.75,.25), col="gray", lty=3)
 abline(h=0, v=0, col="black", lty=1)
 
 # s.arrow(rlq.all_presabs$co,add.plot = T)
 s.arrow2(rlq.all_presabs$co, add.plot = T, labfont = 2, box = F)
 # s.arrow2(rlq.all_presabs$co[1,], colboxes = 'grey60',  collab = 'grey60', arrowcol = 'grey60', add.plot = T, labfont = 2, box = F)
 
 r_var = rownames(rlq.all_presabs$li)
 r_var[4] = "EAST"
 r_var[5] = "WEST"
 r_var[6] = "West ridge"
 r_var[7] = "West talweg"
 r_var[8] = "East ridge"
 r_var[9] = "East talweg"
 
 
 text(rlq.all_presabs$li$Axis1, rlq.all_presabs$li$Axis2, r_var,
      col = c("black", 'grey60','black','black', 'black', 'black', 'black', 'black', 'black'),
      pos = c(3,2,1,3,1,2,3,1,3), font = 2)
 
 arrow1(0,0,rlq.all_presabs$li$Axis1[1], rlq.all_presabs$li$Axis2[1], edge = T, arrowcol = 'black', lty = 2)
 arrow1(0,0,rlq.all_presabs$li$Axis1[2], rlq.all_presabs$li$Axis2[2], edge = T, arrowcol = 'grey60', lty = 2)
 arrow1(0,0,rlq.all_presabs$li$Axis1[3], rlq.all_presabs$li$Axis2[3], edge = T, arrowcol = 'black', lty = 2)
 
 points(rlq.all_presabs$li$Axis1[4:5], rlq.all_presabs$li$Axis2[4:5], col = c('black', 'black'), pch = c(0,0))
 points(rlq.all_presabs$li$Axis1[6:9], rlq.all_presabs$li$Axis2[6:9], col = c('black', 'black'), pch = c(0,0))
 
 
 plot(rlq.all_presabs$mR)
 
 
 

# textbox(tf_co[1,], rlq.all_presabs$co$Comp2[1], row.names(rlq.all_presabs$co)[1], box=TRUE, margin = c(-.04),justify = "c")


plot.new()
textbox(c(-0.1,0.4), .5, c("keep going",rep("and going",10)), justify='c', cex=0.6,
        leading=1, font=4, border="gold", lty=2, lwd=4, margin=0.025)


plot.new()
textbox(c(0,0.3), 1, c("many words","more words","why not?",
                       "keep going",rep("and going",10)), box=TRUE, margin = .3)

#############################################################################################################################################
#############################################################################################################################################

###### TO PLOT ######

#############################################################################################################################################
#############################################################################################################################################
library(pgirmess) 

library(lattice)
library(gtools)
##### init data #####


stat_obs_all = readRDS(file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all_15012018.rds")

md_null_all =  readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all03012018.rds")

compare_all = compare_stats_obs_null(stat_obs_all,md_null_all)


stat_obs_west = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_west03012018.rds")

md_null_west = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west03012018.rds")

compare_west = compare_stats_obs_null(stat_obs_west,md_null_west)


stat_obs_east = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_east03012018.rds")

md_null_east = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east03012018.rds")


compare_east = compare_stats_obs_null(stat_obs_east,md_null_east)



#### plot obs values  all  #####




# compare est - west : on prend quoi? 

# vector colnames : 

names(stat_obs_all$tab_stats_obs)

vecvar = c(3,4,5,6,9,10,11,12,13,14,15,16,17,18,27,28,35,36)


vecvar = c(5,6,13,14,35,36)

vecvar = c(1,2,3,4,27,28)


length(vecvar)

vars=names(stat_obs_all$tab_stats_obs)[vecvar]

vars_title = vars

vars_title =  c("CV A1", "CV A2", "Ranges A1", "Range A2", "CM A1", "CM A2")
vars_title =  c("Ranges A1", "Range A2", "CWV A1", "CWV A2",  "CWM A1", "CWM A2")


X11( width= 10 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F1.png", width = 900, height = 500, pointsize = 14)

layout(matrix(1:12, 2, 6, byrow = T),widths = rep(c(3,1),3) , heights= rep(1,2))

#layout.show(2)
 
twi = log(env_var$TWI)

for(i in vecvar){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  par(mar= c(2.5,2.5,2.5,0))
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_all$tab_stats_obs[,i]), main = vars_title[which(vecvar == i)], type = "n")
  
  # title(main = vars_title[which(vecvar == i)], line = 1, cex = 2)
  points(stat_obs_all$tab_stats_obs[E,i]~twi[E])
  points(stat_obs_all$tab_stats_obs[W,i]~twi[W], pch = 16)
  
  kenw = cor.test(stat_obs_all$tab_stats_obs[W,i],twi[W], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("    O : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  kene = cor.test(stat_obs_all$tab_stats_obs[E,i],twi[E], method="kendall", use="pairwise", exact = F) 
  
  t <-signif(kene$estimate, 2)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  

  
  kt = kruskal.test( stat_obs_all$tab_stats_obs[,i] ~ as.numeric(as.factor(env_var$SIDE) ) )
  
  par(mar= c(2.5,0,2.5,1))
  
  plot(stat_obs_all$tab_stats_obs[,i] ~ as.factor(env_var$SIDE), axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "")
  axis(1, at = c(1,2), labels = c("E", "W"))
  stars <- stars.pval(kt$p.value)
  
  mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.6, adj=1)
  
  
}




dev.off()





######## alternative plot : plot CM and Ranges ( from bernard verdier 2012) ######





df_segA1 = c()
df_segA2 = c()
for (i in 1:40){
mi = apply(stat_obs_all$list_scores_sp_com[[i]],2, min )
ma = apply(stat_obs_all$list_scores_sp_com[[i]],2, max )

t = env_var$TWI[i]
df_segA1 = rbind(df_segA1,c(t, mi[1],t, ma[1]))
df_segA2 = rbind(df_segA2,c(t, mi[2],t, ma[2]))
}

df_seg = list(df_segA1,df_segA2)


X11( width= 10 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F1.png", width = 900, height = 500, pointsize = 14)

####### alter 1 : kenndal cor  ######

layout(matrix(1:4, 2, 2, byrow = T),widths = c(4,1) , heights= rep(1,2))

vars = stat_obs_all$tab_stats_obs[,c(35,36)]
varnames = c("Community Mean A1","Community Mean A2")
for (i in 1:2){

par(mar= c(3.5,2.5,2.5,0))
twi = env_var$TWI

plot(0, 0, xlim = range(twi), ylim = range(df_seg[[i]][,c(2,4)]), main = varnames[i], xlab = "TWI", ylab = "", type = "n", log = "x")

points(vars[W,i] ~ twi[W], pch = 16)

Segments(df_seg[[i]][W], lwd =2, col = "black")

points(vars[E,i] ~ twi[E], pch = 16, col = "grey50")

Segments(df_seg[[i]][E], lwd = 2, col = "grey50")

kenw = cor.test(vars[W,i],twi[W], method="kendall", use="pairwise", exact = F) 

t <- signif(kenw$estimate, 2)
stars <- stars.pval(kenw$p.value)
mtext(text = substitute(paste("    O : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)

kene = cor.test(vars[E,i],twi[E], method="kendall", use="pairwise", exact = F) 

t <-signif(kene$estimate, 2)
stars <- stars.pval(kene$p.value)
mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)


kt = kruskal.test(vars[,i] ~ as.numeric(as.factor(env_var$SIDE) ) )

par(mar= c(2.5,0,2.5,1))

plot(vars[,i] ~ as.factor(env_var$SIDE), axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "", ylim = range(df_seg[[i]][,c(2,4)]))
axis(1, at = c(1,2), labels = c("E", "W"))
stars <- stars.pval(kt$p.value)

mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.6, adj=1)



}


####### 
####### alter 2 : linearreg  ######


layout(matrix(1:4, 2, 2, byrow = T),widths = c(4,1) , heights= rep(1,2))

vars = stat_obs_all$tab_stats_obs[,c(35,36)]
varnames = c("Community Mean (RLQ Axis 1)","Community Mean (RLQ Axis 2)")
##### axe 1 #####

par(mar= c(4,4,1,0))
twi = env_var$TWI
  
  plot(0, 0, xlim = range(twi), ylim = range(df_seg[[1]][,c(2,4)]), xlab = "TWI", ylab = varnames[1], type = "n", log = "x")
  
  points(vars[W,1] ~ twi[W], pch = 16)
  
  Segments(df_seg[[1]][W], lwd =2, col = "black")
  
  points(vars[E,1] ~ twi[E], pch = 16, col = "grey50")
  
  Segments(df_seg[[1]][E], lwd = 2, col = "grey50")

  mw = lm(vars[W,1] ~ twi[W])
  kenw = summary(mw) 
  t <- signif(kenw$r.squared, 2)
  stars <- stars.pval(kenw$coefficients[8])
  mtext(text = substitute(paste("    W : R² = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=1)

  abline(mw, col = "black", untf=T , lty = 2)
  
  me = lm(vars[E,1]~twi[E])
  kene = summary(me)  
  
  t <-signif(kene$r.squared, 2)
  stars <- stars.pval(kene$coefficients[8])
  mtext(text = substitute(paste("     E : R² = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  abline(me, col = "grey50", untf=T , lty = 2)
  
  kt = kruskal.test(vars[,1] ~ as.numeric(as.factor(env_var$SIDE) ) )
  
  par(mar= c(4,0,1,1))
  
  plot(vars[,1] ~ as.factor(env_var$SIDE), axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "", ylim = range(df_seg[[1]][,c(2,4)]))
  axis(1, at = c(1,2), labels = c("E", "W"))
  stars <- stars.pval(kt$p.value)
  
  mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.7, adj=1)
  
  ##### axe 2 #####
  par(mar= c(4,4,1,0))
  twi = env_var$TWI
  
  plot(0, 0, xlim = range(twi), ylim = range(df_seg[[2]][,c(2,4)]), xlab = "TWI", ylab = varnames[2], type = "n", log = "x")
  
  points(vars[W,2] ~ twi[W], pch = 16)
  
  Segments(df_seg[[2]][W], lwd =2, col = "black")
  
  points(vars[E,2] ~ twi[E], pch = 16, col = "grey50")
  
  Segments(df_seg[[2]][E], lwd = 2, col = "grey50")
  
  mw = lm(vars[W,2] ~ twi[W])
  kenw = summary(mw) 
  t <- signif(kenw$r.squared, 2)
  stars <- stars.pval(kenw$coefficients[8])
  mtext(text = "W : R² = NS", side=3, cex=.7, adj=1)
  
  # abline(mw, col = "black", untf=T , lty = 2)
  
  me = lm(vars[E,2]~twi[E])
  kene = summary(me) 
  
  t <-signif(kene$r.squared, 2)
  stars <- stars.pval(kene$coefficients[8])
  mtext(text = substitute(paste("     E : R² = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  abline(me, col = "grey50", untf=T , lty = 2)
  
  kt = kruskal.test(vars[,2] ~ as.numeric(as.factor(env_var$SIDE) ) )
  
  par(mar= c(4,0,1,1))
  
  plot(vars[,2] ~ as.factor(env_var$SIDE), axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "", ylim = range(df_seg[[2]][,c(2,4)]))
  axis(1, at = c(1,2), labels = c("E", "W"))
  stars <- stars.pval(kt$p.value)
  
  mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.7, adj=1)
  




###### plot SES all #####



X11( width= 10 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F2.png", width = 900, height = 500, pointsize = 14)


layout(matrix(1:12, 2, 6, byrow = T),widths = rep(c(3,1),3) , heights= rep(1,2))

SES =  compare_all2$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars]

twi = log(env_var$TWI)
# twi = env_var$TWI

for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,0))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",vars_title[i], sep = " ") )
  points(SES_nona[E,i]~twi[E])
  points(SES_nona[W,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[W,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("    O : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  kene = cor.test(SES_nona[E,i],twi[E], method="kendall", use="pairwise") 
  
  t <-signif(kene$estimate, 3)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
  par(mar= c(2.5,0,2.5,1))
  
  
  kt = kruskal.test( SES_nona[,i] ~ as.numeric(as.factor(env_var$SIDE) ) )
  
  plot(SES_nona[,i] ~ as.factor(env_var$SIDE), axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "")
  axis(1, at = c(1,2), labels = c("East", "West"))
  stars <- stars.pval(kt$p.value)
  
  mtext(text = paste("p =", signif(kt$p.value,digits = 3), stars), side=3, cex=.6, adj=1)
  
}


dev.off()

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
####### TOPO EST OUEST #####
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

##### RLQ west #####


compwest = list()

# with abundances

compwest$spe <- data.frame(community_west, check.names = F)


#   compwest$spe = data.frame(ifelse(community_west == 0, 0, 1)) # transform direct la matrix abundance en pers / abs


compwest$traits = TF_west

# env_var_select = data.frame(env_var_west[,c("TWI", "INS")])

# models du 28/12/2017 avec ( "TWI","INS", "CURV")


# log transform LA & SLA

# compwest$traits$SLA = log(compwest$traits$SLA)

# compwest$traits$LA = log(compwest$traits$LA)

# compwest$traits$LT = 1/(compwest$traits$SLA*compwest$traits$LDMC)




env_var_select = env_var_west[,c( "TWI","INS", "CURV")]    # OK 

compwest$env = env_var_select


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.west <- dudi.coa(compwest$spe, scannf = FALSE)

acpR.west <- dudi.hillsmith(as.data.frame(compwest$env), row.w = afcL.west$lw,
                            scannf = FALSE)


acpQ.west <- dudi.pca(compwest$traits, row.w = afcL.west$cw,
                      scannf = FALSE)


s.arrow(acpQ.west$co)

#### RLQ ####

rlq.west <- rlq(acpR.west, afcL.west, acpQ.west,
                scannf = FALSE)


100*rlq.west$eig/sum(rlq.west$eig)


##### 4corner #####

nrepet <- 9999

four.comb.west <- fourthcorner(compwest$env, compwest$spe,
                               compwest$traits, modeltype = 6, p.adjust.method.G = "none", 
                               p.adjust.method.D = "none", nrepet = nrepet)

plot(four.comb.west, alpha = 0.05, stat="D2")

four.comb.west

Srlq <- fourthcorner2(compwest$env, compwest$spe, compwest$traits,
                      modeltype = 6,  p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ


###### plot RLQ west abund #####

percent_axis = signif(100*rlq.west$eig/sum(rlq.west$eig), digits = 4)

par(mar = c(4,4,1,1))
plot(0, 0, xlim = c(-.50, .55), ylim = c(-.15,.22), type = "n", 
     xlab = paste("RLQ axis 1 (",percent_axis[1],"%)",sep = ""), 
     ylab = paste("RLQ axis 2 (",percent_axis[2],"%)",sep = ""), xaxt='n', yaxt='n')

axis(1, seq(-.5,.75,.25))
axis(2, seq(-.4,.4,.1))

abline(h=seq(-.4,.4,.1), v=seq(-.5,.75,.25), col="gray", lty=3)
abline(h=seq(0), v=seq(0), col="black", lty=1)

# s.arrow(rlq.all_presabs$co,add.plot = T)
s.arrow2(rlq.west$co[c(2,4),], add.plot = T, labfont = 2, box = F)
s.arrow2(rlq.west$co[c(1,3),], colboxes = 'gray60',  collab = 'gray60', arrowcol = 'gray60', add.plot = T, labfont = 2, box = F)

r_var = rownames(rlq.west$li)


text(rlq.west$li$Axis1, rlq.west$li$Axis2, r_var, col = c("black", 'gray60','gray60'), pos = c(3,3,2), font = 2)

arrow1(0,0,rlq.west$li$Axis1[1], rlq.west$li$Axis2[1], edge = T, arrowcol = 'black', lty = 2)
arrow1(0,0,rlq.west$li$Axis1[2], rlq.west$li$Axis2[2], edge = T, arrowcol = 'gray60', lty = 2)
arrow1(0,0,rlq.west$li$Axis1[3], rlq.west$li$Axis2[3], edge = T, arrowcol = 'gray60', lty = 2)


# points(rlq.all_presabs$li$Axis1, rlq.all_presabs$li$Axis2, col = c("gray50", 'gray80','gray50','black', 'black'), pch = c(1,1,1,0,0))




##### RLQ east #####


compeast = list()

# with abundances

compeast$spe <- data.frame(community_east, check.names = F)


#   compeast$spe = data.frame(ifelse(community_east == 0, 0, 1)) # transform direct la matrix abundance en pers / abs


compeast$traits = TF_east

# log transform LA & SLA

# compeast$traits$SLA = log(compeast$traits$SLA)

# compeast$traits$LA = log(compeast$traits$LA)

# compeast$traits$LT = 1/(compeast$traits$SLA*compeast$traits$LDMC)



# env_var_select = data.frame(env_var_east[,c("TWI", "CURV")])

# models du 28/12/2017 avec ( "TWI","INS", "CURV")

env_var_select = env_var_east[,c( "TWI","INS", "CURV")]    # OK 

compeast$env = env_var_select


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.east <- dudi.coa(compeast$spe, scannf = FALSE)

acpR.east <- dudi.hillsmith(as.data.frame(compeast$env), row.w = afcL.east$lw,
                            scannf = FALSE)


acpQ.east <- dudi.pca(compeast$traits, row.w = afcL.east$cw,
                      scannf = FALSE)

s.arrow(acpQ.east$co)

#### RLQ ####

rlq.east <- rlq(acpR.east, afcL.east, acpQ.east,
                scannf = FALSE)


100*rlq.east$eig/sum(rlq.east$eig)


##### 4corner #####

nrepet <- 9999

four.comb.east <- fourthcorner(compeast$env, compeast$spe,
                               compeast$traits, modeltype = 6, p.adjust.method.G = "none", p.adjust.method.D = "none", nrepet = nrepet)
                               

plot(four.comb.east, alpha = 0.05, stat="D2")

four.comb.east


Srlq <- fourthcorner2(compeast$env, compeast$spe, compeast$traits,
                       p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ

testrlq.east <- randtest(rlq.east, modeltype = 6, nrepet = nrepet)

plot(four.comb.east, x.rlq = rlq.east, alpha = 0.05,
     stat = "D2", type = "biplot")

###### plot RLQ east abund #####

percent_axis = signif(100*rlq.east$eig/sum(rlq.east$eig), digits = 4)

par(mar = c(4,4,1,1))
plot(0, 0, xlim = c(-.5, .5), ylim = c(-.12,.15), type = "n", 
     xlab = paste("RLQ axis 1 (",percent_axis[1],"%)",sep = ""), 
     ylab = paste("RLQ axis 2 (",percent_axis[2],"%)",sep = ""), xaxt='n', yaxt='n')

axis(1, seq(-.5,.75,.25))
axis(2, seq(-.4,.4,.1))

abline(h=seq(-.4,.4,.1), v=seq(-.5,.75,.25), col="gray", lty=3)
abline(h=seq(0), v=seq(0), col="black", lty=1)

# s.arrow(rlq.all_presabs$co,add.plot = T)
s.arrow2(rlq.east$co[c(1),], add.plot = T, labfont = 2, box = F)
s.arrow2(rlq.east$co[c(2:4),], colboxes = 'gray60',  collab = 'gray60', arrowcol = 'gray60', add.plot = T, labfont = 2, box = F)
s.arrow2(rlq.east$co[c(5),], colboxes = 'gray60',  collab = 'gray60', arrowcol = 'gray60', add.plot = T, labfont = 2, box = F)

r_var = rownames(rlq.east$li)


text(rlq.east$li$Axis1, rlq.east$li$Axis2, r_var, col = c("black", 'gray60','black'), pos = c(2,3,4), font = 2)

arrow1(0,0,rlq.east$li$Axis1[1], rlq.east$li$Axis2[1], edge = T, arrowcol = 'black', lty = 2)
arrow1(0,0,rlq.east$li$Axis1[2], rlq.east$li$Axis2[2], edge = T, arrowcol = 'gray60', lty = 2)
arrow1(0,0,rlq.east$li$Axis1[3], rlq.east$li$Axis2[3], edge = T, arrowcol = 'black', lty = 2)


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
####### alternative = on ne prend que SES pour l'axe 1 de la RLQ car c'est le seul significatif ######
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
####### WEST ######
names(stat_obs_west$tab_stats_obs)

vecvar_w = c(13,27,3)

SES =  compare_west$SES_stat_obs

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

# SES_nona = SES_nona[names(SES_nona) %in% vars_w]

SES_nona = SES_nona[,c(2,3,1)]

names_var_ew = names(SES_nona)

vars_title_w =  c( "Ranges Axis 1", "CWM Axis 1", "CWV Axis 1")
twi = env_var$TWI


df_seg_w 



layout(matrix(1:3, 1, 3, byrow = T),widths = rep(c(1),3) , heights= 1)

par(mar= c(4,4,2,1))

for (i in 1:ncol(SES_nona)){
  plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,i]), type = "n", xlab = "TWI", ylab = paste("SES",vars_title_w[i], sep = " ") ,
       log = 'x', cex.lab = 1.2)
  points(SES_nona[,i]~twi[W], pch = 16)
  
  df_seg_w = cbind(twi[W], rep(0,length(twi[W])), twi[W], SES_nona[,i])
  Segments(df_seg_w)
  
  kenw = cor.test(SES_nona[,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}






####### EAST ######


SES =  compare_east$SES_stat_obs

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% names_var_ew]

SES_nona = SES_nona[,c(2,3,1)]

vars_title_e =  c( "Ranges Axis 1", "CWM Axis 1", "CWV Axis 1")
twi = env_var$TWI



layout(matrix(1:3, 1, 3, byrow = T),widths = rep(c(1),3) , heights= 1)

par(mar= c(4,4,2,1))

for (i in 1:ncol(SES_nona)){
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", xlab = "TWI", ylab = paste("SES",vars_title_e[i], sep = " ") ,
       log = 'x', cex.lab = 1.2)
  points(SES_nona[,i]~twi[E], pch = 16)
  
  df_seg_e = cbind(twi[E], rep(0,length(twi[E])), twi[E], SES_nona[,i])
  Segments(df_seg_e)
  
  kene = cor.test(SES_nona[,i],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kene$estimate, 3)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}



##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

##### alternative 2 : tout à partir de la même RLQ avec LA et SLA log-transformés #####

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


# data 


compall = readRDS( file = "/home/thesardfou/Documents/thèse/papiers/structure_fonctionnelle_VS_précipitatation&topographie/data/data_for_RLQ.rds")

stat_obs_all = readRDS(file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all_15012018logsla&la.rds")

md_null_all2 =  readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all2_15012018logsla&la.rds")

compare_all2 = compare_stats_obs_null(stat_obs_all,md_null_all)

md_null_west2 = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west2_15012018logsla&la.rds")
md_null_east2 = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east2_15012018logsla&la.rds")


compare_west = compare_stats_obs_null(stat_obs_west,md_null_west)


stat_obs_east = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_east03012018.rds")

md_null_east = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east03012018.rds")

compare_all2 = compare_stats_obs_null(stat_obs_all,md_null_all2)

compare_all3 = compare_stats_obs_null(stat_obs_all,null_model_all3)

compare_west2 = compare_stats_obs_null(stat_obs_all,md_null_west2, sitlist = 1:20)


compare_east2 = compare_stats_obs_null(stat_obs_all,md_null_east2, sitlist = 21:40)



# log transform LA & SLA

compall$traits$SLA = log(compall$traits$SLA)

compall$traits$LA = log(compall$traits$LA)

compall$traits$LT = 1/(compall$traits$SLA*compall$traits$LDMC)


env_var_select = compall$env[,c( "TWI","INS", "CURV")]
# env_var_select = compall$env[,c( "TWI","INS", "CURV", "SIDE", "CLASS")]

#  env_var_select = compall$env[,c( "TWI", "SIDE")]

# env_var_select = compall$env[,c("SIDE", "CLASS")]

compall$covar = compall$env[, "SIDE"]


compall$env = env_var_select


# basic rlq

afcL.all <- dudi.coa(compall$spe, scannf = FALSE)

acpR.all <- dudi.hillsmith(compall$env, row.w = afcL.all$lw,
                           scannf = FALSE)

acpQ.all <- dudi.pca(compall$traits, row.w = afcL.all$cw,
                     scannf = FALSE)

rlq.all <- rlq(acpR.all, afcL.all, acpQ.all,
               scannf = FALSE)

plot(rlq.all )

Srlq <- fourthcorner2(compall$env, compall$spe, compall$traits,
                      modeltype = 6, p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ

testRaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testRaxes.comb.all, stat = "D")

testQaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testQaxes.comb.all, stat = "D")


wrlq.all = wca(rlq.all, fac = compall$covar, scannf = FALSE,nf=2)

brlq.all = bca(rlq.all, fac = compall$covar, scannf = FALSE,nf=2)

wrlq.all
brlq.all
plot(wrlq.all)
plot(brlq.all)

testRaxes.comb.wrlq.all <- fourthcorner.rlq(wrlq.all, modeltype = 6,
                                            typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                            p.adjust.method.D = "none")
print(testRaxes.comb.wrlq.all, stat = "D")

testQaxes.comb.wrlq.all <- fourthcorner.rlq(wrlq.all, modeltype = 6,
                                            typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                            p.adjust.method.D = "none")
print(testQaxes.comb.wrlq.all, stat = "D")









#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

###### alternative à partir de la RLQ globale  #####

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

stat_obs_all2 = readRDS(file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all_15012018logsla&la.rds")

md_null_all2 = readRDS(file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all2_15012018logsla&la.rds")

md_null_west2 = readRDS(file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west2_15012018logsla&la.rds")
md_null_east2 = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east2_15012018logsla&la.rds")


compare_all2 = compare_stats_obs_null(stat_obs_all2,md_null_all2)

compare_west2 = compare_stats_obs_null(stat_obs_all2,md_null_west2, sitlist = 1:20)


compare_east2 = compare_stats_obs_null(stat_obs_all2,md_null_east2, sitlist = 21:40)


##### ranges est et west #####

se= stat_obs_all2$rlq.abund$lQ[rownames(stat_obs_all2$rlq.abund$lQ) %in% colnames(community_east),]
sw = stat_obs_all2$rlq.abund$lQ[rownames(stat_obs_all2$rlq.abund$lQ) %in% colnames(community_west),]
range(sw$AxcQ1)
range(se$AxcQ1)
################################################################################################################################################################################################################

###### PLOT SES from rlq all abundances log la SLA ######

################################################################################################################################################################################################################




SES_nona = compare_all2$SES_stat_obs[,apply(is.na(compare_all2$SES_stat_obs),2,sum) == 0]

colnames(SES_nona)

vecvar = c(1,2,24,25,11,12,9,10,22)

length(vecvar)
colnames(SES_nona)[vecvar]
colnames(SES_nona)[i]
sesnames = c("Ranges Axis 1" , "Ranges Axis 2" , "CWM Axis 1" , "CWM Axis 2" , "Kurtosis Axis 1", "Kurtosis Axis 2", "Skewness Axis 1" , "Skewness Axis 2", "Fdis")

facforest = factor(c(rep("DMF",20), rep("HF",20)))

# layout(matrix(c(1:8,9,9), 5, 2, byrow = T),widths = rep(1,2) , heights= rep(1,5))

layout(matrix(c(rep(1,5),2,4,6,8,10,3,5,7,9,10), 5, 3, byrow = F),widths = c(1,4,4) , heights= rep(1,5))
par(mar= c(0,0,0,0))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

text(x = 0.5, y = 0.5,"Strandard Effect Size", xpd=TRUE , srt=+90, font = 2, cex = 1.3)

for(i in 1:(length(vecvar)-1)){
  
  kt = kruskal.test( SES_nona[,vecvar[i]] ~ facforest )
  
  par(mar= c(2,2,2.7,1))
  
  plot.new()
  plot.window(xlim=c(.5,2.5), ylim = (range(SES_nona[,vecvar[i]])+c(-.5,.5)))
  
  points(as.numeric(facforest), SES_nona[,vecvar[i]], pch = 1)
  box()
  axis(2)
  axis(1, at = c(1,2), labels = c("DMF", "HF"))
  title(sesnames[i], line = 1.2)
  stars <- stars.pval(kt$p.value)
  
  mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.7, adj=1)
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}

par(mar= c(2,7.5,2.7,7.5))

kt = kruskal.test( SES_nona[,vecvar[9]] ~ facforest )

plot.new()
plot.window(xlim=c(.5,2.5), ylim = (range(SES_nona[,vecvar[9]])+c(-.5,.5)))

points(as.numeric(facforest), SES_nona[,vecvar[9]], pch = 1)
box()
axis(2)
axis(1, at = c(1,2), labels = c("DMF", "HF"))
title(sesnames[9], line = 1)
stars <- stars.pval(kt$p.value)

mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.7, adj=1)
abline(h=0)

abline(h=1.96,lty = 3)
abline(h=-1.96,lty = 3)



################################################################################################################################################################################################################

###### PLOT SES from rlq all abundances log la SLA ===> West & Est from RLQ ALL, permutation intra-forest type ######

################################################################################################################################################################################################################



###### plot SES west #####

SES =  compare_west2$SES_stat_obs

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

# choisir un vecteur vars_w pour selectionner les stats à afficher

vecvar = c(33,34,24,25,3,4)

length(vecvar_w)

names(SES_nona)

vars_w = names(compare_west2$SES_stat_obs)[vecvar]

vars_title_w = vars_w

vars_title_w =  c( "Ranges Axis 1", "Range Axis 2", "CWM Axis 1", "CWM Axis 2", "CWV Axis 1", "CWV Axis 2")


#SES_nona = SES_nona[names(SES_nona) %in% vars_w]

twi = env_var$TWI

layout(matrix(1:6, 3, 2, byrow = T),widths = rep(c(1),2) , heights= rep(c(1),3))

par(mar= c(4,4,2,1))

# i in 1:ncol(SES_nona) # pour tout voir

for (i in 1:length(vecvar)){
 # pour tout voir
  # plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,i]), type = "n", xlab = "TWI", ylab = paste("SES",colnames(SES_nona)[i], sep = " ") , log = 'x', cex.lab = 1.2)
  # pour les noms au propre
   plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,vecvar[i]]), type = "n", xlab = "TWI", ylab = paste("SES",vars_title_w[i], sep = " "),
        log = 'x', cex.lab = 1.2, font.lab = 2)
  points(SES_nona[,vecvar[i]]~twi[W], pch = 16)
  
  df_seg_w = cbind(twi[W], rep(0,length(twi[W])), twi[W], SES_nona[,vecvar[i]])
  Segments(df_seg_w)
  
  kenw = cor.test(SES_nona[,vecvar[i]],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
  
  
  
}





###### plot SES east #####

SES =  compare_east2$SES_stat_obs

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]


layout(matrix(1:6, 3, 2, byrow = T),widths = rep(c(1),2) , heights= rep(c(1),3))

par(mar= c(4,4,2,1))

# i in 1:ncol(SES_nona) # pour tout voir

for (i in 1:length(vecvar)){
  # pour tout voir
  # plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", xlab = "TWI", ylab = paste("SES",colnames(SES_nona)[i], sep = " ") , log = 'x', cex.lab = 1.2)
  # pour les noms au propre
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,vecvar[i]]), type = "n", xlab = "TWI", ylab = paste("SES",vars_title_w[i], sep = " "),
       log = 'x', cex.lab = 1.2, font.lab = 2)
  points(SES_nona[,vecvar[i]]~twi[E], pch = 16)
  
  df_seg_w = cbind(twi[E], rep(0,length(twi[E])), twi[E], SES_nona[,vecvar[i]])
  Segments(df_seg_w)
  
  kenw = cor.test(SES_nona[,vecvar[i]],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}

























































##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

##### autres plots #####

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################





###### est avec curvature #####
 twi = env_var$CURV


layout(matrix(1:3, 1, 3, byrow = T),widths = rep(c(1),3) , heights= 1)

par(mar= c(4,4,2,1))

for (i in 1:ncol(SES_nona)){
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", xlab = "CURV", ylab = paste("SES",vars_title_e[i], sep = " ") ,
        cex.lab = 1.2)
  points(SES_nona[,i]~twi[E], pch = 16)
  
  df_seg_e = cbind(twi[E], rep(0,length(twi[E])), twi[E], SES_nona[,i])
  Segments(df_seg_e)
  
  kene = cor.test(SES_nona[,i],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kene$estimate, 3)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.7, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}




##### FILTRE WEST #####

#### plot obs values  west  #####




# compare est - west : on prend quoi? 

# vector colnames : 

names(stat_obs_west$tab_stats_obs)

vecvar_w = c(3,4,5,6,9,10,11,12,13,14,15,16,17,18,27,28,35,36)


vecvar_w = c(3,13,27)


length(vecvar_w)

vars_w=names(stat_obs_west$tab_stats_obs)[vecvar_w]

vars_title_w = vars_w

vars_title_w =  c("CWV Axis 1", "CWV Axis 2", "Ranges Axis 1", "Range Axis 2", "CWM Axis 1", "CWM Axis 2")

library(lattice)
library(gtools)

X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F1w.png", width = 900, height = 500, pointsize = 14)

layout(matrix(1:6, 2, 3, byrow = F),widths = rep(c(1),3) , heights= rep(1,2))

#layout.show(2)

twi = log(env_var$TWI)

for(i in vecvar_w){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  

  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_west$tab_stats_obs[,i]), main = vars_title_w[which(vecvar_w == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_w == i)], line = 1, cex = 2)
  points(stat_obs_west$tab_stats_obs[,i]~twi[W], pch = 16)
  
  kenw = cor.test(stat_obs_west$tab_stats_obs[,i],twi[W], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  

}




dev.off()


###### plot SES west #####



X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F2_w.png", width = 900, height = 500, pointsize = 14)


layout(matrix(1:6, 2, 3, byrow = T),widths = rep(c(1),3) , heights= rep(1,2))

SES =  compare_west$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars_w]

twi = env_var$TWI
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",vars_title_w[i], sep = " ") )
  points(SES_nona[,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}


dev.off()


##### FILTRE EAST #####

#### plot obs values  east  #####



# vector colnames : 

names(stat_obs_east$tab_stats_obs)

vecvar_e = c(3,4,5,6,9,10,11,12,13,14,15,16,17,18,27,28,35,36)


vecvar_e = c(5,6,13,14,35,36)


length(vecvar_e)

vars_e=names(stat_obs_east$tab_stats_obs)[vecvar_e]

vars_title_e = vars_e

vars_title_e =  c("CV A1", "CV A2", "Ranges A1", "Range A2", "CM A1", "CM A2")

library(lattice)
library(gtools)

X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F1e.png", width = 900, height = 500, pointsize = 14)

layout(matrix(1:6, 2, 3, byrow = T),widths = rep(c(1),3) , heights= rep(1,2))

#layout.show(2)

twi = log(env_var$TWI)
# twi = env_var$TWI

for(i in vecvar_e){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_east$tab_stats_obs[,i]), main = vars_title_e[which(vecvar_e == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_e == i)], line = 1, cex = 2)
  points(stat_obs_east$tab_stats_obs[,i]~twi[E], pch = 16)
  
  kenw = cor.test(stat_obs_east$tab_stats_obs[,i],twi[E], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}




dev.off()


###### plot SES west #####



X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F2_e.png", width = 900, height = 500, pointsize = 14)


layout(matrix(1:6, 2, 3, byrow = T),widths = rep(c(1),3) , heights= rep(1,2))

SES =  compare_east$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars_w]

twi = log(env_var$TWI)
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",vars_title_e[i], sep = " ") )
  points(SES_nona[,i]~twi[E], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}


dev.off()














##### shape WEST #####

#### plot obs values  west  #####




# vector colnames : 

names(stat_obs_west$tab_stats_obs)

vecvar_w = c(3,4,5,6,9,10,11,12,13,14,15,16,17,18,27,28,35,36)


vecvar_w = c(3,4,9,10,11,12,23,25,27,28)


length(vecvar_w)

vars_w=names(stat_obs_west$tab_stats_obs)[vecvar_w]

vars_title_w = vars_w

vars_title_w =  c("CWV A1", "CWV A2", "Skewness A1", "Skewness A2", "Kurtosis A1", "Kurtosis A2", "FEve", "FDis", "CWM A1", "CWM A2")

library(lattice)
library(gtools)

X11( width= 4 , height= 10, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F3w.png", width = 400, height = 1000, pointsize = 14)

layout(matrix(1:10, 5, 2, byrow = T),widths = rep(c(1),5) , heights= rep(1,2))

#layout.show(2)

twi = log(env_var$TWI)

for(i in vecvar_w){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_west$tab_stats_obs[,i]), main = vars_title_w[which(vecvar_w == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_w == i)], line = 1, cex = 2)
  points(stat_obs_west$tab_stats_obs[,i]~twi[W], pch = 16)
  
  kenw = cor.test(stat_obs_west$tab_stats_obs[,i],twi[W], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}




dev.off()


###### plot SES west #####



X11( width= 4 , height= 10, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F4w.png", width = 400, height = 1000, pointsize = 14)


layout(matrix(1:10, 5, 2, byrow = T),widths = rep(c(1),5) , heights= rep(1,2))

SES =  compare_west$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars_w]

twi = log(env_var$TWI)
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",vars_title_w[i], sep = " ") )
  points(SES_nona[,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}


dev.off()




##### shape EAST #####

#### plot obs values  east  #####



# vector colnames : 

names(stat_obs_east$tab_stats_obs)

vecvar_e = c(3,4,5,6,9,10,11,12,13,14,15,16,17,18,27,28,35,36)

vecvar_e = c(3,4,9,10,11,12,23,25,27,28)


length(vecvar_e)

vars_e=names(stat_obs_west$tab_stats_obs)[vecvar_e]


vars_title_e =  c("CWV A1", "CWV A2", "Skewness A1", "Skewness A2", "Kurtosis A1", "Kurtosis A2", "FEve", "FDis", "CWM A1", "CWM A2")

library(lattice)
library(gtools)

X11( width= 4 , height= 10, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F3e.png", width = 400, height = 1000, pointsize = 14)

layout(matrix(1:10, 5, 2, byrow = T),widths = rep(c(1),5) , heights= rep(1,2))

#layout.show(2)

twi = log(env_var$TWI)

for(i in vecvar_e){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_east$tab_stats_obs[,i]), main = vars_title_e[which(vecvar_e == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_e == i)], line = 1, cex = 2)
  points(stat_obs_east$tab_stats_obs[,i]~twi[E], pch = 16)
  
  kenw = cor.test(stat_obs_east$tab_stats_obs[,i],twi[E], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}




dev.off()


###### plot SES west #####



X11( width= 4 , height= 10, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F4e.png", width = 400, height = 1000, pointsize = 14)

layout(matrix(1:10, 5, 2, byrow = T),widths = rep(c(1),5) , heights= rep(1,2))

SES =  compare_east$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars_e]

twi = log(env_var$TWI)
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",vars_title_e[i], sep = " ") )
  points(SES_nona[,i]~twi[E], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}


dev.off()


























