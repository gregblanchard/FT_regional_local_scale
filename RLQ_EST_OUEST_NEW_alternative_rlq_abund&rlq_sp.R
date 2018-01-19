######################################################################################################################################################################################

###### RLQ MODEL NULL : TRAITS = POSITIONS ON RLQ AXES 

######################################################################################################################################################################################



##### initialise data #####





library(lattice)
library(gtools)


library(vegan)
library(ade4)

#####  load data ##### 

# occurences

FBA <- read.csv("/home/thesardfou/Documents/data/tabs/Fiche terrain FBA_saisie_10_10_17.csv", sep = ",")

# env var


# tab modif avev points modif dans ""/home/thesardfou/Documents/GIS/echantillonage_kone"
# tab_extract3 <- read.csv('/home/thesardfou/Documents/data/tab_extract3-all.csv')

# shp points kone dans le meme fichier
# tab_extract3 <- read.csv("/home/thesardfou/Documents/GIS/echantillonage_kone/points_all_save/tab_extract3-all.csv")

# REFERENCE ! (pts .shp ref : /home/thesardfou/Documents/GIS/echantillonage_kone/points_all_save/savefindec2017_REFERENCE/oldREFERENCE/points_done_all.shp)

# tab_extract3 <- read.csv("/home/thesardfou/Documents/GIS/echantillonage_kone/points_all_save/tab_extract3-all_findec2017.csv")

# NEW REFERENCE ! (pts .shp ref : /home/thesardfou/Documents/GIS/echantillonage_kone/points_all_save/savefindec2017_REFERENCE/points_done_all.shp)


tab_extract3 <- read.csv('/home/thesardfou/Documents/data/tab_extract3-all_fin2017.csv')


# TF

FT_all_kone = read.csv( '/home/thesardfou/Documents/data/FT_all_kone_111017.csv')


##### make clean ####


#enlève les lignes avec les morts
FBA<-subset(FBA,!FBA$strate==0)
FBA<-subset(FBA,!FBA$Taxon=="")


# which(FBA$Taxon == "Diospyros" )



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
# sites



sites <- unique(FBA$id_pts)
nsites <- length(sites)


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

dim(TF_kone_temp_nona)






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

community_nona = data.frame(community_all, check.names = F)

compall = list()

# with abundances

compall$spe <- data.frame(community_nona, check.names = F)

#   compall$spe = data.frame(ifelse(community_nona == 0, 0, 1)) # transform direct la matrix abundance en pers / abs



compall$traits = TF_kone_com_nona

# compall$traits$SLA = log(compall$traits$SLA)

# compall$traits$LA =log(compall$traits$LA)

# compall$traits$LDMC =log(compall$traits$LDMC)

# variables pour les models du 28/12/2017  : ("PREC", "TWI","INS", "CURV")

# env_var_select = env_var[,c("PMA", "TWI","INS", "CURV")]


# env_var_select = env_var[,c( "TWI","INS", "CURV", "TWIxPREC", "SIDE")]


# model du 03012018 ( "TWI","INS", "CURV", "SIDE")

# model du 10012018 ( "SIDE", "CLASS" )

# 
env_var_select = env_var[,c( "TWI","INS", "CURV", "SIDE")]

####

# env_var_select = env_var[,c("SIDE", "CLASS")]

#  env_var_select = data.frame(SIDE = env_var[,c("SIDE")], nul = rep(0,40))

# test 2 variables twi

#  twiw = c(env_var$TWI[W],rep(mean(env_var$TWI[W]),20))
#  twie = c(rep(mean(env_var$TWI[E]),20),env_var$TWI[E])

#  env_var_select = data.frame(PREC = env_var[,c("SIDE")], twiw , twie)



# env_var_select = env_var[,c("INS","PREC","CURV", "CLASS")]

compall$env = env_var_select


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.all <- dudi.coa(compall$spe, scannf = FALSE)

acpR.all <- dudi.hillsmith(compall$env, row.w = afcL.all$lw,
                           scannf = FALSE)


acpQ.all <- dudi.pca(compall$traits, row.w = afcL.all$cw,
                     scannf = FALSE)

#### RLQ ####

rlq.all <- rlq(acpR.all, afcL.all, acpQ.all,
               scannf = FALSE)


100*rlq.all$eig/sum(rlq.all$eig)

summary(rlq.all)

plot(rlq.all)


#######################################################################################################################
##### FOURTH CORNER ALL DATA : la relation entre traits en environement est-elle significative ??  #####
#######################################################################################################################


nrepet <- 9999

four.comb.all <- fourthcorner(compall$env, compall$spe,
                              compall$traits, modeltype = 6, p.adjust.method.G = "none", 
                              p.adjust.method.D = "none", nrepet = nrepet) 


four.comb.all
four.comb.all.adj <- p.adjust.4thcorner(four.comb.all,
                                        p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")

four.comb.all.adj
plot(four.comb.all.adj, alpha = 0.05, stat="D2")


plot(four.comb.all, alpha = 0.05, stat="D2")


plot(four.comb.all, x.rlq = rlq.all, alpha = 0.05,
     stat = "D2", type = "biplot")

Srlq <- fourthcorner2(compall$env, compall$spe, compall$traits,
                      modeltype = 6, p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ

plot(Srlq, stat = "G")




testRaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")

print(testRaxes.comb.all, stat = "D")


testQaxes.comb.all <- fourthcorner.rlq(rlq.all, modeltype = 6,
                                       typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")

print(testQaxes.comb.all, stat = "D")

testrlq.all <- randtest(rlq.all, modeltype = 6, nrepet = nrepet)

testrlq.all
plot(testrlq.all)


par(mfrow = c(1, 2))
plot(testQaxes.comb.all, alpha = 0.05, type = "table",
     stat = "D2")
plot(testRaxes.comb.all, alpha = 0.05, type = "table",
     stat = "D2")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(testQaxes.comb.all, alpha = 0.05, type = "biplot",
     stat = "D2")
plot(testRaxes.comb.all, alpha = 0.05, type = "biplot",
     stat = "D2")
par(mfrow = c(1, 1))








#######################################################################################################################
##### NULL MODEL  = coordordinates on RLQ = ps traits ==> compute stats ==> swap community matrix (null model) ==> compare (SES)  #####
#######################################################################################################################

######   stat obs all #####
stat_obs_all = tf_stats_from_rlq4(comp = compall ,class = env_var$CLASS, class2 = env_var$SIDE)

# saveRDS(stat_obs_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all10012018.rds")

##### model null all #####

md_null_all = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, class = env_var$CLASS, 
                        class2 = env_var$SIDE, tab_stats_obs = stat_obs_all$tab_stats_obs, compare = F)



# saveRDS(md_null_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all10012018.rds")


#      md_null_all$list_mean_stats_classe2_null = list_mean_stats_classe2_null
#     md_null_all = readRDS(md_null_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all28122017.rds")

#      md_null_all = readRDS(md_null_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all21122017.rds")

##### compare null & obs #####


compare_all = compare_stats_obs_null(stat_obs_all,md_null_all)


##### compare classe SIDE ####


##### plot RLQ ####


plot(rlq.all)
plot(stat_obs_all$mean_stats_classe$CWM.AxcQ2~stat_obs_all$mean_stats_classe$CWM.AxcQ1)

s.arrow(rlq.all$co)
points(stat_obs_all$mean_stats_classe$CWM.AxcQ2~stat_obs_all$mean_stats_classe$CWM.AxcQ1)


#######################################################################################################################
##### test représentation  de groupes fonctionels  (from kleyer JVS 2012 )  #####
#######################################################################################################################





s.label(rlq.all$lQ, clabel = 0)

par(mar = c(0.1, 0.1, 0.1, 0.1))
library(maptools)
pointLabel(rlq.all$lQ,row.names(rlq.all$lQ), cex=0.7)


# pour virer le dendrocnide qui fait chier
# rlq.all$lQ = rlq.all$lQ[!rownames(rlq.all$lQ) == "Dendrocnide latifolia",]



# FD_test = dbFD(rlq.all$lQ, calc.FGR = T, clust.type = "kmeans",km.crit = c("calinski"), km.iter = 100)


FD_test = dbFD(x=rlq.all$lQ,corr = "cailliez",calc.CWM = TRUE,calc.FGR = T,clust.type = "kmeans",km.sup.gr=10) 

# ou : 


# hc2 <- hclust(dist(rlq.all$lQ), method = "ward")
# spe.group2 <- as.factor(cutree(hc2, k = 4))
# plot(hc2)



spe.group2 <- as.factor(FD_test$spfgr)
# pout changer le nom des facteurs #
# levels(spe.group2) <- c("C","B","D","A")

# spe.group2 <- factor(spe.group2, levels=c("A","B","C","D"))

s.class(rlq.all$lQ, spe.group2, col= 1:nlevels(spe.group2))


s.arrow(rlq.all$c1, add.plot = T,clab=0.8)



pointLabel(rlq.all$lQ,row.names(rlq.all$lQ), cex=0.7)
#pointLabel(rlq.all$lQ,row.names(rlq.all$lQ), cex=1.2)

s.arrow(rlq.all$l1, add.plot = T,clab=0.8)

pointLabel(stat_obs_all$mean_stats_classe$CWM.AxcQ2~stat_obs_all$mean_stats_classe$CWM.AxcQ1, levels(env_var$CLASS))


hc2 <- hclust(dist(rlq.all$lQ), method = "ward")
plot(hc2)


# saveRDS(compare_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/compare_all21122017.rds")


#######################################################################################################################
##### test ANOVA sur les axes de la RLQ #####
#######################################################################################################################


df_aov = cbind(compare_all$SES_stat_obs, SIDE = env_var$SIDE, CLASS = env_var$CLASS)

names(df_aov)

fit <- aov(CWM.AxcQ1 ~ SIDE, data=df_aov)

fit <- aov(CWM.AxcQ2 ~ CLASS, data=df_aov)

summary(fit)

#### aov 2 facteurs ####

fit2 <- aov(CWM.AxcQ1 ~ SIDE + SIDE/CLASS, data=df_aov)
fit2 <- aov(CWM.AxcQ1 ~ SIDE + Error(CLASS), data=df_aov)


summary(fit2)
plot(fit2)
library(TukeyC)

tuk = TukeyC(df_aov,
             model = 'CWM.AxcQ1 ~ SIDE/CLASS',
             which = 'CLASS',
             fl1=1,
             sig.level = 0.05)

summary(tuk)

library(lme4)
data.nest.lme <- lmer(CWM.AxcQ1 ~ SIDE+(1|CLASS), df_aov, REML=TRUE)
anova(data.nest.lme)

#######################################################################################################################
##### test ANOVA non parametrique : sur les axes de la RLQ #####
#######################################################################################################################
df_aov = cbind(compare_all$SES_stat_obs, SIDE = env_var$SIDE, CLASS = env_var$CLASS)

names(df_aov)

fit <- kruskal.test(CWM.AxcQ1 ~ SIDE, data=df_aov)
fit 
fit <- kruskal.test(CWM.AxcQ1 ~ CLASS, data=df_aov)
fit

pairwise.wilcox.test(df_aov$CWM.AxcQ1, df_aov$CLASS,
                     p.adjust.method = "BH")

#######################################################################################################################
##### test kendall corelation (non parametrique) : sur les axes de la RLQ #####
#######################################################################################################################


cor.test(compare_all$SES_stat_obs$CWM.AxcQ1 ,  env_var$PREC, method="kendall") 
df_lm = cbind(compare_all$SES_stat_obs, env_var)




#######################################################################################################################
##### dissimilaryty funct vs floristic ====> voir script " BC " #####
#######################################################################################################################
#######################################################################################################################
##### rareté uniqueness (violle et al) #####
#######################################################################################################################


rlq.all$lQ
mindist =function(x){a =sort(x)
return(a[2])} 

dist_sp_pool = as.matrix(dist(rlq.all$lQ, method = "euclidean"))

mindist_sp_pool = apply(dist_sp_pool,2,mindist)

abund_sp_pool = aggregate(FBA_nona$Taxon, by= list(FBA_nona$Taxon), length)

plot(mindist_sp_pool~abund_sp_pool$x, log = "x")

md = lm(mindist_sp_pool~abund_sp_pool$x)
abline(md)
summary(md)

# west
dist_sp_pool = as.matrix(dist(rlq.west$lQ, method = "euclidean"))

mindist_sp_pool = apply(dist_sp_pool,2,mindist)

FBA_west = FBA_nona[FBA_nona$id_pts %in% env_var$id[W],]

abund_sp_pool = aggregate(FBA_west$Taxon, by= list(FBA_west$Taxon), length)

plot(mindist_sp_pool~abund_sp_pool$x, log = "x")

md = lm(mindist_sp_pool~abund_sp_pool$x)
abline(md)
summary(md)

# east
dist_sp_pool = as.matrix(dist(rlq.east$lQ, method = "euclidean"))

mindist_sp_pool = apply(dist_sp_pool,2,mindist)

FBA_west = FBA_nona[FBA_nona$id_pts %in% env_var$id[E],]

abund_sp_pool = aggregate(FBA_west$Taxon, by= list(FBA_west$Taxon), length)

plot(mindist_sp_pool~abund_sp_pool$x, log = "x")

md = lm(mindist_sp_pool~abund_sp_pool$x)
abline(md)
summary(md)

#######################################################################################################################
##### rareté distainctiveness (violle et al) #####
#######################################################################################################################

list_fdistinc = list()
list_abund_com =list()
for (i in 1:length(stat_obs_all$list_scores_sp_com)){
  
  tmp = as.matrix(dist(stat_obs_all$list_scores_sp_com[[i]], method = "euclidean"))
  tmp[tmp == 0] = NA
  
  list_fdistinc[[i]] = apply(tmp,1,mean, na.rm=TRUE)
  list_abund_com[[i]] = community_all[i,community_all[i,]>0]
}


for (i in 1:40){
  plot(list_fdistinc[[i]]~list_abund_com[[i]])
  md = lm(list_fdistinc[[i]]~list_abund_com[[i]])
  print(summary(md))
}


#######################################################################################################################
##### plot SES #####
#######################################################################################################################
TFnum = data.frame(lapply(TF_kone_com_nona,as.numeric))


functional.beta.multi(com01,TFnum)

dim(com01[1:10,])

dim(TFnum)
##### plot ####


SES =  compare_all3$SES_stat_obs

plot(SES$Ranges_Ax1_abund~env_var$PREC)



SES_nona = SES[,apply(is.na(SES),2,sum) == 0]
twi = log(env_var$TWI)
#     twi = env_var$TWI
#### plot SES ####
par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(SES_nona[,i]), type = "n", main = names(SES_nona)[i])
  points(SES_nona[E,i]~twi[E])
  points(SES_nona[W,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[W,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("    O : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  kene =cor.test(SES_nona[E,i],twi[E], method="kendall", use="pairwise") 
  
  t <-signif(kene$estimate, 3)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}

#######################################################################################################################
##### plot stat obs #####
#######################################################################################################################
par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(stat_obs_all$tab_stats_obs)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_all$tab_stats_obs[,i]), type = "n", main = names(stat_obs_all$tab_stats_obs)[i])
  points(stat_obs_all$tab_stats_obs[E,i]~twi[E])
  points(stat_obs_all$tab_stats_obs[W,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(stat_obs_all$tab_stats_obs[W,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("    O : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  kene =cor.test(stat_obs_all$tab_stats_obs[E,i],twi[E], method="kendall", use="pairwise") 
  
  t <-signif(kene$estimate, 3)
  stars <- stars.pval(kene$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  

}


######   community canop & sous-bois #####

community_all_canop = community_all[]
community_all_canop[] = 0
community_all_sbois = community_all_canop
for(i in 1:length(unique(FBA_nona$id_pts))){
  
  all_tmp_canop = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate > 2  , ]
  
  sp_canop_com_tmp = unique(all_tmp_canop$Taxon)
  
  community_all_canop[i, colnames(community_all) %in% sp_canop_com_tmp] = community_all[i, colnames(community_all) %in% sp_canop_com_tmp]
  
  all_tmp_sbois = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate < 3  , ]
  
  sp_sbois_com_tmp = unique(all_tmp_sbois $Taxon)
  
  community_all_sbois[i, colnames(community_all) %in% sp_sbois_com_tmp] = community_all[i, colnames(community_all) %in% sp_sbois_com_tmp]
  
  
}




######   stat obs all canop #####
stat_obs_canop = tf_stats_from_rlq3(orlq = rlq.all , com = community_nona_canop, 
                                    class = env_var$CLASS)



##### model null canop #####

md_null_canop = tf_stats_from_rlq_permut_species(orlq = rlq.all, com = community_nona_canop, iter = 9, class = env_var$CLASS, 
                                                 tab_stats_obs = stat_obs_canop$tab_stats_obs, compare = F)


######   stat obs all sous bois #####
stat_obs_sbois = tf_stats_from_rlq3(orlq = rlq.all , com = community_nona_sbois, 
                                    class = env_var$CLASS)


##### model null sous bois #####

md_null_sbois = tf_stats_from_rlq_permut_species(orlq = rlq.all, com = community_nona_sbois, iter = 9, class = env_var$CLASS)


md_null_sbois = tf_stats_from_rlq_permut_species(orlq = rlq.all, com = community_nona_sbois, iter = 9, class = env_var$CLASS, 
                                                 tab_stats_obs = stat_obs_sbois$tab_stats_obs, compare = F)






###### relations topo ~ Ht & STRATIFICATION #####

nb_ind_sbois_all =  apply(community_all_sbois,1,sum)


htmax = aggregate(HT~id_pts, data=FBA_nona, max)

htmax = htmax[order(as.numeric(as.character(htmax$id_pts))),]


plot(nb_ind_sbois_all[W]~htmax$HT[W])
summary(lm(nb_ind_sbois_all[W]~htmax$HT[W]))


plot(nb_ind_sbois_all[E]~htmax$HT[E])
summary(lm(nb_ind_sbois_all[E]~htmax$HT[E]))


plot(nb_ind_sbois_all[W]~env_var$TWI[W])
summary(lm(nb_ind_sbois_all[W]~env_var$TWI[W]))

plot(nb_ind_sbois_all[E]~env_var$TWI[E])
summary(lm(nb_ind_sbois_all[E]~env_var$TWI[E]))



plot(htmax$HT[W]~env_var$TWI[W])
summary(lm(htmax$HT[W]~env_var$TWI[W]))



plot(htmax$HT[E]~env_var$TWI[E])
summary(lm(htmax$HT[E]~env_var$TWI[E]))










########################################################################################################################
########################################################################################################################

##### Null models on RLQ synthetic traits (positions on axes) for East and West communities #####

########################################################################################################################
########################################################################################################################


########################################################################################################################
##### west #####
########################################################################################################################

compwest = list()

# with abundances

compwest$spe <- data.frame(community_west, check.names = F)


#   compwest$spe = data.frame(ifelse(community_west == 0, 0, 1)) # transform direct la matrix abundance en pers / abs


compwest$traits = TF_west

# env_var_select = data.frame(env_var_west[,c("TWI", "INS")])

# models du 28/12/2017 avec ( "TWI","INS", "CURV")
# models du 03012018 avec ( "TWI","INS", "CURV")


env_var_select = env_var_west[,c( "TWI","INS", "CURV")]    # OK 

env_var_select = env_var_west[,c( "TWI", "CURV")]    # OK 


compwest$env = data.frame(env_var_select, check.names = F)


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.west <- dudi.coa(compwest$spe, scannf = FALSE)

acpR.west <- dudi.hillsmith(as.data.frame(compwest$env), row.w = afcL.west$lw,
                            scannf = FALSE)


acpQ.west <- dudi.pca(compwest$traits, row.w = afcL.west$cw,
                      scannf = FALSE)

#### RLQ ####

rlq.west <- rlq(acpR.west, afcL.west, acpQ.west,
                scannf = FALSE)


100*rlq.west$eig/sum(rlq.west$eig)

summary(rlq.west)

plot(rlq.west)
plot(rlq.west, xax = 1, yax = 2)

##### FOURTH CORNER WEST : la relation entre traits en environement est-elle significative ??  #####



nrepet <- 9999

four.comb.west <- fourthcorner(compwest$env, compwest$spe,
                               compwest$traits, modeltype = 6, p.adjust.method.G = "none", 
                               p.adjust.method.D = "none", nrepet = nrepet)

four.comb.west
four.comb.west.adj <- p.adjust.4thcorner(four.comb.west,
                                         p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")


plot(four.comb.west.adj, alpha = 0.05, stat="D2")


plot(four.comb.west, alpha = 0.05, stat="D2")


plot(four.comb.west, x.rlq = rlq.west, alpha = 0.05,
     stat = "D2", type = "biplot")



testQaxes.comb.west <- fourthcorner.rlq(rlq.west, modeltype = 6,
                                        typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                        p.adjust.method.D = "none")

print(testQaxes.comb.west, stat = "D")

testRaxes.comb.west <- fourthcorner.rlq(rlq.west, modeltype = 6,
                                        typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                        p.adjust.method.D = "none")
print(testRaxes.comb.west, stat = "D")


par(mfrow = c(1, 2))
plot(testQaxes.comb.west, alpha = 0.05, type = "table",
     stat = "D2")
plot(testRaxes.comb.west, alpha = 0.05, type = "table",
     stat = "D2")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(x.rlq = rlq.all, alpha = 0.05, type = "biplot",
     stat = "D2")
plot(testRaxes.comb.west, alpha = 0.05, type = "biplot",
     stat = "D2")
par(mfrow = c(1, 1))

Srlq <- fourthcorner2(compwest$env, compwest$spe,
                      compwest$traits,
                      modeltype = 6, p.adjust.method.G = "none", nrepet = nrepet)

Srlq$trRLQ

######   stat obs all #####
stat_obs_west = tf_stats_from_rlq4(comp = compwest, 
                                   class = env_var_west$CLASS)

stat_obs_west$tab_stats_obs
# saveRDS(stat_obs_west, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_west03012018.rds")

##### model null west from rlq west #####

md_null_west = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_west, iter = 999, class = env_var_west$CLASS, 
                                                compare = F)




 
# md_null_west = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west27122017.rds")

# md_null_west = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west28122017.rds")



library(gtools)


twi = log(env_var_west$TWI)
#     twi = env_var$TWI
#### plot Stats obs WEST  ####
par(mfrow=c(6,6),mar= c(2,2,2,2))
for(i in 1:ncol(stat_obs_west$tab_stats_obs)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(stat_obs_west$tab_stats_obs[,i]~ twi, xlim = range(twi), ylim = range(stat_obs_west$tab_stats_obs[,i]), main = names(stat_obs_west$tab_stats_obs)[i])
  
  
  kenw = cor.test(stat_obs_west$tab_stats_obs[,i],twi, method="kendall", use="pairwise") 
  
  t <-signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("     W : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  
  
  
}

SES_nona = compare_west$SES_stat_obs

SES_nona = SES_nona[,apply(is.na(SES_nona),2,sum) == 0]

par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(stat_obs_west$tab_stats_obs)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(SES_nona[,i]~ twi, xlim = range(twi), ylim = range(SES_nona[,i]), main = names(SES_nona)[i])
  
  
  kenw = cor.test(SES_nona[,i],twi, method="kendall", use="pairwise") 
  
  t <-signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("     W : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
  
}










######   community canop & sous-bois #####

community_west_canop = community_west[]
community_west_canop[] = 0
community_west_sbois = community_west_canop
for(i in 1:length(unique(FBA_nona$id_pts)[W])){
  
  west_tmp_canop = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate > 2  , ]
  
  sp_canop_com_tmp = unique(west_tmp_canop$Taxon)
  
  community_west_canop[i, colnames(community_west) %in% sp_canop_com_tmp] = community_west[i, colnames(community_west) %in% sp_canop_com_tmp]
  
  west_tmp_sbois = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate < 3  , ]
  
  sp_sbois_com_tmp = unique(west_tmp_sbois $Taxon)
  
  community_west_sbois[i, colnames(community_west) %in% sp_sbois_com_tmp] = community_west[i, colnames(community_west) %in% sp_sbois_com_tmp]
  
  
}


nb_ind_sbois_west =  apply(community_west_sbois,1,sum)



######   stat obs west canop #####
stat_obs_west_canop = tf_stats_from_rlq3(orlq = rlq.west , com = community_west_canop, 
                                         class = env_var$CLASS)



##### model null west canop #####

md_null_west_canop = tf_stats_from_rlq_permut_species(orlq = rlq.west, com = community_west_canop, iter = 9, class = env_var_west$CLASS, 
                                                      tab_stats_obs = stat_obs_west_canop$tab_stats_obs, compare = F)
##### compare west canop #####

compare_west_canop  = compare_stats_obs_null(stat_obs_west_canop,md_null_west_canop)


######   stat obs west sous bois #####
stat_obs_west_sbois = tf_stats_from_rlq3(orlq = rlq.west , com = community_west_sbois, 
                                         class = env_var_west$CLASS)


##### model null west sous bois #####


md_null_west_sbois = tf_stats_from_rlq_permut_species(orlq = rlq.west, com = community_west_sbois, iter = 9, class = env_var_west$CLASS, 
                                                      tab_stats_obs = stat_obs_west_sbois$tab_stats_obs, compare = F)

##### compare west sous bois #####


compare_west_sbois  = compare_stats_obs_null(stat_obs_west_sbois,md_null_west_sbois)




#### plot SES WEST   : compare ssbois et canop ####
dim(SES_west_sbois)
dim(SES_west_canop)

names(SES_west_sbois)
names(SES_west_canop)



SES_west = compare_west$SES_stat_obs[,apply(is.na(compare_west$SES_stat_obs),2,sum) == 0]
SES_west_canop = compare_west_canop$SES_stat_obs[,apply(!is.na(compare_west_canop$SES_stat_obs),2,sum) > 1]
SES_west_sbois = compare_west_sbois$SES_stat_obs[,apply(!is.na(compare_west_sbois$SES_stat_obs),2,sum) > 1]

# SES_west_sbois = compare_west_sbois$SES_stat_obs

twi = log(env_var_west$TWI)
#     twi = env_var$TWI
par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(SES_west[,i]), type = "n", main = names(SES_west)[i])
  #points(SES_west[,i]~ twi)
  
  stat_tmp = names(SES_west)[i]
  
  points(SES_west_canop[,stat_tmp]~twi)
  
  
  kenc = cor.test(SES_west_canop[,stat_tmp],twi, method="kendall", use="pairwise") 
  
  t <- signif(kenc$estimate, 3)
  stars <- stars.pval(kenc$p.value)
  mtext(text = substitute(paste("    C : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  if(stat_tmp %in% names(SES_west_sbois)){
    points(SES_west_sbois[,stat_tmp]~twi, pch = 16)
    
    
    kensb =cor.test(SES_west_sbois[,stat_tmp],twi, method="kendall", use="pairwise") 
    
    t <-signif(kensb$estimate, 3)
    stars <- stars.pval(kensb$p.value)
    mtext(text = substitute(paste("     U : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
    
  }
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}



#### plot SES WEST   : compare all (canop + ssboi) et canop (seul) ####
dim(SES_west_sbois)
dim(SES_west_canop)

names(SES_west_sbois)
names(SES_west_canop)



SES_west = compare_west$SES_stat_obs[,apply(is.na(compare_west$SES_stat_obs),2,sum) == 0]
SES_west_canop = compare_west_canop$SES_stat_obs[,apply(!is.na(compare_west_canop$SES_stat_obs),2,sum) > 1]
SES_west_sbois = compare_west_sbois$SES_stat_obs[,apply(!is.na(compare_west_sbois$SES_stat_obs),2,sum) > 1]

# SES_west_sbois = compare_west_sbois$SES_stat_obs

twi = log(env_var_west$TWI)
#     twi = env_var$TWI
par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(compare_west$SES_stat_obs)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(  compare_west$SES_stat_obs[,i]), type = "n", main = names(SES_west)[i])
  #points(SES_west[,i]~ twi)
  
  stat_tmp = names(SES_west)[i]
  
  points(SES_west[,stat_tmp]~twi)
  
  
  kenc = cor.test( compare_west$SES_stat_obs[,stat_tmp],twi, method="kendall", use="pairwise") 
  
  t <- signif(kenc$estimate, 3)
  stars <- stars.pval(kenc$p.value)
  mtext(text = substitute(paste("    A : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  if(stat_tmp %in% names(SES_west_canop)){
    points(SES_west_canop[,stat_tmp]~twi, pch = 16)
    
    
    kensb =cor.test(SES_west_canop[,stat_tmp],twi, method="kendall", use="pairwise") 
    
    t <-signif(kensb$estimate, 3)
    stars <- stars.pval(kensb$p.value)
    mtext(text = substitute(paste("     C : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
    
  }
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}










########################################################################################################################
##### est #####
########################################################################################################################

community_east = data.frame(community_east, check.names = F)

compeast = list()

# with abundances

compeast$spe <- data.frame(community_east, check.names = F)

#   compeast$spe = data.frame(ifelse(community_east == 0, 0, 1)) # transform direct la matrix abundance en pers / abs


compeast$traits = TF_east

# models du 28/12/2017 et 03012018 avec ( "TWI","INS", "CURV")

env_var_select = env_var_east[,c("TWI","INS", "CURV")]

# env_var_select = env_var[,c("INS","PREC","TWI")]

compeast$env = env_var_select


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.east <- dudi.coa(compeast$spe, scannf = FALSE)

acpR.east <- dudi.hillsmith(compeast$env, row.w = afcL.east$lw,
                            scannf = FALSE)


acpQ.east <- dudi.pca(compeast$traits, row.w = afcL.east$cw,
                      scannf = FALSE)

#### RLQ ####

rlq.east <- rlq(acpR.east, afcL.east, acpQ.east,
                scannf = FALSE)


100*rlq.east$eig/sum(rlq.east$eig)

summary(rlq.east)

plot(rlq.east)
plot(rlq.east, xax = 1, yax = 2)



##### FOURTH CORNER WEST : la relation entre traits en environement est-elle significative ??  #####



nrepet <- 9999

four.comb.east <- fourthcorner(compeast$env, compeast$spe,
                               compeast$traits, modeltype = 6, p.adjust.method.G = "none", 
                               p.adjust.method.D = "none", nrepet = nrepet)

four.comb.east
four.comb.east.adj <- p.adjust.4thcorner(four.comb.east,
                                         p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")


plot(four.comb.east.adj, alpha = 0.05, stat="D2")


plot(four.comb.east, alpha = 0.05, stat="D2")


plot(four.comb.east, x.rlq = rlq.east, alpha = 0.05,
     stat = "D2", type = "biplot")





testQaxes.comb.east <- fourthcorner.rlq(rlq.east, modeltype = 6,
                                        typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                        p.adjust.method.D = "none")

print(testQaxes.comb.east, stat = "D")

testRaxes.comb.east <- fourthcorner.rlq(rlq.east, modeltype = 6,
                                        typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "none",
                                        p.adjust.method.D = "none")
print(testRaxes.comb.east, stat = "D")


par(mfrow = c(1, 2))
plot(testQaxes.comb.east, alpha = 0.05, type = "table",
     stat = "D2")
plot(testRaxes.comb.east, alpha = 0.05, type = "table",
     stat = "D2")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(testQaxes.comb.east, alpha = 0.05, type = "biplot",
     stat = "D2")
plot(testRaxes.comb.east, alpha = 0.05, type = "biplot",
     stat = "D2")
par(mfrow = c(1, 1))









######   stat obs all #####
stat_obs_east = tf_stats_from_rlq4(comp = compeast, 
                                   class = env_var$CLASS[21:40])

# saveRDS(stat_obs_east, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_east03012018.rds")

##### model null all #####

md_null_east = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_east, iter = 999, class = env_var$CLASS[21:40], 
                                                tab_stats_obs = stat_obs_east$tab_stats_obs, compare = F)



compare_east = compare_stats_obs_null(stat_obs_east,md_null_east)



# saveRDS(md_null_east, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east03012018.rds")



#   md_null_east = readRDS( file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east28122017.rds")


##### plot east stat obs & SES #####



twi = log(env_var_east$TWI)
#     twi = env_var$TWI
#### plot Stats obs WEST  ####
par(mfrow=c(6,6),mar= c(2,2,2,2))
for(i in 1:ncol(stat_obs_east$tab_stats_obs)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(stat_obs_east$tab_stats_obs[,i]~ twi, xlim = range(twi), ylim = range(stat_obs_east$tab_stats_obs[,i]), main = names(stat_obs_east$tab_stats_obs)[i])
  
  
  kenw = cor.test(stat_obs_east$tab_stats_obs[,i],twi, method="kendall", use="pairwise") 
  
  t <-signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  
  
  
}

SES_nona = compare_east$SES_stat_obs

SES_nona = SES_nona[,apply(is.na(SES_nona),2,sum) == 0]

par(mfrow=c(6,6),mar= c(2,2,2,2))
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(SES_nona[,i]~ twi, xlim = range(twi), ylim = range(SES_nona[,i]), main = names(SES_nona)[i])
  
  
  kenw = cor.test(SES_nona[,i],twi, method="kendall", use="pairwise") 
  
  t <-signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste("     E : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
  
}





######   community canop & sous-bois #####

community_east_canop = community_east[]
community_east_canop[] = 0
community_east_sbois = community_east_canop
for(i in 1:length(unique(FBA_nona$id_pts)[W])){
  
  east_tmp_canop = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate > 2  , ]
  
  sp_canop_com_tmp = unique(east_tmp_canop$Taxon)
  
  community_east_canop[i, colnames(community_east) %in% sp_canop_com_tmp] = community_east[i, colnames(community_east) %in% sp_canop_com_tmp]
  
  east_tmp_sbois = FBA_nona[FBA_nona$id_pts == unique(FBA_nona$id_pts)[i] & FBA_nona$strate < 3  , ]
  
  sp_sbois_com_tmp = unique(east_tmp_sbois $Taxon)
  
  community_east_sbois[i, colnames(community_east) %in% sp_sbois_com_tmp] = community_east[i, colnames(community_east) %in% sp_sbois_com_tmp]
  
  
}


nb_ind_sbois_east =  apply(community_east_sbois,1,sum)



######   stat obs east canop #####
stat_obs_east_canop = tf_stats_from_rlq3(orlq = rlq.east , com = community_east_canop, 
                                         class = env_var$CLASS)



##### model null east canop #####

md_null_east_canop = tf_stats_from_rlq_permut_species(orlq = rlq.east, com = community_east_canop, iter = 999, class = env_var_east$CLASS, 
                                                      tab_stats_obs = stat_obs_east_canop$tab_stats_obs, compare = F)
##### compare east canop #####

compare_east_canop  = compare_stats_obs_null(stat_obs_east_canop,md_null_east_canop)


######   stat obs east sous bois #####
stat_obs_east_sbois = tf_stats_from_rlq3(orlq = rlq.east , com = community_east_sbois, 
                                         class = env_var_east$CLASS)


##### model null east sous bois #####


md_null_east_sbois = tf_stats_from_rlq_permut_species(orlq = rlq.east, com = community_east_sbois, iter = 9, class = env_var_east$CLASS, 
                                                      tab_stats_obs = stat_obs_east_sbois$tab_stats_obs, compare = F)

##### compare east sous bois #####


compare_east_sbois  = compare_stats_obs_null(stat_obs_east_sbois,md_null_east_sbois)




#### plot SES east   : compare ssbois et canop ####
dim(SES_east_sbois)
dim(SES_east_canop)

names(SES_east_sbois)
names(SES_east_canop)



SES_east = compare_east$SES_stat_obs[,apply(is.na(compare_east$SES_stat_obs),2,sum) == 0]
SES_east_canop = compare_east_canop$SES_stat_obs[,apply(!is.na(compare_east_canop$SES_stat_obs),2,sum) > 1]
SES_east_sbois = compare_east_sbois$SES_stat_obs[,apply(!is.na(compare_east_sbois$SES_stat_obs),2,sum) > 1]

# SES_east_sbois = compare_east_sbois$SES_stat_obs

twi = log(env_var_east$TWI)
#     twi = env_var$TWI
par(mfrow=c(6,6),mar= c(2,2,2,2))
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(SES_east[,i]), type = "n", main = names(SES_east)[i])
  #points(SES_east[,i]~ twi)
  
  stat_tmp = names(SES_east)[i]
  
  points(SES_east_canop[,stat_tmp]~twi)
  
  
  kenc = cor.test(SES_east_canop[,stat_tmp],twi, method="kendall", use="pairwise") 
  
  t <- signif(kenc$estimate, 3)
  stars <- stars.pval(kenc$p.value)
  mtext(text = substitute(paste("    C : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  if(stat_tmp %in% names(SES_east_sbois)){
    points(SES_east_sbois[,stat_tmp]~twi, pch = 16)
    
    
    kensb =cor.test(SES_east_sbois[,stat_tmp],twi, method="kendall", use="pairwise") 
    
    t <-signif(kensb$estimate, 3)
    stars <- stars.pval(kensb$p.value)
    mtext(text = substitute(paste("     U : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
    
  }
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}



#### plot SES east   : compare all (canop + ssboi) et canop (seul) ####
dim(SES_east_sbois)
dim(SES_east_canop)

names(SES_east_sbois)
names(SES_east_canop)



SES_east = compare_east$SES_stat_obs[,apply(is.na(compare_east$SES_stat_obs),2,sum) == 0]
SES_east_canop = compare_east_canop$SES_stat_obs[,apply(!is.na(compare_east_canop$SES_stat_obs),2,sum) > 1]
SES_east_sbois = compare_east_sbois$SES_stat_obs[,apply(!is.na(compare_east_sbois$SES_stat_obs),2,sum) > 1]

# SES_east_sbois = compare_east_sbois$SES_stat_obs

twi = log(env_var_east$TWI)
#     twi = env_var$TWI
par(mfrow=c(5,6),mar= c(2,2,2,2))
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  
  
  plot(0, 0, xlim = range(twi), ylim = range(SES_east[,i]), type = "n", main = names(SES_east)[i])
  #points(SES_east[,i]~ twi)
  
  stat_tmp = names(SES_east)[i]
  
  points(SES_east[,stat_tmp]~twi)
  
  
  kenc = cor.test(SES_east[,stat_tmp],twi, method="kendall", use="pairwise") 
  
  t <- signif(kenc$estimate, 3)
  stars <- stars.pval(kenc$p.value)
  mtext(text = substitute(paste("    A : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  if(stat_tmp %in% names(SES_east_canop)){
    points(SES_east_canop[,stat_tmp]~twi, pch = 16)
    
    
    kensb =cor.test(SES_east_canop[,stat_tmp],twi, method="kendall", use="pairwise") 
    
    t <-signif(kensb$estimate, 3)
    stars <- stars.pval(kensb$p.value)
    mtext(text = substitute(paste("     C : ", tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=1)
    
  }
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
  
}





##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

##### model null all de base (script épuré) #####

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 




community_nona = data.frame(community_all, check.names = F)

compall = list()

# with abundances

compall$spe <- data.frame(community_nona, check.names = F)



compall$traits = TF_kone_com_nona


env_var_select = env_var[,c( "TWI","INS", "CURV", "SIDE")]


compall$env = env_var_select




stat_obs_all = tf_stats_from_rlq4(comp = compall ,class = env_var$CLASS, class2 = env_var$SIDE)

# saveRDS(stat_obs_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all_15012018.rds")

md_null_all = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, class = env_var$CLASS, 
                                                class2 = env_var$SIDE, tab_stats_obs = stat_obs_all$tab_stats_obs, compare = F)


# saveRDS(md_null_all, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all_15012018.rds")


# est
community_west = data.frame(community_west, check.names = F)

compwest = list()


compwest$spe <- data.frame(community_west, check.names = F)


compwest$traits = TF_west


env_var_select = env_var_west[,c("TWI","INS", "CURV")]

compwest$env = env_var_select



######   stat obs all #####
stat_obs_west = tf_stats_from_rlq4(comp = compwest, 
                                   class = env_var$CLASS[21:40])

# saveRDS(stat_obs_west, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_west15012018.rds")

##### model null all #####

md_null_west = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_west, iter = 999, class = env_var$CLASS[21:40], 
                                                 tab_stats_obs = stat_obs_west$tab_stats_obs, compare = F)

# saveRDS(md_null_west, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west15012018.rds")


compare_west = compare_stats_obs_null(stat_obs_west,md_null_west)





# est
community_east = data.frame(community_east, check.names = F)

compeast = list()


compeast$spe <- data.frame(community_east, check.names = F)


compeast$traits = TF_east


env_var_select = env_var_east[,c("TWI","INS", "CURV")]

compeast$env = env_var_select



######   stat obs all #####
stat_obs_east = tf_stats_from_rlq4(comp = compeast, 
                                   class = env_var$CLASS[21:40])

#  saveRDS(stat_obs_east, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_east15012018.rds")

##### model null all #####

md_null_east = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_east, iter = 999, class = env_var$CLASS[21:40], 
                                                 tab_stats_obs = stat_obs_east$tab_stats_obs, compare = F)

 # saveRDS(md_null_east, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east15012018.rds")
 

compare_east = compare_stats_obs_null(stat_obs_east,md_null_east)






##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

##### model null all west & east from rlq all with LA & SLA log transformed #####

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 






compall = list()

# with abundances

compall$spe <- data.frame(community_all, check.names = F)


compall$traits = TF_kone_com_nona

# log transform LA & SLA

compall$traits$SLA = log(compall$traits$SLA)

compall$traits$LA =log(compall$traits$LA)





env_var_select = env_var[,c( "TWI","INS", "CURV", "SIDE")]
# env_var_select = env_var[,c( "TWI","INS", "CURV", "SIDE", "CLASS")]

#  env_var_select = env_var[,c( "TWI", "SIDE")]

# env_var_select = env_var[,c("SIDE", "CLASS")]


compall$env = env_var_select


stat_obs_all2 = tf_stats_from_rlq4(comp = compall ,class = env_var$CLASS, class2 = env_var$SIDE)


# saveRDS(stat_obs_all2, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/stat_obs_all_15012018logsla&la.rds")


md_null_all2 = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, class = env_var$CLASS, 
                                                class2 = env_var$SIDE, tab_stats_obs = stat_obs_all$tab_stats_obs, compare = F)


#  saveRDS(md_null_all2, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_all2_15012018logsla&la.rds")

##### model null all from rlq all with LA & SLA log transformed with swap of species WITH RESPECT TO TRAITS #####



null_model_all3 = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, method = "spswap", recompute.rlq.sp = F, 
                                                    class = env_var_west$CLASS, 
                                                    compare = F, class2 = env_var$SIDE)

 #  saveRDS(null_model_all3, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/null_model_all3_15012018logsla&la.rds")


##### null model east & west from rlq all #####

md_null_west2 = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, class = env_var_west$CLASS, 
                                                  compare = F, sitlist = 1:20)



md_null_east2 = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 999, class = env_var_east$CLASS, 
                                                  compare = F, sitlist = 21:40)


debug(tf_stats_from_rlq_permut_species2)

# saveRDS(md_null_west2, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_west2_15012018logsla&la.rds")
 # saveRDS(md_null_east2, file = "/home/thesardfou/Documents/data/Model_null_RLQ_stats_com TF/compare_est_ouest/md_null_east2_15012018logsla&la.rds")


 

 
 
compare_all2 = compare_stats_obs_null(stat_obs_all,md_null_all2)

compare_all3 = compare_stats_obs_null(stat_obs_all,null_model_all3)

compare_west2 = compare_stats_obs_null(stat_obs_all,md_null_west2, sitlist = 1:20)


compare_east2 = compare_stats_obs_null(stat_obs_all,md_null_east2, sitlist = 21:40)



md_null_east2b = tf_stats_from_rlq_permut_species2(stat_obs = stat_obs_all, iter = 3, class = env_var_west$CLASS, 
                                                   compare = F, sitlist = 21:40)



compare_east2b = compare_stats_obs_null(stat_obs_all,md_null_east2b, sitlist = 21:40)

compare_east2b$SES_stat_obs

plot(compare_west2$SES_stat_obs)


##### plot #####
#### all #####


layout(matrix(1:6, 2, 3, byrow = F),widths = rep(c(1),3) , heights= rep(1,2))

#layout.show(2)
vecvar_w = 1:ncol(stat_obs_all$tab_stats_obs)
vars_title_w = colnames(stat_obs_all$tab_stats_obs)
twi = env_var$TWI

for(i in vecvar_w){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_all$tab_stats_obs[,i]), main = vars_title_w[which(vecvar_w == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_w == i)], line = 1, cex = 2)
  points(stat_obs_all$tab_stats_obs[,i]~twi[], pch = 16)
  
  kenw = cor.test(stat_obs_all$tab_stats_obs[,i],twi[], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}




for(i in vecvar_w){
  

kt = kruskal.test( stat_obs_all$tab_stats_obs[,i] ~ as.numeric(as.factor(env_var$SIDE) ) )

par(mar= c(4,4,2.5,1))

plot(stat_obs_all$tab_stats_obs[,i] ~ as.factor(env_var$SIDE), axes=T, frame.plot=TRUE, xlab = "", ylab = "",
     main = vars_title_w[which(vecvar_w == i)])
axis(1, at = c(1,2), labels = c("E", "W"))
stars <- stars.pval(kt$p.value)

mtext(text = paste("p =", signif(kt$p.value,digits = 2), stars), side=3, cex=.6, adj=1)


}

################################################################################################################################################################################################################
################################################################################################################################################################################################################

###### PLOT SES from rlq all abundances log la SLA ===> West & Est from RLQ ALL, permutation intra-forest type ######

################################################################################################################################################################################################################
################################################################################################################################################################################################################



#### west #####

#X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F1w.png", width = 900, height = 500, pointsize = 14)

layout(matrix(1:6, 2, 3, byrow = F),widths = rep(c(1),3) , heights= rep(1,2))

#layout.show(2)
vecvar_w = 1:ncol(stat_obs_all$tab_stats_obs)
vars_title_w = colnames(stat_obs_all$tab_stats_obs)
twi = log(env_var$TWI)

for(i in vecvar_w){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_all$tab_stats_obs[W,i]), main = vars_title_w[which(vecvar_w == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_w == i)], line = 1, cex = 2)
  points(stat_obs_all$tab_stats_obs[W,i]~twi[W], pch = 16)
  
  kenw = cor.test(stat_obs_all$tab_stats_obs[W,i],twi[W], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}


#### est #####
#layout.show(2)
vecvar_w = 1:ncol(stat_obs_all$tab_stats_obs)
vars_title_w = colnames(stat_obs_all$tab_stats_obs)
twi = log(env_var$TWI)

for(i in vecvar_w){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi), ylim = range(stat_obs_all$tab_stats_obs[E,i]), main = vars_title_w[which(vecvar_w == i)], type = "n")
  
  # title(main = vars_title[which(vecvar_w == i)], line = 1, cex = 2)
  points(stat_obs_all$tab_stats_obs[E,i]~twi[E], pch = 16)
  
  kenw = cor.test(stat_obs_all$tab_stats_obs[E,i],twi[E], method="kendall", use="pairwise", exact = F) 
  
  t <- signif(kenw$estimate, 2)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
}



###### plot SES west #####



X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F2_w.png", width = 900, height = 500, pointsize = 14)


layout(matrix(1:32, 4, 8, byrow = T),widths = rep(c(1),8) , heights= rep(1,4))

SES =  compare_west2$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]


twi = env_var$TWI
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[W]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",names(SES_nona)[i], sep = " ") )
  points(SES_nona[,i]~twi[W], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[W], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}




###### plot SES east #####



X11( width= 8 , height= 6, type="cairo")

# png(filename="/home/thesardfou/Documents/images/R_export/figures_P1/F2_w.png", width = 900, height = 500, pointsize = 14)


layout(matrix(1:32, 4, 8, byrow = T),widths = rep(c(1),8) , heights= rep(1,4))

SES =  compare_east2$SES_stat_obs

names(SES_nona)

SES_nona = SES[,apply(is.na(SES),2,sum) == 0]

SES_nona = SES_nona[names(SES_nona) %in% vars_w]

twi = env_var$TWI
for(i in 1:ncol(SES_nona)){
  
  #  xyplot(SES_nona[,i]~env_var$TWI, groups =  env_var$SIDE)
  par(mar= c(2.5,2.5,2.5,1))
  
  
  plot(0, 0, xlim = range(twi[E]), ylim = range(SES_nona[,i]), type = "n", main = paste("SES",names(SES_nona)[i], sep = " ") )
  points(SES_nona[,i]~twi[E], pch = 16)
  
  
  
  kenw = cor.test(SES_nona[,i],twi[E], method="kendall", use="pairwise") 
  
  t <- signif(kenw$estimate, 3)
  stars <- stars.pval(kenw$p.value)
  mtext(text = substitute(paste( tau, " = ", x, y, sep=""), list(x=t, y=stars)), side=3, cex=.6, adj=0)
  
  
  abline(h=0)
  
  abline(h=1.96,lty = 3)
  abline(h=-1.96,lty = 3)
  
}














##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

##### PARTIAL RLQ all with LA & SLA log transformed #####

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 






# data 


compall = readRDS( file = "/home/thesardfou/Documents/thèse/papiers/structure_fonctionnelle_VS_précipitatation&topographie/data/data_for_RLQ.rds")


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

# saveRDS(compall, file = "/home/thesardfou/Documents/thèse/papiers/structure_fonctionnelle_VS_précipitatation&topographie/data/data_for_RLQ.rds")


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


######################################################################################################

# test 4corner with partial rlq 

######################################################################################################


lm(compall$env ~ compall$covar)









