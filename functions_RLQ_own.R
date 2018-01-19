
################################################################################################################################################
################################################################################################################################################

####### compute multiple FT stats for communities from RLQ ==> RLQ from abundance matrix & rlq from pres/abs matrix ######
####### null model : compute multiple FT stats for communities from RLQ with iterations on species permutation ######

###### avec script R RLQ_EST_OUEST_NEW

################################################################################################################################################
################################################################################################################################################

tf_stats_from_rlq <- function(TF, env, com, class = NULL, fba = NULL, iter = NULL, null_model = F){
  
  # compute multiple FT stats for communities from RLQ ==> RLQ from abundance matrix & rlq from pres/abs matrix 
  # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
  

  library(ade4)
  library(FD)
  library(moments)
  library(utils)
  
  

# obs datas

compFBA = list()

# with abundances

compFBA$spe <- data.frame(com, check.names = F)

compFBA$traits = TF

env_var_select = env


compFBA$env = env_var_select


######  perform the separate analyses of each table = CoA & PCA ##### 


afcL.FBA <- dudi.coa(compFBA$spe, scannf = FALSE)

acpR.FBA <- dudi.hillsmith(compFBA$env, row.w = afcL.FBA$lw,
                           scannf = FALSE)


acpQ.FBA <- dudi.pca(compFBA$traits, row.w = afcL.FBA$cw,
                     scannf = FALSE)

#### RLQ ####

rlq.FBA <- rlq(acpR.FBA, afcL.FBA, acpQ.FBA,
               scannf = FALSE)


rlq.FBA$lQ


sites_nona = unique(FBA_nona$id_pts)
list_scores_abund_com = list()
Ranges_Ax1 = c()
Ranges_Ax2 = c()
CWV_Ax1 = c()
CWV_Ax2 = c()
FDis_sp = c()
FDis_abund = c()
skewness_Ax1_abund = c()
skewness_Ax2_abund = c()
kurtosis_Ax1_abund = c()
kurtosis_Ax2_abund = c()
df_FD = data.frame()

for (i in 1:nrow(compFBA$spe)){
  
  # species de la com avec les abondances
  
  sp_com_abund = rep(names(compFBA$spe[i,]), compFBA$spe[i,]) 
  
  # spécies (occurences) de la com
  
  com_tmp = rlq.FBA$lQ[rownames(rlq.FBA$lQ) %in%  sp_com_abund,]
  
  com_abund_tmp = rlq.FBA$lQ[as.character(sp_com_abund),]
  
  list_scores_abund_com[[i]] = com_abund_tmp
  
  Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
  Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
  
  CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
  CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
  
  FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
  
  skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
  skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
  
  kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
  kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
  
  
  invisible(capture.output(FD <- data.frame(dbFD(com_abund_tmp, stand.x =T )))) # from Laliberté et al. 2010
  
  df_FD = rbind(df_FD, FD)
  
}


names(df_FD)[1] = "nbind"


##### sp canop & sous bois #####

  #### with occurences ####
  
  
  compFBA_sp = compFBA
  
  compFBA_sp$spe = data.frame(ifelse(compFBA$spe == 0, 0, 1), check.names = F) # transform direct la matrix abundance en pers / abs
  
  
  
  ######  perform the separate analyses of each table = CoA & PCA ##### 
  
  
  afcL.FBA_sp <- dudi.coa(compFBA_sp$spe, scannf = FALSE)
  
  acpR.FBA_sp <- dudi.hillsmith(compFBA_sp$env, row.w = afcL.FBA_sp$lw,
                                scannf = FALSE)
  
  
  acpQ.FBA_sp <- dudi.pca(compFBA_sp$traits, row.w = afcL.FBA_sp$cw,
                          scannf = FALSE)
  
  #### RLQ ####
  
  rlq.FBA_sp <- rlq(acpR.FBA_sp, afcL.FBA_sp, acpQ.FBA_sp,
                    scannf = FALSE)
  
  
  list_scores_sp_com = list()
  CV_Ax1 = c()
  CV_Ax2 = c()
  FDis_sp = c()
  skewness_Ax1_sp = c()
  skewness_Ax2_sp = c()
  kurtosis_Ax1_sp = c()
  kurtosis_Ax2_sp = c()
  df_FD_sp = data.frame()
  for (i in 1:nrow(compFBA_sp$spe)){
    
    # species de la com 
    
    sp_com = rep(names(compFBA_sp$spe[i,]), compFBA_sp$spe[i,]) 
    
    
    com_tmp = rlq.FBA_sp$lQ[sp_com,]
    
    
    scores_sp_com = list()
    
    
    list_scores_sp_com[[i]] = com_tmp
    
    
    CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
    CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
    
    FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
    
    skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
    skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
    
    kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
    kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD_sp <- data.frame(dbFD(com_tmp, stand.x =T )))) # from Laliberté et al. 2010
    
    df_FD_sp = rbind(df_FD_sp, FD_sp)
  }
  
  
  names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
  
  df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
  
  # dataframe with stat per coms
  
  
  
  tab_stats_obs = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                        skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                        skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
  
  tab_stats_obs = data.frame(tab_stats_obs)
  
  
  ##### stats for com classes #####
  
  if(is.null(class) == F){
    mean_stats_classe = list()
    for (i in 1:length(levels(class))){ # i in classes
      
      coms_tmp = c(1:40)[ class == levels(class)[i] ]
      
      mean_stats_classe[[levels(class)[i]]] = apply(tab_stats_obs[coms_tmp,-ncol(tab_stats_obs)],2,mean)
      
      
    }
    
  }
  

  
  
  
  
  
  
  if(null_model == T){
    
    # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
    
    print(paste("NULL MODEL WITH ", iter, " ITERATIONS"))
    
    list_tab_stats_null = list()
    
    list_mean_stats_classe_null = list()
    
    list_perm_com = list()
    
    
    pb <- txtProgressBar(min = 0, max =  iter, style = 3)
    for (k in 1:iter){
      
      # permute sp
      sp_tmp = sample(rownames(rlq.FBA$lQ))
      
      perm_coord_tmp =  rlq.FBA$lQ
      
      rownames(perm_coord_tmp) = sp_tmp
      
      list_perm_com[[k]] = perm_coord_tmp
      
      ##### compute stats  ####
      
      
      
      sites_nona = unique(FBA_nona$id_pts)
      list_scores_sp_com = list()
      Ranges_Ax1 = c()
      Ranges_Ax2 = c()
      CWV_Ax1 = c()
      CWV_Ax2 = c()
      FDis_sp = c()
      FDis_abund = c()
      skewness_Ax1_abund = c()
      skewness_Ax2_abund = c()
      kurtosis_Ax1_abund = c()
      kurtosis_Ax2_abund = c()
      df_FD = data.frame()
      
      for (i in 1:nrow(compFBA$spe)){
        
        # species de la com avec les abondances
        
        sp_com_abund = rep(names(compFBA$spe[i,]), compFBA$spe[i,]) 
        
        # spécies (occurences) de la com
        
        com_tmp =  perm_coord_tmp[rownames( perm_coord_tmp) %in%  sp_com_abund,]
        
        com_abund_tmp =  perm_coord_tmp[as.character(sp_com_abund),]
        
        list_scores_sp_com = com_abund_tmp
        
        Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
        Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
        
        CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
        CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
        
        FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
        
        skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
        skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
        
        kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
        kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
        
        
        # from Laliberté et al. 2010
        invisible(capture.output(FD <- data.frame(dbFD(com_abund_tmp, stand.x =T , stand.FRic = T))))
        
        df_FD = rbind(df_FD, FD)
        
      }
      
      
      names(df_FD)[1] = "nbind"
      
      
      #### with occurences ####
      
      perm_coord_tmp =  rlq.FBA_sp$lQ
      
      rownames(perm_coord_tmp) = sp_tmp
      
      #### compute stats ####
      
      
      list_scores_sp_com_tmp = list()
      CV_Ax1 = c()
      CV_Ax2 = c()
      FDis_sp = c()
      skewness_Ax1_sp = c()
      skewness_Ax2_sp = c()
      kurtosis_Ax1_sp = c()
      kurtosis_Ax2_sp = c()
      df_FD_sp = data.frame()
      for (i in 1:nrow(compFBA_sp$spe)){
        
        # species de la com 
        
        sp_com = rep(names(compFBA_sp$spe[i,]), compFBA_sp$spe[i,]) 
        
        
        com_tmp = perm_coord_tmp[sp_com,]
        
        
        scores_sp_com = list()
        
        
        list_scores_sp_com_tmp[[i]] = com_tmp
        
        
        CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
        CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
        
        FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
        
        skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
        skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
        
        kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
        kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
        
        # from Laliberté et al. 2010
        invisible(capture.output(FD_sp <- data.frame(dbFD(com_tmp, stand.x =T , stand.FRic = T))))
        
        
        df_FD_sp = rbind(df_FD_sp, FD_sp)
      }
      
      
      names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
      
      paste(names(df_FD_sp),"_sp")
      
      df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
      
      # dataframe with stat per coms
      
      
      
      tab_stats_null_temp = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                                  skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                                  skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
      
      tab_stats_null_temp = data.frame(tab_stats_null_temp)
      
      
      ##### stats for com classes #####
      
      
      if(is.null(class) == F){
        mean_stats_classe_null_temp = list()
        
        for (i in 1:length(levels(class))){ # i in classes
          
          coms_tmp = c(1:40)[ class == levels(class)[i] ]
          
          mean_stats_classe_null_temp[[levels(class)[i]]] = apply(tab_stats_null_temp[coms_tmp,-ncol(tab_stats_null_temp)],2,mean)
          
          
        }
        
        
        list_mean_stats_classe_null[[k]] = mean_stats_classe_null_temp
        
      }  
      
      #### stok stats for each iteration ####
      list_tab_stats_null[[k]] =  tab_stats_null_temp
      
      
      
      
      # update progress bar
      setTxtProgressBar(pb, k)
      
    }
    close(pb)
    

    
    return(list(orlq = rlq.FBA, tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com, 
                list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                list_perm_com = list_perm_com))
  
  }else{
    
     return(list(orlq = rlq.FBA, tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com ))
  }
  


}



################################################################################################################################################
################################################################################################################################################


####### compute multiple FT stats for communities from RLQ  ==> only RLQ from abundancec matrix ######
####### null model : compute multiple FT stats for communities from RLQ with iterations on species permutation ######

################################################################################################################################################
################################################################################################################################################

tf_stats_from_rlq2 <- function(TF, env, com, class = NULL, fba = NULL, iter = NULL, null_model = F){
  
  # compute multiple FT stats for communities from RLQ  ==> only RLQ from abundancec matrix
  # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation
  
  library(ade4)
  library(FD)
  library(moments)
  library(utils)
  
  
  
  # obs datas
  
  compFBA = list()
  
  # with abundances
  
  compFBA$spe <- data.frame(com, check.names = F)
  
  compFBA$traits = TF
  
  env_var_select = env
  
  
  compFBA$env = env_var_select
  
  
  ######  perform the separate analyses of each table = CoA & PCA ##### 
  
  
  afcL.FBA <- dudi.coa(compFBA$spe, scannf = FALSE)
  
  acpR.FBA <- dudi.hillsmith(compFBA$env, row.w = afcL.FBA$lw,
                             scannf = FALSE)
  
  
  acpQ.FBA <- dudi.pca(compFBA$traits, row.w = afcL.FBA$cw,
                       scannf = FALSE)
  
  #### RLQ ####
  
  rlq.FBA <- rlq(acpR.FBA, afcL.FBA, acpQ.FBA,
                 scannf = FALSE)
  
  

  sites_nona = unique(FBA_nona$id_pts)
  list_scores_abund_com = list()
  Ranges_Ax1 = c()
  Ranges_Ax2 = c()
  CWV_Ax1 = c()
  CWV_Ax2 = c()
  FDis_sp = c()
  FDis_abund = c()
  skewness_Ax1_abund = c()
  skewness_Ax2_abund = c()
  kurtosis_Ax1_abund = c()
  kurtosis_Ax2_abund = c()
  df_FD = data.frame()
  
  for (i in 1:nrow(compFBA$spe)){
    
    # species de la com avec les abondances
    
    sp_com_abund = rep(names(compFBA$spe[i,]), compFBA$spe[i,]) 
    
    # spécies (occurences) de la com
    
    com_tmp = rlq.FBA$lQ[rownames(rlq.FBA$lQ) %in%  sp_com_abund,]
    
    com_abund_tmp = rlq.FBA$lQ[as.character(sp_com_abund),]
    
    list_scores_abund_com[[i]] = com_abund_tmp
    
    Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
    Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
    
    CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
    CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
    
    FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
    
    skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
    skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
    
    kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
    kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
    
    
    invisible(capture.output(FD <- data.frame(dbFD(com_abund_tmp, stand.x = F)))) # from Laliberté et al. 2010
    
    df_FD = rbind(df_FD, FD)
    
  }
  
  
  names(df_FD)[1] = "nbind"
  
  
  
  #### with occurences BUT FROM RLQ COMPUTED ONLY WITH ABUNDANCES MATRIX ####
  
  
  
  #### RLQ ####
  
  rlq.FBA_sp 
  
  list_scores_sp_com = list()
  CV_Ax1 = c()
  CV_Ax2 = c()
  FDis_sp = c()
  skewness_Ax1_sp = c()
  skewness_Ax2_sp = c()
  kurtosis_Ax1_sp = c()
  kurtosis_Ax2_sp = c()
  df_FD_sp = data.frame()
  for (i in 1:nrow(compFBA_sp$spe)){
    
    # species de la com 
    
    sp_com = rep(names(compFBA_sp$spe[i,]), compFBA_sp$spe[i,]) 
    
    
    com_tmp = rlq.FBA$lQ[sp_com,]
    
    
    scores_sp_com = list()
    
    
    list_scores_sp_com[[i]] = com_tmp
    
    
    CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
    CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
    
    FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
    
    skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
    skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
    
    kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
    kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD_sp <- data.frame(dbFD(com_tmp, stand.x = F)))) # from Laliberté et al. 2010
    
    df_FD_sp = rbind(df_FD_sp, FD_sp)
  }
  
  
  names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
  
  df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
  
  # dataframe with stat per coms
  
  
  
  tab_stats_obs = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                        skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                        skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
  
  tab_stats_obs = data.frame(tab_stats_obs)
  
  
  ##### stats for com classes #####
  
  if(is.null(class) == F){
    mean_stats_classe = list()
    for (i in 1:length(levels(class))){ # i in classes
      
      coms_tmp = c(1:40)[ class == levels(class)[i] ]
      
      mean_stats_classe[[levels(class)[i]]] = apply(tab_stats_obs[coms_tmp,-ncol(tab_stats_obs)],2,mean)
      
      
    }
    
  }
  
  
  
  
  
  
  
  
  if(null_model == T){
    
    # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
    
    print(paste("NULL MODEL WITH ", iter, " ITERATIONS"))
    
    list_tab_stats_null = list()
    
    list_mean_stats_classe_null = list()
    
    list_perm_com = list()
    
    
    pb <- txtProgressBar(min = 0, max =  iter, style = 3)
    for (k in 1:iter){
      
      # permute sp
      sp_tmp = sample(rownames(rlq.FBA$lQ))
      
      perm_coord_tmp =  rlq.FBA$lQ
      
      rownames(perm_coord_tmp) = sp_tmp
      
      list_perm_com[[k]] = perm_coord_tmp
      
      ##### compute stats  ####
      
      
      
      sites_nona = unique(FBA_nona$id_pts)
      list_scores_sp_com = list()
      Ranges_Ax1 = c()
      Ranges_Ax2 = c()
      CWV_Ax1 = c()
      CWV_Ax2 = c()
      FDis_sp = c()
      FDis_abund = c()
      skewness_Ax1_abund = c()
      skewness_Ax2_abund = c()
      kurtosis_Ax1_abund = c()
      kurtosis_Ax2_abund = c()
      df_FD = data.frame()
      
      for (i in 1:nrow(compFBA$spe)){
        
        # species de la com avec les abondances
        
        sp_com_abund = rep(names(compFBA$spe[i,]), compFBA$spe[i,]) 
        
        # spécies (occurences) de la com
        
        com_tmp =  perm_coord_tmp[rownames( perm_coord_tmp) %in%  sp_com_abund,]
        
        com_abund_tmp =  perm_coord_tmp[as.character(sp_com_abund),]
        
        list_scores_sp_com = com_abund_tmp
        
        Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
        Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
        
        CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
        CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
        
        FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
        
        skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
        skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
        
        kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
        kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
        
        
        # from Laliberté et al. 2010
        invisible(capture.output(FD <- data.frame(dbFD(com_abund_tmp, stand.x = F))))
        
        df_FD = rbind(df_FD, FD)
        
      }
      
      
      names(df_FD)[1] = "nbind"
      
      
      #### with occurences BUT POSITIONS ON RLQ COMPUTED FROM ABUNDANCE MATRIX ONLY  ####
      
   
      #### compute stats ####
      
      
      list_scores_sp_com_tmp = list()
      CV_Ax1 = c()
      CV_Ax2 = c()
      FDis_sp = c()
      skewness_Ax1_sp = c()
      skewness_Ax2_sp = c()
      kurtosis_Ax1_sp = c()
      kurtosis_Ax2_sp = c()
      df_FD_sp = data.frame()
      for (i in 1:nrow(compFBA_sp$spe)){
        
        # species de la com 
        
        sp_com = rep(names(compFBA_sp$spe[i,]), compFBA_sp$spe[i,]) 
        
        
        com_tmp = perm_coord_tmp[sp_com,]
        
        
        scores_sp_com = list()
        
        
        list_scores_sp_com_tmp[[i]] = com_tmp
        
        
        CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
        CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
        
        FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
        
        skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
        skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
        
        kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
        kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
        
        # from Laliberté et al. 2010
        invisible(capture.output(FD_sp <- data.frame(dbFD(com_tmp, stand.x = F))))
        
        
        df_FD_sp = rbind(df_FD_sp, FD_sp)
      }
      
      
      names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
      
      paste(names(df_FD_sp),"_sp")
      
      df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
      
      # dataframe with stat per coms
      
      
      
      tab_stats_null_temp = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                                  skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                                  skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
      
      tab_stats_null_temp = data.frame(tab_stats_null_temp)
      
      
      ##### stats for com classes #####
      
      
      if(is.null(class) == F){
        mean_stats_classe_null_temp = list()
        
        for (i in 1:length(levels(class))){ # i in classes
          
          coms_tmp = c(1:40)[ class == levels(class)[i] ]
          
          mean_stats_classe_null_temp[[levels(class)[i]]] = apply(tab_stats_null_temp[coms_tmp,-ncol(tab_stats_null_temp)],2,mean)
          
          
        }
        
        
        list_mean_stats_classe_null[[k]] = mean_stats_classe_null_temp
        
      }  
      
      #### stok stats for each iteration ####
      list_tab_stats_null[[k]] =  tab_stats_null_temp
      
      
      
      
      # update progress bar
      setTxtProgressBar(pb, k)
      
    }
    close(pb)
    
    
    
    return(list(orlq = rlq.FBA, tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com, 
                list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                list_perm_com = list_perm_com))
    
  }else{
    
    return(list(orlq = rlq.FBA, tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com ))
  }
  
  
  
}



################################################################################################################################################
################################################################################################################################################


####### compute multiple FT stats for communities from RLQ object already existing ==> only RLQ from abundancec matrix ######
###### no error if com is empty #####

################################################################################################################################################
################################################################################################################################################




tf_stats_from_rlq3 <- function(orlq, com, class = NULL, class2 = NULL, fba = NULL, iter = NULL, null_model = F){
  
  
  # compute multiple FT stats for communities from RLQ object already existing ==> only RLQ from abundancec matrix
  library(FD)
  library(moments)
  library(utils)
  
  rbind.all.columns <- function(x, y) {
    
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    
    x[, c(as.character(y.diff))] <- NA
    
    y[, c(as.character(x.diff))] <- NA
    
    return(rbind(x, y))
  }
  
  com = data.frame(com, check.names = F)
  
  mean_stats_classe = NULL
  mean_stats_classe2 = NULL

  #### RLQ ####
  
  rlq.FBA <- orlq
  
  
  
  sites_nona = unique(fba$id_pts)
  list_scores_abund_com = list()
  Ranges_Ax1 = c()
  Ranges_Ax2 = c()
  CWV_Ax1 = c()
  CWV_Ax2 = c()
  FDis_sp = c()
  FDis_abund = c()
  skewness_Ax1_abund = c()
  skewness_Ax2_abund = c()
  kurtosis_Ax1_abund = c()
  kurtosis_Ax2_abund = c()
  df_FD = data.frame()
  
  for (i in 1:nrow(com)){
    
    # species de la com avec les abondances
    
    sp_com_abund = rep(names(com[i,]), com[i,]) 
    

    
    # spécies (occurences) de la com
    
    com_tmp = rlq.FBA$lQ[rownames(rlq.FBA$lQ) %in%  sp_com_abund,]
    
    com_abund_tmp = rlq.FBA$lQ[as.character(sp_com_abund),]
    
    list_scores_abund_com[[i]] = com_abund_tmp
    
    Ranges_Ax1 = c(Ranges_Ax1, ( max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) / max(rlq.FBA$lQ$AxcQ1) - min(rlq.FBA$lQ$AxcQ1) )
    Ranges_Ax2 = c(Ranges_Ax2, ( max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) ) / max(rlq.FBA$lQ$AxcQ2) - min(rlq.FBA$lQ$AxcQ2) )
    
    CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
    CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
    
    FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
    
    skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
    skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
    
    kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
    kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD <- tryCatch({
      
      data.frame(dbFD(com_abund_tmp, stand.x = F))
    }, error=function(e){ a = t(data.frame(rep("NA",10)))
    colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
    rownames(a) =  "Community1"
    return(a)
    }
    )
    )) # from Laliberté et al. 2010
    
    
    if(dim(df_FD)[1] == 0){
      
      df_FD = FD
    }else{
      df_FD = rbind.all.columns(df_FD, FD)
    }
    
  }
  
  
  names(df_FD)[1] = "nbind"
  
  
  
  #### with occurences BUT FROM RLQ COMPUTED ONLY WITH ABUNDANCES MATRIX ####
  
  
  
  #### RLQ ####
  

  
  list_scores_sp_com = list()
  CV_Ax1 = c()
  CV_Ax2 = c()
  FDis_sp = c()
  skewness_Ax1_sp = c()
  skewness_Ax2_sp = c()
  kurtosis_Ax1_sp = c()
  kurtosis_Ax2_sp = c()
  df_FD_sp = data.frame()
  for (i in 1:nrow(com)){
    
    # species de la com 
    
    sp_com = names(com)[ com[i,] > 0] 
    
    
    com_tmp = rlq.FBA$lQ[sp_com,]
    
    
    scores_sp_com = list()
    
    
    list_scores_sp_com[[i]] = com_tmp
    
    
    CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
    CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
    
    FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
    
    skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
    skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
    
    kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
    kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
    
    
   
    invisible(capture.output(FD_sp <- tryCatch({
                               
                               data.frame(dbFD(com_tmp, stand.x = F))
                          }, error=function(e){ a = t(data.frame(rep("NA",10)))
                          colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
                          rownames(a) =  "Community1"
                          return(a)
                          }
                          )
      )) # from Laliberté et al. 2010
    
    if(dim(df_FD_sp)[1] == 0){
      
      df_FD_sp = FD_sp
    }else{
      df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
    }
    
    
  }
  
  
  names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
  
  df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
  
  # dataframe with stat per coms
  
  
  
  tab_stats_obs = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                        skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                        skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
  
  tab_stats_obs = data.frame(tab_stats_obs)
  
  tab_stats_obs[tab_stats_obs == -Inf] <- NA
 
  tab_stats_obs[tab_stats_obs=="<NA>"] = NA 
  tab_stats_obs = as.data.frame(lapply(tab_stats_obs, as.numeric))
  ##### stats for com classes #####
  
  if(is.null(class) == F){
    mean_stats_classe = list()
    for (i in 1:length(levels(class))){ # i in classes
      
      coms_tmp = c(1:length(class))[ class == levels(class)[i] ]
      
      mean_stats_classe[[levels(class)[i]]] = apply(tab_stats_obs[coms_tmp,],2,mean, na.rm = T)
      
      
    }
    mean_stats_classe = data.frame(t(data.frame(mean_stats_classe)))
  }
  
 
  
  
  if(is.null(class2) == F){
    mean_stats_classe2 = list()
    for (i in 1:length(levels(class2))){ # i in classes
      
      coms_tmp = c(1:length(class2))[ class2 == levels(class2)[i] ]
      
      mean_stats_classe2[[levels(class2)[i]]] = apply(tab_stats_obs[coms_tmp,],2,mean, na.rm = T)
      
      
    }
    mean_stats_classe2 = data.frame(t(data.frame(mean_stats_classe2)))
    
  }
  
  
  
  
  
  
  if(null_model == T){
    
    # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
    
    print(paste("NULL MODEL WITH ", iter, " ITERATIONS"))
    
    list_tab_stats_null = list()
    
    list_mean_stats_classe_null = list()
    
    list_perm_com = list()
    
    
    pb <- txtProgressBar(min = 0, max =  iter, style = 3)
    for (k in 1:iter){
      
      # permute sp
      sp_tmp = sample(rownames(rlq.FBA$lQ))
      
      perm_coord_tmp =  rlq.FBA$lQ
      
      rownames(perm_coord_tmp) = sp_tmp
      
      list_perm_com[[k]] = perm_coord_tmp
      
      ##### compute stats  ####
      
      
      
      sites_nona = unique(FBA_nona$id_pts)
      list_scores_sp_com = list()
      Ranges_Ax1 = c()
      Ranges_Ax2 = c()
      CWV_Ax1 = c()
      CWV_Ax2 = c()
      FDis_sp = c()
      FDis_abund = c()
      skewness_Ax1_abund = c()
      skewness_Ax2_abund = c()
      kurtosis_Ax1_abund = c()
      kurtosis_Ax2_abund = c()
      df_FD = data.frame()
      
      for (i in 1:nrow(compFBA$spe)){
        
        # species de la com avec les abondances
        
        sp_com_abund = rep(names(com[i,]), com[i,]) 
        
        # spécies (occurences) de la com
        
        com_tmp =  perm_coord_tmp[rownames( perm_coord_tmp) %in%  sp_com_abund,]
        
        com_abund_tmp =  perm_coord_tmp[as.character(sp_com_abund),]
        
        list_scores_sp_com = com_abund_tmp
        
        Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
        Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
        
        CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
        CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
        
        FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
        
        skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
        skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
        
        kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
        kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
        
        
        
        invisible(capture.output(FD <- tryCatch({
          
          data.frame(dbFD(com_abund_tmp, stand.x = F))
        }, error=function(e){ a = t(data.frame(rep("NA",10)))
        colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
        rownames(a) =  "Community1"
        return(a)
        }
        )
        )) # from Laliberté et al. 2010
        
        
        if(dim(df_FD)[1] == 0){
          
          df_FD = FD
        }else{
          df_FD = rbind.all.columns(df_FD, FD)
        }
        
      }
        
      
      
      
      names(df_FD)[1] = "nbind"
      
      
      #### with occurences BUT POSITIONS ON RLQ COMPUTED FROM ABUNDANCE MATRIX ONLY  ####
      
      
      #### compute stats ####
      
      
      list_scores_sp_com_tmp = list()
      CV_Ax1 = c()
      CV_Ax2 = c()
      FDis_sp = c()
      skewness_Ax1_sp = c()
      skewness_Ax2_sp = c()
      kurtosis_Ax1_sp = c()
      kurtosis_Ax2_sp = c()
      df_FD_sp = data.frame()
      for (i in 1:nrow(compFBA_sp$spe)){
        
        # species de la com 
        
        sp_com = rep(names(com[i,]), com[i,]) 
        
        
        com_tmp = perm_coord_tmp[sp_com,]
        
        
        scores_sp_com = list()
        
        
        list_scores_sp_com_tmp[[i]] = com_tmp
        
        
        CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
        CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
        
        FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
        
        skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
        skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
        
        kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
        kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
        
        # from Laliberté et al. 2010
        invisible(capture.output(FD_sp <- tryCatch({
          
          data.frame(dbFD(com_tmp, stand.x = F))
        }, error=function(e){ a = t(data.frame(rep("NA",10)))
        colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
        rownames(a) =  "Community1"
        return(a)
        }
        )
        )) # from Laliberté et al. 2010
        
        if(dim(df_FD_sp)[1] == 0){
          
          df_FD_sp = FD_sp
        }else{
          df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
        }
        
        
        
        df_FD_sp = rbind(df_FD_sp, FD_sp)
      }
      
      
      names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
      
      paste(names(df_FD_sp),"_sp")
      
      df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
      
      # dataframe with stat per coms
      
      
      
      tab_stats_null_temp = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                                  skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                                  skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
      
      tab_stats_null_temp = data.frame(tab_stats_null_temp)
      
      tab_stats_null_temp[tab_stats_null_temp == -Inf] <- NA
      tab_stats_null_temp[tab_stats_null_temp=="<NA>"] = NA 
      
      ##### stats for com classes #####
      
      
      if(is.null(class) == F){
        mean_stats_classe_null_temp = list()
        
        for (i in 1:length(levels(class))){ # i in classes
          
          coms_tmp = c(1:40)[ class == levels(class)[i] ]
          
          mean_stats_classe_null_temp[[levels(class)[i]]] = apply(tab_stats_null_temp[coms_tmp,-ncol(tab_stats_null_temp)],2,mean)
          
          
        }
        
        mean_stats_classe_null_temp = t(data.frame(mean_stats_classe_null_temp))
        list_mean_stats_classe_null[[k]] = mean_stats_classe_null_temp
        
      }  
      
      #### stok stats for each iteration ####
      list_tab_stats_null[[k]] =  tab_stats_null_temp
      
      
      
      
      # update progress bar
      setTxtProgressBar(pb, k)
      
    }
    close(pb)
    
    
    
    return(list( tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com, 
                list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                list_perm_com = list_perm_com))
    
  }else{
    
    return(list( tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, mean_stats_classe2 = mean_stats_classe2,
                list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com ))
  }
  
  
  
}


################################################################################################################################################
################################################################################################################################################


####### compute multiple FT stats for communities from RLQ object already existing ==> RLQ from abundancec matrix & RLQ from pres/abs matrix ######
###### no error if com is empty #####

################################################################################################################################################
################################################################################################################################################



tf_stats_from_rlq4 <- function(comp, class = NULL, class2 = NULL){
  
  
  # compute multiple FT stats for communities from RLQ object already existing ==> RLQ from abundancec matrix & RLQ from pres/abs matrix
  
  # need "comp" object with RLQ tabs :  $spe (community matrix), $traits (traits/sp matrix), $env (env/site matrix)  
  
  library(FD)
  library(moments)
  library(utils)
  library(ade4)
  
  rbind.all.columns <- function(x, y) {
    
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    
    x[, c(as.character(y.diff))] <- NA
    
    y[, c(as.character(x.diff))] <- NA
    
    return(rbind(x, y))
  }
  
  com = data.frame(comp$spe, check.names = F)
  
  mean_stats_classe = NULL
  mean_stats_classe2 = NULL
  
  #### RLQ ####
  afcL.abund <- dudi.coa(comp$spe, scannf = FALSE)
  
  acpR.abund <- dudi.hillsmith(comp$env, row.w = afcL.abund$lw,
                             scannf = FALSE)
  
  
  acpQ.abund <- dudi.pca(comp$traits, row.w = afcL.abund$cw,
                       scannf = FALSE)
  
  
  
  rlq.abund <- rlq(acpR.abund, afcL.abund, acpQ.abund,
                 scannf = FALSE)

  
  ###### stats ######
  
  list_scores_abund_com = list()
  Ranges_Ax1_abund = c()
  Ranges_Ax2_abund = c()
  CWV_Ax1 = c()
  CWV_Ax2 = c()
  FDis_sp = c()
  FDis_abund = c()
  skewness_Ax1_abund = c()
  skewness_Ax2_abund = c()
  kurtosis_Ax1_abund = c()
  kurtosis_Ax2_abund = c()
  Ranges_Ax1_abund_nst = c() 
  Ranges_Ax2_abund_nst = c()
  df_FD = data.frame()
  
  for (i in 1:nrow(com)){
    
    # species de la com avec les abondances
    
    sp_com_abund = rep(names(com[i,]), com[i,]) 
    
    
    
    # spécies (occurences) de la com
    
    com_tmp = rlq.abund$lQ[rownames(rlq.abund$lQ) %in%  sp_com_abund,]
    
    com_abund_tmp = rlq.abund$lQ[as.character(sp_com_abund),]
    
    list_scores_abund_com[[i]] = com_abund_tmp
    
    Ranges_Ax1_abund_nst = c(Ranges_Ax1_abund_nst,  max(com_abund_tmp$AxcQ1) - min(com_abund_tmp$AxcQ1) ) 
    Ranges_Ax2_abund_nst = c(Ranges_Ax2_abund_nst,  max(com_abund_tmp$AxcQ2) - min(com_abund_tmp$AxcQ2) ) 
    
    
    Ranges_Ax1_abund = c(Ranges_Ax1_abund, ( max(com_abund_tmp$AxcQ1) - min(com_abund_tmp$AxcQ1) ) / (max(rlq.abund$lQ$AxcQ1) - min(rlq.abund$lQ$AxcQ1)) )
    Ranges_Ax2_abund = c(Ranges_Ax2_abund, ( max(com_abund_tmp$AxcQ2) - min(com_abund_tmp$AxcQ2) ) / (max(rlq.abund$lQ$AxcQ2) - min(rlq.abund$lQ$AxcQ2)) )
    
    CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
    CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
    
    FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
    
    skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
    skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
    
    kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
    kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD <- tryCatch({
      
      data.frame(dbFD(com_abund_tmp, stand.x = F))
    }, error=function(e){ a = t(data.frame(rep("NA",10)))
    colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
    rownames(a) =  "Community1"
    return(a)
    }
    )
    )) # from Laliberté et al. 2010
    
    
    if(dim(df_FD)[1] == 0){
      
      df_FD = FD
    }else{
      df_FD = rbind.all.columns(df_FD, FD)
    }
    
  }
  
  
  names(df_FD)[1] = "nbind"
  
  
  
  #### with occurences FROM RLQ COMPUTED  WITH occurences  MATRIX ####
  
  
  
  #### RLQ  sp ####
  
  
  
  comp$spesp = data.frame(ifelse(comp$spe == 0, 0, 1), check.names = F) # transform direct la matrix abundance en pers / abs
  com = data.frame(comp$spesp, check.names = F)
  
  afcL.sp <- dudi.coa(comp$spesp, scannf = FALSE)
  
  acpR.sp <- dudi.hillsmith(comp$env, row.w = afcL.sp$lw,
                         scannf = FALSE)
  
  
  acpQ.sp <- dudi.pca(comp$traits, row.w = afcL.sp$cw,
                   scannf = FALSE)
  
  
  
  rlq.sp <- rlq(acpR.sp, afcL.sp, acpQ.sp,
              scannf = FALSE)
  
  
  
  list_scores_sp_com = list()
  CV_Ax1 = c()
  CV_Ax2 = c()
  FDis_sp = c()
  Ranges_Ax1_sp = c()
  Ranges_Ax2_sp = c()
  Ranges_Ax1_sp_nst = c() 
  Ranges_Ax2_sp_nst = c()
  skewness_Ax1_sp = c()
  skewness_Ax2_sp = c()
  kurtosis_Ax1_sp = c()
  kurtosis_Ax2_sp = c()
  df_FD_sp = data.frame()
  for (i in 1:nrow(com)){
    
    # species de la com 
    
    sp_com = names(com)[ com[i,] > 0] 
    
    
    com_tmp = rlq.sp$lQ[sp_com,]
    
    
    scores_sp_com = list()
    
    
    list_scores_sp_com[[i]] = com_tmp
    
    
    CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
    CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
    
    FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
    
    Ranges_Ax1_sp_nst = c(Ranges_Ax1_sp_nst,  max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) 
    Ranges_Ax2_sp_nst = c(Ranges_Ax2_sp_nst,  max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) )
    
    Ranges_Ax1_sp = c(Ranges_Ax1_sp, ( max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) / (max(rlq.sp$lQ$AxcQ1) - min(rlq.sp$lQ$AxcQ1) ))
    Ranges_Ax2_sp = c(Ranges_Ax2_sp, ( max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) ) / (max(rlq.sp$lQ$AxcQ2) - min(rlq.sp$lQ$AxcQ2) ))
    
    
    skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
    skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
    
    kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
    kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD_sp <- tryCatch({
      
      data.frame(dbFD(com_tmp, stand.x = F))
    }, error=function(e){ a = t(data.frame(rep("NA",10)))
    colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
    rownames(a) =  "Community1"
    return(a)
    }
    )
    )) # from Laliberté et al. 2010
    
    if(dim(df_FD_sp)[1] == 0){
      
      df_FD_sp = FD_sp
    }else{
      df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
    }
    
    
  }
  
  
  names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
  
  df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
  
  # dataframe with stat per coms
  
  
  
  tab_stats_obs = cbind(Ranges_Ax1_abund, Ranges_Ax2_abund, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                        skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                        Ranges_Ax1_sp, Ranges_Ax2_sp, skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp,
                        Ranges_Ax1_abund_nst,Ranges_Ax2_abund_nst, Ranges_Ax1_sp_nst, Ranges_Ax2_sp_nst)
  
  tab_stats_obs = data.frame(tab_stats_obs)
  
  tab_stats_obs[tab_stats_obs == -Inf] <- NA
  
  tab_stats_obs[tab_stats_obs=="<NA>"] = NA 
  tab_stats_obs = as.data.frame(lapply(tab_stats_obs, as.numeric))
  ##### stats for com classes #####
  
  if(is.null(class) == F){
    mean_stats_classe = list()
    for (i in 1:length(levels(class))){ # i in classes
      
      coms_tmp = c(1:length(class))[ class == levels(class)[i] ]
      
      mean_stats_classe[[levels(class)[i]]] = apply(tab_stats_obs[coms_tmp,],2,mean, na.rm = T)
      
      
    }
    mean_stats_classe = data.frame(t(data.frame(mean_stats_classe)))
  }
  
  
  
  
  if(is.null(class2) == F){
    mean_stats_classe2 = list()
    for (i in 1:length(levels(class2))){ # i in classes
      
      coms_tmp = c(1:length(class2))[ class2 == levels(class2)[i] ]
      
      mean_stats_classe2[[levels(class2)[i]]] = apply(tab_stats_obs[coms_tmp,],2,mean, na.rm = T)
      
      
    }
    mean_stats_classe2 = data.frame(t(data.frame(mean_stats_classe2)))
    
  }
  
  

    return(list( tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
                 mean_stats_classe2 = mean_stats_classe2,
                 list_scores_abund_com = list_scores_abund_com, 
                 list_scores_sp_com = list_scores_sp_com , rlq.abund =  rlq.abund, rlq.sp = rlq.sp,
                 species = names(com), commat = comp$spe))
  
  
}

#  a =tf_stats_from_rlq4(comp = compall ,class = env_var$CLASS, class2 = env_var$SIDE)


################################################################################################################################################
################################################################################################################################################

################## null model : compute multiple FT stats for communities from RLQ with iterations on species permutation ##################

################################################################################################################################################
################################################################################################################################################

#orlq = rlq.FBA
#com = community_nona_sbois
#iter = 3
#class = env_var$CLASS

################################################################################################################################################
###### à utiliser avec un objet issu de tf_stats_from_rlq3 ######
################################################################################################################################################


tf_stats_from_rlq_permut_species <- function(orlq, com, iter, class = NULL, class2 = NULL, fba = NULL, tab_stats_obs = NULL,
                                             mean_stats_classe = NULL, compare = F){

  # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
  
  
  list_mean_stats_classe_null = NULL
  list_mean_stats_classe2_null = NULL
  
  rbind.all.columns <- function(x, y) {
    
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    
    x[, c(as.character(y.diff))] <- NA
    
    y[, c(as.character(x.diff))] <- NA
    
    return(rbind(x, y))
  }
  
  
  com = data.frame(com, check.names = F)
  
  
  list_tab_stats_null = list()
  
  list_mean_stats_classe_null = list()
  list_mean_stats_classe2_null = list()
  
  
  list_perm_com = list()
  
  print(paste("Null Model : ",iter," iterations"))
  pb <- txtProgressBar(min = 0, max =  iter, style = 3)
  for (k in 1:iter){
    
    # permute sp
    sp_tmp = sample(rownames(orlq$lQ))
    
    perm_coord_tmp =  orlq$lQ
    
    rownames(perm_coord_tmp) = sp_tmp
    
    list_perm_com[[k]] = perm_coord_tmp
    
    ##### compute stats  ####
    
    
    
    sites_nona = unique(FBA_nona$id_pts)
    list_scores_sp_com = list()
    Ranges_Ax1 = c()
    Ranges_Ax2 = c()
    CWV_Ax1 = c()
    CWV_Ax2 = c()
    FDis_sp = c()
    FDis_abund = c()
    skewness_Ax1_abund = c()
    skewness_Ax2_abund = c()
    kurtosis_Ax1_abund = c()
    kurtosis_Ax2_abund = c()
    df_FD = data.frame()
    
    for (i in 1:nrow(com)){
      
      # species de la com avec les abondances
      
      sp_com_abund = rep(names(com[i,]), com[i,]) 
      
      # spécies (occurences) de la com
      
      com_tmp =  perm_coord_tmp[rownames( perm_coord_tmp) %in%  sp_com_abund,]
      
      com_abund_tmp =  perm_coord_tmp[as.character(sp_com_abund),]
      
      list_scores_sp_com = com_abund_tmp
      
      Ranges_Ax1 = c(Ranges_Ax1, ( max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) / max(orlq$lQ$AxcQ1) - min(orlq$lQ$AxcQ1) )
      Ranges_Ax2 = c(Ranges_Ax2, ( max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) ) / max(orlq$lQ$AxcQ2) - min(orlq$lQ$AxcQ2) )
      
      CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
      CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
      
      FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
      
      skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
      skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
      
      kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
      kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
      
      
      # from Laliberté et al. 2010
      
      invisible(capture.output(FD <- tryCatch({
        
        data.frame(dbFD(com_abund_tmp, stand.x = F))
      }, error=function(e){ a = t(data.frame(rep("NA",10)))
      colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
      rownames(a) =  "Community1"
      return(a)
      }
      )
      )) # from Laliberté et al. 2010
      
    
      
      if(dim(df_FD)[1] == 0){
        
        df_FD = FD
      }else{
        df_FD = rbind.all.columns(df_FD, FD)
      }
      
      
    }
    
    
    names(df_FD)[1] = "nbind"
    
    
    #### with occurences ####
    

    #### compute stats ####
    
    
    list_scores_sp_com_tmp = list()
    CV_Ax1 = c()
    CV_Ax2 = c()
    FDis_sp = c()
    skewness_Ax1_sp = c()
    skewness_Ax2_sp = c()
    kurtosis_Ax1_sp = c()
    kurtosis_Ax2_sp = c()
    df_FD_sp = data.frame()
    for (i in 1:nrow(com)){
      
      # species de la com 
      
      sp_com = names(com)[com[i,]>0] 
      
 
      com_tmp = perm_coord_tmp[sp_com,]
      
      
      scores_sp_com = list()
      
      
      list_scores_sp_com_tmp[[i]] = com_tmp
      
      
      CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
      CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
      
      FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
      
      skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
      skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
      
      kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
      kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
      
      # from Laliberté et al. 2010
      invisible(capture.output(FD_sp <- tryCatch({
        
        data.frame(dbFD(com_tmp, stand.x = F))
      }, error=function(e){ a = t(data.frame(rep("NA",10)))
      colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
      rownames(a) =  "Community1"
      return(a)
      }
      )
      )) # from Laliberté et al. 2010
      
      if(dim(df_FD_sp)[1] == 0){
        
        df_FD_sp = FD_sp
      }else{
        df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
      }
      
      
    }
    
    
    names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
    
    paste(names(df_FD_sp),"_sp")
    
    df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
    
    # dataframe with stat per coms
    
    
    
    tab_stats_null_temp = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                                skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                                skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
    
    tab_stats_null_temp = data.frame(tab_stats_null_temp)
    
    tab_stats_null_temp[tab_stats_null_temp == -Inf] <- NA
    
    tab_stats_null_temp[tab_stats_null_temp=="<NA>"] = NA 
    
    tab_stats_null_temp = as.data.frame(lapply(tab_stats_null_temp, as.numeric))
    ##### stats for com classes #####
    
    
    if(is.null(class) == F){
      mean_stats_classe_null_temp = list()
      
      for (i in 1:length(levels(class))){ # i in classes
        
        coms_tmp = c(1:length(class))[ class == levels(class)[i] ]
        
        mean_stats_classe_null_temp[[levels(class)[i]]] = apply(tab_stats_null_temp[coms_tmp,],2,mean, na.rm = T)
        
        
      }
      mean_stats_classe_null_temp = t(data.frame(mean_stats_classe_null_temp))
      list_mean_stats_classe_null[[k]] = mean_stats_classe_null_temp
      
    }  
    
    if(is.null(class2) == F){
      mean_stats_classe2_null_temp = list()
    
      for (i in 1:length(levels(class2))){ # i in classes
        
        coms_tmp = c(1:length(class2))[ class2 == levels(class2)[i] ]
        
        mean_stats_classe2_null_temp[[levels(class2)[i]]] = apply(tab_stats_null_temp[coms_tmp,],2,mean, na.rm = T)
        
        
      }
      mean_stats_classe2_null_temp = t(data.frame(mean_stats_classe2_null_temp))
      list_mean_stats_classe2_null[[k]] = mean_stats_classe2_null_temp
      
    }
    
    
    #### stok stats for each iteration ####
    list_tab_stats_null[[k]] =  tab_stats_null_temp
    
    
   
    
  
    
    # update progress bar
    setTxtProgressBar(pb, k)
    
  }
  close(pb)
  
  
  
  if(compare == F){
    return(list(list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                list_mean_stats_classe2_null = list_mean_stats_classe2_null, list_perm_com = list_perm_com))
    
  }
  
  if(compare == T){
    
    #### compute SES and significativity of obs stats ####
    
    if(is.null(tab_stats_obs) == T){
      print("tab_stats_obs argument empty, tab_stats_obs needed!")
    }else{
        
      
      
      
      # SES
      # turn stats null to a 3d array 
      
      list2ary = function(input.list){  #input a list of lists
        
        rows.cols <- dim(input.list[[1]])
        
        sheets <- length(input.list)
        
        output.ary <- array(as.numeric(unlist(input.list)), dim = c(rows.cols, sheets))
        
        colnames(output.ary) <- colnames(input.list[[1]])
        
        row.names(output.ary) <- row.names(input.list[[1]])
        
        return(output.ary)    # output as a 3-D array
        
      }
      
      
      
      
      array_stats_null =  list2ary(list_tab_stats_null)
      
      mean_stats_null = apply(array_stats_null,c(1,2), mean)
      
      sd_stats_null = apply(array_stats_null,c(1,2), sd)
      
      
      SES_stat_obs = ( tab_stats_obs - mean_stats_null) / sd_stats_null
      
      # significativity
      
      print(paste("compare null model with ",iter," iterations and stats obs"))
      pb <- txtProgressBar(min = 0, max =  iter, style = 3)
      
      tab_stat_obs_is_greater = tab_stats_obs
      tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
      tab_stat_obs_is_smaller = tab_stat_obs_is_greater
      tab_stat_obs_is_different = tab_stat_obs_is_greater
      for (i in 1:length( list_tab_stats_null)){
        
        null = list_tab_stats_null[[i]]
        tmp = tab_stats_obs > data.frame( null) 
        
        tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
        
        tmp2 = tab_stats_obs < data.frame( null) 
        
        tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
        
        tmp3 = ( tab_stats_obs < data.frame( null) ) | ( tab_stats_obs > data.frame( null) )
        
        tab_stat_obs_is_different = tab_stat_obs_is_different + tmp3
        
        setTxtProgressBar(pb, i)
      }
      

      
      tab_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
      tab_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
      tab_stat_obs_is_different_pvalues = 1 - tab_stat_obs_is_different/(iter+1)
      
      
      
    # compare mean stats of classes
      
      if(!is.null(class) == T){
        
        
        array_stats_null =  list2ary(list_mean_stats_classe_null)
        
        mean_stats_null = apply(array_stats_null,c(1,2), mean)
        
        sd_stats_null = apply(array_stats_null,c(1,2), sd)
        
        
        SES_mean_stat_obs = ( tab_stats_obs - mean_stats_null) / sd_stats_null
        
        # significativity
        
        print(paste("compare null model with ",iter," iterations and mean stats obs for classes"))
        pb <- txtProgressBar(min = 0, max =  iter, style = 3)
        
        tab_stat_obs_is_greater = mean_stats_classe
        tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
        tab_stat_obs_is_smaller = tab_stat_obs_is_greater
        tab_stat_obs_is_different = tab_stat_obs_is_greater
        for (i in 1:length( list_mean_stats_classe_null)){
          

          null = list_mean_stats_classe_null[[i]]
          tmp = tab_stats_obs > data.frame( null) 
          
          tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
          
          tmp2 = tab_stats_obs < data.frame( null) 
          
          tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
          
          tmp3 = ( tab_stats_obs < data.frame( null) ) | ( tab_stats_obs > data.frame( null) )
          
          tab_stat_obs_is_different = tab_stat_obs_is_different + tmp3
          
          setTxtProgressBar(pb, i)
        }
        
        
        
        tab_mean_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
        tab_mean_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
        tab_mean_stat_obs_is_different_pvalues = 1 - tab_stat_obs_is_different/(iter+1)
        
      }
      
     
      
      
      
      return(list(list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                  list_perm_com = list_perm_com, 
                  SES_stat_obs = SES_stat_obs, 
                  tab_stat_obs_is_greater_pvalues = tab_stat_obs_is_greater_pvalues, 
                  tab_stat_obs_is_smaller_pvalues = tab_stat_obs_is_smaller_pvalues,
                  tab_stat_obs_is_different_pvalues = tab_stat_obs_is_different_pvalues, 
                  SES_mean_stat_obs = SES_mean_stat_obs, 
                  tab_mean_stat_obs_is_greater_pvalues = tab_mean_stat_obs_is_greater_pvalues, 
                  tab_mean_stat_obs_is_smaller_pvalues = tab_mean_stat_obs_is_smaller_pvalues,
                  tab_mean_stat_obs_is_different_pvalues = tab_mean_stat_obs_is_different_pvalues))
      

      
      }
    
    
    close(pb)
  }
  
  
}






################################################################################################################################################
###### à utiliser avec un objet "stat_obs" issu de tf_stats_from_rlq4 ######
################################################################################################################################################

stat_obs = stat_obs_all

sitlist = 1:20

tf_stats_from_rlq_permut_species2 <- function(stat_obs, iter, method = "traitswap", recompute.rlq.sp = T, class = NULL, class2 = NULL, 
                                              fba = NULL, tab_stats_obs = NULL,
                                             mean_stats_classe = NULL, compare = F, sitlist = NULL ){
  
  ##### null model : compute multiple FT stats for communities from RLQ with iterations on species permutation :
  ##### method "traitswap" permute synthetic traits (rlq coordinates)
  ##### method "spswap" permute species in the community matrix with respect to global sp abondances & nb of individuals per plot 
  ###### à utiliser avec un objet "stat_obs" issu de tf_stats_from_rlq4 ######
  
  
  library(vegan)
  library(FD)
  
  
  list_mean_stats_classe_null = NULL
  list_mean_stats_classe2_null = NULL
  
  rbind.all.columns <- function(x, y) {
    
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    
    x[, c(as.character(y.diff))] <- NA
    
    y[, c(as.character(x.diff))] <- NA
    
    return(rbind(x, y))
  }
  
  
  species = stat_obs$species
  
  com = data.frame(stat_obs$commat, check.names = F)
  
  if(!is.null(sitlist)){
    
    com = com[sitlist,]
    
    com = com[,!apply(com,2,sum)==0]
    
    species = colnames(com)
    
    stat_obs$rlq.abund$lQ = stat_obs$rlq.abund$lQ[rownames(stat_obs$rlq.abund$lQ) %in% colnames(com),]
    stat_obs$rlq.sp$lQ = stat_obs$rlq.sp$lQ[rownames(stat_obs$rlq.sp$lQ) %in% colnames(com),]
    
  }
  
  
  list_tab_stats_null = list()
  
  list_mean_stats_classe_null = list()
  list_mean_stats_classe2_null = list()
  list_perm_com_abund = list()
  list_perm_com_sp = list()

  
  perm_coord_tmp_abund =  stat_obs$rlq.abund$lQ
 
  
  print(paste("Null Model : ",iter," iterations"))
  pb <- txtProgressBar(min = 0, max =  iter, style = 3)
  for (k in 1:iter){
    
    
    
    if(method == "spswap"){
      
      perm = permatfull(com,fixedmar = "both", times = 1)
      
      com = data.frame(perm$perm, check.names = F)
      
      list_perm_com_abund[[k]] = com
      
      sp_tmp = species
      
    }else if(method == "traitswap"){
      
      # permute sp
      
      sp_tmp = sample(species)
      
      # merge with rlq coordinates
      
      rownames(perm_coord_tmp_abund) = sp_tmp
      
      list_perm_com_abund[[k]] = perm_coord_tmp_abund
    }else{ stop("method must be traitswap or spswap")}
    

    
    ##### compute stats  ####
    
    
    orlq = stat_obs$rlq.abund
    list_scores_sp_com = list()
    Ranges_Ax1_abund = c()
    Ranges_Ax2_abund = c()
    CWV_Ax1 = c()
    CWV_Ax2 = c()
    FDis_sp = c()
    FDis_abund = c()
    skewness_Ax1_abund = c()
    skewness_Ax2_abund = c()
    kurtosis_Ax1_abund = c()
    kurtosis_Ax2_abund = c()
    Ranges_Ax1_abund_nst = c() 
    Ranges_Ax2_abund_nst = c() 
    
    df_FD = data.frame()
    
    for (i in 1:nrow(com)){
      
      # species de la com avec les abondances
      
      sp_com_abund = rep(names(com[i,]), com[i,]) 
      
      # spécies (occurences) de la com
      
      com_tmp =  perm_coord_tmp_abund[rownames( perm_coord_tmp_abund) %in%  sp_com_abund,]
      
      com_abund_tmp =  perm_coord_tmp_abund[as.character(sp_com_abund),]
      
      list_scores_sp_com = com_abund_tmp
      
      Ranges_Ax1_abund_nst = c(Ranges_Ax1_abund_nst,  max(com_abund_tmp$AxcQ1) - min(com_abund_tmp$AxcQ1) ) 
      Ranges_Ax2_abund_nst = c(Ranges_Ax2_abund_nst,  max(com_abund_tmp$AxcQ2) - min(com_abund_tmp$AxcQ2) ) 
      
      
      Ranges_Ax1_abund = c(Ranges_Ax1_abund, ( max(com_abund_tmp$AxcQ1) - min(com_abund_tmp$AxcQ1) ) / (max(orlq$lQ$AxcQ1) - min(orlq$lQ$AxcQ1)) )
      Ranges_Ax2_abund = c(Ranges_Ax2_abund, ( max(com_abund_tmp$AxcQ2) - min(com_abund_tmp$AxcQ2) ) / (max(orlq$lQ$AxcQ2) - min(orlq$lQ$AxcQ2)) )
      
      CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
      CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
      
      FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
      
      skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
      skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
      
      kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
      kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
      
      
      # from Laliberté et al. 2010
      
      invisible(capture.output(FD <- tryCatch({
        
        data.frame(dbFD(com_abund_tmp, stand.x = F))
      }, error=function(e){ a = t(data.frame(rep("NA",10)))
      colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
      rownames(a) =  "Community1"
      return(a)
      }
      )
      )) # from Laliberté et al. 2010
      
      
      
      if(dim(df_FD)[1] == 0){
        
        df_FD = FD
      }else{
        df_FD = rbind.all.columns(df_FD, FD)
      }
      
      
    }
    
    
    names(df_FD)[1] = "nbind"
    
    
    #### with occurences ####
    # perm_coord_tmp_sp
    
    # merge with rlq coordinates
    
    if(recompute.rlq.sp == T){
      perm_coord_tmp_sp =  stat_obs$rlq.sp$lQ
      
      rownames(perm_coord_tmp_sp) = sp_tmp
      
      list_perm_com_sp[[k]] = perm_coord_tmp_sp
    }else{
      
      perm_coord_tmp_sp =  stat_obs$rlq.abund$lQ
      
    }
    
    
    
    
    
    
    
    #### compute stats ####
    
    orlq = stat_obs$rlq.sp
    
    list_scores_sp_com_tmp = list()
    CV_Ax1 = c()
    CV_Ax2 = c()
    FDis_sp = c()
    Ranges_Ax1_sp_nst = c() 
    Ranges_Ax2_sp_nst = c() 
    Ranges_Ax1_sp = c()
    Ranges_Ax2_sp = c()
    skewness_Ax1_sp = c()
    skewness_Ax2_sp = c()
    kurtosis_Ax1_sp = c()
    kurtosis_Ax2_sp = c()
    df_FD_sp = data.frame()
    for (i in 1:nrow(com)){
      
      # species de la com 
      
      sp_com = names(com)[com[i,]>0] 
      
      
      com_tmp = perm_coord_tmp_sp[sp_com,]
      
      
      scores_sp_com = list()
      
      
      list_scores_sp_com_tmp[[i]] = com_tmp
      
      Ranges_Ax1_sp_nst = c(Ranges_Ax1_sp_nst,  max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) 
      Ranges_Ax2_sp_nst = c(Ranges_Ax2_sp_nst,  max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) ) 
      
      Ranges_Ax1_sp = c(Ranges_Ax1_sp, ( max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1) ) / (max(orlq$lQ$AxcQ1) - min(orlq$lQ$AxcQ1)) )
      Ranges_Ax2_sp = c(Ranges_Ax2_sp, ( max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2) ) / (max(orlq$lQ$AxcQ2) - min(orlq$lQ$AxcQ2)) )
      
      CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
      CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
      
      FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
      
      skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
      skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
      
      kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
      kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
      
      # from Laliberté et al. 2010
      invisible(capture.output(FD_sp <- tryCatch({
        
        data.frame(dbFD(com_tmp, stand.x = F))
      }, error=function(e){ a = t(data.frame(rep("NA",10)))
      colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
      rownames(a) =  "Community1"
      return(a)
      }
      )
      )) # from Laliberté et al. 2010
      
      if(dim(df_FD_sp)[1] == 0){
        
        df_FD_sp = FD_sp
      }else{
        df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
      }
      
      
    }
    
    
    names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
    
    paste(names(df_FD_sp),"_sp")
    
    df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
    
    # dataframe with stat per coms
    
    
    
    tab_stats_null_temp = cbind(Ranges_Ax1_abund, Ranges_Ax2_abund, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                                skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                                Ranges_Ax1_sp, Ranges_Ax2_sp, skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp, 
                                Ranges_Ax1_abund_nst, Ranges_Ax2_abund_nst, Ranges_Ax1_sp_nst, Ranges_Ax2_sp_nst)
    
    tab_stats_null_temp = data.frame(tab_stats_null_temp)
    
    tab_stats_null_temp[tab_stats_null_temp == -Inf] <- NA
    
    tab_stats_null_temp[tab_stats_null_temp=="<NA>"] = NA 
    
    tab_stats_null_temp = as.data.frame(lapply(tab_stats_null_temp, as.numeric))
    ##### stats for com classes #####
    
    
    if(is.null(class) == F){
      mean_stats_classe_null_temp = list()
      
      for (i in 1:length(levels(class))){ # i in classes
        
        coms_tmp = c(1:length(class))[ class == levels(class)[i] ]
        
        mean_stats_classe_null_temp[[levels(class)[i]]] = apply(tab_stats_null_temp[coms_tmp,],2,mean, na.rm = T)
        
        
      }
      mean_stats_classe_null_temp = t(data.frame(mean_stats_classe_null_temp))
      list_mean_stats_classe_null[[k]] = mean_stats_classe_null_temp
      
    }  
    
    if(is.null(class2) == F){
      mean_stats_classe2_null_temp = list()
      
      for (i in 1:length(levels(class2))){ # i in classes
        
        coms_tmp = c(1:length(class2))[ class2 == levels(class2)[i] ]
        
        mean_stats_classe2_null_temp[[levels(class2)[i]]] = apply(tab_stats_null_temp[coms_tmp,],2,mean, na.rm = T)
        
        
      }
      mean_stats_classe2_null_temp = t(data.frame(mean_stats_classe2_null_temp))
      list_mean_stats_classe2_null[[k]] = mean_stats_classe2_null_temp
      
    }
    
    
    #### stok stats for each iteration ####
    list_tab_stats_null[[k]] =  tab_stats_null_temp
    
    
    
    
    
    
    # update progress bar
    setTxtProgressBar(pb, k)
    
  }
  close(pb)
  
  
  
  if(compare == F){
    return(list(list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                list_mean_stats_classe2_null = list_mean_stats_classe2_null, 
                list_perm_com_abund = list_perm_com_abund, list_perm_com_sp = list_perm_com_sp, iterations = iter, method = method))
    
  }
  
  if(compare == T){
    
    #### compute SES and significativity of obs stats ####
    
    if(is.null(tab_stats_obs) == T){
      print("tab_stats_obs argument empty, tab_stats_obs needed!")
    }else{
      
      
      
      
      # SES
      # turn stats null to a 3d array 
      
      list2ary = function(input.list){  #input a list of lists
        
        rows.cols <- dim(input.list[[1]])
        
        sheets <- length(input.list)
        
        output.ary <- array(as.numeric(unlist(input.list)), dim = c(rows.cols, sheets))
        
        colnames(output.ary) <- colnames(input.list[[1]])
        
        row.names(output.ary) <- row.names(input.list[[1]])
        
        return(output.ary)    # output as a 3-D array
        
      }
      
      
      
      
      array_stats_null =  list2ary(list_tab_stats_null)
      
      mean_stats_null = apply(array_stats_null,c(1,2), mean)
      
      sd_stats_null = apply(array_stats_null,c(1,2), sd)
      
      
      SES_stat_obs = ( tab_stats_obs - mean_stats_null) / sd_stats_null
      
      # significativity
      
      print(paste("compare null model with ",iter," iterations and stats obs"))
      pb <- txtProgressBar(min = 0, max =  iter, style = 3)
      
      tab_stat_obs_is_greater = tab_stats_obs
      tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
      tab_stat_obs_is_smaller = tab_stat_obs_is_greater
      tab_stat_obs_is_different = tab_stat_obs_is_greater
      for (i in 1:length( list_tab_stats_null)){
        
        null = list_tab_stats_null[[i]]
        tmp = tab_stats_obs > data.frame( null) 
        
        tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
        
        tmp2 = tab_stats_obs < data.frame( null) 
        
        tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
        
        tmp3 = ( tab_stats_obs < data.frame( null) ) | ( tab_stats_obs > data.frame( null) )
        
        tab_stat_obs_is_different = tab_stat_obs_is_different + tmp3
        
        setTxtProgressBar(pb, i)
      }
      
      
      
      tab_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
      tab_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
      tab_stat_obs_is_different_pvalues = 1 - tab_stat_obs_is_different/(iter+1)
      
      
      
      # compare mean stats of classes
      
      if(!is.null(class) == T){
        
        
        array_stats_null =  list2ary(list_mean_stats_classe_null)
        
        mean_stats_null = apply(array_stats_null,c(1,2), mean)
        
        sd_stats_null = apply(array_stats_null,c(1,2), sd)
        
        
        SES_mean_stat_obs = ( tab_stats_obs - mean_stats_null) / sd_stats_null
        
        # significativity
        
        print(paste("compare null model with ",iter," iterations and mean stats obs for classes"))
        pb <- txtProgressBar(min = 0, max =  iter, style = 3)
        
        tab_stat_obs_is_greater = mean_stats_classe
        tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
        tab_stat_obs_is_smaller = tab_stat_obs_is_greater
        tab_stat_obs_is_different = tab_stat_obs_is_greater
        for (i in 1:length( list_mean_stats_classe_null)){
          
          
          null = list_mean_stats_classe_null[[i]]
          tmp = tab_stats_obs > data.frame( null) 
          
          tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
          
          tmp2 = tab_stats_obs < data.frame( null) 
          
          tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
          
          tmp3 = ( tab_stats_obs < data.frame( null) ) | ( tab_stats_obs > data.frame( null) )
          
          tab_stat_obs_is_different = tab_stat_obs_is_different + tmp3
          
          setTxtProgressBar(pb, i)
        }
        
        
        
        tab_mean_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
        tab_mean_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
        tab_mean_stat_obs_is_different_pvalues = 1 - tab_stat_obs_is_different/(iter+1)
        
      }
      
      
      
      
      
      return(list(list_tab_stats_null = list_tab_stats_null, list_mean_stats_classe_null = list_mean_stats_classe_null, 
                  list_perm_com = list_perm_com, 
                  SES_stat_obs = SES_stat_obs, 
                  tab_stat_obs_is_greater_pvalues = tab_stat_obs_is_greater_pvalues, 
                  tab_stat_obs_is_smaller_pvalues = tab_stat_obs_is_smaller_pvalues,
                  tab_stat_obs_is_different_pvalues = tab_stat_obs_is_different_pvalues, 
                  SES_mean_stat_obs = SES_mean_stat_obs, 
                  tab_mean_stat_obs_is_greater_pvalues = tab_mean_stat_obs_is_greater_pvalues, 
                  tab_mean_stat_obs_is_smaller_pvalues = tab_mean_stat_obs_is_smaller_pvalues,
                  tab_mean_stat_obs_is_different_pvalues = tab_mean_stat_obs_is_different_pvalues))
      
      
      
    }
    
    
    close(pb)
  }
  
  
}





################################################################################################################################################
################################################################################################################################################


#### NULL MODEL 3 : permut traits in community = permute (traits) abundances  #####

################################################################################################################################################
################################################################################################################################################

list_coms_traits = stat_obs_west$list_scores_abund_com



tf_permut_intracom = function(list_coms_traits, iter = 99){
  
  
  #### NULL MODEL  : permut traits in community = permute (traits) abundances  #####
  
  list_perm_com = list()
  list_tab_stats_null = list()
  for (k in 1:iter){
  
  tmp_list_prem_abund_com = list()
  skewness_Ax1_abund = c()
  skewness_Ax2_abund = c()
  kurtosis_Ax1_abund = c()
  kurtosis_Ax2_abund = c()
  CWV_Ax1 = c()
  CWV_Ax2 = c()
  for (i in 1:length(list_coms_traits)){
    
    com_tmp = list_coms_traits[[i]]
    sp_perm = sample(rownames(com_tmp))
    perm_com_temp = com_tmp
    rownames(perm_com_temp) = sp_perm
    
    
    skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(perm_com_temp$AxcQ1))
    skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(perm_com_temp$AxcQ2))
    
    kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(perm_com_temp$AxcQ1))
    kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(perm_com_temp$AxcQ2))
    
    CWV_Ax1 = c(CWV_Ax1, var(perm_com_temp$AxcQ1))
    CWV_Ax2 = c(CWV_Ax2, var(perm_com_temp$AxcQ2))
    
    tmp_list_prem_abund_com[[i]] = perm_com_temp
    
    
    
  }
  
  tab_stats_null_temp = cbind(CWV_Ax1, CWV_Ax2, skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund)
  
  list_tab_stats_null[[k]] = tab_stats_null_temp
  
  list_perm_com[[k]] = tmp_list_prem_abund_com
}
  
 
  
  return(list(list_tab_stats_null = list_tab_stats_null,
              list_perm_com = list_perm_com, iterations = iter))
             
}

a = tf_permut_intracom(list_coms_traits, iter = 2)
a
################################################################################################################################################
################################################################################################################################################


#### just compare stats obs and null #####

################################################################################################################################################
################################################################################################################################################

obs = stat_obs_all

#obs = stat_obs_sbois 
 

#      obs = stat_obs_all 
#      null = md_null_all 


#      obs = a 
#      null = b 


compare_stats_obs_null <- function(obs, null, sitlist = NULL, meanstat = F){
  

  tab_stats_obs = obs$tab_stats_obs
  
  list_tab_stats_null = null$list_tab_stats_null
  
  iter = length(list_tab_stats_null)

  if (!is.null(sitlist)){
    
    tab_stats_obs = tab_stats_obs[sitlist,]
    
  }
  
  
  #### compute SES and significativity of obs stats ####
  

    
    
    # SES
    # turn stats null to a 3d array 
    
    list2ary = function(input.list){  #input a list of lists
      
      rows.cols <- dim(input.list[[1]])
      
      sheets <- length(input.list)
      
      output.ary <- array(as.numeric(unlist(input.list)), dim = c(rows.cols, sheets))
      
      colnames(output.ary) <- colnames(input.list[[1]])
      
      row.names(output.ary) <- row.names(input.list[[1]])
      
      return(output.ary)    # output as a 3-D array
      
    }
    
    
    
    
    array_stats_null =  list2ary(list_tab_stats_null)
    
    mean_stats_null = apply(array_stats_null,c(1,2), mean, na.rm = T)
    
    sd_stats_null = apply(array_stats_null,c(1,2), sd, na.rm = T)
    
    
    SES_stat_obs = ( tab_stats_obs - mean_stats_null) / sd_stats_null
    
    # significativity
    
    print(paste("compare null model with ",iter," iterations and stats obs"))
    pb <- txtProgressBar(min = 0, max =  iter, style = 3)
    
    tab_stat_obs_is_greater = tab_stats_obs
    tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
    tab_stat_obs_is_smaller = tab_stat_obs_is_greater
    for (i in 1:iter){
      
      # difference en siprimant la derniere colonne : colone classe donc pas numeric
      snull = list_tab_stats_null[[i]]
      tmp = tab_stats_obs > data.frame(snull) 
      
    #  snull[!is.na(snull)] = 0
      
      tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
      
      tmp2 = tab_stats_obs < data.frame(snull) 
      
      tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
      
      setTxtProgressBar(pb, i)
    }
    

    
    tab_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
    tab_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
    tab_stat_obs_is_different_pvalues = pmin(tab_stat_obs_is_greater_pvalues,tab_stat_obs_is_smaller_pvalues)
    
    
    
    # compare mean stats of classes
    
    if(meanstat == T){
      
       if( ( !is.null(null$list_mean_stats_classe_null) & !is.null(obs$mean_stats_classe) ) == T ){
      
          mean_stats_classe = obs$mean_stats_classe
          
          list_mean_stats_classe_null = null$list_mean_stats_classe_null
          
          array_stats_null =  list2ary(list_mean_stats_classe_null)
          
          mean_stats_null = apply(array_stats_null,c(1,2), mean, na.rm = T)
          
          sd_stats_null = apply(array_stats_null,c(1,2), sd, na.rm = T)
          
          
          SES_mean_stat_obs = ( mean_stats_classe - mean_stats_null) / sd_stats_null
          SES_mean_stat_obs = data.frame(SES_mean_stat_obs)
          # significativity
          
          print(paste("compare null model with ",iter," iterations and mean stats obs for classes"))
          pb <- txtProgressBar(min = 0, max =  iter, style = 3)
          
          tab_stat_obs_is_greater = mean_stats_classe
          tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
          tab_stat_obs_is_smaller = tab_stat_obs_is_greater
          tab_stat_obs_is_different = tab_stat_obs_is_greater
          for (i in 1:iter){
            
            # difference en siprimant la derniere colonne : colone classe donc pas numeric
            
            snull = list_mean_stats_classe_null[[i]]
            tmp = mean_stats_classe > data.frame( snull) 
            #  snull[!is.na(snull)] = 0
            
            tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
            
            tmp2 = mean_stats_classe < data.frame( snull) 
            
            tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
            
            setTxtProgressBar(pb, i)
          }
          
          
          
          tab_mean_stat_obs_is_greater_pvalues = 1 - tab_stat_obs_is_greater/(iter+1)
          tab_mean_stat_obs_is_smaller_pvalues = 1 - tab_stat_obs_is_smaller/(iter+1)
          tab_mean_stat_obs_is_different_pvalues = pmin(tab_mean_stat_obs_is_greater_pvalues, tab_mean_stat_obs_is_smaller_pvalues)
          
       
          
       }
      
        if( ( !is.null(null$list_mean_stats_classe2_null) & !is.null(obs$mean_stats_classe2) ) == T ){
      
          mean_stats_classe = obs$mean_stats_classe2
          
          list_mean_stats_classe_null = null$list_mean_stats_classe2_null
          
          array_stats_null =  list2ary(list_mean_stats_classe_null)
          
          mean_stats_null = apply(array_stats_null,c(1,2), mean, na.rm = T)
          
          sd_stats_null = apply(array_stats_null,c(1,2), sd, na.rm = T)
          
          
          SES_mean_stat_obs2 = ( mean_stats_classe - mean_stats_null) / sd_stats_null
          SES_mean_stat_obs2 = data.frame(SES_mean_stat_obs2)
          # significativity
          
          print(paste("compare null model with ",iter," iterations and mean stats obs for classes 2"))
          pb <- txtProgressBar(min = 0, max =  iter, style = 3)
          
          tab_stat_obs_is_greater = mean_stats_classe
          tab_stat_obs_is_greater = tab_stat_obs_is_greater[] = 0
          tab_stat_obs_is_smaller = tab_stat_obs_is_greater
          tab_stat_obs_is_different = tab_stat_obs_is_greater
          for (i in 1:iter){
            
            # difference en siprimant la derniere colonne : colone classe donc pas numeric
            
            
            snull = list_mean_stats_classe_null[[i]]
            tmp = mean_stats_classe > data.frame( snull) 
            #  snull[!is.na(snull)] = 0
            
            tab_stat_obs_is_greater = tab_stat_obs_is_greater + tmp
            
            tmp2 = mean_stats_classe < data.frame( snull) 
            
            tab_stat_obs_is_smaller = tab_stat_obs_is_smaller + tmp2
            
            setTxtProgressBar(pb, i)
          }
          
          
          
          tab_mean_stat_obs_is_greater_pvalues2 = 1 - tab_stat_obs_is_greater/(iter+1)
          tab_mean_stat_obs_is_smaller_pvalues2 = 1 - tab_stat_obs_is_smaller/(iter+1)
          tab_mean_stat_obs_is_different_pvalues2 = pmin(tab_mean_stat_obs_is_greater_pvalues2, tab_mean_stat_obs_is_smaller_pvalues2)
          
          
          return(list(SES_stat_obs = SES_stat_obs, 
                      tab_stat_obs_is_greater_pvalues = tab_stat_obs_is_greater_pvalues, 
                      tab_stat_obs_is_smaller_pvalues = tab_stat_obs_is_smaller_pvalues,
                      tab_stat_obs_is_different_pvalues = tab_stat_obs_is_different_pvalues, 
                      SES_mean_stat_obs = SES_mean_stat_obs, 
                      tab_mean_stat_obs_is_greater_pvalues = tab_mean_stat_obs_is_greater_pvalues, 
                      tab_mean_stat_obs_is_smaller_pvalues = tab_mean_stat_obs_is_smaller_pvalues,
                      tab_mean_stat_obs_is_different_pvalues = tab_mean_stat_obs_is_different_pvalues,
                      SES_mean_stat_obs2 = SES_mean_stat_obs2, 
                      tab_mean_stat_obs_is_greater_pvalues2 = tab_mean_stat_obs_is_greater_pvalues2, 
                      tab_mean_stat_obs_is_smaller_pvalues2 = tab_mean_stat_obs_is_smaller_pvalues2,
                      tab_mean_stat_obs_is_different_pvalues2 = tab_mean_stat_obs_is_different_pvalues2))
          
        } else if ( ( !is.null(null$list_mean_stats_classe_null) & !is.null(obs$mean_stats_classe) ) == T ) {
            return(list(SES_stat_obs = SES_stat_obs, 
                      tab_stat_obs_is_greater_pvalues = tab_stat_obs_is_greater_pvalues, 
                      tab_stat_obs_is_smaller_pvalues = tab_stat_obs_is_smaller_pvalues,
                      tab_stat_obs_is_different_pvalues = tab_stat_obs_is_different_pvalues, 
                      SES_mean_stat_obs = SES_mean_stat_obs, 
                      tab_mean_stat_obs_is_greater_pvalues = tab_mean_stat_obs_is_greater_pvalues, 
                      tab_mean_stat_obs_is_smaller_pvalues = tab_mean_stat_obs_is_smaller_pvalues,
                      tab_mean_stat_obs_is_different_pvalues = tab_mean_stat_obs_is_different_pvalues))
          
        } 
      
      } else {
        return(list(SES_stat_obs = SES_stat_obs, 
                  tab_stat_obs_is_greater_pvalues = tab_stat_obs_is_greater_pvalues, 
                  tab_stat_obs_is_smaller_pvalues = tab_stat_obs_is_smaller_pvalues,
                  tab_stat_obs_is_different_pvalues = tab_stat_obs_is_different_pvalues))
      
    }
    
    
    
    
    

    # update progress bar
    
    
  
  close(pb)
}





################################################################################################################################################
################################################################################################################################################


# test pour oecosimu

################################################################################################################################################
################################################################################################################################################






tf_stats_for_oecosimu <- function(x){
  
  # compute multiple FT stats for communities from RLQ ==> RLQ from abundance matrix & rlq from pres/abs matrix 
  # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
  
  
  library(ade4)
  library(FD)
  library(moments)
  library(utils)
  
  
  
  # obs datas
  
  compFBA = list()
  
  # with abundances
  
  compFBA$spe <- data.frame(x, check.names = F)
  
  compFBA$traits = TF_kone_com_nona
  
  env_var_select = env_var_keep
  
  
  compFBA$env = env_var_select
  
  
  ######  perform the separate analyses of each table = CoA & PCA ##### 
  
  
  afcL.FBA <- dudi.coa(compFBA$spe, scannf = FALSE)
  
  acpR.FBA <- dudi.hillsmith(compFBA$env, row.w = afcL.FBA$lw,
                             scannf = FALSE)
  
  
  acpQ.FBA <- dudi.pca(compFBA$traits, row.w = afcL.FBA$cw,
                       scannf = FALSE)
  
  #### RLQ ####
  
  rlq.FBA <- rlq(acpR.FBA, afcL.FBA, acpQ.FBA,
                 scannf = FALSE)
  
  
  
  
  sites_nona = unique(fba$id_pts)
  list_scores_abund_com = list()
  Ranges_Ax1 = c()
  Ranges_Ax2 = c()
  CWV_Ax1 = c()
  CWV_Ax2 = c()
  FDis_sp = c()
  FDis_abund = c()
  skewness_Ax1_abund = c()
  skewness_Ax2_abund = c()
  kurtosis_Ax1_abund = c()
  kurtosis_Ax2_abund = c()
  df_FD = data.frame()
  
  for (i in 1:nrow(com)){
    
    # species de la com avec les abondances
    
    sp_com_abund = rep(names(com[i,]), com[i,]) 
    
    
    
    # spécies (occurences) de la com
    
    com_tmp = rlq.FBA$lQ[rownames(rlq.FBA$lQ) %in%  sp_com_abund,]
    
    com_abund_tmp = rlq.FBA$lQ[as.character(sp_com_abund),]
    
    list_scores_abund_com[[i]] = com_abund_tmp
    
    Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))
    Ranges_Ax2 = c(Ranges_Ax2, max(com_tmp$AxcQ2) - min(com_tmp$AxcQ2))
    
    CWV_Ax1 = c(CWV_Ax1, var(com_abund_tmp$AxcQ1))
    CWV_Ax2 = c(CWV_Ax2, var(com_abund_tmp$AxcQ2))
    
    FDis_abund = c(FDis_abund, mean(dist(com_abund_tmp)))
    
    skewness_Ax1_abund = c(skewness_Ax1_abund, skewness(com_abund_tmp$AxcQ1))
    skewness_Ax2_abund = c(skewness_Ax2_abund, skewness(com_abund_tmp$AxcQ2))
    
    kurtosis_Ax1_abund = c(kurtosis_Ax1_abund, kurtosis(com_abund_tmp$AxcQ1))
    kurtosis_Ax2_abund = c(kurtosis_Ax2_abund, kurtosis(com_abund_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD <- tryCatch({
      
      data.frame(dbFD(com_abund_tmp, stand.x = F))
    }, error=function(e){ a = t(data.frame(rep("NA",10)))
    colnames(a) =   c("nbsp", "sing.sp","FRic", "qual.FRic", "FEve", "FDis",  "RaoQ",  "CWM.AxcQ1", "CWM.AxcQ2",  "FDiv")
    rownames(a) =  "Community1"
    return(a)
    }
    )
    )) # from Laliberté et al. 2010
    
    
    if(dim(df_FD)[1] == 0){
      
      df_FD = FD
    }else{
      df_FD = rbind.all.columns(df_FD, FD)
    }
    
  }
  
  
  names(df_FD)[1] = "nbind"
  
  
  
  #### with occurences BUT FROM RLQ COMPUTED ONLY WITH ABUNDANCES MATRIX ####
  
  
  
  #### RLQ ####
  
  
  
  list_scores_sp_com = list()
  CV_Ax1 = c()
  CV_Ax2 = c()
  FDis_sp = c()
  skewness_Ax1_sp = c()
  skewness_Ax2_sp = c()
  kurtosis_Ax1_sp = c()
  kurtosis_Ax2_sp = c()
  df_FD_sp = data.frame()
  for (i in 1:nrow(com)){
    
    # species de la com 
    
    sp_com = names(com)[ com[i,] > 0] 
    
    
    com_tmp = rlq.FBA$lQ[sp_com,]
    
    
    scores_sp_com = list()
    
    
    list_scores_sp_com[[i]] = com_tmp
    
    
    CV_Ax1 = c(CV_Ax1, var(com_tmp$AxcQ1))
    CV_Ax2 = c(CV_Ax2, var(com_tmp$AxcQ2))
    
    FDis_sp = c(FDis_sp, mean(dist(com_tmp)))
    
    skewness_Ax1_sp = c(skewness_Ax1_sp, skewness(com_tmp$AxcQ1))
    skewness_Ax2_sp = c(skewness_Ax2_sp,skewness(com_tmp$AxcQ2))
    
    kurtosis_Ax1_sp = c(kurtosis_Ax1_sp, kurtosis(com_tmp$AxcQ1))
    kurtosis_Ax2_sp = c(kurtosis_Ax2_sp,kurtosis(com_tmp$AxcQ2))
    
    
    
    invisible(capture.output(FD_sp <- tryCatch({
      
      data.frame(dbFD(com_tmp, stand.x = F))
    }, error=function(e){ a = t(data.frame(rep("NA",10)))
    colnames(a) =    c("nbsp", "sing.sp","FRic", "qual.FRic","FEve","FDiv","FDis", "RaoQ","CWM.AxcQ1","CWM.AxcQ2")
    rownames(a) =  "Community1"
    return(a)
    }
    )
    )) # from Laliberté et al. 2010
    
    if(dim(df_FD_sp)[1] == 0){
      
      df_FD_sp = FD_sp
    }else{
      df_FD_sp = rbind.all.columns(df_FD_sp, FD_sp)
    }
    
    
  }
  
  
  
  names(df_FD_sp) =  paste(names(df_FD_sp),"_sp", sep = "")
  
  df_FD_sp = df_FD_sp[,3:dim(df_FD_sp)[2]]  # supprime les column redondantes avec df_FD
  
  # dataframe with stat per coms
  
  
  
  tab_stats_obs = cbind(Ranges_Ax1, Ranges_Ax2, CWV_Ax1, CWV_Ax2, CV_Ax1, CV_Ax2, FDis_sp, FDis_abund,
                        skewness_Ax1_abund, skewness_Ax2_abund, kurtosis_Ax1_abund,  kurtosis_Ax2_abund,
                        skewness_Ax1_sp, skewness_Ax2_sp, kurtosis_Ax1_sp, kurtosis_Ax2_sp, df_FD, df_FD_sp)
  
  tab_stats_obs = data.frame(tab_stats_obs)
  
  tab_stats_obs[tab_stats_obs == -Inf] <- NA
  
  tab_stats_obs[tab_stats_obs=="<NA>"] = NA 
  tab_stats_obs = as.data.frame(lapply(tab_stats_obs, as.numeric))

  
  return(list( stat = tab_stats_obs$Ranges_Ax1, tab_stats_obs = tab_stats_obs, mean_stats_classe = mean_stats_classe, 
               list_scores_abund_com = list_scores_abund_com, list_scores_sp_com = list_scores_sp_com ))
}



tf_stats_for_oecosimu <- function(x){
  
  # compute multiple FT stats for communities from RLQ ==> RLQ from abundance matrix & rlq from pres/abs matrix 
  # null model : compute multiple FT stats for communities from RLQ with iterations on species permutation 
  
  
  library(ade4)
  library(FD)
  library(moments)
  library(utils)
  
  com = x
  
  # obs datas
  
  compFBA = list()
  
  # with abundances
  
  compFBA$spe <- data.frame(x, check.names = F)
  
  compFBA$traits = TF_kone_com_nona
  
  env_var_select = env_var_keep
  
  
  compFBA$env = env_var_select
  
  
  ######  perform the separate analyses of each table = CoA & PCA ##### 
  
  
  afcL.FBA <- dudi.coa(compFBA$spe, scannf = FALSE)
  
  acpR.FBA <- dudi.hillsmith(compFBA$env, row.w = afcL.FBA$lw,
                             scannf = FALSE)
  
  
  acpQ.FBA <- dudi.pca(compFBA$traits, row.w = afcL.FBA$cw,
                       scannf = FALSE)
  
  #### RLQ ####
  
  rlq.FBA <- rlq(acpR.FBA, afcL.FBA, acpQ.FBA,
                 scannf = FALSE)
  

  
 
  Ranges_Ax1 = c()

  
  for (i in 1:nrow(com)){
    
    # species de la com avec les abondances
    
    sp_com_abund = rep(names(com[i,]), com[i,]) 
    
    
    
    # spécies (occurences) de la com
    
    com_tmp = rlq.FBA$lQ[rownames(rlq.FBA$lQ) %in%  sp_com_abund,]

    
    Ranges_Ax1 = c(Ranges_Ax1, max(com_tmp$AxcQ1) - min(com_tmp$AxcQ1))

    
  }
  
  
  names(df_FD)[1] = "nbind"
  

  
  return(stat = tab_stats_obs$Ranges_Ax1)
}













