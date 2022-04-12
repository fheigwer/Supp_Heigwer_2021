library(broom)
library(dbplyr)
library(FitAR)
library(Hmisc)
library(limma)
library(lubridate)
library(MASS)
library(parallel)
library(preprocessCore)
library(reshape2)
library(rlang)
library(RPostgreSQL)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

call_interactions<-function (x, TP, TargetNeg, QueryNeg, eps = 1e-04, maxiter = 100, na.rm = TRUE){
  if (missing(TP)) {
    MTP = 1
    TP = rep(1, dim(x)[1])
  }
  else {
    MTP = max(TP)
  }
  z <- x
  nr <- nrow(z)
  nc <- ncol(z)
  t <- 0
  r <- numeric(nr)
  c <- matrix(numeric(MTP * nc), nrow = MTP)
  oldsum <- 0
  for (iter in 1L:maxiter) {
    rdelta <- apply(z, 1L, median, na.rm = na.rm)
    z <- z - matrix(rdelta, nrow = nr, ncol = nc)
    r <- r + rdelta
    cdelta <- apply(z, 2L, function(s) {
      tapply(s, TP, median, na.rm = na.rm)
    })
    if (MTP == 1) {
      z <- z - t(matrix(cdelta, nrow = nc, ncol = nr))
    }
    else {
      z <- z - apply(cdelta, 2, function(s) {
        s[TP]
      })
    }
    c <- c + cdelta
    if (missing(QueryNeg)) {
      delta = median(c, na.rm = na.rm)
      if (!is.finite(delta)) {
        delta = 0
      }
    }
    else {
      delta <- mean(c[, QueryNeg], na.rm = TRUE)
      if (!is.finite(delta)) {
        delta = median(c, na.rm = na.rm)
        if (!is.finite(delta)) {
          delta = 0
        }
      }
    }
    c <- c - delta
    t <- t + delta
    if (missing(TargetNeg)) {
      delta = median(r, na.rm = na.rm)
      if (!is.finite(delta)) {
        delta = 0
      }
    }
    else {
      delta <- mean(r[TargetNeg], na.rm = TRUE)
      if (!is.finite(delta)) {
        delta = median(r, na.rm = na.rm)
        if (!is.finite(delta)) {
          delta = 0
        }
      }
    }
    r <- r - delta
    t <- t + delta
    newsum <- sum(abs(z), na.rm = na.rm)
    converged <- newsum == 0 || abs(newsum - oldsum) < eps * 
      newsum
    if (converged) 
      break
    oldsum <- newsum
  }
  if (!converged) {
    warning(gettextf("maineffects() did not converge in %d iterations", 
                     maxiter), domain = NA)
  }
  ans <- list(neg = t, targetMainEffect = r, queryMainEffect = t(c), 
              pi = z)
  ans
}

#call_int<-function(curr_feat_string){

curr_feat_string<-args[1]
  
  test_db <- src_postgres(dbname = "incell2000_test",
                          host = "b110-sc2sn01",
                          user = "florianH",
                          password = "x!Kl9R_p7XZYjLhg")
  
  
  successful_screens <- 
    tbl(test_db,'screenings_screening') %>% 
    filter(completed==T,qc_general==T) %>%
    collect() %>% 
    separate(plate_barcode,c("sid","screen")) %>% .$screen %>% unique()
  
  curr_feat_sym<-sym(curr_feat_string)
  print(curr_feat_string)
  
  syndata_1 <-
    tbl(test_db,'D1086_trimmed_mean') %>% 
    filter(screen %in% c(successful_screens,"S64")) %>%
    dplyr::select(screen,plate,well,field,value=curr_feat_string) %>% # change feature name here
    distinct() %>%  
    collect(n=Inf) 
  syndata_2 <- 
    tbl(test_db,'D1086_InCell6000_trimmed_mean') %>% 
    filter(screen %in% c(successful_screens,"S64")) %>%
    dplyr::select(screen,plate,well,field,value=curr_feat_string) %>% # change feature name here
    distinct() %>%
    collect(n=Inf) 
  
  syndata<-backup<-bind_rows(list("2200"=syndata_1,"InCell6000"=syndata_2),.id = "microscope")
  
  syndata %<>%
    mutate(well=gsub(well,pattern="_",replacement = "" ,perl = T)) %>% 
    mutate(well=gsub(well,pattern="([A-Z])0(\\d)",replacement = "\\1\\2" ,perl = T)) %>% 
    mutate(screen=gsub(screen,pattern="^S(\\d{2}$)",replacement = "S0\\1" ,perl = T)) %>% 
    mutate(plate=gsub(plate,pattern="Crtl",replacement = "CTRL" ,perl = T)) %>%
    complete(screen,plate,well,field)
  
  save(syndata,file=paste0("/data/heigwer/SYNGENE_interactions/syngene_data_",curr_feat_string,".RData"))
  
 # load(paste0("/data/heigwer/SYNGENE_interactions/syngene_data_",curr_feat_string,".RData"))
  
  syndata<-syndata %>%
    mutate(value=if_else(plate=="CTRL2" & screen %in% c("S111","S079","S190","S191","S192"),NA_real_,value),
           screen=factor(screen),
           plate=factor(plate),
           well=factor(well),
           field=factor(field)) %>%
    mutate(value=if_else(value==0,NA_real_,value))
  #####################################################################################################################################
  #Let's aggregate the cell number per well as the sum of all fields
  #####################################################################################################################################
  
  syndata %<>%
    group_by(screen,plate,well) %>%
    summarise(value=mean(value,na.rm=T)) %>%
    ungroup() 
  
  #normalization function that normalizes each plate to its median 
  #except the control plates who get normalized by the average of all sample plates
  norm_to_median<-function(x){
    screenval <- x %>% pull(screen) %>% unique()
    val <- ref_vals %>% filter(screen==screenval) %>% pull(m)
    if(filter(x,kind=="sample") %>% n_distinct() > 1){
      return(x %>% mutate(value=value/median(value[content=="sample"],na.rm=T),normed=1))
    }else{
      return(x %>% mutate(value=value/val,normed=0))
    }
  }
  #z-transform each data point by the median of the negative controls
  #of the control plates and the mad of the rluc controls of each sample plate
  norm_to_ctrlplates<-function(x){
    m <- 
      x %>% 
      filter(kind=="ctrl"#,
             #content=="rluc"
      ) %>% 
      summarise(m=median(value,na.rm=T))
    s <- 
      x %>% 
      filter(kind=="sample",content=="rluc") %>% 
      summarise(s=mad(value,na.rm=T))
    
    x %>% mutate(value=(value-m$m)/s$s)
  }
  #the former function introduced some artifact such that the variance per screen has ben artificially skewed
  #so now we z-transform each data point by the median of the control-plates
  #and scale the result by the mad of the rluc controls of the entire screened set
  norm_to_ctrlplates_2<-function(x){
    m <- 
      x %>% 
      filter(kind=="ctrl"#,
             #content=="rluc"
      ) %>% 
      summarise(m=median(value,na.rm=T)) %>%
      pull(m)
    x %>% mutate(value=(value-m)/ref_mad)
  }
  
  #collect the library annotations
  syngene_annotation_db <- tbl(test_db,"D1086_annotation") %>% collect(n=Inf)
  
  #merge raw data with annotation
  syndata %<>% left_join(syngene_annotation_db) %>% dplyr::select(screen,plate,well,content,value) 
  
  #digest average query effect of the respective feature on these screen
  ref_vals <- 
    syndata %>%  
    mutate(kind=if_else(grepl("CTRL",plate),"ctrl","sample"),
           value=FitAR::glog(value,a = quantile(value,0.03,na.rm=T))) %>%
    filter(kind=="sample") %>%
    group_by(screen) %>%
    summarise(m=median(value[content=="sample"],na.rm=T)) %>%
    ungroup() 
  
  #we norm each plate by its median of sample well if it is a sample plate,
  #else we normlize the plate by the average query effect
  feature_wise <-
    syndata %>%
    mutate(kind=if_else(grepl("CTRL",plate),"ctrl","sample"),
           rawvalue=value,
           value=FitAR::glog(value,a = quantile(value,0.03,na.rm=T))
    ) %>%
    group_by(screen,plate) %>%
    do(norm_to_median(.)) %>%
    ungroup()
  
  #we digest the median absolute deviation of the negative controls over all sample plates 
  #that were screened to estimate the expected effect variances that could be observed by chance
  ref_mad <-
    feature_wise %>% 
    filter(kind=="sample",content=="rluc") %>%
    summarise(ref_mad=mad(value,na.rm = T)) %>%
    pull(ref_mad)
  
  feature_wise %<>%
    group_by(screen) %>%
    do(norm_to_ctrlplates_2(.)) %>%
    mutate(design=ifelse(grepl("^20",plate),2,1)) %>% 
    extract(plate,"tmp","(\\d{2}$)",remove = F) %>% 
    mutate(targetid=paste(tmp,well,sep = "_")) %>% 
    ungroup() %>% 
    distinct()
  
  save(feature_wise,file=paste0("/data/heigwer/SYNGENE_interactions/normalized_data_",curr_feat_string,".RData"))
  
  annotated_data_to_plot <- feature_wise %>%
    ungroup() %>%
    dplyr::select(value,rawvalue,plate,screen,well,content) %>%
    #filter(content=="sample") %>%
    tidyr::extract(well,c("row","column"),"(\\w)(\\d+)",remove = F) %>%
    mutate(column=gsub("^(\\d)$","0\\1",column)) %>%
    mutate(
      row=factor(row,levels = sort(unique(row),decreasing = T)),
      column=factor(column)
    ) %>%
    mutate(row1=as.numeric(row),column1=as.numeric(column)) %>% 
    filter(!is.na(value)) %>%
    ungroup()
  
  
  #collect the query annotations
  syngene_query_annotation <- 
    tbl(test_db,"screenings_screening") %>% 
    collect(n=Inf) %>%
    tidyr::extract(plate_barcode,c("project","screen"),regex = "(D1086|D0186)_(.+)")
  
  annotated_data_to_plot %<>% 
    left_join(syngene_annotation_db) %>%
    dplyr::select(-id) %>%
    left_join(syngene_query_annotation) %>%
    tidyr::extract(hd3_id,c("design","gene_id"),regex = "(HD3.)(.+)") %>%
    tidyr::extract(plate,c("plt"),regex = "(\\w\\w)$",remove = F) %>%
    unite(targetid,plt,well,remove = F) %>%
    mutate(kind=if_else(plate %in% c("CTRL1","CTRL2"),"control","sample"))
  
  quantnorm<-function(x){
    y<-x %>% 
      dplyr::select(value,plate,targetid,screen) %>% 
      spread(screen,value) 
    
    z<-y %>%
      .[,3:4] %>%
      as.matrix() %>%
      normalize.quantiles() %>%
      as.data.frame() %>%
      tbl_df()
    
    y[,3:4]<-z
    
    q<-y %>%
      gather(screen,value,-plate,-targetid) %>%
      pull(value)
    
    if(length(q)==nrow(x)){
      #print(x)
    }else{
      print(x)
    }
    
    q
  }
  
  outlierremo<-function(x){
    y<-x %>% 
      dplyr::select(value,plate,targetid,screen) %>% 
      spread(screen,value) 
    
    mdl<-rlm(x=y[,3][[1]],y=y[,4][[1]])
    
    y$w<-if_else(mdl$w<1/5*mean(mdl$w,na.rm=T),"outlier","normal")
    
    y <- y %>% dplyr::select(plate,targetid,w)
    
    q <- x %>% left_join(y) %>% pull(w)
    
    if(length(q)==nrow(x)){
      #print(x)
    }else{
      print(x)
    }
    q
  }
  
  
  query_main<-annotated_data_to_plot %>% 
    filter(kind!="control") %>% 
    group_by(query_name) %>%
    summarise(query_main=median(value,na.rm=T))
  
  annotated_data_to_plot %<>%
    left_join(query_main) %>%
    filter(kind!="control") %>% 
    group_by(well,screen) %>%
    mutate(normvalue=value-median(value,na.rm=T)) %>%
    ungroup() %>%
    mutate(value=normvalue+query_main) 
  
  annotated_data_to_plot_woctrl_wellnorm <- annotated_data_to_plot %>% 
    group_by(query_name,plate,targetid) %>%
    filter(n()==2) %>%
    group_by(query_name) %>%
    do(mutate(.,value=quantnorm(.))) %>%
    do(mutate(.,outlier=outlierremo(.))) %>%
    ungroup() %>%
    mutate(value=if_else(outlier=="outlier",NA_real_,value))

  
  queries_in_replicate <-
    annotated_data_to_plot_woctrl_wellnorm %>% 
    group_by(screen,query_name) %>% 
    summarise(n=n()) %>% 
    group_by(query_name) %>% 
    filter(n()==2) %>% 
    pull(query_name) %>% 
    unique()
  
  filtered_data<-annotated_data_to_plot_woctrl_wellnorm %>%
    filter(!grepl("CTRL",plate),
           content=="sample",
           query_name %in% queries_in_replicate,
           !(query_name %in% c("Snr1","asp","Dlg5","SREBP","Hrs"))
    ) %>%
    dplyr::select(fbgn,design,gene_id,targetid,screen,query_name,value)
  
  save(filtered_data,file=paste0("/data/heigwer/SYNGENE_interactions/filtered_data_",curr_feat_string,".RData"))
  
  target_annotation<-annotated_data_to_plot_woctrl_wellnorm %>%
    filter(!grepl("CTRL",plate), content=="sample") %>%
    dplyr::select(targetid,fbgn,gene_id,gene_symbol) %>%
    distinct() %>% 
    group_by(targetid) %>% 
    do(head(.,n=1)) %>%
    ungroup()
  
  int_data <- filtered_data %>% 
    unite(query_name,query_name,screen) %>%
    spread(query_name,value)
  
  confounder <- filtered_data %>%
    group_by(targetid,query_name) %>%
    group_by(query_name) %>%
    mutate(query_main=median(value,na.rm = T)) %>%
    group_by(targetid) %>%
    mutate(target_main=median(value,na.rm = T)) %>%
    ungroup() %>%
    dplyr::select(-value,-screen,-fbgn,-design,-gene_id) %>%
    distinct() 
  
  interaction <- call_interactions(int_data[,5:ncol(int_data)]) #medpolish(int_data[,3:ncol(int_data)],eps=1e-4,maxiter=100,na.rm=T)$quanvalue
  
  interaction <- cbind.data.frame(int_data[,1:4],interaction$pi) %>% 
    gather(query_name,value,-fbgn,-design,-gene_id,-targetid) %>% 
    tbl_df() %>% 
    separate(query_name,c("query_name","screen"),sep = "_")
  
  interaction_stats <- 
    interaction %>% 
    group_by(query_name,targetid) %>% 
    mutate(repl=c(1:4)) %>%
    dplyr::select(targetid,query_name,value,repl) %>%
    drop_na() %>%
    ungroup()
  
  interaction_stats<-interaction_stats %>% left_join(confounder) %>% left_join(target_annotation)
  
  interaction_stats_biased<-
    interaction_stats %>% 
    group_by(targetid,query_name) %>% 
    summarise(value=mean(value,na.rm=T),query_main=query_main[1],target_main=target_main[1],gene_symbol=gene_symbol[1])
  
  fun<-function(x){
    x$combo=x$target_main+x$query_main
    mdl<- rlm(value~combo,data = x)
    x$correction=predict(mdl,x)
    x$correctedpi=x$value-x$correction
    return(x)
  }
  
  interaction_stats %<>% group_by(query_name) %>% do(
    fun(.)
  ) %>%
    dplyr::select(-value,-correction,-combo) %>%
    separate(targetid,c("plate","well"),remove = "F",sep="_") %>%
    group_by(well,query_name) %>%
    mutate(correctedpi=correctedpi-median(correctedpi,na.rm=T)) %>%
    ungroup() %>%
    spread(repl,correctedpi) %>%
    ungroup() %>% 
    mutate(
      pval=eBayes(lmFit(.[,10:13]))$p.value,
      mpi=rowMeans(.[,10:13],na.rm = T)) %>%
    group_by(query_name) %>%
    mutate(mpi=mpi/sd(mpi,na.rm=T)) %>%
    ungroup()
  
  interaction_stats$fdr <- p.adjust(interaction_stats$pval,method = "BH")
  
  save(interaction_stats,file=paste0("/data/heigwer/SYNGENE_interactions/statistically_tested_pi_",curr_feat_string,".RData"))
  
 # rm(list = as.list(ls()))
#}




