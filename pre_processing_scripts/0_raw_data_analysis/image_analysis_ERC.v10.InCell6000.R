#!/usr/bin/R
library(EBImage)
library(FNN)
library(geometry)
library(stringr)
library(e1071)
library(dplyr)
library(RPostgreSQL)

test_db <- src_postgres(
  dbname = "incell2000_test",
  host = "b110-sc2sn01",
  user = "florianH",
  password = "x!Kl9R_p7XZYjLhg"
)

analyze<-function(image_string){
  start.time <- Sys.time()
  images=unlist(strsplit(x=image_string,split="::")) 
  nuc_image_name=images[1] 
  body_image_name=images[2] 
  mito_image_name=images[3]
  if(file.exists(nuc_image_name) && file.exists(body_image_name) && file.exists(mito_image_name) ){
    
    #badcells=c()
    
    identifier=sub(".*/(.+)_\\w+.tif","\\1",nuc_image_name,perl=TRUE);
    
    TEST_IMAGE_nuclei_raw=readImage(nuc_image_name)#[500:1500,500:1500]
    TEST_IMAGE_nuclei=gblur(normalize(log(TEST_IMAGE_nuclei_raw)),sigma=1)
    
    TEST_IMAGE_tubulin_raw=readImage(mito_image_name)#[500:1500,500:1500]
    TEST_IMAGE_tubulin=gblur(normalize(log(TEST_IMAGE_tubulin_raw)),sigma=1)  
    
    TEST_IMAGE_actin_raw=readImage(body_image_name)#[500:1500,500:1500]
    TEST_IMAGE_actin=gblur(normalize(log(TEST_IMAGE_actin_raw)),sigma=1)
    
    med=median(TEST_IMAGE_nuclei)
    mit=mean(TEST_IMAGE_nuclei) 
    
    actin_binary= TEST_IMAGE_actin > 0.2
    strange_behavior=c(0,0,0)
    names(strange_behavior)=c("dark","fussel","weak_tubulin")
    
    if(med<0.15){
      strange_behavior[1]=1
    }
    if(abs(med-mit)>0.06){
      strange_behavior[2]=1
    }
    #display(TEST_IMAGE_nuclei)
    if(diff(range(TEST_IMAGE_tubulin_raw))>0.1
       &&
       skewness(TEST_IMAGE_tubulin_raw)>1.6){
      if(strange_behavior[1]==1){
        
        nuclei_objects = bwlabel(fillHull(opening(thresh(TEST_IMAGE_nuclei,w=31,h=31,offset=0.065))))      
        
        body_binary=fillHull(opening( thresh(TEST_IMAGE_actin,w=120,h=120,offset = 0.007) + (TEST_IMAGE_tubulin >0.25)))
      }else{
        
        nuclei_objects = bwlabel(fillHull(opening(thresh(TEST_IMAGE_nuclei,w=21,h=21,offset=0.1),makeBrush(5))))
        
        body_binary=opening( (TEST_IMAGE_actin > 0.25) + (TEST_IMAGE_tubulin >0.20) )
      }
    }else{
      nuclei_objects = bwlabel(fillHull(opening(thresh(TEST_IMAGE_nuclei,w=21,h=21,offset=0.1),makeBrush(5))))
      
      body_binary=opening( (TEST_IMAGE_actin > 0.25) )
      strange_behavior[3]=1 
    }
    
    
    if(max(nuclei_objects)>30){    
      DNA_features = computeFeatures(nuclei_objects, ref=TEST_IMAGE_nuclei_raw, 
                                     methods.noref=c("computeFeatures.shape"),
                                     methods.ref=c("computeFeatures.basic", "computeFeatures.moment", "computeFeatures.haralick"),
                                     xname="DNA", refnames="0", properties=FALSE, expandRef=NULL, 
                                     basic.quantiles=c(0.05), haralick.scales=c(1)
      )   
      colnames(DNA_features) = sub(".0.", ".", colnames(DNA_features))    
      
      cell_bodies_objects= propagate(TEST_IMAGE_actin, seeds= nuclei_objects, mask=body_binary,lambda=3e-4)
      rm(body_binary)
      actin_features = computeFeatures(cell_bodies_objects, ref=TEST_IMAGE_actin_raw, 
                                       methods.noref=c("computeFeatures.shape"),
                                       methods.ref=c("computeFeatures.basic", "computeFeatures.moment", "computeFeatures.haralick"),
                                       xname="actin", refnames="0", properties=FALSE, expandRef=NULL, 
                                       basic.quantiles=c(0.05), haralick.scales=c(1)
      )   
      colnames(actin_features) = sub(".0.", ".", colnames(actin_features))     
      
      tubulin_features =  computeFeatures(cell_bodies_objects, ref=TEST_IMAGE_tubulin_raw,
                                          methods.noref = NULL,
                                          methods.ref=c("computeFeatures.basic", "computeFeatures.haralick"),
                                          xname="tubulin", refnames="0", properties=FALSE, expandRef=NULL, 
                                          basic.quantiles=c(0.05), haralick.scales=c(1)
      )
      colnames(tubulin_features) = sub(".0.", ".", colnames(tubulin_features))    
      img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei, red=TEST_IMAGE_tubulin)
      res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
      res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='white')
      writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 10, 8 )
      rm(TEST_IMAGE_actin,TEST_IMAGE_nuclei,TEST_IMAGE_tubulin) 
      
      nearest.neighbours = get.knn(actin_features[, c("actin.m.cx", "actin.m.cy")],k=30)[["nn.dist"]][,c(10,20,30)]
      colnames(nearest.neighbours) = c("dist.10.nn","dist.20.nn","dist.30.nn")
      
      nucleus.displacement = sqrt(rowSums((actin_features[,c("actin.m.cx", "actin.m.cy")] - DNA_features[,c("DNA.m.cx", "DNA.m.cy")])^2))
      names(nucleus.displacement) = c("nuclear.displacement")
      
      cell_features=cbind(
        tubulin_features,
        actin_features,
        DNA_features,
        nearest.neighbours,
        nucleus.displacement,
        rep(strange_behavior[1],times=nrow(nearest.neighbours)),
        rep(strange_behavior[2],times=nrow(nearest.neighbours)),
        rep(strange_behavior[3],times=nrow(nearest.neighbours))
      )
      colnames(cell_features)=c(
        colnames(tubulin_features),
        colnames(actin_features),
        colnames(DNA_features),
        colnames(nearest.neighbours),
        "nuclear.displacement",
        "dark",
        "fussel",
        "weak_tubulin"
      )
      
      barcodes=str_match(identifier,"(ERC_20x_4t_SYNA_(\\w+?)_(\\w+?)_(\\w+?)_(\\d+))")
      
      tmp=colnames(cell_features)
      tmpf=cell_features
      cell_features=
        cbind.data.frame(
          as.data.frame(matrix(barcodes[2:6],nrow=nrow(cell_features),ncol=length(barcodes[2:6]),byrow=TRUE)),
          cell_features
        )
      colnames(cell_features)=c("barcode","plate","screen","well","field",tmp)
      save(cell_features,file=paste(dir,"/",identifier,"_single_cell.RData",sep=""))
      
      #copy_to(test_db, cell_features, name="D1086_InCell6000_single_cell", temporary = FALSE, indexes = list("plate","screen", "well","field"))
      db_insert_into( con = test_db$con, table = "D1086_InCell6000_single_cell", values = cell_features)
      
      trim_mean_data=apply(tmpf, 2, mean,trim=0.01, na.rm=TRUE)
      names(trim_mean_data) <- featnames <- paste0(names(trim_mean_data),".tmean")      
      
            
      sd_data=apply(tmpf, 2, sd, na.rm=TRUE)
      names(sd_data) <- sdnames <- paste0(names(sd_data),".sd")  
      trim_mean_data=c(nrow(tmpf),trim_mean_data,sd_data)
      names(trim_mean_data) <- featnames <- c("cells",featnames,sdnames)
      
      trim_mean_data=c(barcodes[2:6],trim_mean_data)
      dim(trim_mean_data)=c(1,length(trim_mean_data))
      trim_mean_data=as.data.frame(trim_mean_data)
      names(trim_mean_data)=c("barcode","plate","screen","well","field",featnames)
      
      write.table(trim_mean_data,file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
      #copy_to(test_db, trim_mean_data, name="D1086_trimmed_mean", temporary = FALSE, indexes = list("plate","screen", "well","field"))
      db_insert_into( con = test_db$con, table = "D1086_InCell6000_trimmed_mean", values = trim_mean_data)
      
    }else{
      mocknames=c("tubulin.b.mean",
                  "tubulin.b.sd",
                  "tubulin.b.mad",
                  "tubulin.b.q005",
                  "tubulin.h.asm.s1",
                  "tubulin.h.con.s1",
                  "tubulin.h.cor.s1",
                  "tubulin.h.var.s1",
                  "tubulin.h.idm.s1",
                  "tubulin.h.sav.s1",
                  "tubulin.h.sva.s1",
                  "tubulin.h.sen.s1",
                  "tubulin.h.ent.s1",
                  "tubulin.h.dva.s1",
                  "tubulin.h.den.s1",
                  "tubulin.h.f12.s1",
                  "tubulin.h.f13.s1",
                  "actin.s.area",
                  "actin.s.perimeter",
                  "actin.s.radius.mean",
                  "actin.s.radius.sd",
                  "actin.s.radius.min",
                  "actin.s.radius.max",
                  "actin.b.mean",
                  "actin.b.sd",
                  "actin.b.mad",
                  "actin.b.q005",
                  "actin.m.cx",
                  "actin.m.cy",
                  "actin.m.majoraxis",
                  "actin.m.eccentricity",
                  "actin.m.theta",
                  "actin.h.asm.s1",
                  "actin.h.con.s1",
                  "actin.h.cor.s1",
                  "actin.h.var.s1",
                  "actin.h.idm.s1",
                  "actin.h.sav.s1",
                  "actin.h.sva.s1",
                  "actin.h.sen.s1",
                  "actin.h.ent.s1",
                  "actin.h.dva.s1",
                  "actin.h.den.s1",
                  "actin.h.f12.s1",
                  "actin.h.f13.s1",
                  "DNA.s.area",
                  "DNA.s.perimeter",
                  "DNA.s.radius.mean",
                  "DNA.s.radius.sd",
                  "DNA.s.radius.min",
                  "DNA.s.radius.max",
                  "DNA.b.mean",
                  "DNA.b.sd",
                  "DNA.b.mad",
                  "DNA.b.q005",
                  "DNA.m.cx",
                  "DNA.m.cy",
                  "DNA.m.majoraxis",
                  "DNA.m.eccentricity",
                  "DNA.m.theta",
                  "DNA.h.asm.s1",
                  "DNA.h.con.s1",
                  "DNA.h.cor.s1",
                  "DNA.h.var.s1",
                  "DNA.h.idm.s1",
                  "DNA.h.sav.s1",
                  "DNA.h.sva.s1",
                  "DNA.h.sen.s1",
                  "DNA.h.ent.s1",
                  "DNA.h.dva.s1",
                  "DNA.h.den.s1",
                  "DNA.h.f12.s1",
                  "DNA.h.f13.s1",
                  "dist.10.nn",
                  "dist.20.nn",
                  "dist.30.nn",
                  "nuclear.displacement",
                  "dark",
                  "fussel",
                  "weak_tubulin")
      mock=rep(NA_real_,times=length(mocknames))
      barcodes=str_match(identifier,"(ERC_20x_4t_SYNA_(\\w+?)_(\\w+?)_(\\w+?)_(\\d+))")
      #barcodes=str_match("ERC_20x_4t_SYNA_1002_S19_F10_2","(ERC_20x_4t_SYNA_(\\w+?)_(\\w+?)_(\\w+?)_(\\d+))")
      
      mock=c(barcodes[2:6],mock)
      dim(mock)=c(1,length(mock))
      mock=as.data.frame(mock)
      names(mock)=c("barcode","plate","screen","well","field",mocknames)
      
      for(nms in mocknames){
        mock[[nms]]=as.numeric(mock[[nms]])
        mock[[nms]]=NA_real_
      }
      
      save(mock,file=paste(dir,"/",identifier,"_single_cell.RData",sep=""))
      
      #copy_to(test_db, mock, name="D1086_single_cell", temporary = FALSE, indexes = list("plate","screen", "well","field"))
      db_insert_into( con = test_db$con, table = "D1086_InCell6000_single_cell", values = mock)
          
      mock$tmp=NA_real_
      mock=cbind.data.frame(mock,t(rep(NA_real_,times=length(mocknames))))
      names(mock)=c("barcode","plate","screen","well","field",
                    "cells",
                    "tubulin.b.mean.tmean",
                    "tubulin.b.sd.tmean",
                    "tubulin.b.mad.tmean",
                    "tubulin.b.q005.tmean",
                    "tubulin.h.asm.s1.tmean",
                    "tubulin.h.con.s1.tmean",
                    "tubulin.h.cor.s1.tmean",
                    "tubulin.h.var.s1.tmean",
                    "tubulin.h.idm.s1.tmean",
                    "tubulin.h.sav.s1.tmean",
                    "tubulin.h.sva.s1.tmean",
                    "tubulin.h.sen.s1.tmean",
                    "tubulin.h.ent.s1.tmean",
                    "tubulin.h.dva.s1.tmean",
                    "tubulin.h.den.s1.tmean",
                    "tubulin.h.f12.s1.tmean",
                    "tubulin.h.f13.s1.tmean",
                    "actin.s.area.tmean",
                    "actin.s.perimeter.tmean",
                    "actin.s.radius.mean.tmean",
                    "actin.s.radius.sd.tmean",
                    "actin.s.radius.min.tmean",
                    "actin.s.radius.max.tmean",
                    "actin.b.mean.tmean",
                    "actin.b.sd.tmean",
                    "actin.b.mad.tmean",
                    "actin.b.q005.tmean",
                    "actin.m.cx.tmean",
                    "actin.m.cy.tmean",
                    "actin.m.majoraxis.tmean",
                    "actin.m.eccentricity.tmean",
                    "actin.m.theta.tmean",
                    "actin.h.asm.s1.tmean",
                    "actin.h.con.s1.tmean",
                    "actin.h.cor.s1.tmean",
                    "actin.h.var.s1.tmean",
                    "actin.h.idm.s1.tmean",
                    "actin.h.sav.s1.tmean",
                    "actin.h.sva.s1.tmean",
                    "actin.h.sen.s1.tmean",
                    "actin.h.ent.s1.tmean",
                    "actin.h.dva.s1.tmean",
                    "actin.h.den.s1.tmean",
                    "actin.h.f12.s1.tmean",
                    "actin.h.f13.s1.tmean",
                    "DNA.s.area.tmean",
                    "DNA.s.perimeter.tmean",
                    "DNA.s.radius.mean.tmean",
                    "DNA.s.radius.sd.tmean",
                    "DNA.s.radius.min.tmean",
                    "DNA.s.radius.max.tmean",
                    "DNA.b.mean.tmean",
                    "DNA.b.sd.tmean",
                    "DNA.b.mad.tmean",
                    "DNA.b.q005.tmean",
                    "DNA.m.cx.tmean",
                    "DNA.m.cy.tmean",
                    "DNA.m.majoraxis.tmean",
                    "DNA.m.eccentricity.tmean",
                    "DNA.m.theta.tmean",
                    "DNA.h.asm.s1.tmean",
                    "DNA.h.con.s1.tmean",
                    "DNA.h.cor.s1.tmean",
                    "DNA.h.var.s1.tmean",
                    "DNA.h.idm.s1.tmean",
                    "DNA.h.sav.s1.tmean",
                    "DNA.h.sva.s1.tmean",
                    "DNA.h.sen.s1.tmean",
                    "DNA.h.ent.s1.tmean",
                    "DNA.h.dva.s1.tmean",
                    "DNA.h.den.s1.tmean",
                    "DNA.h.f12.s1.tmean",
                    "DNA.h.f13.s1.tmean",
                    "dist.10.nn.tmean",
                    "dist.20.nn.tmean",
                    "dist.30.nn.tmean",
                    "nuclear.displacement.tmean",
                    "dark.tmean",
                    "fussel.tmean",
                    "weak_tubulin.tmean",
                    "tubulin.b.mean.sd",
                    "tubulin.b.sd.sd",
                    "tubulin.b.mad.sd",
                    "tubulin.b.q005.sd",
                    "tubulin.h.asm.s1.sd",
                    "tubulin.h.con.s1.sd",
                    "tubulin.h.cor.s1.sd",
                    "tubulin.h.var.s1.sd",
                    "tubulin.h.idm.s1.sd",
                    "tubulin.h.sav.s1.sd",
                    "tubulin.h.sva.s1.sd",
                    "tubulin.h.sen.s1.sd",
                    "tubulin.h.ent.s1.sd",
                    "tubulin.h.dva.s1.sd",
                    "tubulin.h.den.s1.sd",
                    "tubulin.h.f12.s1.sd",
                    "tubulin.h.f13.s1.sd",
                    "actin.s.area.sd",
                    "actin.s.perimeter.sd",
                    "actin.s.radius.mean.sd",
                    "actin.s.radius.sd.sd",
                    "actin.s.radius.min.sd",
                    "actin.s.radius.max.sd",
                    "actin.b.mean.sd",
                    "actin.b.sd.sd",
                    "actin.b.mad.sd",
                    "actin.b.q005.sd",
                    "actin.m.cx.sd",
                    "actin.m.cy.sd",
                    "actin.m.majoraxis.sd",
                    "actin.m.eccentricity.sd",
                    "actin.m.theta.sd",
                    "actin.h.asm.s1.sd",
                    "actin.h.con.s1.sd",
                    "actin.h.cor.s1.sd",
                    "actin.h.var.s1.sd",
                    "actin.h.idm.s1.sd",
                    "actin.h.sav.s1.sd",
                    "actin.h.sva.s1.sd",
                    "actin.h.sen.s1.sd",
                    "actin.h.ent.s1.sd",
                    "actin.h.dva.s1.sd",
                    "actin.h.den.s1.sd",
                    "actin.h.f12.s1.sd",
                    "actin.h.f13.s1.sd",
                    "DNA.s.area.sd",
                    "DNA.s.perimeter.sd",
                    "DNA.s.radius.mean.sd",
                    "DNA.s.radius.sd.sd",
                    "DNA.s.radius.min.sd",
                    "DNA.s.radius.max.sd",
                    "DNA.b.mean.sd",
                    "DNA.b.sd.sd",
                    "DNA.b.mad.sd",
                    "DNA.b.q005.sd",
                    "DNA.m.cx.sd",
                    "DNA.m.cy.sd",
                    "DNA.m.majoraxis.sd",
                    "DNA.m.eccentricity.sd",
                    "DNA.m.theta.sd",
                    "DNA.h.asm.s1.sd",
                    "DNA.h.con.s1.sd",
                    "DNA.h.cor.s1.sd",
                    "DNA.h.var.s1.sd",
                    "DNA.h.idm.s1.sd",
                    "DNA.h.sav.s1.sd",
                    "DNA.h.sva.s1.sd",
                    "DNA.h.sen.s1.sd",
                    "DNA.h.ent.s1.sd",
                    "DNA.h.dva.s1.sd",
                    "DNA.h.den.s1.sd",
                    "DNA.h.f12.s1.sd",
                    "DNA.h.f13.s1.sd",
                    "dist.10.nn.sd",
                    "dist.20.nn.sd",
                    "dist.30.nn.sd",
                    "nuclear.displacement.sd",
                    "dark.sd",
                    "fussel.sd",
                    "weak_tubulin.sd")
      db_insert_into( con = test_db$con, table = "D1086_InCell6000_trimmed_mean", values = mock)
      #copy_to(test_db, mock, name="D1086_trimmed_mean", temporary = FALSE, indexes = list("plate","screen", "well","field"))
      
      for(nms in names(mock)){
        if(is.na(mock[[nms]])){
          mock[[nms]]=0
        }
      }
      
      write.table(mock,file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
      
      img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei, red=TEST_IMAGE_tubulin)
      res = paintObjects(nuclei_objects, img, opac=c(0.3),thick=0.3, col='white')
      writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 10, 8 )
      
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

args=commandArgs(trailingOnly = TRUE)
options("scipen"=100, "digits"=4,warn=-1,MulticoreParam=quote(MulticoreParam(workers=6)))
dir=args[length(args)];
args=as.list(args[-length(args)])
result=lapply(as.list(args),analyze)


