library(dbplyr)
library(RPostgreSQL)
library(dplyr)


test_db <- src_postgres(dbname = "incell2000_test",
                        host = "b110-sc2sn01",
                        user = "florianH",
                        password = "x!Kl9R_p7XZYjLhg")

feature_list<-tbl(test_db,'D1086_trimmed_mean') %>% collect(n=1) %>% names() %>% .[c(6:83,87:163)] %>% .[!grepl("cx|cy|theta",.)]

for(i in feature_list){
  system(command = paste0("echo 'Rscript /data/heigwer/SYNGENE_interactions/interaction_feature_analysis_v4.R ",i,"' | qsub -l walltime=20:00:00 -l nodes=b110-sc2cn01.inet.dkfz-heidelberg.de:ppn=2"))
  print(paste0("echo 'Rscript /data/heigwer/SYNGENE_interactions/interaction_feature_analysis_v4.R ",i,"' | qsub -l walltime=20:00:00 -l nodes=b110-sc2cn01.inet.dkfz-heidelberg.de"))
}
