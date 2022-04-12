
library(tidyverse)
library(enrichR)
library(httr)
library(rjson)
library(mgsa)

gene_grps<- read_delim("objects/gene_group_data_fb_2019_02_mod.tsv",delim = "\t") %>% rename(gene_symbol=Group_member_FB_gene_symbol)

conversion<-read_delim("objects/SYNGENE_target_id_to_modern_symbol.txt",delim = "\t") %>% distinct() %>% rename(gene_symbol=current_symbol)

annotated_genes <- conversion %>% left_join(gene_grps) %>% drop_na()

gene_sets<- annotated_genes %>% select(FB_group_name,gene_symbol) %>% group_by(FB_group_name) %>% filter(n()>1) %>% ungroup() %>% distinct() %>% split(.$FB_group_name,drop = T) %>% lapply(function(x){x %>% pull(gene_symbol)})

# run once and save results
#go_grps<- read_delim("objects/go_term_bio_process_gene_symbol_filtered.txt",delim = "\t")
#go_grps_ref<-go_grps %>% gather(go_term,value,-gene_symbol) %>% filter(gene_symbol %in% unique(conversion$gene_symbol)) %>% group_by(go_term) %>% filter(value!=0) %>% filter(n()>=3) %>% ungroup()  %>% split(.$go_term,drop = T) %>% lapply(function(x){x %>% pull(gene_symbol)})
#saveRDS(go_grps_ref,"objects/go_term_bio_process_formattedformgsa.rds")

go_grps_ref <- readRDS("objects/go_term_bio_process_formattedformgsa.rds")

listEnrichrDbs<-function ()
{
  dfSAF <- options()$stringsAsFactors
  options(stringsAsFactors = FALSE)
  dbs <- GET(url = "http://amp.pharm.mssm.edu/FlyEnrichr/datasetStatistics")$content
  dbs <- intToUtf8(dbs)
  dbs <- fromJSON(dbs)
  dbs <- lapply(dbs$statistics, function(x) do.call(cbind.data.frame,
                                                    x))
  dbs <- do.call(rbind.data.frame, dbs)
  options(stringsAsFactors = dfSAF)
  dbs
}

enrichr<-function (genes, databases = NULL)
{
  cat("Uploading data to Enrichr... ")
  if (is.vector(genes)) {
    temp <- POST(url = "http://amp.pharm.mssm.edu/FlyEnrichr/enrich",
                 body = list(list = paste(genes, collapse = "\n")))
  }
  else if (is.data.frame(genes)) {
    temp <- POST(url = "http://amp.pharm.mssm.edu/FlyEnrichr/enrich",
                 body = list(list = paste(paste(genes[, 1], genes[,
                                                                  2], sep = ","), collapse = "\n")))
  }
  else {
    warning("genes must be a vector of gene names or a dataframe with genes and score.")
  }
  GET(url = "http://amp.pharm.mssm.edu/FlyEnrichr/share")
  cat("Done.\n")
  dbs <- as.list(databases)
  dfSAF <- options()$stringsAsFactors
  options(stringsAsFactors = FALSE)
  result <- lapply(dbs, function(x) {
    cat("  Querying ", x, "... ", sep = "")
    r <- GET(url = "http://amp.pharm.mssm.edu/FlyEnrichr/export",
             query = list(file = "API", backgroundType = x))
    r <- intToUtf8(r$content)
    tc <- textConnection(r)
    r <- read.table(tc, sep = "\t", header = TRUE, quote = "")
    close(tc)
    cat("Done.\n")
    return(r)
  })
  options(stringsAsFactors = dfSAF)
  cat("Parsing results... ")
  names(result) <- dbs
  cat("Done.\n")
  return(result)
}


sudden_gain_big <-
  results_data %>%
  as_tibble() %>%
  unite("QueTarg",query_name,target_name,plate,well) %>%
  unite("Population",context,group) %>%
  filter(fdr<=0.1,mpi>0) %>%
  dplyr::select(QueTarg,Population) %>%
  mutate(interacts=1) %>%
  spread(Population,interacts,fill=0) %>%
  mutate(total =  rowSums(across(where(is.numeric)))) %>%
  filter(total==1, Isolated_Elongated == 1) %>% # isolated condensed, isolated elongated isolated big, isolated irregular nucleus
  tidyr::separate(col=QueTarg,into=c("query_name","target_name","plate","well"),sep="_") %>%
  select_if(is.character) %>%
  filter(!grepl("^Rp",target_name))



geneset = sudden_gain_big %>%
  group_by(target_name) %>%
  count() %>%
  arrange(desc(n)) %>%
  head(100) %>%
  pull(target_name)

dbs <- c("Coexpression_Predicted_GO_Biological_Process_2018")#,"GO_Biological_Process_2018")
cat("Enrichment in progress\n")
print(sort(unique(geneset)))
cat("--------------------------------------------\n")
enriched <- enrichr(unique(geneset), dbs)
cat("--------------------------------------------\n")
#cat("Enrichment in Coexpression_Predicted_GO_Biological_Process_2018\n")
print(head(enriched$Coexpression_Predicted_GO_Biological_Process_2018 %>% arrange(desc(Combined.Score))))
cat("--------------------------------------------\n")
cat("Enrichment in GO_Biological_Process_2018\n")
#print(head(enriched$GO_Biological_Process_2018 %>% arrange(desc(Combined.Score))))
cat("--------------------------------------------\n")
cat("Enrichment in Flybase_gene_families\n")
res<-mgsa(unique(geneset),gene_sets)
setsResults(res) %>% as_tibble(rownames = 'gene_family')%>% arrange(desc(estimate)) %>% head(20) %>% print()
cat("Enrichment in Flybase_gene GO_terms\n")
res2<-mgsa(unique(geneset),go_grps_ref)
setsResults(res2) %>% as_tibble(rownames = 'GO_terms') %>% arrange(desc(estimate)) %>% head(10) %>% print()


p1 <- enriched$Coexpression_Predicted_GO_Biological_Process_2018 %>%
  arrange(desc(Combined.Score)) %>%
  head(20) %>%
  ggplot(aes(x=-log10(P.value),y=reorder(Term,-P.value,mean),size=Combined.Score)) +
    geom_point() +
    theme_cowplot() +
    ylab("enriched terms") +
    ggtitle('Enriched Coexpression BioProcess Elongated')

p2 <- setsResults(res) %>%
  arrange(desc(estimate)) %>%
  head(20) %>%
  rownames_to_column("Term") %>%
  ggplot(aes(x=estimate,y=reorder(Term,estimate,mean),size=inStudySet/inPopulation)) +
  geom_point() +
  theme_cowplot() +
  ylab("enriched terms") +
  ggtitle('Enriched BioProcess Elongated')

p3 <- p1 +p2 + plot_layout(nrow = 2,ncol = 1)


ggsave("Genesetenrichment_suddenElongated2.pdf",width = 16,height = 20)
