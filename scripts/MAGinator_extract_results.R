# requires : tools gcc/7.4.0 intel/perflibs/2020_update4 R/4.0.0 

#Rscript --vanilla MAGinator_extract_results.R data/simulated_data/MAGinator/collectionID_order.txt data/simulated_data/MAGinator/signature_genes/screened data/simulated_data/MAGinator/screened_flat


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("use: Rscript --vanilla MAGinator_extract_results.R <collectionID_order.txt> <screened/> <outdir>", call.=FALSE)
}

cluster_map <- args[1] #"collectionID_order.txt"
clust_file_dir <- args[2] #screened
outdir = args[3]

#print(cluster_map, clust_file_dir, outdir)

cluster_map <- read.csv(file=cluster_map, header = FALSE)
clust_files <- list.files(clust_file_dir)

for (clust_file in clust_files) {
  clust_id <- as.integer(strsplit(clust_file, "_")[[1]][2])
  clust_id_b1 <- clust_id+1
  #print(clust_id_b1)
  
  #Find matching clustername
  clust_name <- cluster_map$V1[clust_id_b1]
  #print(clust_name)
  clust_file_full <- paste0(clust_file_dir, "/", clust_file)
  clust_data <- readRDS(clust_file_full)
  
  
  desc_out <- paste0(outdir, "/", clust_name, "_desc.csv")
  kmers_out <- paste0(outdir, "/",clust_name, "_kmers.csv")
  
  write.csv(clust_data$MSE, desc_out, row.names=F)
  write.csv(clust_data$genes, kmers_out, row.names=F)
  
}

