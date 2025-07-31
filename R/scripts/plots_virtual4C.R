plot_virtual_4C <- function(s, my_hiccup_files, resolution){

requires(GenomicRanges, 
         GenomicFeatures, 
         TxDb.Dmelanogaster.UCSC.dm6.ensGene,
         rtracklayer, 
         org.Dm.eg.db, 
         dplyr, 
         tidyr, 
         ggplot2)
 
  
# Subset hiccup results to region of interest containing bait start-end
tryCatch({
  my_hiccup_files[[s]]
}, error = function(e) {
  stop(paste("Error: HICCUP files for sample", s, "not found."))
})
  
bedpe <- my_hiccup_files[[s]] |> 
  filter(  ((((x1-14308589) * (x2-14304455)) < 0 & chr1 == "chr3R")| 
              (((y1-14308589) * (y2-14304455)) < 0 & chr2 == "chr3R")) & 
             fdrBL <= 0.05 & 
             fdrDonut <= 0.05 & 
             fdrH <= 0.05 & 
             fdrV <= 0.05) |> 
  mutate(interaction_distance = (centroid1 - centroid2) / 1000) |> 
  mutate(archight = interaction_distance / max(interaction_distance)) |> 
  arrange(centroid2)

tryCatch({
  if(nrow(bedpe) == 0) {
    stop(paste("\n\nError: No valid HICCUP interactions found for sample", s, ".\n\n"))
  }
}, error = function(e) {
  warning(e$message)
  return(NULL)
})

region_start <- min(c(bedpe$x1, bedpe$y1)) - 50000
region_end <- max(c(bedpe$x2, bedpe$y2)) + 50000


# Extract genes overlapping HICCUP loop coordinates for region of interest

# Define a GRanges object for a specific region
my_regions_of_interest <- GRanges(seqnames = "chr3R",
                                  ranges = IRanges(start = c(bedpe$x1, bedpe$y1), 
                                                   end = c(bedpe$x2, bedpe$y2)))

# Extract all gene ranges from the TxDb object
all_genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Find overlaps between significant chromosomal contact coordinates and the gene ranges
overlaps <- findOverlaps(my_regions_of_interest, all_genes, type = "any", ignore.strand = TRUE)

# Extract the gene IDs of the overlapping genes
overlapping_gene_ids <- all_genes$gene_id[subjectHits(overlaps)]

# Convert flybase IDs into symbols
overlapping_gene_symbols <- select(org.Dm.eg.db, keys=overlapping_gene_ids, columns = c("SYMBOL"), keytype = "FLYBASE")

# Make 2-column df with genes to highlight in plot below
my_interacting_genes <- data.frame(gene = overlapping_gene_symbols$SYMBOL, color = "darkred")




# Import HiC data for region of interest
hicFile <- list.files(path = paste(WD,"juicer_results",s,sep = "/"),
                      pattern = paste0(s,"_inter_30.hic"), 
                      full.names = TRUE, , recursive = TRUE)

hicDataChromRegion <- readHic(file = hicFile,
                              chrom = "3R", assembly = Dm6,
                              chromstart = region_start, chromend = region_end,
                              resolution = resolution, res_scale = "BP", norm = "KR"
)


################################
# Make virtual-4C plot
################################

pdf(file = paste0("./virtual-4C-analysis_files/Plots/virtual4c-",s,"-",resolution,".pdf"), 
    height = 16, width = 16)

  # 1- create page -------------------------------
  pageCreate(width = 16, height = 20, default.units = "inches", ygrid = 0, xgrid = 0, showGuides = FALSE)
  
  # 2- Plot Hi-C data in region -------------------------------
  hicPlot <- plotHicTriangle(
    data = hicDataChromRegion, colorTrans = "linear",
    chrom = "chr3R",
    chromstart = region_start,
    chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = 11, width = 15, height = 8.5,
    just = c("left","bottom"), default.units = "inches",
    zrange = c(1,quantile(hicDataChromRegion$counts, probs=c(.95)))
  )
  
  # 3- Add genome label -------------------------------
  annoGenomeLabel(
    plot = hicPlot, scale = "Mb",
    x = 0.5, y = 11.1
  )
  
  # 4- Plot genes -------------------------------
  plotGenes(
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = 11.2, width = 15, height = 1, 
    just = c("left", "top"), default.units = "inches", stroke = 0.05, fontsize = 6,
    geneHighlights = my_interacting_genes
  )
  
  
  
  # 5- Plot signal track data Mock SuHw overlap -------------------------------
  
  plotText(
    label = "GSE200993 Mock SuHw overlap", fontsize = 8,
    x = 0.5, y = 12.3, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = GSE200993_Mock_SuHw.overlap.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = 12.4, width = 15, height = 0.3, linecolor = "#33ccff", 
    just = c("left", "top"), default.units = "inches"
  )
  
  # 6- Plot signal track data SuHw ChIRP -------------------------------
  plotText(
    label = "SuHw ChIRP", fontsize = 8,
    x = 0.5, y = 12.8, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = SuHw_ChIRP_conc.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = 12.9, width = 15, height = 0.3, linecolor = "#ff9900", 
    just = c("left", "top"), default.units = "inches"
  )
  
  # 7- Plot Virtual 4C -------------------------------
  plotText(
    label = paste("Virtual 4C", s), fontsize = 8,
    x = 0.5, y = 13.3, just = "left",
    default.units = "inches",
  )
  
  plotPairsArches(
    data = bedpe,
    chrom = "chr3R",
    chromstart = region_start,
    chromend = region_end,
    assembly = Dm6,
    fill = "#6600cc",
    linecolor = "#6600cc",
    flip = TRUE,
    archHeight = bedpe$archight,
    alpha = 1,
    x = 0.5, y = 13.4, width = 15, height = 1,
    just = c("left", "top"),
    default.units = "inches", baseline = TRUE, lwd = 0.2
  )
  dev.off()
  
  return()
}