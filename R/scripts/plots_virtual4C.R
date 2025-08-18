plot_virtual_4C <- function(s, my_hiccup_files, resolution){

my_key <- paste0(s,":",resolution)
bedpe <- tibble()

# Import HiC data for region of interest
hicFile <- list.files(path = paste(WD,"juicer_results",s,sep = "/"),
                      pattern = paste0("*inter_30.hic"), 
                      full.names = TRUE, , recursive = TRUE)

print(paste("Importing Hi-C data from file:", hicFile))
hicDataChromRegion <- readHic(file = hicFile,
                              chrom = "3R", assembly = Dm6,
                              resolution = resolution, res_scale = "BP", norm = "KR"
)

# Add contact scores to my_roi_all.gr for plotting ------------------

# Entering observed twice to account for from and to contact coords that are in two diff rows

contacts.df <- hicDataChromRegion |> 
                rename("x1"="3R_A", "y1"="3R_B") |> 
                mutate(x2=x1+resolution-2, y2=y1+resolution-2) |> 
                filter( ( ( (x1-14308801) * (x2-14304362) ) < 0) | 
                        ( ( (y1-14308801) * (y2-14304362) ) < 0)
                        )

# Define a GRanges object for a specific region
my_roi_all.gr <- GRanges(seqnames = "chr3R",
                         ranges = IRanges(start = c(contacts.df$x1, contacts.df$y1), 
                                          end = c(contacts.df$x2, contacts.df$y2)))

mcols(my_roi_all.gr)$counts <- rep(contacts.df$counts, n=2)

# Merge overlapping ranges
my_roi_all_merged.gr <- reduce(my_roi_all.gr)

# Find overlaps between my_regions_of_interest.gr and my_regions_of_interest_merged.gr
overlaps <- findOverlaps(my_roi_all.gr, my_roi_all_merged.gr, type = "any", ignore.strand = TRUE)

# Aggregate metadata for the merged ranges (example: concatenate gene names)
# Create an empty list to store aggregated metadata
aggregated_metadata <- list()
for (i in seq_along(my_roi_all_merged.gr)) {
  # Get original ranges that overlap with the current reduced range
  original_indices <- queryHits(overlaps)[subjectHits(overlaps) == i]
  
  # Concatenate gene names
  aggregated_metadata[[i]] <- sum(my_roi_all.gr$counts[original_indices])
}

# Add aggregated metadata to the reduced (merged) GRanges object
mcols(my_roi_all_merged.gr)$counts_sum <- unlist(aggregated_metadata)

# delete my_roi_all_merged.gr intervals overlapping coords 14304362-14308801 spanning SuHw gene
# my_roi_all_merged.gr <- my_roi_all_merged.gr[(end(my_roi_all_merged.gr) < 14304362 | start(my_roi_all_merged.gr) > 14308801)]

# Convert the GRanges object to a data frame for plotting (bed format)
# If you want to use log10 scale for counts, uncomment the next line

# my_roi_all_merged.df <- as.data.frame(my_roi_all_merged.gr)[, 1:3] |>
#   mutate(score = if_else(my_roi_all_merged.gr$counts_sum == 0, 0, 
#                          log10(mcols(my_roi_all_merged.gr)$counts_sum)), strand = ".") |>
#   filter( ((start - 14308801) * (end - 14304362)) > 0)

my_roi_all_merged.df <- as.data.frame(my_roi_all_merged.gr)[, 1:3] |>
  mutate(score = if_else(my_roi_all_merged.gr$counts_sum == 0, 0, 
                         mcols(my_roi_all_merged.gr)$counts_sum), strand = ".") # |> filter( ((start - 14308801) * (end - 14304362)) > 0)

# Save Virtual4C HiC data as bed file for plotting with other programs
write_csv(x = my_roi_all_merged.df |> mutate(score = 10**score ), file = paste0(WD,"/R/render/virtual-4C-analysis_files/Tables/hic_contacts-",sample,"-",resolution,".bed"), 
          col_names = FALSE)


# Subset hiccup results to region of interest containing bait start-end -----------------------
tryCatch({
  my_hiccup_files[[my_key]]
}, error = function(e) {
  stop(paste("Error: HICCUP files for sample", my_key, "not found."))
})

tryCatch({
  # Fetch contact info from hic file just for bins overlapping region of interest
  # Region of interest is defined as the DpnII coord closest to SuHw gene boundaries.
  bedpe <- my_hiccup_files[[my_key]] |> 
    filter(  ((((x1-14308801) * (x2-14304362)) < 0 & chr1 == "chr3R")| 
                (((y1-14308801) * (y2-14304362)) < 0 & chr2 == "chr3R")) & 
               fdrBL <= 0.05 & 
               fdrDonut <= 0.05 & 
               fdrH <= 0.05 & 
               fdrV <= 0.05) |> 
    mutate(interaction_distance = (centroid1 - centroid2) / 1000) |> 
    mutate(archight = interaction_distance / max(interaction_distance)) |> 
    arrange(centroid2)
  
}, error = function(e) {
  warning(paste0("The following error accurred:",e$message))
  return(NULL)
}, warning = function(w) {
  warning(paste0("The following warning accurred:",w$message))
  return(NULL)
})

if (nrow(bedpe) == 0) {
  warning("***** No significant interactions found in the specified region. *****")
  return(NULL)
}

print(paste("***** Resolution", resolution, "****"))

# Extract significant contacts within regions of interest


# Set coordinates of genomic region to be plotted.
region_start <- min(c(bedpe$x1, bedpe$y1)) - 50000
region_end <- max(c(bedpe$x2, bedpe$y2)) + 50000

# Extract genes overlapping HICCUP loop coordinates for region of interest  ------------------
hiccup_contacts.df <- data.frame(seqnamess = "chr3R",
                                start = unique(c(bedpe$x1, bedpe$y1)),
                                end = unique(c(bedpe$x1, bedpe$y1)) + resolution - 1,
                                score = 10,
                                strand = ".")
hiccup_contacts.df <- hiccup_contacts.df[order(hiccup_contacts.df$start), ] # Sorting intervals, otherwise some are not plotted.

# Define a GRanges object for a specific region
my_regions_of_interest.gr <- GRanges(seqnames = "chr3R",
                                  ranges = IRanges(start = c(bedpe$x1, bedpe$y1), 
                                                   end = c(bedpe$x2, bedpe$y2)))


# Extract all gene ranges from the TxDb object
all_genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Find overlaps between significant chromosomal contact coordinates and the gene ranges
overlaps <- findOverlaps(my_regions_of_interest.gr, all_genes, type = "any", ignore.strand = TRUE)

# Extract the gene IDs of the overlapping genes
overlapping_gene_ids <- all_genes$gene_id[subjectHits(overlaps)]

# Convert flybase IDs into symbols
overlapping_gene_symbols <- AnnotationDbi::select(org.Dm.eg.db, keys=overlapping_gene_ids, columns = c("SYMBOL"), keytype = "FLYBASE") |> unique()

# Make 2-column df with genes to highlight in plot below
my_interacting_genes <- data.frame(gene = overlapping_gene_symbols$SYMBOL, color = ifelse(overlapping_gene_symbols$SYMBOL == "su(Hw)", "blue","darkred"))

# move row containing gene su(Hw) to the first row in my_interacting_genes
my_interacting_genes <- my_interacting_genes[order(my_interacting_genes$gene != "su(Hw)"),]


# Search regions of interest (Contacts) that overlap with SuHw_ChIRP_conc.gr ------------------
SuHw_ChIRP_overlap <- findOverlaps(SuHw_ChIRP_conc.gr, my_regions_of_interest.gr, type = "any", ignore.strand = TRUE, maxgap = 0)
SuHw_ChIRP_conc.overlap.gr <- SuHw_ChIRP_conc.gr[from(SuHw_ChIRP_overlap) |> unique(),]
SuHw_ChIRP_conc.overlap.df <- as.data.frame(SuHw_ChIRP_conc.overlap.gr)[,1:3] |> mutate(score=10, strand = '.')

my_regions_of_interest_overlap.gr <- my_regions_of_interest.gr[to(SuHw_ChIRP_overlap) |> unique(),]
my_regions_of_interest_overlap.start <- start(my_regions_of_interest_overlap.gr) |> unique()
bedpe <- bedpe |> mutate(overlap = if_else((x1 %in% my_regions_of_interest_overlap.start & 
                                              ((x1-(14308589+resolution)) * (x2-(14304455-resolution))) > 0) | 
                                            (y1 %in% my_regions_of_interest_overlap.start & 
                                               ((y1-(14308589+resolution)) * (y2-(14304455-resolution))) > 0) , TRUE, FALSE))


# Subet hicDataChromRegion and my_roi_all_merged.df to region of interest
hicDataChromRegion <- hicDataChromRegion[hicDataChromRegion$`3R_A` >= region_start & hicDataChromRegion$`3R_A` <= region_end & 
                                           hicDataChromRegion$`3R_B` >= region_start & hicDataChromRegion$`3R_B` <= region_end,]

my_roi_all_merged.df <- my_roi_all_merged.df[my_roi_all_merged.df$start >= region_start & 
                                              my_roi_all_merged.df$end <= region_end,]

################################
# Make virtual-4C plot
################################

dir.create(paste0("~/data/ELISSA_LEI/TK_208/R/render/virtual-4C-analysis_files/Plots/"), showWarnings = FALSE, recursive = TRUE)

  pdf(file = paste0("~/data/ELISSA_LEI/TK_208/R/render/virtual-4C-analysis_files/Plots/virtual4c-",s,"-",resolution,".pdf"), height = 8, width = 16)
  #png(file = paste0("./virtual-4C-analysis_files/Plots/virtual4c-",s,"-",resolution,".png"), height = 4600, width = 4600, res=300)

  y_base = -0.5
  
  # 1- create page -------------------------------
  pageCreate(width = 16, height = 8, default.units = "inches", ygrid = 0, xgrid = 0, showGuides = FALSE)
  
  # 2- Plot Hi-C matrix data in region -------------------------------
  
  # hicPlot <- plotHicTriangle(
  #   data = hicDataChromRegion, colorTrans = "linear",
  #   chrom = "chr3R",
  #   chromstart = region_start,
  #   chromend = region_end,
  #   assembly = Dm6,
  #   x = 0.5, y = y_base + 11, width = 15, height = 8.5,
  #   just = c("left","bottom"), default.units = "inches",
  #   zrange = c(1,quantile(hicDataChromRegion$counts, probs=c(.95)))
  # )
  # 
  
  # 3- Plot all Virtual 4C contacts, including N.S.  -------------------------------
  plotText(
    label = paste("Virtual 4C Su(Hw), HiC resolution =", resolution,"bp"), fontsize = 8,
    x = 0.5, y = y_base + 1, just = "left",
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
    flip = FALSE,
    archHeight = bedpe$archight,
    alpha = 1,
    x = 0.5, y = y_base + 1, width = 15, height = 1,
    just = c("left", "top"),
    default.units = "inches", baseline = FALSE, lwd = 0.2
  )
  
  plotSignal( # histogram of contact scores
    data = my_roi_all_merged.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 2.0, width = 15, height = 1, linecolor = "#6600cc", fill = "#6600cc", ymax = 0.001,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  plotSignal( # contact sites
    data = hiccup_contacts.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 3.1, width = 15, height = 0.3, linecolor = "#6600cc", fill = "#6600cc", 
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  
  # 4- Plot signal track data Mock SuHw overlap -------------------------------
  
  plotText(
    label = "GSE200993 Mock SuHw overlap", fontsize = 8,
    x = 0.5, y = y_base + 3.5, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = GSE200993_Mock_SuHw.overlap.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 3.6, width = 15, height = 0.3, linecolor = "#33ccff", fill = "#33ccff", 
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 5- Plot signal track data SuHw ChIRP -------------------------------
  plotText(
    label = "SuHw ChIRP", fontsize = 8,
    x = 0.5, y = y_base + 4.0, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = SuHw_ChIRP_conc.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 4.1, width = 15, height = 0.3, linecolor = "#ff9900", fill = "#ff9900", ymax = 0.5,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # Overlap with regions of interest
  plotSignal(
    data = SuHw_ChIRP_conc.overlap.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 4.1, width = 15, height = 0.3, linecolor = "blue3",
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8, baseline = FALSE,
  )
  
  # Add extra orange line to mask blue baseline fro overlap. baseline = FALSE parameter does not work.
  plotSignal(
    data = SuHw_ChIRP_conc.overlap.df |> mutate(score = 0),
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 4.1, width = 15, height = 0.3, linecolor = "#ff9900",
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 6- Plot signal track data pool even ChIRP -------------------------------
  plotText(
    label = "Pool even ChIRP", fontsize = 8,
    x = 0.5, y = y_base + 4.5, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = Pool_even.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 4.6, width = 15, height = 0.3, linecolor = "#ff9900", fill = "#ff9900", ymax = 0.5,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 7- Plot signal track data pool odd ChIRP -------------------------------
  plotText(
    label = "Pool odd ChIRP", fontsize = 8,
    x = 0.5, y = y_base + 5.0, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = Pool_odd.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 5.1, width = 15, height = 0.3, linecolor = "#ff9900", fill = "#ff9900", ymax = 0.5,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 8- Plot signal track data pool3 ChIRP -------------------------------
  plotText(
    label = "Pool3 ChIRP", fontsize = 8,
    x = 0.5, y = y_base + 5.5, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = Pool3.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 5.6, width = 15, height = 0.3, linecolor = "#ff9900", fill = "#ff9900", ymax = 0.5,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 9- Plot signal track data pool4 ChIRP -------------------------------
  plotText(
    label = "Pool4 ChIRP", fontsize = 8,
    x = 0.5, y = y_base + 6.0, just = "left",
    default.units = "inches",
  )
  
  plotSignal(
    data = Pool4.df,
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 6.1, width = 15, height = 0.3, linecolor = "#ff9900", fill = "#ff9900", ymax = 0.5,
    just = c("left", "top"), default.units = "inches", scale = TRUE, fontsize = 8,
  )
  
  # 10- Plot genes -------------------------------
  plotGenes(
    chrom = "chr3R", chromstart = region_start, chromend = region_end,
    assembly = Dm6,
    x = 0.5, y = y_base + 6.5, width = 15, height = 1, 
    just = c("left", "top"), default.units = "inches", stroke = 0.05, fontsize = 6,
    geneHighlights = my_interacting_genes, geneOrder = my_interacting_genes$gene
  )
  
  # 11- Add genome label -------------------------------
  annoGenomeLabel(
    plot = hicPlot, scale = "Mb",
    x = 0.5, y = y_base + 7.6
  )
   dev.off()
  
  return(bedpe)
}
