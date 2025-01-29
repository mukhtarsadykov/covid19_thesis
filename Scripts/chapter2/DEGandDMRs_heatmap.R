# Read the data
degs <- read.csv("~/Desktop/tmp.csv")

# Define the gene symbols to filter
gene_symbols <- sort(c("PDCD1", "IRF7", "PTGES", "IFITM1", "LITAF", "SPATA2", "STAP2", 
                  "S1PR5", "MAPK3", "CFAP92", "SAXO3", "EML3", "ARPC5L", "FLII", 
                  "SEPTIN9", "PKP3", "EPB41L3", "AHNAK", "AHNAK2", "PPL", "DNM1", 
                  "LIMK2", "TPM4", "BIN2", "NEIL1", "PARP10", "SLC7A5", "ENG", 
                  "ELF3", "EPHA2", "HSPB9", "BRD3", "EPAS1", "PTK6", "PRKCZ", 
                  "HIC1", "ARHGAP22", "ARHGAP45"))

# Filter the data for the relevant GeneSymbols
filtered_genes <- degs %>%
  filter(geneSymbol %in% gene_symbols) %>%
  select(geneSymbol, GeneID)%>%
  distinct()%>%
  arrange(geneSymbol) 


# show two heatmaps DEGs and DMRs

### DNAm ###

# from DMRicher coverage 1, 11 samples
setwd("~/Desktop/phd/Thesis/ch2/dnam/")
dmrs <- read.table("dmricher/DMRs/DMR_individual_smoothed_methylation.txt", header = T, sep = "\t")
#rownames <- dmrs$GeneID
rownames <- dmrs$geneSymbol

#rownames_anotation <- cbind(dmrs$GeneID, dmrs$geneSymbol)
#rownames_anotation <- data.frame(rownames_anotation)
#names(rownames_anotation) <- c("GeneID", "geneSymbol")
#row.names(rownames_anotation) <- rownames_anotation$GeneID

dmrs <- dmrs[,20:length(dmrs)]
names(dmrs) <- c("ICU_01",	"ICU_02",	"ICU_03",	"ICU_04",	"NON_01",	"NON_02",	"NON_03",	"NON_04",	"NON_05",	"NON_06",	"NON_07")
row.names(dmrs) <- rownames

study <- read.csv("~/Desktop/phd/Thesis/ch2/meta_11.csv", row.names=1) # dim is 80 x 19

study <- study[1:11,]

study$AltName <- c("ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02",
                   "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")

save_pheatmap_pdf <- function(x, filename, width=9, height=9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Reodered study according to heatmap clustering
#study <- read.csv("~/Desktop/rnaseq/last/raw/80s_new_metadata/meta_80s.csv", row.names = 1)

mycol <- colorRampPalette(c("blue","white","red"))(20)

#annotation info
annotation <- data.frame(study$AltName, study$Outcome,
                         study$Age, study$Hospital, study$ICU, study$Gender, study$Comorbidities, row.names = 1)
colnames(annotation) <- gsub("study.", "", colnames(annotation))
rownames(annotation) <- colnames(lcpm) # check out the row names of annotation
##colors
Age_col <- as.character(wes_palette("IsleofDogs2",  
                                    40, 
                                    type = "continuous"))
names(Age_col) <- unique(study$Age)
Hospital_col <- as.character(wes_palette("Chevalier1",  
                                         length(unique(study$Hospital)), 
                                         type = "continuous"))
names(Hospital_col) <- unique(study$Hospital)
libsize_col <- as.character(wes_palette("IsleofDogs2",  
                                        76, 
                                        type = "continuous"))
names(libsize_col) <- unique(study$lib.size)

ann_colors = list(
  ICU = c(Yes="#972D15", No ="#D8B70B"),
  Comorbidities = c(none="#E1BD6D", one="#9FA682", two="#274150", three="#F2300F"),
  #Nationality = c(Saudi="#F1BB7B", Yemeni="#FD6467", Other="#5B1A18"),
  Gender = c(Female = "#D783FF", Male = "#00DED8"),
  Outcome = c(Dead="#972D15", Alive ="#D8B70B"),
  Age = Age_col,
  Hospital = Hospital_col
  #Library = c(low="#972D15", high ="#D8B70B")
  #lib.size = libsize_col
)
annotation <- annotation[, c("ICU", "Outcome")]

#save_pheatmap_pdf(heatmap1, "~/Desktop/phd/Thesis/ch2/dnam/dmricher/Heatmap_annotation_11_ordered.pdf")

# Reorder the indices to match the alphabetical order of filtered_genes$GeneSymbol
ordered_genes <- filtered_genes %>%
  arrange(geneSymbol)


# Get the row indices based on the ordered GeneSymbol
i <- which(row.names(dmrs) %in% ordered_genes$GeneSymbol)

# Reorder dmrs matrix based on the order of filtered_genes$GeneSymbol
dmrs_ordered <- dmrs[match(ordered_genes$geneSymbol, row.names(dmrs)), ]

# Create heatmap with rows in alphabetical order
heatmap1 <- pheatmap(dmrs_ordered, 
                     cluster_cols = F, 
                     cluster_rows = F,
                     scale = "row",
                     color = colorRampPalette(c("blue", "white", "red"))(20),
                     annotation_col = annotation, 
                     labels_col = study$AltName, 
                     annotation_colors = ann_colors,
                     treeheight_row = 0, show_rownames =0
)

### 38genes in DEG ##

deg_counts <- read.csv("~/Desktop/phd/Thesis/ch2/batch_normalized_11_newNames.csv", row.names = 1)
deg_counts <- deg_counts[,1:11]
ordered_genes <- filtered_genes %>%
  arrange(geneSymbol)


# Get the row indices based on the ordered GeneSymbol
i <- which(row.names(deg_counts) %in% ordered_genes$geneSymbol)

# Reorder dmrs matrix based on the order of filtered_genes$GeneSymbol
dmrs_ordered <- deg_counts[match(ordered_genes$GeneID, row.names(deg_counts)), ]

# Create heatmap with rows in alphabetical order
heatmap1 <- pheatmap(dmrs_ordered, 
                     cluster_cols = F, 
                     cluster_rows = F,
                     scale = "row",
                     color = colorRampPalette(c("blue", "white", "red"))(20),
                     annotation_col = annotation, 
                     labels_col = study$AltName, 
                     annotation_colors = ann_colors,
                     treeheight_row = 0, show_rownames =0
)

tt_sig <- read.csv('../topTable_11_all.csv')

tt_sig_38genes <- tt_sig %>%
  filter(SYMBOL %in% gene_symbols) %>%
  select(ENSGID, SYMBOL, ICUvsNON_LFC, ICUvsNON_sig0, ICUvsNON_sig1, ICUvsNON_sig2)%>%
  arrange(SYMBOL) 
